// ***************************************************************************
// BamHttp_p.cpp (c) 2011 Derek Barnett
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 24 July 2013 (DB)
// ---------------------------------------------------------------------------
// Provides reading/writing of BAM files on HTTP server
// ***************************************************************************

#include "api/BamAux.h"
#include "api/internal/io/BamHttp_p.h"
#include "api/internal/io/HttpHeader_p.h"
#include "api/internal/io/TcpSocket_p.h"
using namespace BamTools;
using namespace BamTools::Internal;

#include <cassert>
#include <cctype>
#include <cstdlib>
#include <algorithm>
#include <sstream>
using namespace std;

namespace BamTools {
namespace Internal {

// -----------
// constants
// -----------

static const string HTTP_PORT   = "80";
static const string HTTP_PREFIX = "http://";
static const size_t HTTP_PREFIX_LENGTH = 7;

static const string DOUBLE_NEWLINE = "\n\n";

static const string GET_METHOD   = "GET";
static const string HEAD_METHOD  = "HEAD";
static const string HOST_HEADER  = "Host";
static const string RANGE_HEADER = "Range";
static const string BYTES_PREFIX = "bytes=";
static const string CONTENT_LENGTH_HEADER = "Content-Length";

static const char HOST_SEPARATOR  = '/';
static const char PROXY_SEPARATOR = ':';

// -----------------
// utility methods
// -----------------

static inline
bool endsWith(const string& source, const string& pattern) {
    return ( source.find(pattern) == (source.length() - pattern.length()) );
}

static inline
string toLower(const string& s) {
    string out;
    const size_t sSize = s.size();
    out.reserve(sSize);
    for ( size_t i = 0; i < sSize; ++i )
        out[i] = tolower(s[i]);
    return out;
}

} // namespace Internal
} // namespace BamTools

// ------------------------
// BamHttp implementation
// ------------------------

BamHttp::BamHttp(const string& url)
    : IBamIODevice()
    , m_socket(new TcpSocket)
    , m_port(HTTP_PORT)
    , m_request(0)
    , m_response(0)
    , m_isUrlParsed(false)
    , m_filePosition(-1)
    , m_fileEndPosition(-1)
    , m_rangeEndPosition(-1)
{
    ParseUrl(url);
}

BamHttp::~BamHttp(void) {

    // close connection & clean up
    Close();
    if ( m_socket )
        delete m_socket;
}

void BamHttp::ClearResponse(void) {
    if ( m_response ) {
        delete m_response;
        m_response = 0;
    }
}

void BamHttp::Close(void) {

    // disconnect socket & clear related resources
    DisconnectSocket();

    // reset state
    m_isUrlParsed = false;
    m_filePosition     = -1;
    m_fileEndPosition  = -1;
    m_rangeEndPosition = -1;
    m_mode = IBamIODevice::NotOpen;
}

bool BamHttp::ConnectSocket(void) {

    BT_ASSERT_X(m_socket, "null socket?");

    // any state checks, etc?
    if ( !m_socket->ConnectToHost(m_hostname, m_port, m_mode) ) {
        SetErrorString("BamHttp::ConnectSocket", m_socket->GetErrorString());
        return false;
    }

    // return success
    return true;
}

void BamHttp::DisconnectSocket(void) {

    // disconnect socket & clean up
    m_socket->DisconnectFromHost();
    ClearResponse();
    if ( m_request )  {
        delete m_request;
        m_request = 0;
    }
}

bool BamHttp::EnsureSocketConnection(void) {
    if ( m_socket->IsConnected() )
        return true;
    return ConnectSocket();
}

bool BamHttp::IsOpen(void) const {
    return IBamIODevice::IsOpen() && m_isUrlParsed;
}

bool BamHttp::IsRandomAccess(void) const {
    return true;
}

bool BamHttp::Open(const IBamIODevice::OpenMode mode) {

    // BamHttp only supports read-only access
    if ( mode != IBamIODevice::ReadOnly ) {
        SetErrorString("BamHttp::Open", "writing on this device is not supported");
        return false;
    }
    m_mode = mode;

    // attempt connection to socket
    if ( !ConnectSocket() ) {
        SetErrorString("BamHttp::Open", m_socket->GetErrorString());
        return false;
    }

    // initialize our file positions
    m_filePosition     = 0;
    m_fileEndPosition  = 0;
    m_rangeEndPosition = 0;

    // attempt to send initial request (just 'HEAD' to check connection)
    if ( !SendHeadRequest() ) {
        SetErrorString("BamHttp::Open", m_socket->GetErrorString());
        return false;
    }

    // clear response from HEAD request, not needed
    ClearResponse();

    // return success
    return true;
}

void BamHttp::ParseUrl(const string& url) {

    // clear flag to start
    m_isUrlParsed = false;

    // make sure url starts with "http://", case-insensitive
    string tempUrl(url);
    toLower(tempUrl);
    const size_t prefixFound = tempUrl.find(HTTP_PREFIX);
    if ( prefixFound == string::npos )
        return;

    // find end of host name portion (first '/' hit after the prefix)
    const size_t firstSlashFound = tempUrl.find(HOST_SEPARATOR, HTTP_PREFIX_LENGTH);
    if ( firstSlashFound == string::npos ) {
        ;  // no slash found... no filename given along with host?
    }

    // fetch hostname (check for proxy port)
    string hostname = tempUrl.substr(HTTP_PREFIX_LENGTH, (firstSlashFound - HTTP_PREFIX_LENGTH));
    const size_t colonFound = hostname.find(PROXY_SEPARATOR);
    if ( colonFound != string::npos ) {
        ; // TODO: handle proxy port (later, just skip for now)
    } else {
        m_hostname = hostname;
        m_port = HTTP_PORT;
    }

    // store remainder of URL as filename (must be non-empty)
    string filename = tempUrl.substr(firstSlashFound);
    if ( filename.empty() )
        return;
    m_filename = filename;

    // set parsed OK flag
    m_isUrlParsed = true;
}

int64_t BamHttp::Read(char* data, const unsigned int numBytes) {

    // if BamHttp not in a valid state
    if ( !IsOpen() )
        return -1;

    int64_t numBytesReadSoFar = 0;
    while ( numBytesReadSoFar < numBytes ) {

        const size_t remaining = static_cast<size_t>( numBytes - numBytesReadSoFar );

        // if we're not holding a valid GET reponse, get one
        if ( m_response == 0 ) {
            if ( !SendGetRequest(remaining) )
                return -1;
        }
        BT_ASSERT_X(m_response, "null HTTP response");

        // check response status code
        const int statusCode = m_response->GetStatusCode();

        // if we receieved full file contents in response
        if ( statusCode == 200 ) {

            // try to read 'remaining' bytes from socket
            const int64_t socketBytesRead = ReadFromSocket(data+numBytesReadSoFar, remaining);

            // if error
            if ( socketBytesRead < 0 ) {
                SetErrorString("BamHttp::Read", m_socket->GetErrorString());
                return -1;
            }

            // EOF
            else if ( socketBytesRead == 0 )
                return numBytesReadSoFar;

            // update counters
            numBytesReadSoFar += socketBytesRead;
            m_filePosition    += socketBytesRead;

        }

        // else if we received a range of bytes in response
        else if ( statusCode == 206 ) {

            // if we've exhausted the last request
            if ( m_filePosition == m_rangeEndPosition ) {
                if ( !SendGetRequest(remaining) )
                    return -1;
            }

            else {

                // try to read 'remaining' bytes from socket
                const int64_t socketBytesRead = ReadFromSocket(data+numBytesReadSoFar, remaining);

                // if error
                if ( socketBytesRead < 0 ) {
                    SetErrorString("BamHttp::Read", m_socket->GetErrorString());
                    return -1;
                }

                // maybe EOF
                else if ( socketBytesRead == 0 ) {

                    // if we know we're not at end position, fire off a new request
                    if ( m_fileEndPosition > 0 && m_filePosition < m_fileEndPosition ) {
                        if ( !SendGetRequest() )
                            return -1;
                    } else
                        return numBytesReadSoFar;
                }

                // update counters
                numBytesReadSoFar += socketBytesRead;
                m_filePosition    += socketBytesRead;
            }
        }


        // else some other HTTP status
        else {
            SetErrorString("BamHttp::Read", "unsupported status code in response");
            return -1;
        }
    }

    // return actual number of bytes read
    return numBytesReadSoFar;
}

int64_t BamHttp::ReadFromSocket(char* data, const unsigned int maxNumBytes) {
    return m_socket->Read(data, maxNumBytes);
}

bool BamHttp::ReceiveResponse(void) {

    // fetch header, up until double new line
    string responseHeader;
    do {

        // make sure we can read a line
        if ( !m_socket->WaitForReadLine() )
            return false;

        // read line & append to full header
        const string headerLine = m_socket->ReadLine();
        responseHeader += headerLine;

    } while ( !endsWith(responseHeader, DOUBLE_NEWLINE) );

    // sanity check
    if ( responseHeader.empty() ) {
        SetErrorString("BamHttp::ReceiveResponse", "empty HTTP response");
        Close();
        return false;
    }

    // create response from header text
    m_response = new HttpResponseHeader(responseHeader);
    if ( !m_response->IsValid() ) {
        SetErrorString("BamHttp::ReceiveResponse", "could not parse HTTP response");
        Close();
        return false;
    }

    // if we get here, success
    return true;
}

bool BamHttp::Seek(const int64_t& position, const int origin) {

    // if HTTP device not in a valid state
    if ( !IsOpen() ) {
        SetErrorString("BamHttp::Seek", "cannot seek on unopen connection");
        return false;
    }

    // reset the connection
    DisconnectSocket();
    if ( !ConnectSocket() ) {
        SetErrorString("BamHttp::Seek", m_socket->GetErrorString());
        return false;
    }

    // udpate file position
    switch ( origin ) {
        case SEEK_CUR : m_filePosition += position; break;
        case SEEK_SET : m_filePosition  = position; break;
        default :
            SetErrorString("BamHttp::Seek", "unsupported seek origin");
            return false;
    }

    // return success
    return true;
}

bool BamHttp::SendGetRequest(const size_t numBytes) {

    // clear previous data
    ClearResponse();
    if ( m_request )
        delete m_request;
    m_socket->ClearBuffer();

    // make sure we're connected
    if ( !EnsureSocketConnection() )
        return false;

    // create range string
    const int64_t endPosition = m_filePosition + std::max(static_cast<size_t>(0x10000), numBytes);
    stringstream range("");
    range << BYTES_PREFIX << m_filePosition << '-' << endPosition;

    // create request
    m_request = new HttpRequestHeader(GET_METHOD, m_filename);
    m_request->SetField(HOST_HEADER,  m_hostname);
    m_request->SetField(RANGE_HEADER, range.str());

    // send request
    const string requestHeader = m_request->ToString();
    const size_t headerSize    = requestHeader.size();
    if ( WriteToSocket(requestHeader.c_str(), headerSize) != headerSize ) {
        SetErrorString("BamHttp::SendHeadRequest", m_socket->GetErrorString());
        return false;
    }

    // ensure clean buffer
    m_socket->ClearBuffer();

    // wait for response
    if ( !ReceiveResponse() ) {
        SetErrorString("BamHttp::SendGetRequest", m_socket->GetErrorString());
        Close();
        return false;
    }
    BT_ASSERT_X(m_response, "BamHttp::SendGetRequest : null HttpResponse");
    BT_ASSERT_X(m_response->IsValid(), "BamHttp::SendGetRequest : invalid HttpResponse");

    // check response status code
    const int statusCode = m_response->GetStatusCode();
    switch ( statusCode ) {

        // ranged response, as requested
        case 206 :
            // get content length if available
            if ( m_response->ContainsKey(CONTENT_LENGTH_HEADER) ) {
                const string contentLengthString = m_response->GetValue(CONTENT_LENGTH_HEADER);
                m_rangeEndPosition = m_filePosition + atoi( contentLengthString.c_str() );
            }
            return true;

        // full contents, not range
        case 200 :
        {
            // skip up to current file position
            RaiiBuffer tmp(0x8000);
            int64_t numBytesRead = 0;
            while ( numBytesRead < m_filePosition ) {

                // read data from response
                const int64_t remaining = m_filePosition - numBytesRead;
                const size_t bytesToRead = static_cast<size_t>( (remaining > 0x8000) ? 0x8000 : remaining );
                const int64_t socketBytesRead = ReadFromSocket(tmp.Buffer, bytesToRead);

                // if error
                if ( socketBytesRead < 0 ) {
                    SetErrorString("BamHttp::SendGetRequest", m_socket->GetErrorString());
                    Close();
                    return false;
                }

                // else if EOF
                else if ( socketBytesRead == 0 && m_socket->BufferBytesAvailable() == 0 )
                    break;

                // update byte counter
                numBytesRead += socketBytesRead;
            }

            // return success
            return ( numBytesRead == m_filePosition);
        }

        // any other status codes
        default:
            break;
    }

    // fail on unexpected status code
    SetErrorString("BamHttp::SendGetRequest", "unsupported status code in response");
    Close();
    return false;
}

bool BamHttp::SendHeadRequest(void) {

    // ensure clean slate
    ClearResponse();
    if ( m_request )
        delete m_request;
    m_socket->ClearBuffer();

    // make sure we're connected
    if ( !EnsureSocketConnection() )
        return false;

    // create request
    m_request = new HttpRequestHeader(HEAD_METHOD, m_filename);
    m_request->SetField(HOST_HEADER, m_hostname);

    // send request
    const string requestHeader = m_request->ToString();
    const size_t headerSize    = requestHeader.size();
    if ( WriteToSocket(requestHeader.c_str(), headerSize) != headerSize ) {
        SetErrorString("BamHttp::SendHeadRequest", m_socket->GetErrorString());
        return false;
    }

    m_socket->ClearBuffer();

    // wait for response from server
    if ( !ReceiveResponse() ) {
        SetErrorString("BamHttp::SendHeadRequest", m_socket->GetErrorString());
        Close();
        return false;
    }
    BT_ASSERT_X(m_response, "BamHttp::SendHeadRequest : null HttpResponse");
    BT_ASSERT_X(m_response->IsValid(), "BamHttp::SendHeadRequest : invalid HttpResponse");

    // get content length if available
    if ( m_response->ContainsKey(CONTENT_LENGTH_HEADER) ) {
        const string contentLengthString = m_response->GetValue(CONTENT_LENGTH_HEADER);
        m_fileEndPosition = atoi( contentLengthString.c_str() ) - 1;
    }

    // return whether we found any errors
    return m_socket->GetError() == TcpSocket::NoError;
}

int64_t BamHttp::Tell(void) const {
    return ( IsOpen() ? m_filePosition : -1 );
}

int64_t BamHttp::Write(const char* data, const unsigned int numBytes) {
    (void)data;
    (void)numBytes;
    BT_ASSERT_X(false, "BamHttp::Write : write-mode not supported on this device");
    SetErrorString("BamHttp::Write", "write-mode not supported on this device");
    return -1;
}

int64_t BamHttp::WriteToSocket(const char* data, const unsigned int numBytes) {
    if ( !m_socket->IsConnected() )
        return -1;
    m_socket->ClearBuffer();
    return m_socket->Write(data, numBytes);
}
