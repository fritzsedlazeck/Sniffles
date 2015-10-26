// ***************************************************************************
// bamtools_split.cpp (c) 2010 Derek Barnett, Erik Garrison
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 24 July 2013 (DB)
// ---------------------------------------------------------------------------
// Splits a BAM file on user-specified property, creating a new BAM output
// file for each value found
// ***************************************************************************

#include "bamtools_split.h"

#include <api/BamConstants.h>
#include <api/BamReader.h>
#include <api/BamWriter.h>
#include <utils/bamtools_options.h>
#include <utils/bamtools_variant.h>
using namespace BamTools;

#include <ctime>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>
using namespace std;

namespace BamTools {
  
// string constants
static const string SPLIT_MAPPED_TOKEN    = ".MAPPED";
static const string SPLIT_UNMAPPED_TOKEN  = ".UNMAPPED";
static const string SPLIT_PAIRED_TOKEN    = ".PAIRED_END";
static const string SPLIT_SINGLE_TOKEN    = ".SINGLE_END";
static const string SPLIT_REFERENCE_TOKEN = ".REF_";
static const string SPLIT_TAG_TOKEN       = ".TAG_";

string GetTimestampString(void) {

    // get human readable timestamp
    time_t currentTime;
    time(&currentTime);
    stringstream timeStream("");
    timeStream << ctime(&currentTime);

    // convert whitespace to '_'
    string timeString = timeStream.str();
    size_t found = timeString.find(" ");
    while (found != string::npos) {
        timeString.replace(found, 1, "_");
        found = timeString.find(" ", found+1);
    }
    return timeString;
}

// remove copy of filename without extension
// (so /path/to/file.txt becomes /path/to/file )
string RemoveFilenameExtension(const string& filename) {
    size_t found = filename.rfind(".");
    return filename.substr(0, found);
}
    
} // namespace BamTools

// ---------------------------------------------
// SplitSettings implementation

struct SplitTool::SplitSettings {

    // flags
    bool HasInputFilename;
    bool HasCustomOutputStub;
    bool HasCustomRefPrefix;
    bool HasCustomTagPrefix;
    bool IsSplittingMapped;
    bool IsSplittingPaired;
    bool IsSplittingReference;
    bool IsSplittingTag;
    
    // string args
    string CustomOutputStub;
    string CustomRefPrefix;
    string CustomTagPrefix;
    string InputFilename;
    string TagToSplit;
    
    // constructor
    SplitSettings(void)
        : HasInputFilename(false)
        , HasCustomOutputStub(false)
        , HasCustomRefPrefix(false)
        , HasCustomTagPrefix(false)
        , IsSplittingMapped(false)
        , IsSplittingPaired(false)
        , IsSplittingReference(false)
        , IsSplittingTag(false)
        , CustomOutputStub("")
        , CustomRefPrefix("")
        , CustomTagPrefix("")
        , InputFilename(Options::StandardIn())
        , TagToSplit("")
    { } 
};  

// ---------------------------------------------
// SplitToolPrivate declaration

class SplitTool::SplitToolPrivate {
      
    // ctor & dtor
    public:
        SplitToolPrivate(SplitTool::SplitSettings* settings)
            : m_settings(settings)
        { }

        ~SplitToolPrivate(void) {
            m_reader.Close();
        }
        
    // 'public' interface
    public:
        bool Run(void);
        
    // internal methods
    private:
        // close & delete BamWriters in map
        template<typename T>
        void CloseWriters(map<T, BamWriter*>& writers);
        // calculate output stub based on IO args given
        void DetermineOutputFilenameStub(void);
        // open our BamReader
        bool OpenReader(void);
        // split alignments in BAM file based on isMapped property
        bool SplitMapped(void);
        // split alignments in BAM file based on isPaired property
        bool SplitPaired(void);
        // split alignments in BAM file based on refID property
        bool SplitReference(void);
        // finds first alignment and calls corresponding SplitTagImpl<> 
        // depending on tag type
        bool SplitTag(void);
        // templated split tag implementation 
        // handle the various types that are possible for tags
        template<typename T>
        bool SplitTagImpl(BamAlignment& al);    
        
    // data members
    private:
        SplitTool::SplitSettings* m_settings;
        string m_outputFilenameStub;
        BamReader m_reader;
        string m_header;
        RefVector m_references;
};

void SplitTool::SplitToolPrivate::DetermineOutputFilenameStub(void) {
  
    // if user supplied output filename stub, use that
    if ( m_settings->HasCustomOutputStub ) 
        m_outputFilenameStub = m_settings->CustomOutputStub;
    
    // else if user supplied input BAM filename, use that (minus ".bam" extension) as stub
    else if ( m_settings->HasInputFilename )
        m_outputFilenameStub = RemoveFilenameExtension(m_settings->InputFilename);
        
    // otherwise, user did not specify -stub, and input is coming from STDIN
    // generate stub from timestamp
    else m_outputFilenameStub = GetTimestampString();      
}

bool SplitTool::SplitToolPrivate::OpenReader(void) {

    // attempt to open BAM file
    if ( !m_reader.Open(m_settings->InputFilename) ) {
        cerr << "bamtools split ERROR: could not open BAM file: " << m_settings->InputFilename << endl;
        return false;
    }

    // save file 'metadata' & return success
    m_header     = m_reader.GetHeaderText();
    m_references = m_reader.GetReferenceData();
    return true;
}

bool SplitTool::SplitToolPrivate::Run(void) {
  
    // determine output stub
    DetermineOutputFilenameStub();

    // open up BamReader
    if ( !OpenReader() )
        return false;
    
    // determine split type from settings
    if ( m_settings->IsSplittingMapped )    return SplitMapped();
    if ( m_settings->IsSplittingPaired )    return SplitPaired();
    if ( m_settings->IsSplittingReference ) return SplitReference();
    if ( m_settings->IsSplittingTag )       return SplitTag();

    // if we get here, no property was specified 
    cerr << "bamtools split ERROR: no property given to split on... " << endl
         << "Please use -mapped, -paired, -reference, or -tag TAG to specifiy desired split behavior." << endl;
    return false;
}    

bool SplitTool::SplitToolPrivate::SplitMapped(void) {
    
    // set up splitting data structure
    map<bool, BamWriter*> outputFiles;
    map<bool, BamWriter*>::iterator writerIter;
    
    // iterate through alignments
    BamAlignment al;
    BamWriter* writer;
    bool isCurrentAlignmentMapped;
    while ( m_reader.GetNextAlignment(al) ) {
      
        // see if bool value exists
        isCurrentAlignmentMapped = al.IsMapped();
        writerIter = outputFiles.find(isCurrentAlignmentMapped);
          
        // if no writer associated with this value
        if ( writerIter == outputFiles.end() ) {
        
            // open new BamWriter
            const string outputFilename = m_outputFilenameStub + ( isCurrentAlignmentMapped
                                                                  ? SPLIT_MAPPED_TOKEN
                                                                  : SPLIT_UNMAPPED_TOKEN ) + ".bam";
            writer = new BamWriter;
            if ( !writer->Open(outputFilename, m_header, m_references) ) {
                cerr << "bamtools split ERROR: could not open " << outputFilename
                     << " for writing." << endl;
                return false;
            }
          
            // store in map
            outputFiles.insert( make_pair(isCurrentAlignmentMapped, writer) );
        } 
        
        // else grab corresponding writer
        else writer = (*writerIter).second;
        
        // store alignment in proper BAM output file 
        if ( writer )
            writer->SaveAlignment(al);
    }
    
    // clean up BamWriters 
    CloseWriters(outputFiles);
    
    // return success
    return true;
}

bool SplitTool::SplitToolPrivate::SplitPaired(void) {
  
    // set up splitting data structure
    map<bool, BamWriter*> outputFiles;
    map<bool, BamWriter*>::iterator writerIter;
    
    // iterate through alignments
    BamAlignment al;
    BamWriter* writer;
    bool isCurrentAlignmentPaired;
    while ( m_reader.GetNextAlignment(al) ) {
      
        // see if bool value exists
        isCurrentAlignmentPaired = al.IsPaired();
        writerIter = outputFiles.find(isCurrentAlignmentPaired);
          
        // if no writer associated with this value
        if ( writerIter == outputFiles.end() ) {
        
            // open new BamWriter
            const string outputFilename = m_outputFilenameStub + ( isCurrentAlignmentPaired
                                                                  ? SPLIT_PAIRED_TOKEN
                                                                  : SPLIT_SINGLE_TOKEN ) + ".bam";
            writer = new BamWriter;
            if ( !writer->Open(outputFilename, m_header, m_references) ) {
                cerr << "bamtool split ERROR: could not open " << outputFilename
                     << " for writing." << endl;
                return false;
            }
          
            // store in map
            outputFiles.insert( make_pair(isCurrentAlignmentPaired, writer) );
        } 
        
        // else grab corresponding writer
        else writer = (*writerIter).second;
        
        // store alignment in proper BAM output file 
        if ( writer ) 
            writer->SaveAlignment(al);
    }
    
    // clean up BamWriters 
    CloseWriters(outputFiles);
    
    // return success
    return true;  
}

bool SplitTool::SplitToolPrivate::SplitReference(void) {
  
    // set up splitting data structure
    map<int32_t, BamWriter*> outputFiles;
    map<int32_t, BamWriter*>::iterator writerIter;
    
    // determine reference prefix
    string refPrefix = SPLIT_REFERENCE_TOKEN;
    if ( m_settings->HasCustomRefPrefix )
        refPrefix = m_settings->CustomRefPrefix;

    // make sure prefix starts with '.'
    const size_t dotFound = refPrefix.find('.');
    if ( dotFound != 0 )
        refPrefix = string(".") + refPrefix;

    // iterate through alignments
    BamAlignment al;
    BamWriter* writer;
    int32_t currentRefId;
    while ( m_reader.GetNextAlignment(al) ) {
      
        // see if bool value exists
        currentRefId = al.RefID;
        writerIter = outputFiles.find(currentRefId);
          
        // if no writer associated with this value
        if ( writerIter == outputFiles.end() ) {
        
            // fetch reference name for ID
            string refName;
            if ( currentRefId == -1 )
                refName = "unmapped";
            else
                refName = m_references.at(currentRefId).RefName;

            // construct new output filename
            const string outputFilename = m_outputFilenameStub + refPrefix + refName + ".bam";

            // open new BamWriter
            writer = new BamWriter;
            if ( !writer->Open(outputFilename, m_header, m_references) ) {
                cerr << "bamtools split ERROR: could not open " << outputFilename
                     << " for writing." << endl;
                return false;
            }

            // store in map
            outputFiles.insert( make_pair(currentRefId, writer) );
        } 
        
        // else grab corresponding writer
        else writer = (*writerIter).second;
        
        // store alignment in proper BAM output file 
        if ( writer ) 
            writer->SaveAlignment(al);
    }
    
    // clean up BamWriters 
    CloseWriters(outputFiles);
    
    // return success
    return true;
}

// finds first alignment and calls corresponding SplitTagImpl<>() depending on tag type
bool SplitTool::SplitToolPrivate::SplitTag(void) {  
  
    // iterate through alignments, until we hit TAG
    BamAlignment al;
    while ( m_reader.GetNextAlignment(al) ) {
      
        // look for tag in this alignment and get tag type
        char tagType(0);
        if ( !al.GetTagType(m_settings->TagToSplit, tagType) )
            continue;
        
        // request split method based on tag type
        // pass it the current alignment found
        switch ( tagType ) {
          
            case (Constants::BAM_TAG_TYPE_INT8)  :
            case (Constants::BAM_TAG_TYPE_INT16) :
            case (Constants::BAM_TAG_TYPE_INT32) :
                return SplitTagImpl<int32_t>(al);
                
            case (Constants::BAM_TAG_TYPE_UINT8)  :
            case (Constants::BAM_TAG_TYPE_UINT16) :
            case (Constants::BAM_TAG_TYPE_UINT32) :
                return SplitTagImpl<uint32_t>(al);
              
            case (Constants::BAM_TAG_TYPE_FLOAT)  :
                return SplitTagImpl<float>(al);
            
            case (Constants::BAM_TAG_TYPE_ASCII)  :
            case (Constants::BAM_TAG_TYPE_STRING) :
            case (Constants::BAM_TAG_TYPE_HEX)    :
                return SplitTagImpl<string>(al);

            case (Constants::BAM_TAG_TYPE_ARRAY) :
                cerr << "bamtools split ERROR: array tag types are not supported" << endl;
                return false;
          
            default:
                cerr << "bamtools split ERROR: unknown tag type encountered: " << tagType << endl;
                return false;
        }
    }
    
    // tag not found, but that's not an error - return success
    return true;
}

// --------------------------------------------------------------------------------
// template method implementation
// *Technical Note* - use of template methods declared & defined in ".cpp" file
//                    goes against normal practices, but works here because these
//                    are purely internal (no one can call from outside this file)

// close BamWriters & delete pointers
template<typename T>
void SplitTool::SplitToolPrivate::CloseWriters(map<T, BamWriter*>& writers) {
  
    typedef map<T, BamWriter*> WriterMap;
    typedef typename WriterMap::iterator WriterMapIterator;
  
    // iterate over writers
    WriterMapIterator writerIter = writers.begin();
    WriterMapIterator writerEnd  = writers.end();
    for ( ; writerIter != writerEnd; ++writerIter ) {
        BamWriter* writer = (*writerIter).second;
        if ( writer == 0 ) continue;

        // close BamWriter
        writer->Close();

        // destroy BamWriter
        delete writer;
        writer = 0;
    }

    // clear the container (destroying the items doesn't remove them)
    writers.clear();
}

// handle the various types that are possible for tags
template<typename T>
bool SplitTool::SplitToolPrivate::SplitTagImpl(BamAlignment& al) {
  
    typedef T TagValueType;
    typedef map<TagValueType, BamWriter*> WriterMap;
    typedef typename WriterMap::iterator WriterMapIterator;
  
    // set up splitting data structure
    WriterMap outputFiles;
    WriterMapIterator writerIter;

    // determine tag prefix
    string tagPrefix = SPLIT_TAG_TOKEN;
    if ( m_settings->HasCustomTagPrefix )
        tagPrefix = m_settings->CustomTagPrefix;

    // make sure prefix starts with '.'
    const size_t dotFound = tagPrefix.find('.');
    if ( dotFound != 0 )
        tagPrefix = string(".") + tagPrefix;

    // local variables
    const string tag = m_settings->TagToSplit;
    BamWriter* writer;
    stringstream outputFilenameStream("");
    TagValueType currentValue;
    
    // retrieve first alignment tag value
    if ( al.GetTag(tag, currentValue) ) {
      
        // open new BamWriter, save first alignment
        outputFilenameStream << m_outputFilenameStub << tagPrefix << tag << "_" << currentValue << ".bam";
        writer = new BamWriter;
        if ( !writer->Open(outputFilenameStream.str(), m_header, m_references) ) {
            cerr << "bamtools split ERROR: could not open " << outputFilenameStream.str()
                 << " for writing." << endl;
            return false;
        }
        writer->SaveAlignment(al);
        
        // store in map
        outputFiles.insert( make_pair(currentValue, writer) );
        
        // reset stream
        outputFilenameStream.str("");
    }
    
    // iterate through remaining alignments
    while ( m_reader.GetNextAlignment(al) ) {
      
        // skip if this alignment doesn't have TAG 
        if ( !al.GetTag(tag, currentValue) ) continue;
        
        // look up tag value in map
        writerIter = outputFiles.find(currentValue);
          
        // if no writer associated with this value
        if ( writerIter == outputFiles.end() ) {
        
            // open new BamWriter
            outputFilenameStream << m_outputFilenameStub << tagPrefix << tag << "_" << currentValue << ".bam";
            writer = new BamWriter;
            if ( !writer->Open(outputFilenameStream.str(), m_header, m_references) ) {
                cerr << "bamtool split ERROR: could not open " << outputFilenameStream.str()
                     << " for writing." << endl;
                return false;
            }

            // store in map
            outputFiles.insert( make_pair(currentValue, writer) );
            
            // reset stream
            outputFilenameStream.str("");
        } 
        
        // else grab corresponding writer
        else writer = (*writerIter).second;
        
        // store alignment in proper BAM output file 
        if ( writer ) 
            writer->SaveAlignment(al);
    }
    
    // clean up BamWriters  
    CloseWriters(outputFiles);
    
    // return success
    return true;  
}

// ---------------------------------------------
// SplitTool implementation

SplitTool::SplitTool(void)
    : AbstractTool()
    , m_settings(new SplitSettings)
    , m_impl(0)
{
    // set program details
    const string name = "bamtools split";
    const string description = "splits a BAM file on user-specified property, creating a new BAM output file for each value found";
    const string args = "[-in <filename>] [-stub <filename stub>] < -mapped | -paired | -reference [-refPrefix <prefix>] | -tag <TAG> > ";
    Options::SetProgramInfo(name, description, args);
    
    // set up options 
    OptionGroup* IO_Opts = Options::CreateOptionGroup("Input & Output");
    Options::AddValueOption("-in",   "BAM filename",  "the input BAM file",  "", m_settings->HasInputFilename,  m_settings->InputFilename,  IO_Opts, Options::StandardIn());
    Options::AddValueOption("-refPrefix", "string", "custom prefix for splitting by references. Currently files end with REF_<refName>.bam. This option allows you to replace \"REF_\" with a prefix of your choosing.", "",
                            m_settings->HasCustomRefPrefix, m_settings->CustomRefPrefix, IO_Opts);
    Options::AddValueOption("-tagPrefix", "string", "custom prefix for splitting by tags. Current files end with TAG_<tagname>_<tagvalue>.bam. This option allows you to replace \"TAG_\" with a prefix of your choosing.", "",
                            m_settings->HasCustomTagPrefix, m_settings->CustomTagPrefix, IO_Opts);
    Options::AddValueOption("-stub", "filename stub", "prefix stub for output BAM files (default behavior is to use input filename, without .bam extension, as stub). If input is stdin and no stub provided, a timestamp is generated as the stub.", "",
                            m_settings->HasCustomOutputStub, m_settings->CustomOutputStub, IO_Opts);
    
    OptionGroup* SplitOpts = Options::CreateOptionGroup("Split Options");
    Options::AddOption("-mapped",    "split mapped/unmapped alignments",       m_settings->IsSplittingMapped,    SplitOpts);
    Options::AddOption("-paired",    "split single-end/paired-end alignments", m_settings->IsSplittingPaired,    SplitOpts);
    Options::AddOption("-reference", "split alignments by reference",          m_settings->IsSplittingReference, SplitOpts);
    Options::AddValueOption("-tag", "tag name", "splits alignments based on all values of TAG encountered (i.e. -tag RG creates a BAM file for each read group in original BAM file)", "", 
                            m_settings->IsSplittingTag, m_settings->TagToSplit, SplitOpts);
}

SplitTool::~SplitTool(void) {
    
    delete m_settings;
    m_settings = 0;
    
    delete m_impl;
    m_impl = 0;
}

int SplitTool::Help(void) {
    Options::DisplayHelp();
    return 0;
}

int SplitTool::Run(int argc, char* argv[]) {
  
    // parse command line arguments
    Options::Parse(argc, argv, 1);
    
    // initialize SplitTool with settings
    m_impl = new SplitToolPrivate(m_settings);
    
    // run SplitTool, return success/fail
    if ( m_impl->Run() ) 
        return 0;
    else 
        return 1;
}
