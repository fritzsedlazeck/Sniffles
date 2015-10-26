// ***************************************************************************
// bamtools_convert.cpp (c) 2010 Derek Barnett, Erik Garrison
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 10 December 2012
// ---------------------------------------------------------------------------
// Converts between BAM and a number of other formats
// ***************************************************************************

#include "bamtools_convert.h"

#include <api/BamConstants.h>
#include <api/BamMultiReader.h>
#include <utils/bamtools_fasta.h>
#include <utils/bamtools_options.h>
#include <utils/bamtools_pileup_engine.h>
#include <utils/bamtools_utilities.h>
using namespace BamTools;

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
using namespace std;
  
namespace BamTools { 
  
// ---------------------------------------------
// ConvertTool constants

// supported conversion format command-line names
static const string FORMAT_BED    = "bed";
static const string FORMAT_FASTA  = "fasta";
static const string FORMAT_FASTQ  = "fastq";
static const string FORMAT_JSON   = "json";
static const string FORMAT_SAM    = "sam";
static const string FORMAT_PILEUP = "pileup";
static const string FORMAT_YAML   = "yaml";

// other constants
static const unsigned int FASTA_LINE_MAX = 50;

// ---------------------------------------------
// ConvertPileupFormatVisitor declaration

class ConvertPileupFormatVisitor : public PileupVisitor {

    // ctor & dtor
    public:
        ConvertPileupFormatVisitor(const RefVector& references,
                                   const string& fastaFilename,
                                   const bool isPrintingMapQualities,
                                   ostream* out);
        ~ConvertPileupFormatVisitor(void);

    // PileupVisitor interface implementation
    public:
        void Visit(const PileupPosition& pileupData);

    // data members
    private:
        Fasta     m_fasta;
        bool      m_hasFasta;
        bool      m_isPrintingMapQualities;
        ostream*  m_out;
        RefVector m_references;
};
    
} // namespace BamTools
  
// ---------------------------------------------
// ConvertSettings implementation

struct ConvertTool::ConvertSettings {

    // flag
    bool HasInput;
    bool HasInputFilelist;
    bool HasOutput;
    bool HasFormat;
    bool HasRegion;

    // pileup flags
    bool HasFastaFilename;
    bool IsOmittingSamHeader;
    bool IsPrintingPileupMapQualities;
    
    // options
    vector<string> InputFiles;
    string InputFilelist;
    string OutputFilename;
    string Format;
    string Region;
    
    // pileup options
    string FastaFilename;

    // constructor
    ConvertSettings(void)
        : HasInput(false)
        , HasInputFilelist(false)
        , HasOutput(false)
        , HasFormat(false)
        , HasRegion(false)
        , HasFastaFilename(false)
        , IsOmittingSamHeader(false)
        , IsPrintingPileupMapQualities(false)
        , OutputFilename(Options::StandardOut())
        , FastaFilename("")
    { } 
};    

// ---------------------------------------------
// ConvertToolPrivate implementation  
  
struct ConvertTool::ConvertToolPrivate {
  
    // ctor & dtor
    public:
        ConvertToolPrivate(ConvertTool::ConvertSettings* settings)
            : m_settings(settings)
            , m_out(cout.rdbuf())
        { }

        ~ConvertToolPrivate(void) { }
    
    // interface
    public:
        bool Run(void);
        
    // internal methods
    private:
        void PrintBed(const BamAlignment& a);
        void PrintFasta(const BamAlignment& a);
        void PrintFastq(const BamAlignment& a);
        void PrintJson(const BamAlignment& a);
        void PrintSam(const BamAlignment& a);
        void PrintYaml(const BamAlignment& a);
        
        // special case - uses the PileupEngine
        bool RunPileupConversion(BamMultiReader* reader);
        
    // data members
    private: 
        ConvertTool::ConvertSettings* m_settings;
        RefVector m_references;
        ostream m_out;
};

bool ConvertTool::ConvertToolPrivate::Run(void) {
 
    // ------------------------------------
    // initialize conversion input/output
        
    // set to default input if none provided
    if ( !m_settings->HasInput && !m_settings->HasInputFilelist )
        m_settings->InputFiles.push_back(Options::StandardIn());
    
    // add files in the filelist to the input file list
    if ( m_settings->HasInputFilelist ) {

        ifstream filelist(m_settings->InputFilelist.c_str(), ios::in);
        if ( !filelist.is_open() ) {
            cerr << "bamtools convert ERROR: could not open input BAM file list... Aborting." << endl;
            return false;
        }

        string line;
        while ( getline(filelist, line) )
            m_settings->InputFiles.push_back(line);
    }

    // open input files
    BamMultiReader reader;
    if ( !reader.Open(m_settings->InputFiles) ) {
        cerr << "bamtools convert ERROR: could not open input BAM file(s)... Aborting." << endl;
        return false;
    }

    // if input is not stdin & a region is provided, look for index files
    if ( m_settings->HasInput && m_settings->HasRegion ) {
        if ( !reader.LocateIndexes() ) {
            cerr << "bamtools convert ERROR: could not locate index file(s)... Aborting." << endl;
            return false;
        }
    }

    // retrieve reference data
    m_references = reader.GetReferenceData();

    // set region if specified
    BamRegion region;
    if ( m_settings->HasRegion ) {
        if ( Utilities::ParseRegionString(m_settings->Region, reader, region) ) {

            if ( reader.HasIndexes() ) {
                if ( !reader.SetRegion(region) ) {
                    cerr << "bamtools convert ERROR: set region failed. Check that REGION describes a valid range" << endl;
                    reader.Close();
                    return false;
                }
            }

        } else {
            cerr << "bamtools convert ERROR: could not parse REGION: " << m_settings->Region << endl;
            cerr << "Check that REGION is in valid format (see documentation) and that the coordinates are valid"
                 << endl;
            reader.Close();
            return false;
        }
    }
        
    // if output file given
    ofstream outFile;
    if ( m_settings->HasOutput ) {
      
        // open output file stream
        outFile.open(m_settings->OutputFilename.c_str());
        if ( !outFile ) {
            cerr << "bamtools convert ERROR: could not open " << m_settings->OutputFilename
                 << " for output" << endl;
            return false; 
        }
        
        // set m_out to file's streambuf
        m_out.rdbuf(outFile.rdbuf()); 
    }
    
    // -------------------------------------
    // do conversion based on format
    
     bool convertedOk = true;
    
    // pileup is special case
    // conversion not done per alignment, like the other formats
    if ( m_settings->Format == FORMAT_PILEUP )
        convertedOk = RunPileupConversion(&reader);
    
    // all other formats
    else {
    
        bool formatError = false;
        
        // set function pointer to proper conversion method
        void (BamTools::ConvertTool::ConvertToolPrivate::*pFunction)(const BamAlignment&) = 0;
        if      ( m_settings->Format == FORMAT_BED )   pFunction = &BamTools::ConvertTool::ConvertToolPrivate::PrintBed;
        else if ( m_settings->Format == FORMAT_FASTA ) pFunction = &BamTools::ConvertTool::ConvertToolPrivate::PrintFasta;
        else if ( m_settings->Format == FORMAT_FASTQ ) pFunction = &BamTools::ConvertTool::ConvertToolPrivate::PrintFastq;
        else if ( m_settings->Format == FORMAT_JSON )  pFunction = &BamTools::ConvertTool::ConvertToolPrivate::PrintJson;
        else if ( m_settings->Format == FORMAT_SAM )   pFunction = &BamTools::ConvertTool::ConvertToolPrivate::PrintSam;
        else if ( m_settings->Format == FORMAT_YAML )  pFunction = &BamTools::ConvertTool::ConvertToolPrivate::PrintYaml;
        else { 
            cerr << "bamtools convert ERROR: unrecognized format: " << m_settings->Format << endl;
            cerr << "Please see documentation for list of supported formats " << endl;
            formatError = true;
            convertedOk = false;
        }
        
        // if format selected ok
        if ( !formatError ) {
        
            // if SAM format & not omitting header, print SAM header first
            if ( (m_settings->Format == FORMAT_SAM) && !m_settings->IsOmittingSamHeader ) 
                m_out << reader.GetHeaderText();
            
            // iterate through file, doing conversion
            BamAlignment a;
            while ( reader.GetNextAlignment(a) )
                (this->*pFunction)(a);
            
            // set flag for successful conversion
            convertedOk = true;
        }
    }
    
    // ------------------------
    // clean up & exit
    reader.Close();
    if ( m_settings->HasOutput )
        outFile.close();
    return convertedOk;   
}

// ----------------------------------------------------------
// Conversion/output methods
// ----------------------------------------------------------

void ConvertTool::ConvertToolPrivate::PrintBed(const BamAlignment& a) { 
  
    // tab-delimited, 0-based half-open 
    // (e.g. a 50-base read aligned to pos 10 could have BED coordinates (10, 60) instead of BAM coordinates (10, 59) )
    // <chromName> <chromStart> <chromEnd> <readName> <score> <strand>

    m_out << m_references.at(a.RefID).RefName << "\t"
          << a.Position << "\t"
          << a.GetEndPosition() << "\t"
          << a.Name << "\t"
          << a.MapQuality << "\t"
          << (a.IsReverseStrand() ? "-" : "+") << endl;
}

// print BamAlignment in FASTA format
// N.B. - uses QueryBases NOT AlignedBases
void ConvertTool::ConvertToolPrivate::PrintFasta(const BamAlignment& a) { 
    
    // >BamAlignment.Name
    // BamAlignment.QueryBases (up to FASTA_LINE_MAX bases per line)
    // ...
    //
    // N.B. - QueryBases are reverse-complemented if aligned to reverse strand
  
    // print header
    m_out << ">" << a.Name << endl;
    
    // handle reverse strand alignment - bases 
    string sequence = a.QueryBases;
    if ( a.IsReverseStrand() )
        Utilities::ReverseComplement(sequence);
    
    // if sequence fits on single line
    if ( sequence.length() <= FASTA_LINE_MAX )
        m_out << sequence << endl;
    
    // else split over multiple lines
    else {
      
        size_t position = 0;
        size_t seqLength = sequence.length(); // handle reverse strand alignment - bases & qualitiesth();
        
        // write subsequences to each line
        while ( position < (seqLength - FASTA_LINE_MAX) ) {
            m_out << sequence.substr(position, FASTA_LINE_MAX) << endl;
            position += FASTA_LINE_MAX;
        }
        
        // write final subsequence
        m_out << sequence.substr(position) << endl;
    }
}

// print BamAlignment in FASTQ format
// N.B. - uses QueryBases NOT AlignedBases
void ConvertTool::ConvertToolPrivate::PrintFastq(const BamAlignment& a) { 
  
    // @BamAlignment.Name
    // BamAlignment.QueryBases
    // +
    // BamAlignment.Qualities
    //
    // N.B. - QueryBases are reverse-complemented (& Qualities reversed) if aligned to reverse strand .
    //        Name is appended "/1" or "/2" if paired-end, to reflect which mate this entry is.
  
    // handle paired-end alignments
    string name = a.Name;
    if ( a.IsPaired() )
        name.append( (a.IsFirstMate() ? "/1" : "/2") );
  
    // handle reverse strand alignment - bases & qualities
    string qualities = a.Qualities;
    string sequence  = a.QueryBases;
    if ( a.IsReverseStrand() ) {
        Utilities::Reverse(qualities);
        Utilities::ReverseComplement(sequence);
    }
  
    // write to output stream
    m_out << "@" << name << endl
          << sequence    << endl
          << "+"         << endl
          << qualities   << endl;
}

// print BamAlignment in JSON format
void ConvertTool::ConvertToolPrivate::PrintJson(const BamAlignment& a) {
  
    // write name & alignment flag
    m_out << "{\"name\":\"" << a.Name << "\",\"alignmentFlag\":\"" << a.AlignmentFlag << "\",";
    
    // write reference name
    if ( (a.RefID >= 0) && (a.RefID < (int)m_references.size()) ) 
        m_out << "\"reference\":\"" << m_references[a.RefID].RefName << "\",";
    
    // write position & map quality
    m_out << "\"position\":" << a.Position+1 << ",\"mapQuality\":" << a.MapQuality << ",";
    
    // write CIGAR
    const vector<CigarOp>& cigarData = a.CigarData;
    if ( !cigarData.empty() ) {
        m_out << "\"cigar\":[";
        vector<CigarOp>::const_iterator cigarBegin = cigarData.begin();
        vector<CigarOp>::const_iterator cigarIter  = cigarBegin;
        vector<CigarOp>::const_iterator cigarEnd   = cigarData.end();
        for ( ; cigarIter != cigarEnd; ++cigarIter ) {
            const CigarOp& op = (*cigarIter);
            if (cigarIter != cigarBegin)
                m_out << ",";
            m_out << "\"" << op.Length << op.Type << "\"";
        }
        m_out << "],";
    }
    
    // write mate reference name, mate position, & insert size
    if ( a.IsPaired() && (a.MateRefID >= 0) && (a.MateRefID < (int)m_references.size()) ) {
        m_out << "\"mate\":{"
              << "\"reference\":\"" << m_references[a.MateRefID].RefName << "\","
              << "\"position\":" << a.MatePosition+1
              << ",\"insertSize\":" << a.InsertSize << "},";
    }
    
    // write sequence
    if ( !a.QueryBases.empty() ) 
        m_out << "\"queryBases\":\"" << a.QueryBases << "\",";
    
    // write qualities
    if ( !a.Qualities.empty() && a.Qualities.at(0) != (char)0xFF ) {
        string::const_iterator s = a.Qualities.begin();
        m_out << "\"qualities\":[" << static_cast<short>(*s) - 33;
        ++s;
        for ( ; s != a.Qualities.end(); ++s )
            m_out << "," << static_cast<short>(*s) - 33;
        m_out << "],";
    }
    
    // write alignment's source BAM file
    m_out << "\"filename\":\"" << a.Filename << "\",";

    // write tag data
    const char* tagData = a.TagData.c_str();
    const size_t tagDataLength = a.TagData.length();
    size_t index = 0;
    if ( index < tagDataLength ) {

        m_out << "\"tags\":{";
        
        while ( index < tagDataLength ) {

            if ( index > 0 )
                m_out << ",";
            
            // write tag name
            m_out << "\"" << a.TagData.substr(index, 2) << "\":";
            index += 2;
            
            // get data type
            char type = a.TagData.at(index);
            ++index;
            switch ( type ) {
                case (Constants::BAM_TAG_TYPE_ASCII) :
                    m_out << "\"" << tagData[index] << "\"";
                    ++index; 
                    break;
                
                case (Constants::BAM_TAG_TYPE_INT8) :
                    // force value into integer-type (instead of char value)
                    m_out << static_cast<int16_t>(tagData[index]);
                    ++index;
                    break;

                case (Constants::BAM_TAG_TYPE_UINT8) :
                    // force value into integer-type (instead of char value)
                    m_out << static_cast<uint16_t>(tagData[index]);
                    ++index; 
                    break;
                
                case (Constants::BAM_TAG_TYPE_INT16) :
                    m_out << BamTools::UnpackSignedShort(&tagData[index]);
                    index += sizeof(int16_t);
                    break;

                case (Constants::BAM_TAG_TYPE_UINT16) :
                    m_out << BamTools::UnpackUnsignedShort(&tagData[index]);
                    index += sizeof(uint16_t);
                    break;
                    
                case (Constants::BAM_TAG_TYPE_INT32) :
                    m_out << BamTools::UnpackSignedInt(&tagData[index]);
                    index += sizeof(int32_t);
                    break;

                case (Constants::BAM_TAG_TYPE_UINT32) :
                    m_out << BamTools::UnpackUnsignedInt(&tagData[index]);
                    index += sizeof(uint32_t);
                    break;

                case (Constants::BAM_TAG_TYPE_FLOAT) :
                    m_out << BamTools::UnpackFloat(&tagData[index]);
                    index += sizeof(float);
                    break;
                
                case (Constants::BAM_TAG_TYPE_HEX)    :
                case (Constants::BAM_TAG_TYPE_STRING) :
                    m_out << "\""; 
                    while (tagData[index]) {
                        if (tagData[index] == '\"')
                            m_out << "\\\""; // escape for json
                        else
                            m_out << tagData[index];
                        ++index;
                    }
                    m_out << "\""; 
                    ++index; 
                    break;      
            }
            
            if ( tagData[index] == '\0') 
                break;
        }

        m_out << "}";
    }

    m_out << "}" << endl;
}

// print BamAlignment in SAM format
void ConvertTool::ConvertToolPrivate::PrintSam(const BamAlignment& a) {
  
    // tab-delimited
    // <QNAME> <FLAG> <RNAME> <POS> <MAPQ> <CIGAR> <MRNM> <MPOS> <ISIZE> <SEQ> <QUAL> [ <TAG>:<VTYPE>:<VALUE> [...] ]
  
    // write name & alignment flag
   m_out << a.Name << "\t" << a.AlignmentFlag << "\t";

    // write reference name
    if ( (a.RefID >= 0) && (a.RefID < (int)m_references.size()) ) 
        m_out << m_references[a.RefID].RefName << "\t";
    else 
        m_out << "*\t";
    
    // write position & map quality
    m_out << a.Position+1 << "\t" << a.MapQuality << "\t";
    
    // write CIGAR
    const vector<CigarOp>& cigarData = a.CigarData;
    if ( cigarData.empty() ) m_out << "*\t";
    else {
        vector<CigarOp>::const_iterator cigarIter = cigarData.begin();
        vector<CigarOp>::const_iterator cigarEnd  = cigarData.end();
        for ( ; cigarIter != cigarEnd; ++cigarIter ) {
            const CigarOp& op = (*cigarIter);
            m_out << op.Length << op.Type;
        }
        m_out << "\t";
    }
    
    // write mate reference name, mate position, & insert size
    if ( a.IsPaired() && (a.MateRefID >= 0) && (a.MateRefID < (int)m_references.size()) ) {
        if ( a.MateRefID == a.RefID )
            m_out << "=\t";
        else
           m_out << m_references[a.MateRefID].RefName << "\t";
        m_out << a.MatePosition+1 << "\t" << a.InsertSize << "\t";
    } 
    else
        m_out << "*\t0\t0\t";
    
    // write sequence
    if ( a.QueryBases.empty() )
        m_out << "*\t";
    else
        m_out << a.QueryBases << "\t";
    
    // write qualities
    if ( a.Qualities.empty() || (a.Qualities.at(0) == (char)0xFF) )
        m_out << "*";
    else
        m_out << a.Qualities;
    
    // write tag data
    const char* tagData = a.TagData.c_str();
    const size_t tagDataLength = a.TagData.length();
    
    size_t index = 0;
    while ( index < tagDataLength ) {

        // write tag name   
        string tagName = a.TagData.substr(index, 2);
        m_out << "\t" << tagName << ":";
        index += 2;
        
        // get data type
        char type = a.TagData.at(index);
        ++index;
        switch ( type ) {
            case (Constants::BAM_TAG_TYPE_ASCII) :
                m_out << "A:" << tagData[index];
                ++index;
                break;

            case (Constants::BAM_TAG_TYPE_INT8) :
                // force value into integer-type (instead of char value)
                m_out << "i:" << static_cast<int16_t>(tagData[index]);
                ++index;
                break;

            case (Constants::BAM_TAG_TYPE_UINT8) :
                // force value into integer-type (instead of char value)
                m_out << "i:" << static_cast<uint16_t>(tagData[index]);
                ++index;
                break;

            case (Constants::BAM_TAG_TYPE_INT16) :
                m_out << "i:" << BamTools::UnpackSignedShort(&tagData[index]);
                index += sizeof(int16_t);
                break;

            case (Constants::BAM_TAG_TYPE_UINT16) :
                m_out << "i:" << BamTools::UnpackUnsignedShort(&tagData[index]);
                index += sizeof(uint16_t);
                break;

            case (Constants::BAM_TAG_TYPE_INT32) :
                m_out << "i:" << BamTools::UnpackSignedInt(&tagData[index]);
                index += sizeof(int32_t);
                break;

            case (Constants::BAM_TAG_TYPE_UINT32) :
                m_out << "i:" << BamTools::UnpackUnsignedInt(&tagData[index]);
                index += sizeof(uint32_t);
                break;

            case (Constants::BAM_TAG_TYPE_FLOAT) :
                m_out << "f:" << BamTools::UnpackFloat(&tagData[index]);
                index += sizeof(float);
                break;

            case (Constants::BAM_TAG_TYPE_HEX)    : // fall-through
            case (Constants::BAM_TAG_TYPE_STRING) :
                m_out << type << ":";
                while (tagData[index]) {
                    m_out << tagData[index];
                    ++index;
                }
                ++index;
                break;
        }

        if ( tagData[index] == '\0' )
            break;
    }

    m_out << endl;
}

// Print BamAlignment in YAML format
void ConvertTool::ConvertToolPrivate::PrintYaml(const BamAlignment& a) {

    // write alignment name
    m_out << "---" << endl;
    m_out << a.Name << ":" << endl;

    // write alignment data
    m_out << "   " << "AlndBases: "     << a.AlignedBases << endl;
    m_out << "   " << "Qualities: "     << a.Qualities << endl;
    m_out << "   " << "Name: "          << a.Name << endl;
    m_out << "   " << "Length: "        << a.Length << endl;
    m_out << "   " << "TagData: "       << a.TagData << endl;
    m_out << "   " << "RefID: "         << a.RefID << endl;
    m_out << "   " << "RefName: "       << m_references[a.RefID].RefName << endl;
    m_out << "   " << "Position: "      << a.Position << endl;
    m_out << "   " << "Bin: "           << a.Bin << endl;
    m_out << "   " << "MapQuality: "    << a.MapQuality << endl;
    m_out << "   " << "AlignmentFlag: " << a.AlignmentFlag << endl;
    m_out << "   " << "MateRefID: "     << a.MateRefID << endl;
    m_out << "   " << "MatePosition: "  << a.MatePosition << endl;
    m_out << "   " << "InsertSize: "    << a.InsertSize << endl;
    m_out << "   " << "Filename: "      << a.Filename << endl;

    // write Cigar data
    const vector<CigarOp>& cigarData = a.CigarData;
    if ( !cigarData.empty() ) {
        m_out << "   " <<  "Cigar: ";
        vector<CigarOp>::const_iterator cigarBegin = cigarData.begin();
        vector<CigarOp>::const_iterator cigarIter  = cigarBegin;
        vector<CigarOp>::const_iterator cigarEnd   = cigarData.end();
        for ( ; cigarIter != cigarEnd; ++cigarIter ) {
            const CigarOp& op = (*cigarIter);
            m_out << op.Length << op.Type;
        }
        m_out << endl;
    }
}

bool ConvertTool::ConvertToolPrivate::RunPileupConversion(BamMultiReader* reader) {
  
    // check for valid BamMultiReader
    if ( reader == 0 ) return false;
  
    // set up our pileup format 'visitor'
    ConvertPileupFormatVisitor* v = new ConvertPileupFormatVisitor(m_references, 
                                                                   m_settings->FastaFilename,
                                                                   m_settings->IsPrintingPileupMapQualities, 
                                                                   &m_out);

    // set up PileupEngine
    PileupEngine pileup;
    pileup.AddVisitor(v);
    
    // iterate through data
    BamAlignment al;
    while ( reader->GetNextAlignment(al) )
        pileup.AddAlignment(al);
    pileup.Flush();
    
    // clean up
    delete v;
    v = 0;
    
    // return success
    return true;
}       

// ---------------------------------------------
// ConvertTool implementation

ConvertTool::ConvertTool(void)
    : AbstractTool()
    , m_settings(new ConvertSettings)
    , m_impl(0)
{
    // set program details
    Options::SetProgramInfo("bamtools convert", "converts BAM to a number of other formats",
                            "-format <FORMAT> [-in <filename> -in <filename> ... | -list <filelist>] [-out <filename>] [-region <REGION>] [format-specific options]");
    
    // set up options 
    OptionGroup* IO_Opts = Options::CreateOptionGroup("Input & Output");
    Options::AddValueOption("-in",     "BAM filename", "the input BAM file(s)", "", m_settings->HasInput,   m_settings->InputFiles,     IO_Opts, Options::StandardIn());
    Options::AddValueOption("-list",   "filename", "the input BAM file list, one line per file", "", m_settings->HasInputFilelist,  m_settings->InputFilelist, IO_Opts);
    Options::AddValueOption("-out",    "BAM filename", "the output BAM file",   "", m_settings->HasOutput,  m_settings->OutputFilename, IO_Opts, Options::StandardOut());
    Options::AddValueOption("-format", "FORMAT", "the output file format - see README for recognized formats", "", m_settings->HasFormat, m_settings->Format, IO_Opts);
    Options::AddValueOption("-region", "REGION", "genomic region. Index file is recommended for better performance, and is used automatically if it exists. See \'bamtools help index\' for more details on creating one", "", m_settings->HasRegion, m_settings->Region, IO_Opts);
    
    OptionGroup* PileupOpts = Options::CreateOptionGroup("Pileup Options");
    Options::AddValueOption("-fasta", "FASTA filename", "FASTA reference file", "", m_settings->HasFastaFilename, m_settings->FastaFilename, PileupOpts);
    Options::AddOption("-mapqual", "print the mapping qualities", m_settings->IsPrintingPileupMapQualities, PileupOpts);
    
    OptionGroup* SamOpts = Options::CreateOptionGroup("SAM Options");
    Options::AddOption("-noheader", "omit the SAM header from output", m_settings->IsOmittingSamHeader, SamOpts);
}

ConvertTool::~ConvertTool(void) {

    delete m_settings;
    m_settings = 0;
    
    delete m_impl;
    m_impl = 0;
}

int ConvertTool::Help(void) {
    Options::DisplayHelp();
    return 0;
}

int ConvertTool::Run(int argc, char* argv[]) {
  
    // parse command line arguments
    Options::Parse(argc, argv, 1);
    
    // initialize ConvertTool with settings
    m_impl = new ConvertToolPrivate(m_settings);
    
    // run ConvertTool, return success/fail
    if ( m_impl->Run() ) 
        return 0;
    else 
        return 1;
}

// ---------------------------------------------
// ConvertPileupFormatVisitor implementation

ConvertPileupFormatVisitor::ConvertPileupFormatVisitor(const RefVector& references, 
                                                       const string& fastaFilename,
                                                       const bool isPrintingMapQualities,
                                                       ostream* out)
    : PileupVisitor()
    , m_hasFasta(false)
    , m_isPrintingMapQualities(isPrintingMapQualities)
    , m_out(out)
    , m_references(references)
{ 
    // set up Fasta reader if file is provided
    if ( !fastaFilename.empty() ) {
      
        // check for FASTA index
        string indexFilename = "";
        if ( Utilities::FileExists(fastaFilename + ".fai") ) 
            indexFilename = fastaFilename + ".fai";
      
        // open FASTA file
        if ( m_fasta.Open(fastaFilename, indexFilename) ) 
            m_hasFasta = true;
    }
}

ConvertPileupFormatVisitor::~ConvertPileupFormatVisitor(void) { 
    // be sure to close Fasta reader
    if ( m_hasFasta ) {
        m_fasta.Close();
        m_hasFasta = false;
    }
}

void ConvertPileupFormatVisitor::Visit(const PileupPosition& pileupData ) {
  
    // skip if no alignments at this position
    if ( pileupData.PileupAlignments.empty() ) return;
  
    // retrieve reference name
    const string& referenceName = m_references[pileupData.RefId].RefName;
    const int& position   = pileupData.Position;
    
    // retrieve reference base from FASTA file, if one provided; otherwise default to 'N'
    char referenceBase('N');
    if ( m_hasFasta && (pileupData.Position < m_references[pileupData.RefId].RefLength) ) {
        if ( !m_fasta.GetBase(pileupData.RefId, pileupData.Position, referenceBase ) ) {
            cerr << "bamtools convert ERROR: pileup conversion - could not read reference base from FASTA file" << endl;
            return;
        }
    }
    
    // get count of alleles at this position
    const int numberAlleles = pileupData.PileupAlignments.size();
    
    // -----------------------------------------------------------
    // build strings based on alleles at this positionInAlignment
    
    stringstream bases("");
    stringstream baseQualities("");
    stringstream mapQualities("");
    
    // iterate over alignments at this pileup position
    vector<PileupAlignment>::const_iterator pileupIter = pileupData.PileupAlignments.begin();
    vector<PileupAlignment>::const_iterator pileupEnd  = pileupData.PileupAlignments.end();
    for ( ; pileupIter != pileupEnd; ++pileupIter ) {
        const PileupAlignment pa = (*pileupIter);
        const BamAlignment& ba = pa.Alignment;
        
        // if beginning of read segment
        if ( pa.IsSegmentBegin )
            bases << '^' << ( ((int)ba.MapQuality > 93) ? (char)126 : (char)((int)ba.MapQuality+33) );
        
        // if current base is not a DELETION
        if ( !pa.IsCurrentDeletion ) {
          
            // get base at current position
            char base = ba.QueryBases.at(pa.PositionInAlignment);
            
            // if base matches reference
            if ( base == '=' || 
                 toupper(base) == toupper(referenceBase) ||
                 tolower(base) == tolower(referenceBase) ) 
            {
                base = ( ba.IsReverseStrand() ? ',' : '.' );
            }
            
            // mismatches reference
            else base = ( ba.IsReverseStrand() ? tolower(base) : toupper(base) );
            
            // store base
            bases << base;
          
            // if next position contains insertion
            if ( pa.IsNextInsertion ) {
                bases << '+' << pa.InsertionLength;
                for (int i = 1; i <= pa.InsertionLength; ++i) {
                    char insertedBase = (char)ba.QueryBases.at(pa.PositionInAlignment + i);
                    bases << (ba.IsReverseStrand() ? (char)tolower(insertedBase) : (char)toupper(insertedBase) );
                }
            }
            
            // if next position contains DELETION
            else if ( pa.IsNextDeletion ) {
                bases << '-' << pa.DeletionLength;
                for (int i = 1; i <= pa.DeletionLength; ++i) {
                    char deletedBase('N');
                    if ( m_hasFasta && (pileupData.Position+i < m_references[pileupData.RefId].RefLength) ) {
                        if ( !m_fasta.GetBase(pileupData.RefId, pileupData.Position+i, deletedBase ) ) {
                            cerr << "bamtools convert ERROR: pileup conversion - could not read reference base from FASTA file" << endl;
                            return;
                        }
                    }
                    bases << (ba.IsReverseStrand() ? (char)tolower(deletedBase) : (char)toupper(deletedBase) );
                }
            }
        }
        
        // otherwise, DELETION
        else bases << '*';
        
        // if end of read segment
        if ( pa.IsSegmentEnd )
            bases << '$';
        
        // store current base quality
        baseQualities << ba.Qualities.at(pa.PositionInAlignment);
        
        // save alignment map quality if desired
        if ( m_isPrintingMapQualities )
            mapQualities << ( ((int)ba.MapQuality > 93) ? (char)126 : (char)((int)ba.MapQuality+33) );
    }
    
    // ----------------------
    // print results 
    
    // tab-delimited
    // <refName> <1-based pos> <refBase> <numberAlleles> <bases> <qualities> [mapQuals]
    
    const string TAB = "\t";
    *m_out << referenceName       << TAB 
           << position + 1        << TAB 
           << referenceBase       << TAB 
           << numberAlleles       << TAB 
           << bases.str()         << TAB 
           << baseQualities.str() << TAB
           << mapQualities.str()  << endl;
}
