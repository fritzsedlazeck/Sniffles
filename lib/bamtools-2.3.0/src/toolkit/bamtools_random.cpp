// ***************************************************************************
// bamtools_random.cpp (c) 2010 Derek Barnett, Erik Garrison
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 24 July 2013 (DB)
// ---------------------------------------------------------------------------
// Grab a random subset of alignments (testing tool)
// ***************************************************************************

#include "bamtools_random.h"

#include <api/BamMultiReader.h>
#include <api/BamWriter.h>
#include <utils/bamtools_options.h>
#include <utils/bamtools_utilities.h>
using namespace BamTools;

#include <ctime>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
using namespace std;
  
namespace BamTools {
  
// define constants
const unsigned int RANDOM_MAX_ALIGNMENT_COUNT = 10000;

// utility methods for RandomTool
int getRandomInt(const int& lowerBound, const int& upperBound) {
    const int range = (upperBound - lowerBound) + 1;
    return ( lowerBound + (int)(range * (double)rand()/((double)RAND_MAX + 1)) );
}
    
} // namespace BamTools
  
// ---------------------------------------------  
// RandomSettings implementation

struct RandomTool::RandomSettings {

    // flags
    bool HasAlignmentCount;
    bool HasInput;
    bool HasInputFilelist;
    bool HasOutput;
    bool HasRandomNumberSeed;
    bool HasRegion;
    bool IsForceCompression;

    // parameters
    unsigned int AlignmentCount;
    vector<string> InputFiles;
    string InputFilelist;
    string OutputFilename;
    unsigned int RandomNumberSeed;
    string Region;
    
    // constructor
    RandomSettings(void)
        : HasAlignmentCount(false)
        , HasInput(false)
        , HasInputFilelist(false)
        , HasOutput(false)
        , HasRandomNumberSeed(false)
        , HasRegion(false)
        , IsForceCompression(false)
        , AlignmentCount(RANDOM_MAX_ALIGNMENT_COUNT)
        , OutputFilename(Options::StandardOut())
        , RandomNumberSeed(0)
    { }  
};  

// ---------------------------------------------
// RandomToolPrivate implementation

struct RandomTool::RandomToolPrivate {

    // ctor & dtor
    public:
        RandomToolPrivate(RandomTool::RandomSettings* settings)
            : m_settings(settings)
        { }

        ~RandomToolPrivate(void) { }

    // interface
    public:
        bool Run(void);

    // data members
    private:
        RandomTool::RandomSettings* m_settings;
};

bool RandomTool::RandomToolPrivate::Run(void) {

    // set to default stdin if no input files provided
    if ( !m_settings->HasInput && !m_settings->HasInputFilelist )
        m_settings->InputFiles.push_back(Options::StandardIn());

    // add files in the filelist to the input file list
    if ( m_settings->HasInputFilelist ) {

        ifstream filelist(m_settings->InputFilelist.c_str(), ios::in);
        if ( !filelist.is_open() ) {
            cerr << "bamtools random ERROR: could not open input BAM file list... Aborting." << endl;
            return false;
        }

        string line;
        while ( getline(filelist, line) )
            m_settings->InputFiles.push_back(line);
    }

    // open our reader
    BamMultiReader reader;
    if ( !reader.Open(m_settings->InputFiles) ) {
        cerr << "bamtools random ERROR: could not open input BAM file(s)... Aborting." << endl;
        return false;
    }

    // look up index files for all BAM files
    reader.LocateIndexes();

    // make sure index data is available
    if ( !reader.HasIndexes() ) {
        cerr << "bamtools random ERROR: could not load index data for all input BAM file(s)... Aborting." << endl;
        reader.Close();
        return false;
    }

    // get BamReader metadata
    const string headerText = reader.GetHeaderText();
    const RefVector references = reader.GetReferenceData();
    if ( references.empty() ) {
        cerr << "bamtools random ERROR: no reference data available... Aborting." << endl;
        reader.Close();
        return false;
    }

    // determine compression mode for BamWriter
    bool writeUncompressed = ( m_settings->OutputFilename == Options::StandardOut() &&
                              !m_settings->IsForceCompression );
    BamWriter::CompressionMode compressionMode = BamWriter::Compressed;
    if ( writeUncompressed ) compressionMode = BamWriter::Uncompressed;

    // open BamWriter
    BamWriter writer;
    writer.SetCompressionMode(compressionMode);
    if ( !writer.Open(m_settings->OutputFilename, headerText, references) ) {
        cerr << "bamtools random ERROR: could not open " << m_settings->OutputFilename
             << " for writing... Aborting." << endl;
        reader.Close();
        return false;
    }

    // if user specified a REGION constraint, attempt to parse REGION string
    BamRegion region;
    if ( m_settings->HasRegion && !Utilities::ParseRegionString(m_settings->Region, reader, region) ) {
        cerr << "bamtools random ERROR: could not parse REGION: " << m_settings->Region << endl;
        cerr << "Check that REGION is in valid format (see documentation) and that the coordinates are valid"
             << endl;
        reader.Close();
        writer.Close();
        return false;
    }

    // seed our random number generator
    if ( m_settings->HasRandomNumberSeed )
        srand( m_settings->RandomNumberSeed );
    else
        srand( time(NULL) );

    // grab random alignments
    BamAlignment al;
    unsigned int i = 0;
    while ( i < m_settings->AlignmentCount ) {

        int randomRefId    = 0;
        int randomPosition = 0;

        // use REGION constraints to select random refId & position
        if ( m_settings->HasRegion ) {

            // select a random refId
            randomRefId = getRandomInt(region.LeftRefID, region.RightRefID);

            // select a random position based on randomRefId
            const int lowerBoundPosition = ( (randomRefId == region.LeftRefID)
                                             ? region.LeftPosition
                                             : 0 );
            const int upperBoundPosition = ( (randomRefId == region.RightRefID)
                                             ? region.RightPosition
                                             : (references.at(randomRefId).RefLength - 1) );
            randomPosition = getRandomInt(lowerBoundPosition, upperBoundPosition);
        }

        // otherwise select from all possible random refId & position
        else {

            // select random refId
            randomRefId = getRandomInt(0, (int)references.size() - 1);

            // select random position based on randomRefId
            const int lowerBoundPosition = 0;
            const int upperBoundPosition = references.at(randomRefId).RefLength - 1;
            randomPosition = getRandomInt(lowerBoundPosition, upperBoundPosition);
        }

        // if jump & read successful, save first alignment that overlaps random refId & position
        if ( reader.Jump(randomRefId, randomPosition) ) {
            while ( reader.GetNextAlignmentCore(al) ) {
                if ( al.RefID == randomRefId && al.Position >= randomPosition ) {
                    writer.SaveAlignment(al);
                    ++i;
                    break;
                }
            }
        }
    }

    // cleanup & exit
    reader.Close();
    writer.Close();
    return true;
}

// ---------------------------------------------
// RandomTool implementation

RandomTool::RandomTool(void) 
    : AbstractTool()
    , m_settings(new RandomSettings)
    , m_impl(0)
{ 
    // set program details
    Options::SetProgramInfo("bamtools random", "grab a random subset of alignments",
                            "[-in <filename> -in <filename> ... | -list <filelist>] [-out <filename>] [-forceCompression] [-n] [-region <REGION>]");
    
    // set up options 
    OptionGroup* IO_Opts = Options::CreateOptionGroup("Input & Output");
    Options::AddValueOption("-in",     "BAM filename", "the input BAM file",                         "", m_settings->HasInput,          m_settings->InputFiles,     IO_Opts, Options::StandardIn());
    Options::AddValueOption("-list",   "filename",     "the input BAM file list, one line per file", "", m_settings->HasInputFilelist,  m_settings->InputFilelist,  IO_Opts);
    Options::AddValueOption("-out",    "BAM filename", "the output BAM file",                        "", m_settings->HasOutput,         m_settings->OutputFilename, IO_Opts, Options::StandardOut());
    Options::AddValueOption("-region", "REGION",       "only pull random alignments from within this genomic region. Index file is recommended for better performance, and is used automatically if it exists. See \'bamtools help index\' for more details on creating one", "", m_settings->HasRegion, m_settings->Region, IO_Opts);
    Options::AddOption("-forceCompression", "if results are sent to stdout (like when piping to another tool), default behavior is to leave output uncompressed. Use this flag to override and force compression", m_settings->IsForceCompression, IO_Opts);
    
    OptionGroup* SettingsOpts = Options::CreateOptionGroup("Settings");
    Options::AddValueOption("-n", "count", "number of alignments to grab. Note - no duplicate checking is performed", "",
                            m_settings->HasAlignmentCount, m_settings->AlignmentCount, SettingsOpts, RANDOM_MAX_ALIGNMENT_COUNT);
    Options::AddValueOption("-seed", "unsigned integer", "random number generator seed (for repeatable results). Current time is used if no seed value is provided.", "",
                            m_settings->HasRandomNumberSeed, m_settings->RandomNumberSeed, SettingsOpts);
}

RandomTool::~RandomTool(void) { 

    delete m_settings;
    m_settings = 0;

    delete m_impl;
    m_impl = 0;
}

int RandomTool::Help(void) { 
    Options::DisplayHelp();
    return 0;
} 

int RandomTool::Run(int argc, char* argv[]) { 

    // parse command line arguments
    Options::Parse(argc, argv, 1);

    // initialize RandomTool with settings
    m_impl = new RandomToolPrivate(m_settings);

    // run RandomTool, return success/fail
    if ( m_impl->Run() )
        return 0;
    else
        return 1;
}
