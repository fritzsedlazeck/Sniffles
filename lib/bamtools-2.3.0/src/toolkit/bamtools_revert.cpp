// ***************************************************************************
// bamtools_revert.cpp (c) 2010 Derek Barnett, Alistair Ward
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 7 April 2011
// ---------------------------------------------------------------------------
// Removes duplicate marks and restores original base qualities
// ***************************************************************************

#include "bamtools_revert.h"

#include <api/BamReader.h>
#include <api/BamWriter.h>
#include <utils/bamtools_options.h>
#include <utils/bamtools_utilities.h>
using namespace BamTools;

#include <iostream>
#include <string>
using namespace std;

namespace BamTools {

static const string OQ_TAG = "OQ";

} // namespace BamTools;

// ---------------------------------------------
// RevertSettings implementation

struct RevertTool::RevertSettings {

    // flags
    bool HasInput;
    bool HasOutput;
    bool IsForceCompression;
    bool IsKeepDuplicateFlag;
    bool IsKeepQualities;

    // filenames
    string InputFilename;
    string OutputFilename;
    
    // constructor
    RevertSettings(void)
        : HasInput(false)
        , HasOutput(false)
        , IsForceCompression(false)
        , IsKeepDuplicateFlag(false)
        , IsKeepQualities(false)
        , InputFilename(Options::StandardIn())
        , OutputFilename(Options::StandardOut())
    { }
};  

// ---------------------------------------------
// RevertToolPrivate implementation

struct RevertTool::RevertToolPrivate {
  
    // ctor & dtor
    public:
        RevertToolPrivate(RevertTool::RevertSettings* settings)
            : m_settings(settings)
        { }
        ~RevertToolPrivate(void) { }
  
    // 'public' interface
    public:
        bool Run(void);
        
    // internal methods
    private:
        void RevertAlignment(BamAlignment& al);
        
    // data members
    private:
        RevertTool::RevertSettings* m_settings;
};

// 'reverts' a BAM alignment
// default behavior (for now) is:
//   1 - replace Qualities with OQ contents
//   2 - clear IsDuplicate flag
// can override default behavior using command line options
void RevertTool::RevertToolPrivate::RevertAlignment(BamAlignment& al) {

    // replace Qualities with OQ contents, if requested
    if ( !m_settings->IsKeepQualities ) {
        string originalQualities;
        if ( al.GetTag(OQ_TAG, originalQualities) ) {
            al.Qualities = originalQualities;
            al.RemoveTag(OQ_TAG);
        }
    }

    // clear duplicate flag, if requested
    if ( !m_settings->IsKeepDuplicateFlag )
        al.SetIsDuplicate(false);
}

bool RevertTool::RevertToolPrivate::Run(void) {
  
    // opens the BAM file without checking for indexes
    BamReader reader;
    if ( !reader.Open(m_settings->InputFilename) ) {
        cerr << "bamtools revert ERROR: could not open " << m_settings->InputFilename
             << " for reading... Aborting." << endl;
        return false;
    }

    // get BAM file metadata
    const string& headerText = reader.GetHeaderText();
    const RefVector& references = reader.GetReferenceData();
    
    // determine compression mode for BamWriter
    bool writeUncompressed = ( m_settings->OutputFilename == Options::StandardOut() &&
                              !m_settings->IsForceCompression );
    BamWriter::CompressionMode compressionMode = BamWriter::Compressed;
    if ( writeUncompressed ) compressionMode = BamWriter::Uncompressed;

    // open BamWriter
    BamWriter writer;
    writer.SetCompressionMode(compressionMode);
    if ( !writer.Open(m_settings->OutputFilename, headerText, references) ) {
        cerr << "bamtools revert ERROR: could not open " << m_settings->OutputFilename
             << " for writing... Aborting." << endl;
        reader.Close();
        return false;
    }

    // plow through file, reverting alignments
    BamAlignment al;
    while ( reader.GetNextAlignment(al) ) {
        RevertAlignment(al);
        writer.SaveAlignment(al);
    }
    
    // clean and exit
    reader.Close();
    writer.Close();
    return true; 
}

// ---------------------------------------------
// RevertTool implementation

RevertTool::RevertTool(void)
    : AbstractTool()
    , m_settings(new RevertSettings)
    , m_impl(0)
{
    // set program details
    Options::SetProgramInfo("bamtools revert", "removes duplicate marks and restores original (non-recalibrated) base qualities", "[-in <filename> -in <filename> ...] [-out <filename> | [-forceCompression]] [revertOptions]");
    
    // set up options 
    OptionGroup* IO_Opts = Options::CreateOptionGroup("Input & Output");
    Options::AddValueOption("-in",  "BAM filename", "the input BAM file",  "", m_settings->HasInput,  m_settings->InputFilename,  IO_Opts, Options::StandardIn());
    Options::AddValueOption("-out", "BAM filename", "the output BAM file", "", m_settings->HasOutput, m_settings->OutputFilename, IO_Opts, Options::StandardOut());
    Options::AddOption("-forceCompression", "if results are sent to stdout (like when piping to another tool), default behavior is to leave output uncompressed. Use this flag to override and force compression", m_settings->IsForceCompression, IO_Opts);

    OptionGroup* RevertOpts = Options::CreateOptionGroup("Revert Options");
    Options::AddOption("-keepDuplicate", "keep duplicates marked", m_settings->IsKeepDuplicateFlag, RevertOpts);
    Options::AddOption("-keepQualities", "keep base qualities (do not replace with OQ contents)", m_settings->IsKeepQualities, RevertOpts);
}

RevertTool::~RevertTool(void) {

    delete m_settings;
    m_settings = 0;
    
    delete m_impl;
    m_impl = 0;
}

int RevertTool::Help(void) {
    Options::DisplayHelp();
    return 0;
}

int RevertTool::Run(int argc, char* argv[]) {
  
    // parse command line arguments
    Options::Parse(argc, argv, 1);

    // intialize RevertTool with settings
    m_impl = new RevertToolPrivate(m_settings);
    
    // run RevertTool, return success/fail
    if ( m_impl->Run() )
        return 0;
    else
        return 1;
}
