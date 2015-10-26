// ***************************************************************************
// bamtools.cpp (c) 2010 Derek Barnett, Erik Garrison
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 12 October 2012 (DB)
// ---------------------------------------------------------------------------
// Integrates a number of BamTools functionalities into a single executable.
// ***************************************************************************

#include "bamtools_convert.h"
#include "bamtools_count.h"
#include "bamtools_coverage.h"
#include "bamtools_filter.h"
#include "bamtools_header.h"
#include "bamtools_index.h"
#include "bamtools_merge.h"
#include "bamtools_random.h"
#include "bamtools_resolve.h"
#include "bamtools_revert.h"
#include "bamtools_sort.h"
#include "bamtools_split.h"
#include "bamtools_stats.h"
#include "bamtools_version.h"
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>
using namespace BamTools;
using namespace std;

// bamtools subtool names
static const string CONVERT  = "convert";
static const string COUNT    = "count";
static const string COVERAGE = "coverage";
static const string FILTER   = "filter";
static const string HEADER   = "header";
static const string INDEX    = "index";
static const string MERGE    = "merge";
static const string RANDOM   = "random";
static const string RESOLVE  = "resolve";
static const string REVERT   = "revert";
static const string SORT     = "sort";
static const string SPLIT    = "split";
static const string STATS    = "stats";

// bamtools help/version constants
static const string HELP          = "help";
static const string LONG_HELP     = "--help";
static const string SHORT_HELP    = "-h";
static const string VERSION       = "version";
static const string LONG_VERSION  = "--version";
static const string SHORT_VERSION = "-v";

// determine if string is a help constant
static bool IsHelp(char* str) {
    return ( str == HELP ||
             str == LONG_HELP ||
             str == SHORT_HELP );
}

// determine if string is a version constant
static bool IsVersion(char* str) {
    return ( str == VERSION ||
             str == LONG_VERSION ||
             str == SHORT_VERSION );
}

// subtool factory method
AbstractTool* CreateTool(const string& arg) {
  
    // determine tool type based on arg
    if ( arg == CONVERT )  return new ConvertTool;
    if ( arg == COUNT )    return new CountTool;
    if ( arg == COVERAGE ) return new CoverageTool;
    if ( arg == FILTER )   return new FilterTool;
    if ( arg == HEADER )   return new HeaderTool;
    if ( arg == INDEX )    return new IndexTool;
    if ( arg == MERGE )    return new MergeTool;
    if ( arg == RANDOM )   return new RandomTool;
    if ( arg == RESOLVE )  return new ResolveTool;
    if ( arg == REVERT )   return new RevertTool;
    if ( arg == SORT )     return new SortTool;
    if ( arg == SPLIT )    return new SplitTool;
    if ( arg == STATS )    return new StatsTool;

    // unknown arg
    return 0;
}

// print help info
int Help(int argc, char* argv[]) {
  
    // check for 'bamtools help COMMAND' to print tool-specific help message
    if (argc > 2) {
        
	// determine desired sub-tool
        AbstractTool* tool = CreateTool( argv[2] );

        // if tool known, print its help screen
        if ( tool ) return tool->Help();
    }

    // print general BamTools help message
    cerr << endl;
    cerr << "usage: bamtools [--help] COMMAND [ARGS]" << endl;
    cerr << endl;
    cerr << "Available bamtools commands:" << endl;
    cerr << "\tconvert         Converts between BAM and a number of other formats" << endl;
    cerr << "\tcount           Prints number of alignments in BAM file(s)" << endl;
    cerr << "\tcoverage        Prints coverage statistics from the input BAM file" << endl;    
    cerr << "\tfilter          Filters BAM file(s) by user-specified criteria" << endl;
    cerr << "\theader          Prints BAM header information" << endl;
    cerr << "\tindex           Generates index for BAM file" << endl;
    cerr << "\tmerge           Merge multiple BAM files into single file" << endl;
    cerr << "\trandom          Select random alignments from existing BAM file(s), intended more as a testing tool." << endl;
    cerr << "\tresolve         Resolves paired-end reads (marking the IsProperPair flag as needed)" << endl;
    cerr << "\trevert          Removes duplicate marks and restores original base qualities" << endl;
    cerr << "\tsort            Sorts the BAM file according to some criteria" << endl;
    cerr << "\tsplit           Splits a BAM file on user-specified property, creating a new BAM output file for each value found" << endl;
    cerr << "\tstats           Prints some basic statistics from input BAM file(s)" << endl;
    cerr << endl;
    cerr << "See 'bamtools help COMMAND' for more information on a specific command." << endl;
    cerr << endl;
    return EXIT_SUCCESS;
}

// print version info
int Version(void) {

    stringstream versionStream("");
    versionStream << BAMTOOLS_VERSION_MAJOR << "."
                  << BAMTOOLS_VERSION_MINOR << "."
                  << BAMTOOLS_VERSION_BUILD;

    cout << endl;
    cout << "bamtools " << versionStream.str() << endl;
    cout << "Part of BamTools API and toolkit" << endl;
    cout << "Primary authors: Derek Barnett, Erik Garrison, Michael Stromberg" << endl;
    cout << "(c) 2009-2012 Marth Lab, Biology Dept., Boston College" << endl;
    cout << endl;
    return EXIT_SUCCESS;
}

// toolkit entry point
int main(int argc, char* argv[]) {

    // just 'bamtools'
    if ( (argc == 1) ) return Help(argc, argv);
    
    // 'bamtools help', 'bamtools --help', or 'bamtools -h'
    if ( IsHelp(argv[1]) ) return Help(argc, argv); 
    
    // 'bamtools version', 'bamtools --version', or 'bamtools -v'
    if ( IsVersion(argv[1]) ) return Version(); 
        
    // determine desired sub-tool, run if found
    AbstractTool* tool = CreateTool( argv[1] );
    if ( tool ) return tool->Run(argc, argv);

    // no tool matched, show help
    return Help(argc, argv);
}
