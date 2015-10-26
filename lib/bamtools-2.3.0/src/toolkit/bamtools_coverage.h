// ***************************************************************************
// bamtools_coverage.h (c) 2010 Derek Barnett, Erik Garrison
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 1 August 2010
// ---------------------------------------------------------------------------
// Prints coverage data for a single BAM file 
// ***************************************************************************

#ifndef BAMTOOLS_COVERAGE_H
#define BAMTOOLS_COVERAGE_H

#include "bamtools_tool.h"

namespace BamTools {
  
class CoverageTool : public AbstractTool {
  
    public:
        CoverageTool(void);
        ~CoverageTool(void);
  
    public:
        int Help(void);
        int Run(int argc, char* argv[]); 
        
    private:  
        struct CoverageSettings;
        CoverageSettings* m_settings;
        
        struct CoverageToolPrivate;
        CoverageToolPrivate* m_impl;
};
  
} // namespace BamTools

#endif // BAMTOOLS_COVERAGE_H
