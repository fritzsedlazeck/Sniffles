// ***************************************************************************
// bamtools_stats.h (c) 2010 Derek Barnett, Erik Garrison
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 7 April 2011
// ---------------------------------------------------------------------------
// Prints general statistics for a single BAM file
// ***************************************************************************

#ifndef BAMTOOLS_STATS_H
#define BAMTOOLS_STATS_H

#include "bamtools_tool.h"

namespace BamTools {
  
class StatsTool : public AbstractTool {
  
    public:
        StatsTool(void);
        ~StatsTool(void);
  
    public:
        int Help(void);
        int Run(int argc, char* argv[]); 
        
    private:
        struct StatsSettings;
        StatsSettings* m_settings;
        
        struct StatsToolPrivate;
        StatsToolPrivate* m_impl;
};
  
} // namespace BamTools

#endif // BAMTOOLS_STATS_H
