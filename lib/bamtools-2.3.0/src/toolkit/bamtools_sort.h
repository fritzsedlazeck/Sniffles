// ***************************************************************************
// bamtools_sort.h (c) 2010 Derek Barnett, Erik Garrison
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 7 April 2011 (DB)
// ---------------------------------------------------------------------------
// Sorts a BAM file
// ***************************************************************************

#ifndef BAMTOOLS_SORT_H
#define BAMTOOLS_SORT_H

#include "bamtools_tool.h"

namespace BamTools {
  
class SortTool : public AbstractTool {
  
    public:
        SortTool(void);
        ~SortTool(void);
  
    public:
        int Help(void);
        int Run(int argc, char* argv[]); 
        
    private:
        struct SortSettings;
        SortSettings* m_settings;
        
        struct SortToolPrivate;
        SortToolPrivate* m_impl;
};
  
} // namespace BamTools

#endif // BAMTOOLS_SORT_H
