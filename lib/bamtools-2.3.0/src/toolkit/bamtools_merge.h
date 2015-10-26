// ***************************************************************************
// bamtools_merge.h (c) 2010 Derek Barnett, Erik Garrison
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 7 April 2011
// ---------------------------------------------------------------------------
// Merges multiple BAM files into one
// ***************************************************************************

#ifndef BAMTOOLS_MERGE_H
#define BAMTOOLS_MERGE_H

#include "bamtools_tool.h"

namespace BamTools {
  
class MergeTool : public AbstractTool {
  
    public:
        MergeTool(void);
        ~MergeTool(void);
  
    public:
        int Help(void);
        int Run(int argc, char* argv[]); 
        
    private:
        struct MergeSettings;
        MergeSettings* m_settings;

        struct MergeToolPrivate;
        MergeToolPrivate* m_impl;
};
  
} // namespace BamTools

#endif // BAMTOOLS_MERGE_H
