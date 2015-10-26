// ***************************************************************************
// bamtools_count.h (c) 2010 Derek Barnett, Erik Garrison
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 7 April 2011
// ---------------------------------------------------------------------------
// Prints alignment count for BAM file(s)
// ***************************************************************************

#ifndef BAMTOOLS_COUNT_H
#define BAMTOOLS_COUNT_H

#include "bamtools_tool.h"

namespace BamTools { 
  
class CountTool : public AbstractTool {
  
    public:
        CountTool(void);
        ~CountTool(void);

    public:
        int Help(void);
        int Run(int argc, char* argv[]); 
        
    private: 
        struct CountSettings;
        CountSettings* m_settings;

        struct CountToolPrivate;
        CountToolPrivate* m_impl;
};
  
} // namespace BamTools

#endif // BAMTOOLS_COUNT_H
