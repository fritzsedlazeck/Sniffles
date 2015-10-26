// ***************************************************************************
// bamtools_convert.h (c) 2010 Derek Barnett, Erik Garrison
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 9 July 2010
// ---------------------------------------------------------------------------
// Converts between BAM and a number of other formats
// ***************************************************************************

#ifndef BAMTOOLS_CONVERT_H
#define BAMTOOLS_CONVERT_H

#include "bamtools_tool.h"

namespace BamTools { 
  
class ConvertTool : public AbstractTool {
  
    public:
        ConvertTool(void);
        ~ConvertTool(void);

    public:
        int Help(void);
        int Run(int argc, char* argv[]); 
        
    private: 
        struct ConvertSettings;
        ConvertSettings* m_settings;
        
        struct ConvertToolPrivate;
        ConvertToolPrivate* m_impl;
};
  
} // namespace BamTools

#endif // BAMTOOLS_CONVERT_H
