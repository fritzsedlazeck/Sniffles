// ***************************************************************************
// bamtools_header.h (c) 2010 Derek Barnett, Erik Garrison
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 7 April 2011
// ---------------------------------------------------------------------------
// Prints the SAM-style header from a single BAM file ( or merged header from
// multiple BAM files) to stdout
// ***************************************************************************

#ifndef BAMTOOLS_HEADER_H
#define BAMTOOLS_HEADER_H

#include "bamtools_tool.h"

namespace BamTools {
  
class HeaderTool : public AbstractTool {
  
    public:
        HeaderTool(void);
        ~HeaderTool(void);
  
    public:
        int Help(void);
        int Run(int argc, char* argv[]); 
        
    private:
        struct HeaderSettings;
        HeaderSettings* m_settings;

        struct HeaderToolPrivate;
        HeaderToolPrivate* m_impl;
};
  
} // namespace BamTools

#endif // BAMTOOLS_HEADER_H
