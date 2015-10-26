// ***************************************************************************
// bamtools_tool.h (c) 2010 Derek Barnett, Erik Garrison
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 2 June 2010
// ---------------------------------------------------------------------------
// Base class for all other BamTools sub-tools
// All derived classes must provide Help() and Run() methods
// ***************************************************************************

#ifndef BAMTOOLS_ABSTRACTTOOL_H
#define BAMTOOLS_ABSTRACTTOOL_H

#include <string>

namespace BamTools { 
  
class AbstractTool {
  
    public:
        AbstractTool(void) { }
        virtual ~AbstractTool(void) { }

    public:
        virtual int Help(void) =0;
        virtual int Run(int argc, char* argv[]) =0; 

    // derived classes should also provide:
    // static std::string Description(void);
    // static std::String Name(void);
};
  
} // namespace BamTools

#endif // BAMTOOLS_ABSTRACTTOOL_H
