// ***************************************************************************
// bamtools_resolve.h (c) 2011 Derek Barnett
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 23 June 2011
// ---------------------------------------------------------------------------
// Resolves paired-end reads (marking the IsProperPair flag as needed).
// ***************************************************************************

#ifndef BAMTOOLS_RESOLVE_H
#define BAMTOOLS_RESOLVE_H

#include "bamtools_tool.h"

namespace BamTools {

class ResolveTool : public AbstractTool {

    public:
        ResolveTool(void);
        ~ResolveTool(void);

    public:
        int Help(void);
        int Run(int argc, char* argv[]);

    private:
        struct ResolveSettings;
        ResolveSettings* m_settings;

        struct ResolveToolPrivate;
        ResolveToolPrivate* m_impl;

        struct ReadNamesFileReader;
        struct ReadNamesFileWriter;
        struct StatsFileReader;
        struct StatsFileWriter;
};

} // namespace BamTools

#endif // BAMTOOLS_RESOLVE_H
