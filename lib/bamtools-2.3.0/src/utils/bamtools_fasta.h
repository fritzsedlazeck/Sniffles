// ***************************************************************************
// bamtools_fasta.h (c) 2010 Derek Barnett, Erik Garrison
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 10 October 2011
// ---------------------------------------------------------------------------
// Provides FASTA reading/indexing functionality.
// ***************************************************************************

#ifndef BAMTOOLS_FASTA_H
#define BAMTOOLS_FASTA_H

#include "utils/utils_global.h"
#include <string>

namespace BamTools {

class UTILS_EXPORT Fasta {
  
    // ctor & dtor
    public:
        Fasta(void);
        ~Fasta(void);
        
    // file-handling methods
    public:
        bool Close(void);
        bool Open(const std::string& filename, const std::string& indexFilename = "");
        
    // sequence access methods
    public:
        bool GetBase(const int& refID, const int& position, char& base);
        bool GetSequence(const int& refId, const int& start, const int& stop, std::string& sequence);
        
    // index-handling methods
    public:
        bool CreateIndex(const std::string& indexFilename);

    // internal implementation
    private:
        struct FastaPrivate;
        FastaPrivate* d;
};
  
} // BAMTOOLS_FASTA_H
  
#endif // BAMTOOLS_FASTA_H
