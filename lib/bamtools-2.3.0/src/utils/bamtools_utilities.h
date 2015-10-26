// ***************************************************************************
// bamtools_utilities.h (c) 2010 Derek Barnett, Erik Garrison
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 7 October 2011
// ---------------------------------------------------------------------------
// Provides general utilities used by BamTools sub-tools.
// ***************************************************************************

#ifndef BAMTOOLS_UTILITIES_H
#define BAMTOOLS_UTILITIES_H

#include <api/BamAux.h>
#include <utils/utils_global.h>
#include <string>
#include <vector>

#define BAMTOOLS_ASSERT_UNREACHABLE BT_ASSERT_UNREACHABLE
#define BAMTOOLS_ASSERT_MESSAGE( condition, message ) BT_ASSERT_X( condition, message )

namespace BamTools {

class BamReader;
class BamMultiReader;

class UTILS_EXPORT Utilities {
  
    public: 
        // returns true if 'source' contains 'pattern' or 'c'
        static bool Contains(const std::string& source, const std::string& pattern);
        static bool Contains(const std::string& source, const char c);

        // returns true if 'source' ends with 'pattern' or 'c'
        static bool EndsWith(const std::string& source, const std::string& pattern);
        static bool EndsWith(const std::string& source, const char c);

        // check if a file exists
        static bool FileExists(const std::string& fname);
        
        // Parses a region string, uses reader to do validation (valid ID's, positions), stores in Region struct
        // Returns success (true/false)
        static bool ParseRegionString(const std::string& regionString,
                                      const BamReader& reader,
                                      BamRegion& region);
        // Same as above, but accepts a BamMultiReader
        static bool ParseRegionString(const std::string& regionString,
                                      const BamMultiReader& reader,
                                      BamRegion& region);

        // sequence utilities
        static void Reverse(std::string& sequence);
        static void ReverseComplement(std::string& sequence);

        // split string on delimiter character (or string of allowed delimiters)
        static std::vector<std::string> Split(const std::string& source, const char delim);
        static std::vector<std::string> Split(const std::string& source, const std::string& delims);

        // returns true if 'source' starts with 'pattern' or 'c'
        static bool StartsWith(const std::string& source, const std::string& pattern);
        static bool StartsWith(const std::string &source, const char c);
};

} // namespace BamTools
  
#endif // BAMTOOLS_UTILITIES_H
