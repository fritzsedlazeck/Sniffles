// ***************************************************************************
// bamtools_utilities.cpp (c) 2010 Derek Barnett, Erik Garrison
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 8 October 2011
// ---------------------------------------------------------------------------
// Provides general utilities used by BamTools sub-tools.
// ***************************************************************************

#include <api/BamMultiReader.h>
#include <api/BamReader.h>
#include <utils/bamtools_utilities.h>
using namespace BamTools;

#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
using namespace std;

namespace BamTools {
  
const char REVCOMP_LOOKUP[] = {'T',  0,  'G', 'H',
                                0,   0,  'C', 'D',
                                0,   0,   0,   0,
                               'K', 'N',  0,   0,
                                0,  'Y', 'W', 'A',
                               'A', 'B', 'S', 'X',
                               'R',  0 };
  
} // namespace BamTools 
  
// returns true if 'source' contains 'pattern'
bool Utilities::Contains(const string& source, const string& pattern) {
    return ( source.find(pattern) != string::npos );
}

// returns true if 'source' contains 'c'
bool Utilities::Contains(const std::string &source, const char c) {
    return ( source.find(c) != string::npos );
}

// returns true if 'source' ends with 'pattern'
bool Utilities::EndsWith(const string& source, const string& pattern) {
    return ( source.find(pattern) == (source.length() - pattern.length()) );
}

// returns true if 'source' ends with 'c'
bool Utilities::EndsWith(const std::string& source, const char c) {
    return ( source.find(c) == (source.length() - 1) );
}

// check if a file exists
bool Utilities::FileExists(const string& filename) {
    ifstream f(filename.c_str(), ifstream::in);
    return !f.fail();
}

// Parses a region string, does validation (valid ID's, positions), stores in Region struct
// Returns success (true/false)
bool Utilities::ParseRegionString(const string& regionString,
                                  const BamReader& reader,
                                  BamRegion& region)
{
    // -------------------------------
    // parse region string
  
    // check first for empty string
    if ( regionString.empty() ) 
        return false;   
    
    // non-empty string, look for a colom
    size_t foundFirstColon = regionString.find(':');
    
    // store chrom strings, and numeric positions
    string startChrom;
    string stopChrom;
    int startPos;
    int stopPos;
    
    // no colon found
    // going to use entire contents of requested chromosome 
    // just store entire region string as startChrom name
    // use BamReader methods to check if its valid for current BAM file
    if ( foundFirstColon == string::npos ) {
        startChrom = regionString;
        startPos   = 0;
        stopChrom  = regionString;
        stopPos    = 0;
    }
    
    // colon found, so we at least have some sort of startPos requested
    else {
      
        // store start chrom from beginning to first colon
        startChrom = regionString.substr(0,foundFirstColon);
        
        // look for ".." after the colon
        size_t foundRangeDots = regionString.find("..", foundFirstColon+1);
        
        // no dots found
        // so we have a startPos but no range
        // store contents before colon as startChrom, after as startPos
        if ( foundRangeDots == string::npos ) {
            startPos   = atoi( regionString.substr(foundFirstColon+1).c_str() ); 
            stopChrom  = startChrom;
            stopPos    = -1;
        } 
        
        // ".." found, so we have some sort of range selected
        else {
          
            // store startPos between first colon and range dots ".."
            startPos = atoi( regionString.substr(foundFirstColon+1, foundRangeDots-foundFirstColon-1).c_str() );
          
            // look for second colon
            size_t foundSecondColon = regionString.find(':', foundRangeDots+1);
            
            // no second colon found
            // so we have a "standard" chrom:start..stop input format (on single chrom)
            if ( foundSecondColon == string::npos ) {
                stopChrom  = startChrom;
                stopPos    = atoi( regionString.substr(foundRangeDots+2).c_str() );
            }
            
            // second colon found
            // so we have a range requested across 2 chrom's
            else {
                stopChrom  = regionString.substr(foundRangeDots+2, foundSecondColon-(foundRangeDots+2));
                stopPos    = atoi( regionString.substr(foundSecondColon+1).c_str() );
            }
        }
    }

    // -------------------------------
    // validate reference IDs & genomic positions
    
    const RefVector references = reader.GetReferenceData();
    
    // if startRefID not found, return false
    int startRefID = reader.GetReferenceID(startChrom);
    if ( startRefID == -1 ) return false;
    
    // startPos cannot be greater than or equal to reference length
    const RefData& startReference = references.at(startRefID);
    if ( startPos >= startReference.RefLength ) return false;
    
    // if stopRefID not found, return false
    int stopRefID = reader.GetReferenceID(stopChrom);
    if ( stopRefID == -1 ) return false;
    
    // stopPosition cannot be larger than reference length
    const RefData& stopReference = references.at(stopRefID);
    if ( stopPos > stopReference.RefLength ) return false;
    
    // if no stopPosition specified, set to reference end
    if ( stopPos == -1 ) stopPos = stopReference.RefLength;  
    
    // -------------------------------
    // set up Region struct & return
    
    region.LeftRefID     = startRefID;
    region.LeftPosition  = startPos;
    region.RightRefID    = stopRefID;;
    region.RightPosition = stopPos;
    return true;
}

// Same as ParseRegionString() above, but accepts a BamMultiReader
bool Utilities::ParseRegionString(const string& regionString,
                                  const BamMultiReader& reader,
                                  BamRegion& region)
{
    // -------------------------------
    // parse region string
  
    // check first for empty string
    if ( regionString.empty() ) 
        return false;   
    
    // non-empty string, look for a colom
    size_t foundFirstColon = regionString.find(':');
    
    // store chrom strings, and numeric positions
    string startChrom;
    string stopChrom;
    int startPos;
    int stopPos;
    
    // no colon found
    // going to use entire contents of requested chromosome 
    // just store entire region string as startChrom name
    // use BamReader methods to check if its valid for current BAM file
    if ( foundFirstColon == string::npos ) {
        startChrom = regionString;
        startPos   = 0;
        stopChrom  = regionString;
        stopPos    = -1;
    }
    
    // colon found, so we at least have some sort of startPos requested
    else {
      
        // store start chrom from beginning to first colon
        startChrom = regionString.substr(0,foundFirstColon);
        
        // look for ".." after the colon
        size_t foundRangeDots = regionString.find("..", foundFirstColon+1);
        
        // no dots found
        // so we have a startPos but no range
        // store contents before colon as startChrom, after as startPos
        if ( foundRangeDots == string::npos ) {
            startPos   = atoi( regionString.substr(foundFirstColon+1).c_str() ); 
            stopChrom  = startChrom;
            stopPos    = -1;
        } 
        
        // ".." found, so we have some sort of range selected
        else {
          
            // store startPos between first colon and range dots ".."
            startPos = atoi( regionString.substr(foundFirstColon+1, foundRangeDots-foundFirstColon-1).c_str() );
          
            // look for second colon
            size_t foundSecondColon = regionString.find(':', foundRangeDots+1);
            
            // no second colon found
            // so we have a "standard" chrom:start..stop input format (on single chrom)
            if ( foundSecondColon == string::npos ) {
                stopChrom  = startChrom;
                stopPos    = atoi( regionString.substr(foundRangeDots+2).c_str() );
            }
            
            // second colon found
            // so we have a range requested across 2 chrom's
            else {
                stopChrom  = regionString.substr(foundRangeDots+2, foundSecondColon-(foundRangeDots+2));
                stopPos    = atoi( regionString.substr(foundSecondColon+1).c_str() );
            }
        }
    }

    // -------------------------------
    // validate reference IDs & genomic positions

    const RefVector references = reader.GetReferenceData();

    // if startRefID not found, return false
    int startRefID = reader.GetReferenceID(startChrom);
    if ( startRefID == -1 ) return false;

    // startPos cannot be greater than or equal to reference length
    const RefData& startReference = references.at(startRefID);
    if ( startPos >= startReference.RefLength ) return false;

    // if stopRefID not found, return false
    int stopRefID = reader.GetReferenceID(stopChrom);
    if ( stopRefID == -1 ) return false;

    // stopPosition cannot be larger than reference length
    const RefData& stopReference = references.at(stopRefID);
    if ( stopPos > stopReference.RefLength ) return false;

    // if no stopPosition specified, set to reference end
    if ( stopPos == -1 ) stopPos = stopReference.RefLength;

    // -------------------------------
    // set up Region struct & return

    region.LeftRefID     = startRefID;
    region.LeftPosition  = startPos;
    region.RightRefID    = stopRefID;;
    region.RightPosition = stopPos;
    return true;
}

void Utilities::Reverse(string& sequence) {
    reverse(sequence.begin(), sequence.end());
}

void Utilities::ReverseComplement(string& sequence) {
    
    // do complement, in-place
    size_t seqLength = sequence.length();
    for ( size_t i = 0; i < seqLength; ++i )
        sequence.replace(i, 1, 1, REVCOMP_LOOKUP[(int)sequence.at(i) - 65]);
    
    // reverse it
    Reverse(sequence);
}

vector<string> Utilities::Split(const string& source, const char delim) {

    stringstream ss(source);
    string field;
    vector<string> fields;

    while ( getline(ss, field, delim) )
        fields.push_back(field);
    return fields;
}

vector<string> Utilities::Split(const string& source, const string& delims) {

    vector<string> fields;

    char* tok;
    char* cchars = new char[source.size()+1];
    char* cstr = &cchars[0];
    strcpy(cstr, source.c_str());
    tok = strtok(cstr, delims.c_str());
    while (tok != NULL) {
        fields.push_back(tok);
        tok = strtok(NULL, delims.c_str());
    }

    delete[] cchars;

    return fields;
}

// returns true if 'source' starts with 'pattern'
bool Utilities::StartsWith(const string& source, const string& pattern) {
    return ( source.find(pattern) == 0 );
}

// returns true if 'source' starts with 'c'
bool Utilities::StartsWith(const std::string &source, const char c) {
    return ( source.find(c) == 0 );
}
