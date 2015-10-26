// ***************************************************************************
// bamtools_filter.cpp (c) 2010 Derek Barnett, Erik Garrison
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 3 May 2013
// ---------------------------------------------------------------------------
// Filters BAM file(s) according to some user-specified criteria
// ***************************************************************************

#include "bamtools_filter.h"

#include <api/BamMultiReader.h>
#include <api/BamWriter.h>
#include <utils/bamtools_filter_engine.h>
#include <utils/bamtools_options.h>
#include <utils/bamtools_utilities.h>
using namespace BamTools;

#include <jsoncpp/json.h>
using namespace Json;

#include <cstdio>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
using namespace std;
  
namespace BamTools {
  
// -------------------------------  
// string literal constants  

// property names
const string ALIGNMENTFLAG_PROPERTY       = "alignmentFlag";
const string CIGAR_PROPERTY               = "cigar";
const string INSERTSIZE_PROPERTY          = "insertSize";
const string ISDUPLICATE_PROPERTY         = "isDuplicate";
const string ISFAILEDQC_PROPERTY          = "isFailedQC";
const string ISFIRSTMATE_PROPERTY         = "isFirstMate";
const string ISMAPPED_PROPERTY            = "isMapped";
const string ISMATEMAPPED_PROPERTY        = "isMateMapped";
const string ISMATEREVERSESTRAND_PROPERTY = "isMateReverseStrand";
const string ISPAIRED_PROPERTY            = "isPaired";
const string ISPRIMARYALIGNMENT_PROPERTY  = "isPrimaryAlignment";
const string ISPROPERPAIR_PROPERTY        = "isProperPair";
const string ISREVERSESTRAND_PROPERTY     = "isReverseStrand";
const string ISSECONDMATE_PROPERTY        = "isSecondMate";
const string ISSINGLETON_PROPERTY         = "isSingleton";
const string LENGTH_PROPERTY              = "length";
const string MAPQUALITY_PROPERTY          = "mapQuality";
const string MATEPOSITION_PROPERTY        = "matePosition";
const string MATEREFERENCE_PROPERTY       = "mateReference";
const string NAME_PROPERTY                = "name";
const string POSITION_PROPERTY            = "position";
const string QUERYBASES_PROPERTY          = "queryBases";
const string REFERENCE_PROPERTY           = "reference";
const string TAG_PROPERTY                 = "tag";

// boolalpha
const string TRUE_STR  = "true";
const string FALSE_STR = "false";
    
RefVector filterToolReferences;    
    
struct BamAlignmentChecker {
    bool check(const PropertyFilter& filter, const BamAlignment& al) {
      
        bool keepAlignment = true;
        const PropertyMap& properties = filter.Properties;
        PropertyMap::const_iterator propertyIter = properties.begin();
        PropertyMap::const_iterator propertyEnd  = properties.end();
        for ( ; propertyIter != propertyEnd; ++propertyIter ) {
          
            // check alignment data field depending on propertyName
            const string& propertyName = (*propertyIter).first;
            const PropertyFilterValue& valueFilter = (*propertyIter).second;
            
            if      ( propertyName == ALIGNMENTFLAG_PROPERTY )  keepAlignment &= valueFilter.check(al.AlignmentFlag);
            else if ( propertyName == CIGAR_PROPERTY ) {
                stringstream cigarSs;
                const vector<CigarOp>& cigarData = al.CigarData;
                if ( !cigarData.empty() ) {
                    vector<CigarOp>::const_iterator cigarBegin = cigarData.begin();
                    vector<CigarOp>::const_iterator cigarIter = cigarBegin;
                    vector<CigarOp>::const_iterator cigarEnd  = cigarData.end();
                    for ( ; cigarIter != cigarEnd; ++cigarIter ) {
                        const CigarOp& op = (*cigarIter);
                        cigarSs << op.Length << op.Type;
                    }
                    keepAlignment &= valueFilter.check(cigarSs.str());
                }
            }
            else if ( propertyName == INSERTSIZE_PROPERTY )           keepAlignment &= valueFilter.check(al.InsertSize);
            else if ( propertyName == ISDUPLICATE_PROPERTY )          keepAlignment &= valueFilter.check(al.IsDuplicate());
            else if ( propertyName == ISFAILEDQC_PROPERTY )           keepAlignment &= valueFilter.check(al.IsFailedQC());
            else if ( propertyName == ISFIRSTMATE_PROPERTY )          keepAlignment &= valueFilter.check(al.IsFirstMate());
            else if ( propertyName == ISMAPPED_PROPERTY )             keepAlignment &= valueFilter.check(al.IsMapped());
            else if ( propertyName == ISMATEMAPPED_PROPERTY )         keepAlignment &= valueFilter.check(al.IsMateMapped());
            else if ( propertyName == ISMATEREVERSESTRAND_PROPERTY )  keepAlignment &= valueFilter.check(al.IsMateReverseStrand());
            else if ( propertyName == ISPAIRED_PROPERTY )             keepAlignment &= valueFilter.check(al.IsPaired());
            else if ( propertyName == ISPRIMARYALIGNMENT_PROPERTY )   keepAlignment &= valueFilter.check(al.IsPrimaryAlignment());
            else if ( propertyName == ISPROPERPAIR_PROPERTY )         keepAlignment &= valueFilter.check(al.IsProperPair());
            else if ( propertyName == ISREVERSESTRAND_PROPERTY )      keepAlignment &= valueFilter.check(al.IsReverseStrand());
            else if ( propertyName == ISSECONDMATE_PROPERTY )         keepAlignment &= valueFilter.check(al.IsSecondMate());
            else if ( propertyName == ISSINGLETON_PROPERTY ) {
                const bool isSingleton = al.IsPaired() && al.IsMapped() && !al.IsMateMapped();
                keepAlignment &= valueFilter.check(isSingleton);
            }
            else if ( propertyName == LENGTH_PROPERTY )               keepAlignment &= valueFilter.check(al.Length);
            else if ( propertyName == MAPQUALITY_PROPERTY )           keepAlignment &= valueFilter.check(al.MapQuality);
            else if ( propertyName == MATEPOSITION_PROPERTY )         keepAlignment &= ( al.IsPaired() && al.IsMateMapped() && valueFilter.check(al.MateRefID) );
            else if ( propertyName == MATEREFERENCE_PROPERTY ) {
                if ( !al.IsPaired() || !al.IsMateMapped() ) return false;
                BAMTOOLS_ASSERT_MESSAGE( (al.MateRefID>=0 && (al.MateRefID<(int)filterToolReferences.size())), "Invalid MateRefID");
                const string& refName = filterToolReferences.at(al.MateRefID).RefName;
                keepAlignment &= valueFilter.check(refName);
            }
            else if ( propertyName == NAME_PROPERTY )                 keepAlignment &= valueFilter.check(al.Name);
            else if ( propertyName == POSITION_PROPERTY )             keepAlignment &= valueFilter.check(al.Position);
            else if ( propertyName == QUERYBASES_PROPERTY )           keepAlignment &= valueFilter.check(al.QueryBases);
            else if ( propertyName == REFERENCE_PROPERTY ) {
                BAMTOOLS_ASSERT_MESSAGE( (al.RefID>=0 && (al.RefID<(int)filterToolReferences.size())), "Invalid RefID");
                const string& refName = filterToolReferences.at(al.RefID).RefName;
                keepAlignment &= valueFilter.check(refName);
            }
            else if ( propertyName == TAG_PROPERTY ) keepAlignment &= checkAlignmentTag(valueFilter, al);
            else BAMTOOLS_ASSERT_UNREACHABLE;
            
            // if alignment fails at ANY point, just quit and return false
            if ( !keepAlignment ) return false;
        }
      
        BAMTOOLS_ASSERT_MESSAGE( keepAlignment, "Error in BamAlignmentChecker... keepAlignment should be true here");
        return keepAlignment;
    }
    
    bool checkAlignmentTag(const PropertyFilterValue& valueFilter, const BamAlignment& al) {
     
        // ensure filter contains string data
        Variant entireTagFilter = valueFilter.Value;
        if ( !entireTagFilter.is_type<string>() ) return false;

        // localize string from variant
        const string& entireTagFilterString = entireTagFilter.get<string>();

        // ensure we have at least "XX:x"
        if ( entireTagFilterString.length() < 4 ) return false;

        // get tagName & lookup in alignment
        // if found, set tagType to tag type character
        // if not found, return false
        const string& tagName = entireTagFilterString.substr(0,2);
        char tagType = '\0';
        if ( !al.GetTagType(tagName, tagType) ) return false;

        // remove tagName & ":" from beginning tagFilter
        string tagFilterString = entireTagFilterString.substr(3);

        // switch on tag type to set tag query value & parse filter token
        int8_t   asciiFilterValue,   asciiQueryValue;
        int32_t  intFilterValue,    intQueryValue;
        uint32_t uintFilterValue,   uintQueryValue;
        float    realFilterValue,   realQueryValue;
        string   stringFilterValue, stringQueryValue;

        PropertyFilterValue tagFilter;
        PropertyFilterValue::ValueCompareType compareType;
        bool keepAlignment = false;
        switch (tagType) {

            // ASCII tag type
            case 'A':
                if ( al.GetTag(tagName, asciiQueryValue) ) {
                    if ( FilterEngine<BamAlignmentChecker>::parseToken(tagFilterString, asciiFilterValue, compareType) ) {
                        tagFilter.Value = asciiFilterValue;
                        tagFilter.Type  = compareType;
                        keepAlignment   = tagFilter.check(asciiQueryValue);
                    }
                }
                break;

            // signed int tag type
            case 'c' :
            case 's' :
            case 'i' :
                if ( al.GetTag(tagName, intQueryValue) ) {
                    if ( FilterEngine<BamAlignmentChecker>::parseToken(tagFilterString, intFilterValue, compareType) ) {
                        tagFilter.Value = intFilterValue;
                        tagFilter.Type  = compareType;
                        keepAlignment   = tagFilter.check(intQueryValue);
                    }
                }
                break;

            // unsigned int tag type
            case 'C' :
            case 'S' :
            case 'I' :
                if ( al.GetTag(tagName, uintQueryValue) ) {
                    if ( FilterEngine<BamAlignmentChecker>::parseToken(tagFilterString, uintFilterValue, compareType) ) {
                        tagFilter.Value = uintFilterValue;
                        tagFilter.Type  = compareType;
                        keepAlignment   = tagFilter.check(uintQueryValue);
                    }
                }
                break;

            // 'real' tag type
            case 'f' :
                if ( al.GetTag(tagName, realQueryValue) ) {
                    if ( FilterEngine<BamAlignmentChecker>::parseToken(tagFilterString, realFilterValue, compareType) ) {
                        tagFilter.Value = realFilterValue;
                        tagFilter.Type  = compareType;
                        keepAlignment   = tagFilter.check(realQueryValue);
                    }
                }
                break;

            // string tag type

            case 'Z':
            case 'H':
                if ( al.GetTag(tagName, stringQueryValue) ) {
                    if ( FilterEngine<BamAlignmentChecker>::parseToken(tagFilterString, stringFilterValue, compareType) ) {
                        tagFilter.Value = stringFilterValue;
                        tagFilter.Type  = compareType;
                        keepAlignment   = tagFilter.check(stringQueryValue);
                    }
                }
                break;

            // unknown tag type
            default :
                keepAlignment = false;
        }

        return keepAlignment;
    }
};    
    
} // namespace BamTools
  
// ---------------------------------------------
// FilterSettings implementation

struct FilterTool::FilterSettings {

    // ----------------------------------
    // IO opts

    // flags
    bool HasInput;
    bool HasInputFilelist;
    bool HasOutput;
    bool HasRegion;
    bool HasScript;
    bool IsForceCompression;

    // filenames
    vector<string> InputFiles;
    string InputFilelist;
    string OutputFilename;
    string Region;
    string ScriptFilename;

    // -----------------------------------
    // General filter opts

    // flags
    bool HasAlignmentFlagFilter;
    bool HasInsertSizeFilter;
    bool HasLengthFilter;
    bool HasMapQualityFilter;
    bool HasNameFilter;
    bool HasQueryBasesFilter;
    bool HasTagFilter; //(s)

    // filters
    string AlignmentFlagFilter;
    string InsertSizeFilter;
    string LengthFilter;
    string MapQualityFilter;
    string NameFilter;
    string QueryBasesFilter;
    string TagFilter;  // support multiple ?

    // -----------------------------------
    // AlignmentFlag filter opts

    // flags
    bool HasIsDuplicateFilter;
    bool HasIsFailedQCFilter;
    bool HasIsFirstMateFilter;
    bool HasIsMappedFilter;
    bool HasIsMateMappedFilter;
    bool HasIsMateReverseStrandFilter;
    bool HasIsPairedFilter;
    bool HasIsPrimaryAlignmentFilter;
    bool HasIsProperPairFilter;
    bool HasIsReverseStrandFilter;
    bool HasIsSecondMateFilter;
    bool HasIsSingletonFilter;

    // filters
    string IsDuplicateFilter;
    string IsFailedQCFilter;
    string IsFirstMateFilter;
    string IsMappedFilter;
    string IsMateMappedFilter;
    string IsMateReverseStrandFilter;
    string IsPairedFilter;
    string IsPrimaryAlignmentFilter;
    string IsProperPairFilter;
    string IsReverseStrandFilter;
    string IsSecondMateFilter;
    string IsSingletonFilter;

    // ---------------------------------
    // constructor

    FilterSettings(void)
        : HasInput(false)
        , HasInputFilelist(false)
        , HasOutput(false)
        , HasRegion(false)
        , HasScript(false)
        , IsForceCompression(false)
        , OutputFilename(Options::StandardOut())
        , HasAlignmentFlagFilter(false)
        , HasInsertSizeFilter(false)
        , HasLengthFilter(false)
        , HasMapQualityFilter(false)
        , HasNameFilter(false)
        , HasQueryBasesFilter(false)
        , HasTagFilter(false)
        , HasIsDuplicateFilter(false)
        , HasIsFailedQCFilter(false)
        , HasIsFirstMateFilter(false)
        , HasIsMappedFilter(false)
        , HasIsMateMappedFilter(false)
        , HasIsMateReverseStrandFilter(false)
        , HasIsPairedFilter(false)
        , HasIsPrimaryAlignmentFilter(false)
        , HasIsProperPairFilter(false)
        , HasIsReverseStrandFilter(false)
        , HasIsSecondMateFilter(false)
        , HasIsSingletonFilter(false)
        , IsDuplicateFilter(TRUE_STR)
        , IsFailedQCFilter(TRUE_STR)
        , IsFirstMateFilter(TRUE_STR)
        , IsMappedFilter(TRUE_STR)
        , IsMateMappedFilter(TRUE_STR)
        , IsMateReverseStrandFilter(TRUE_STR)
        , IsPairedFilter(TRUE_STR)
        , IsPrimaryAlignmentFilter(TRUE_STR)
        , IsProperPairFilter(TRUE_STR)
        , IsReverseStrandFilter(TRUE_STR)
        , IsSecondMateFilter(TRUE_STR)
        , IsSingletonFilter(TRUE_STR)
    { }
};

// ---------------------------------------------
// FilterToolPrivate declaration

class FilterTool::FilterToolPrivate {
      
    // ctor & dtor
    public:
        FilterToolPrivate(FilterTool::FilterSettings* settings);
        ~FilterToolPrivate(void);  
        
    // 'public' interface
    public:
        bool Run(void);
        
    // internal methods
    private:
        bool AddPropertyTokensToFilter(const string& filterName, const map<string, string>& propertyTokens);
        bool CheckAlignment(const BamAlignment& al);
        const string GetScriptContents(void);
        void InitProperties(void);
        bool ParseCommandLine(void);
        bool ParseFilterObject(const string& filterName, const Json::Value& filterObject);
        bool ParseScript(void);
        bool SetupFilters(void);
        
    // data members
    private:
        vector<string> m_propertyNames;
        FilterTool::FilterSettings* m_settings;
        FilterEngine<BamAlignmentChecker> m_filterEngine;
};
 
// ---------------------------------------------
// FilterToolPrivate implementation
  
// constructor  
FilterTool::FilterToolPrivate::FilterToolPrivate(FilterTool::FilterSettings* settings)
    : m_settings(settings)
{ }  
  
// destructor
FilterTool::FilterToolPrivate::~FilterToolPrivate(void) { }

bool FilterTool::FilterToolPrivate::AddPropertyTokensToFilter(const string& filterName,
                                                              const map<string,
                                                              string>& propertyTokens)
{
    // dummy temp values for token parsing
    bool boolValue;
    int32_t int32Value;
    uint16_t uint16Value;
    uint32_t uint32Value;
    string stringValue;
    PropertyFilterValue::ValueCompareType type;
  
    // iterate over property token map
    map<string, string>::const_iterator mapIter = propertyTokens.begin();
    map<string, string>::const_iterator mapEnd  = propertyTokens.end();
    for ( ; mapIter != mapEnd; ++mapIter ) {
      
        const string& propertyName = (*mapIter).first;
        const string& token        = (*mapIter).second;
        
        // ------------------------------
        // convert token to value & compare type 
        // then add to filter engine
        
        // bool conversion
        if ( propertyName == ISDUPLICATE_PROPERTY ||
             propertyName == ISFAILEDQC_PROPERTY ||
             propertyName == ISFIRSTMATE_PROPERTY ||
             propertyName == ISMAPPED_PROPERTY ||
             propertyName == ISMATEMAPPED_PROPERTY ||
             propertyName == ISMATEREVERSESTRAND_PROPERTY ||
             propertyName == ISPAIRED_PROPERTY ||
             propertyName == ISPRIMARYALIGNMENT_PROPERTY ||
             propertyName == ISPROPERPAIR_PROPERTY ||
             propertyName == ISREVERSESTRAND_PROPERTY ||
             propertyName == ISSECONDMATE_PROPERTY ||
             propertyName == ISSINGLETON_PROPERTY
           ) 
        {
            FilterEngine<BamAlignmentChecker>::parseToken(token, boolValue, type);
            m_filterEngine.setProperty(filterName, propertyName, boolValue, type);
        }
        
        // int32_t conversion
        else if ( propertyName == INSERTSIZE_PROPERTY ||
                  propertyName == LENGTH_PROPERTY ||
                  propertyName == MATEPOSITION_PROPERTY ||
                  propertyName == POSITION_PROPERTY 
                ) 
        {
            FilterEngine<BamAlignmentChecker>::parseToken(token, int32Value, type);
            m_filterEngine.setProperty(filterName, propertyName, int32Value, type);
        }
        
        // uint16_t conversion
        else if ( propertyName == MAPQUALITY_PROPERTY )
        {
            FilterEngine<BamAlignmentChecker>::parseToken(token, uint16Value, type);
            m_filterEngine.setProperty(filterName, propertyName, uint16Value, type);
        }
        
        // uint32_t conversion
        else if ( propertyName == ALIGNMENTFLAG_PROPERTY )
        {
            FilterEngine<BamAlignmentChecker>::parseToken(token, uint32Value, type);
            m_filterEngine.setProperty(filterName, propertyName, uint32Value, type);
        }
        
        // string conversion
        else if ( propertyName == CIGAR_PROPERTY || 
                  propertyName == MATEREFERENCE_PROPERTY ||
                  propertyName == NAME_PROPERTY ||
                  propertyName == QUERYBASES_PROPERTY ||
                  propertyName == REFERENCE_PROPERTY 
                ) 
        {
            FilterEngine<BamAlignmentChecker>::parseToken(token, stringValue, type);
            m_filterEngine.setProperty(filterName, propertyName, stringValue, type);
        }
      
        else if ( propertyName == TAG_PROPERTY ) {
            // this will be stored directly as the TAG:VALUE token
            // (VALUE may contain compare ops, will be parsed out later)
            m_filterEngine.setProperty(filterName, propertyName, token, PropertyFilterValue::EXACT);
        }
      
        // else unknown property 
        else {
            cerr << "bamtools filter ERROR: unknown property - " << propertyName << endl;
            return false;
        }
    }
    return true;
}

bool FilterTool::FilterToolPrivate::CheckAlignment(const BamAlignment& al) {
    return m_filterEngine.check(al);
}

const string FilterTool::FilterToolPrivate::GetScriptContents(void) {
  
    // open script for reading
    FILE* inFile = fopen(m_settings->ScriptFilename.c_str(), "rb");
    if ( !inFile ) {
        cerr << "bamtools filter ERROR: could not open script: "
             << m_settings->ScriptFilename << " for reading" << endl;
        return string();
    }
    
    // read in entire script contents  
    char buffer[1024];
    ostringstream docStream("");
    while ( true ) {
        
        // peek ahead, make sure there is data available
        char ch = fgetc(inFile);
        ungetc(ch, inFile);
        if( feof(inFile) )
            break;
        
        // read next block of data
        if ( fgets(buffer, 1024, inFile) == 0 ) {
            cerr << "bamtools filter ERROR: could not read script contents" << endl;
            return string();
        }
        
        docStream << buffer;
    }
    
    // close script file
    fclose(inFile);
    
    // import buffer contents to document, return
    return docStream.str();
}

void FilterTool::FilterToolPrivate::InitProperties(void) {
  
    // store property names in vector 
    m_propertyNames.push_back(ALIGNMENTFLAG_PROPERTY);
    m_propertyNames.push_back(CIGAR_PROPERTY);
    m_propertyNames.push_back(INSERTSIZE_PROPERTY);
    m_propertyNames.push_back(ISDUPLICATE_PROPERTY);
    m_propertyNames.push_back(ISFAILEDQC_PROPERTY);
    m_propertyNames.push_back(ISFIRSTMATE_PROPERTY);
    m_propertyNames.push_back(ISMAPPED_PROPERTY);
    m_propertyNames.push_back(ISMATEMAPPED_PROPERTY);
    m_propertyNames.push_back(ISMATEREVERSESTRAND_PROPERTY);
    m_propertyNames.push_back(ISPAIRED_PROPERTY);
    m_propertyNames.push_back(ISPRIMARYALIGNMENT_PROPERTY);
    m_propertyNames.push_back(ISPROPERPAIR_PROPERTY);
    m_propertyNames.push_back(ISREVERSESTRAND_PROPERTY);
    m_propertyNames.push_back(ISSECONDMATE_PROPERTY);
    m_propertyNames.push_back(ISSINGLETON_PROPERTY);
    m_propertyNames.push_back(LENGTH_PROPERTY);
    m_propertyNames.push_back(MAPQUALITY_PROPERTY);
    m_propertyNames.push_back(MATEPOSITION_PROPERTY);
    m_propertyNames.push_back(MATEREFERENCE_PROPERTY);
    m_propertyNames.push_back(NAME_PROPERTY);
    m_propertyNames.push_back(POSITION_PROPERTY);
    m_propertyNames.push_back(QUERYBASES_PROPERTY);
    m_propertyNames.push_back(REFERENCE_PROPERTY);
    m_propertyNames.push_back(TAG_PROPERTY);
    
    // add vector contents to FilterEngine<BamAlignmentChecker>
    vector<string>::const_iterator propertyNameIter = m_propertyNames.begin();
    vector<string>::const_iterator propertyNameEnd  = m_propertyNames.end();
    for ( ; propertyNameIter != propertyNameEnd; ++propertyNameIter )
        m_filterEngine.addProperty((*propertyNameIter));
}

bool FilterTool::FilterToolPrivate::ParseCommandLine(void) {
  
    // add a rule set to filter engine
    const string CMD = "COMMAND_LINE";
    m_filterEngine.addFilter(CMD);

    // map property names to command line args
    map<string, string> propertyTokens;
    if ( m_settings->HasAlignmentFlagFilter )       propertyTokens.insert( make_pair(ALIGNMENTFLAG_PROPERTY,       m_settings->AlignmentFlagFilter) );
    if ( m_settings->HasInsertSizeFilter )          propertyTokens.insert( make_pair(INSERTSIZE_PROPERTY,          m_settings->InsertSizeFilter) );
    if ( m_settings->HasIsDuplicateFilter )         propertyTokens.insert( make_pair(ISDUPLICATE_PROPERTY,         m_settings->IsDuplicateFilter) );
    if ( m_settings->HasIsFailedQCFilter )          propertyTokens.insert( make_pair(ISFAILEDQC_PROPERTY,          m_settings->IsFailedQCFilter) );
    if ( m_settings->HasIsFirstMateFilter )         propertyTokens.insert( make_pair(ISFIRSTMATE_PROPERTY,         m_settings->IsFirstMateFilter) );
    if ( m_settings->HasIsMappedFilter )            propertyTokens.insert( make_pair(ISMAPPED_PROPERTY,            m_settings->IsMappedFilter) );
    if ( m_settings->HasIsMateMappedFilter )        propertyTokens.insert( make_pair(ISMATEMAPPED_PROPERTY,        m_settings->IsMateMappedFilter) );
    if ( m_settings->HasIsMateReverseStrandFilter ) propertyTokens.insert( make_pair(ISMATEREVERSESTRAND_PROPERTY, m_settings->IsMateReverseStrandFilter) );
    if ( m_settings->HasIsPairedFilter )            propertyTokens.insert( make_pair(ISPAIRED_PROPERTY,            m_settings->IsPairedFilter) );
    if ( m_settings->HasIsPrimaryAlignmentFilter )  propertyTokens.insert( make_pair(ISPRIMARYALIGNMENT_PROPERTY,  m_settings->IsPrimaryAlignmentFilter) );
    if ( m_settings->HasIsProperPairFilter )        propertyTokens.insert( make_pair(ISPROPERPAIR_PROPERTY,        m_settings->IsProperPairFilter) );
    if ( m_settings->HasIsReverseStrandFilter )     propertyTokens.insert( make_pair(ISREVERSESTRAND_PROPERTY,     m_settings->IsReverseStrandFilter) );
    if ( m_settings->HasIsSecondMateFilter )        propertyTokens.insert( make_pair(ISSECONDMATE_PROPERTY,        m_settings->IsSecondMateFilter) );
    if ( m_settings->HasIsSingletonFilter )         propertyTokens.insert( make_pair(ISSINGLETON_PROPERTY,         m_settings->IsSingletonFilter) );
    if ( m_settings->HasLengthFilter )              propertyTokens.insert( make_pair(LENGTH_PROPERTY,              m_settings->LengthFilter) );
    if ( m_settings->HasMapQualityFilter )          propertyTokens.insert( make_pair(MAPQUALITY_PROPERTY,          m_settings->MapQualityFilter) );
    if ( m_settings->HasNameFilter )                propertyTokens.insert( make_pair(NAME_PROPERTY,                m_settings->NameFilter) );
    if ( m_settings->HasQueryBasesFilter )          propertyTokens.insert( make_pair(QUERYBASES_PROPERTY,          m_settings->QueryBasesFilter) );
    if ( m_settings->HasTagFilter )                 propertyTokens.insert( make_pair(TAG_PROPERTY,                 m_settings->TagFilter) );
    
    // send add these properties to filter set "COMMAND_LINE"
    return AddPropertyTokensToFilter(CMD, propertyTokens);
}

bool FilterTool::FilterToolPrivate::ParseFilterObject(const string& filterName, const Json::Value& filterObject) {
  
    // filter object parsing variables
    Json::Value null(Json::nullValue);
    Json::Value propertyValue;
    
    // store results
    map<string, string> propertyTokens;
    
    // iterate over known properties
    vector<string>::const_iterator propertyNameIter = m_propertyNames.begin();
    vector<string>::const_iterator propertyNameEnd  = m_propertyNames.end();
    for ( ; propertyNameIter != propertyNameEnd; ++propertyNameIter ) {
        const string& propertyName = (*propertyNameIter);
        
        // if property defined in filter, add to token list
        propertyValue = filterObject.get(propertyName, null);
        if ( propertyValue != null ) 
            propertyTokens.insert( make_pair(propertyName, propertyValue.asString()) );
    }
  
    // add this filter to engin
    m_filterEngine.addFilter(filterName);
  
    // add token list to this filter
    return AddPropertyTokensToFilter(filterName, propertyTokens);
}

bool FilterTool::FilterToolPrivate::ParseScript(void) {
  
    // read in script contents from file
    const string document = GetScriptContents();
    
    // set up JsonCPP reader and attempt to parse script
    Json::Value root;
    Json::Reader reader;
    if ( !reader.parse(document, root) ) {
        // use built-in error reporting mechanism to alert user what was wrong with the script
        cerr  << "bamtools filter ERROR: failed to parse script - see error message(s) below" << endl
              << reader.getFormatedErrorMessages();
        return false;     
    }

    // initialize return status
    bool success = true;
    
    // see if root object contains multiple filters
    const Json::Value filters = root["filters"];
    if ( !filters.isNull() ) {
      
        // iterate over any filters found
        int filterIndex = 0;
        Json::Value::const_iterator filtersIter = filters.begin();
        Json::Value::const_iterator filtersEnd  = filters.end();
        for ( ; filtersIter != filtersEnd; ++filtersIter, ++filterIndex ) {
            Json::Value filter = (*filtersIter);
            
            // convert filter index to string
            string filterName;
            
            // if id tag supplied
            const Json::Value id = filter["id"];
            if ( !id.isNull() ) 
                filterName = id.asString();
            
            // use array index 
            else {
                stringstream convert;
                convert << filterIndex;
                filterName = convert.str();
            }
            
            // create & parse filter 
            success &= ParseFilterObject(filterName, filter);
        }
        
        // see if user defined a "rule" for these filters
        // otherwise, use filter engine's default rule behavior
        string ruleString("");
        const Json::Value rule = root["rule"];
        if ( rule.isString() )
            ruleString = rule.asString();
        m_filterEngine.setRule(ruleString);
          
        // return success/fail
        return success;
    } 
    
    // otherwise, root is the only filter (just contains properties)
    // create & parse filter named "ROOT"
    else success = ParseFilterObject("ROOT", root);
    
    // return success/failure
    return success;
}

bool FilterTool::FilterToolPrivate::Run(void) {
    
    // set to default input if none provided
    if ( !m_settings->HasInput && !m_settings->HasInputFilelist )
        m_settings->InputFiles.push_back(Options::StandardIn());

    // add files in the filelist to the input file list
    if ( m_settings->HasInputFilelist ) {

        ifstream filelist(m_settings->InputFilelist.c_str(), ios::in);
        if ( !filelist.is_open() ) {
            cerr << "bamtools filter ERROR: could not open input BAM file list... Aborting." << endl;
            return false;
        }

        string line;
        while ( getline(filelist, line) )
            m_settings->InputFiles.push_back(line);
    }

    // initialize defined properties & user-specified filters
    // quit if failed
    if ( !SetupFilters() )
        return false;

    // open reader without index
    BamMultiReader reader;
    if ( !reader.Open(m_settings->InputFiles) ) {
        cerr << "bamtools filter ERROR: could not open input files for reading." << endl;
        return false;
    }

    // retrieve reader header & reference data
    const string headerText = reader.GetHeaderText();
    filterToolReferences = reader.GetReferenceData();
    
    // determine compression mode for BamWriter
    bool writeUncompressed = ( m_settings->OutputFilename == Options::StandardOut() &&
                              !m_settings->IsForceCompression );
    BamWriter::CompressionMode compressionMode = BamWriter::Compressed;
    if ( writeUncompressed ) compressionMode = BamWriter::Uncompressed;

    // open BamWriter
    BamWriter writer;
    writer.SetCompressionMode(compressionMode);
    if ( !writer.Open(m_settings->OutputFilename, headerText, filterToolReferences) ) {
        cerr << "bamtools filter ERROR: could not open " << m_settings->OutputFilename << " for writing." << endl;
        reader.Close();
        return false;
    }

    // if no region specified, filter entire file 
    BamAlignment al;
    if ( !m_settings->HasRegion ) {
        while ( reader.GetNextAlignment(al) ) {
            if ( CheckAlignment(al) ) 
                writer.SaveAlignment(al);
        }
    }
    
    // otherwise attempt to use region as constraint
    else {
        
        // if region string parses OK
        BamRegion region;
        if ( Utilities::ParseRegionString(m_settings->Region, reader, region) ) {

            // attempt to find index files
            reader.LocateIndexes();

            // if index data available for all BAM files, we can use SetRegion
            if ( reader.HasIndexes() ) {

                // attempt to use SetRegion(), if failed report error
                if ( !reader.SetRegion(region.LeftRefID, region.LeftPosition, region.RightRefID, region.RightPosition) ) {
                    cerr << "bamtools filter ERROR: set region failed. Check that REGION describes a valid range" << endl;
                    reader.Close();
                    return false;
                } 
              
                // everything checks out, just iterate through specified region, filtering alignments
                while ( reader.GetNextAlignment(al) )
                    if ( CheckAlignment(al) ) 
                        writer.SaveAlignment(al);
                }
            
            // no index data available, we have to iterate through until we
            // find overlapping alignments
            else {
                while ( reader.GetNextAlignment(al) ) {
                    if ( (al.RefID >= region.LeftRefID)  && ((al.Position + al.Length) >= region.LeftPosition) &&
                         (al.RefID <= region.RightRefID) && ( al.Position <= region.RightPosition) ) 
                    {
                        if ( CheckAlignment(al) ) 
                            writer.SaveAlignment(al);
                    }
                }
            }
        } 
        
        // error parsing REGION string
        else {
            cerr << "bamtools filter ERROR: could not parse REGION: " << m_settings->Region << endl;
            cerr << "Check that REGION is in valid format (see documentation) and that the coordinates are valid"
                 << endl;
            reader.Close();
            return false;
        }
    }

    // clean up & exit
    reader.Close();
    writer.Close();
    return true;
}

bool FilterTool::FilterToolPrivate::SetupFilters(void) {
  
    // set up filter engine with supported properties
    InitProperties();
    
    // parse script for filter rules, if given
    if ( m_settings->HasScript )
        return ParseScript();
    
    // otherwise check command line for filters
    else return ParseCommandLine();
}

// ---------------------------------------------
// FilterTool implementation

FilterTool::FilterTool(void)
    : AbstractTool()
    , m_settings(new FilterSettings)
    , m_impl(0)
{
    // ----------------------------------
    // set program details

    const string usage = "[-in <filename> -in <filename> ... | -list <filelist>] "
                         "[-out <filename> | [-forceCompression]] [-region <REGION>] "
                         "[ [-script <filename] | [filterOptions] ]";

    Options::SetProgramInfo("bamtools filter", "filters BAM file(s)", usage );

    // ----------------------------------
    // I/O options

    OptionGroup* IO_Opts = Options::CreateOptionGroup("Input & Output");

    const string inDesc     = "the input BAM file(s)";
    const string listDesc   = "the input BAM file list, one line per file";
    const string outDesc    = "the output BAM file";
    const string regionDesc = "only read data from this genomic region (see documentation for more details)";
    const string scriptDesc = "the filter script file (see documentation for more details)";
    const string forceDesc  = "if results are sent to stdout (like when piping to another tool), "
                              "default behavior is to leave output uncompressed. Use this flag to "
                              "override and force compression";

    Options::AddValueOption("-in",     "BAM filename", inDesc,     "", m_settings->HasInput,  m_settings->InputFiles,     IO_Opts, Options::StandardIn());
    Options::AddValueOption("-list",   "filename",     listDesc,   "", m_settings->HasInputFilelist,  m_settings->InputFilelist, IO_Opts);
    Options::AddValueOption("-out",    "BAM filename", outDesc,    "", m_settings->HasOutput, m_settings->OutputFilename, IO_Opts, Options::StandardOut());
    Options::AddValueOption("-region", "REGION",       regionDesc, "", m_settings->HasRegion, m_settings->Region,         IO_Opts);
    Options::AddValueOption("-script", "filename",     scriptDesc, "", m_settings->HasScript, m_settings->ScriptFilename, IO_Opts);
    Options::AddOption("-forceCompression",forceDesc, m_settings->IsForceCompression, IO_Opts);

    // ----------------------------------
    // general filter options

    OptionGroup* FilterOpts = Options::CreateOptionGroup("General Filters");

    const string flagDesc    = "keep reads with this *exact* alignment flag (for more detailed queries, see below)";
    const string insertDesc  = "keep reads with insert size that matches pattern";
    const string lengthDesc  = "keep reads with length that matches pattern";
    const string mapQualDesc = "keep reads with map quality that matches pattern";
    const string nameDesc    = "keep reads with name that matches pattern";
    const string queryDesc   = "keep reads with motif that matches pattern";
    const string tagDesc     = "keep reads with this key=>value pair";

    Options::AddValueOption("-alignmentFlag", "int",       flagDesc,    "", m_settings->HasAlignmentFlagFilter, m_settings->AlignmentFlagFilter, FilterOpts);
    Options::AddValueOption("-insertSize",    "int",       insertDesc,  "", m_settings->HasInsertSizeFilter,    m_settings->InsertSizeFilter,    FilterOpts);
    Options::AddValueOption("-length",        "int",       lengthDesc,  "", m_settings->HasLengthFilter,        m_settings->LengthFilter,        FilterOpts);
    Options::AddValueOption("-mapQuality",    "[0-255]",   mapQualDesc, "", m_settings->HasMapQualityFilter,    m_settings->MapQualityFilter,    FilterOpts);
    Options::AddValueOption("-name",          "string",    nameDesc,    "", m_settings->HasNameFilter,          m_settings->NameFilter,          FilterOpts);
    Options::AddValueOption("-queryBases",    "string",    queryDesc,   "", m_settings->HasQueryBasesFilter,    m_settings->QueryBasesFilter,    FilterOpts);
    Options::AddValueOption("-tag",           "TAG:VALUE", tagDesc,     "", m_settings->HasTagFilter,           m_settings->TagFilter,           FilterOpts);

    // ----------------------------------
    // alignment flag filter options

    OptionGroup* AlignmentFlagOpts = Options::CreateOptionGroup("Alignment Flag Filters");

    const string boolArg           = "true/false";
    const string isDupDesc         = "keep only alignments that are marked as duplicate?";
    const string isFailQcDesc      = "keep only alignments that failed QC?";
    const string isFirstMateDesc   = "keep only alignments marked as first mate?";
    const string isMappedDesc      = "keep only alignments that were mapped?";
    const string isMateMappedDesc  = "keep only alignments with mates that mapped";
    const string isMateReverseDesc = "keep only alignments with mate on reverese strand?";
    const string isPairedDesc      = "keep only alignments that were sequenced as paired?";
    const string isPrimaryDesc     = "keep only alignments marked as primary?";
    const string isProperPairDesc  = "keep only alignments that passed PE resolution?";
    const string isReverseDesc     = "keep only alignments on reverse strand?";
    const string isSecondMateDesc  = "keep only alignments marked as second mate?";
    const string isSingletonDesc   = "keep only singletons";

    Options::AddValueOption("-isDuplicate",         boolArg, isDupDesc,         "", m_settings->HasIsDuplicateFilter,         m_settings->IsDuplicateFilter,         AlignmentFlagOpts, TRUE_STR);
    Options::AddValueOption("-isFailedQC",          boolArg, isFailQcDesc,      "", m_settings->HasIsFailedQCFilter,          m_settings->IsFailedQCFilter,          AlignmentFlagOpts, TRUE_STR);
    Options::AddValueOption("-isFirstMate",         boolArg, isFirstMateDesc,   "", m_settings->HasIsFirstMateFilter,         m_settings->IsFirstMateFilter,         AlignmentFlagOpts, TRUE_STR);
    Options::AddValueOption("-isMapped",            boolArg, isMappedDesc,      "", m_settings->HasIsMappedFilter,            m_settings->IsMappedFilter,            AlignmentFlagOpts, TRUE_STR);
    Options::AddValueOption("-isMateMapped",        boolArg, isMateMappedDesc,  "", m_settings->HasIsMateMappedFilter,        m_settings->IsMateMappedFilter,        AlignmentFlagOpts, TRUE_STR);
    Options::AddValueOption("-isMateReverseStrand", boolArg, isMateReverseDesc, "", m_settings->HasIsMateReverseStrandFilter, m_settings->IsMateReverseStrandFilter, AlignmentFlagOpts, TRUE_STR);
    Options::AddValueOption("-isPaired",            boolArg, isPairedDesc,      "", m_settings->HasIsPairedFilter,            m_settings->IsPairedFilter,            AlignmentFlagOpts, TRUE_STR);
    Options::AddValueOption("-isPrimaryAlignment",  boolArg, isPrimaryDesc,     "", m_settings->HasIsPrimaryAlignmentFilter,  m_settings->IsPrimaryAlignmentFilter,  AlignmentFlagOpts, TRUE_STR);
    Options::AddValueOption("-isProperPair",        boolArg, isProperPairDesc,  "", m_settings->HasIsProperPairFilter,        m_settings->IsProperPairFilter,        AlignmentFlagOpts, TRUE_STR);
    Options::AddValueOption("-isReverseStrand",     boolArg, isReverseDesc,     "", m_settings->HasIsReverseStrandFilter,     m_settings->IsReverseStrandFilter,     AlignmentFlagOpts, TRUE_STR);
    Options::AddValueOption("-isSecondMate",        boolArg, isSecondMateDesc,  "", m_settings->HasIsSecondMateFilter,        m_settings->IsSecondMateFilter,        AlignmentFlagOpts, TRUE_STR);
    Options::AddValueOption("-isSingleton",         boolArg, isSingletonDesc,   "", m_settings->HasIsSingletonFilter,         m_settings->IsSingletonFilter,         AlignmentFlagOpts, TRUE_STR);
}

FilterTool::~FilterTool(void) {

    delete m_settings;
    m_settings = 0;

    delete m_impl;
    m_impl = 0;
}

int FilterTool::Help(void) {
    Options::DisplayHelp();
    return 0;
}

int FilterTool::Run(int argc, char* argv[]) {

    // parse command line arguments
    Options::Parse(argc, argv, 1);

    // initialize FilterTool with settings
    m_impl = new FilterToolPrivate(m_settings);

    // run FilterTool, return success/fail
    if ( m_impl->Run() )
        return 0;
    else
        return 1;
}
