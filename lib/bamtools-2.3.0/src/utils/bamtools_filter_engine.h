// ***************************************************************************
// bamtools_filter_engine.h (c) 2010 Derek Barnett, Erik Garrison
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 3 May 2013
// ---------------------------------------------------------------------------
// Provides a generic filter engine based on filter-sets of properties,
// with possible "rules" (compound logical expressions) to create more complex
// queries on a data set.
//
// FilterEngine consists, most importantly, of :
//
//     a list of possible properties (each tagged whether it has been 'enabled' as a filter)
//     a map of filterName => propertySet
//     queue for compound rule expression (i.e. "(filter1 AND filter2) OR !filter3" )
//     
// Each propertySet is a list of properties enabled for this particular filter object
//
//     Implemented as a map of propertyNames to propertyFilterValue
//     ( "property1" => pfv1 
//       "property2" => pfv2 
//       "property4" => pfv4
//       etc. )  
//
//     Any properties that are 'possible', via FilterEngine::addProperty(), but not enabled 
//     via FilterEngine::setProperty() (in our example, say "property3"), evaluate to true 
//     for any query.  Meaning that if a property is not set on this filter, we don't care 
//     about it here, so it passes though OK.
//
// A propertyFilterValue contains a value and comparison type
//
//    ( pfv1: Value = 50,    Type = GREATER_THAN_EQUAL
//      pfv2: Value = "foo", Type = STARTS_WITH
//      pfv4: Value = "bar", Type = CONTAINS
//      etc. )  
//
//    This allows for more complex queries (than simple isEqual?) against a variety of data types.
// 
// ***************************************************************************

#ifndef BAMTOOLS_FILTER_ENGINE_H
#define BAMTOOLS_FILTER_ENGINE_H

#include "utils/utils_global.h"
#include "utils/bamtools_filter_properties.h"
#include "utils/bamtools_filter_ruleparser.h"
#include "utils/bamtools_utilities.h"

#include <algorithm>
#include <iostream>
#include <map>
#include <queue>
#include <sstream>
#include <stack>
#include <string>
#include <utility>
#include <vector>

namespace BamTools {

struct UTILS_EXPORT FilterCompareType {
    enum Type { AND = 0
              , NOT
              , OR
    };
};
  
// -----------------------------------------------------------
// FilterEngine
  
template <typename FilterChecker>
class UTILS_EXPORT FilterEngine {
  
    // ctor & dtor
    public:
        FilterEngine(void) 
            : m_ruleString("")
            , m_isRuleQueueGenerated(false)
            , m_defaultCompareType(FilterCompareType::OR)
            , AND_OPERATOR("&")
            , OR_OPERATOR("|")
            , NOT_OPERATOR("!")
        { }
        
        ~FilterEngine(void) { }
  
    // 'filter set' methods
    public:
        // creates a new filter set, returns true if created, false if error or already exists
        bool addFilter(const std::string& filterName);       
        
        // return list of current filter names
        const std::vector<std::string> filterNames(void);    
  
    // 'property' methods
    public:
      
        // add a new known property (& type) to engine
        bool addProperty(const std::string& propertyName);
  
        // sets property filter (value, type) for propertyName, on a particular filter set 
        // setProperty("filter1", "mapQuality", 50, GREATER_THAN_EQUAL)
        template<typename T>
        bool setProperty(const std::string& filterName, 
                         const std::string& propertyName, 
                         const T& value,
                         const PropertyFilterValue::ValueCompareType& type = PropertyFilterValue::EXACT);
        
        // returns list of all properties known by FilterEngine  ( any created using addProperty() )
        const std::vector<std::string> allPropertyNames(void);
        
        // returns list of property names that are 'enabled' ( only those touched by setProperty() )
        const std::vector<std::string> enabledPropertyNames(void);  
  
     // 'rule' methods
    public:   
        
        // sets comparison operator between filters if no rule string given
        // default is to do an OR on each filter
        void setDefaultCompareType(const FilterCompareType::Type& type = FilterCompareType::OR);

        // sets rule string for building expression queue
        // if empty, creates 
        void setRule(const std::string& ruleString = "");
        
    // token parsing (for property filter generation)
    public:
        template<typename T>
        static bool parseToken(const std::string& token, T& value, PropertyFilterValue::ValueCompareType& type);
        
    // query evaluation
    public:
        // returns true if query passes all filters in FilterEngine
        template<typename T>
        bool check(const T& query);

    // internal rule-handling methods
    private:
        void buildDefaultRuleString(void);
        void buildRuleQueue(void);
        template<typename T>
        bool evaluateFilterRules(const T& query);
        
    // data members
    private:
        // all 'filter sets'
        FilterMap m_filters;
        
        // all known properties
        std::vector<Property> m_properties; 
        
        // infix expression of filter-set comparison rules 
        std::string m_ruleString;
        
        // postfix expression of tokens (filterNames) and operators (as strings)
        // if this is empty, uses m_compareType to build default expression queue
        std::queue<std::string> m_ruleQueue;
        
        // flag to test if the rule expression queue has been generated
        bool m_isRuleQueueGenerated;
        
        // 'default' comparison operator between filters if no rule string given
        // if this is changed, m_ruleString is used to build new m_ruleQueue
        FilterCompareType::Type m_defaultCompareType;
        
        // client-specified checking type ( provides method: bool check(PropertyFilter, T object) )
        FilterChecker m_checker;
        
        // token-parsing constants
        static const int NOT_CHAR          = (int)'!';
        static const int EQUAL_CHAR        = (int)'=';
        static const int GREATER_THAN_CHAR = (int)'>';
        static const int LESS_THAN_CHAR    = (int)'<';
        static const int WILDCARD_CHAR     = (int)'*';
        
        // filter evaluation constants
        const std::string AND_OPERATOR;
        const std::string OR_OPERATOR;
        const std::string NOT_OPERATOR;
};

// creates a new filter set, returns true if created, false if error or already exists
template<typename FilterChecker>
inline bool FilterEngine<FilterChecker>::addFilter(const std::string& filterName) {
    return (m_filters.insert(std::make_pair(filterName, PropertyFilter()))).second;
}

// add a new known property & type to engine
template<typename FilterChecker>
inline bool FilterEngine<FilterChecker>::addProperty(const std::string& propertyName) {
    const std::vector<std::string> propertyNames = allPropertyNames();
    bool found = std::binary_search( propertyNames.begin(), propertyNames.end(), propertyName );
    if ( found ) return false;
    m_properties.push_back( Property(propertyName) );
    std::sort( m_properties.begin(), m_properties.end() );
    return true;
}

// returns list of all properties known by FilterEngine 
// ( any that were created using addProperty() )
template<typename FilterChecker>
inline const std::vector<std::string> FilterEngine<FilterChecker>::allPropertyNames(void) {
    // set up stringlist
    std::vector<std::string> names;
    names.reserve(m_properties.size());
    // iterate through all properties, appending to stringlist
    std::vector<Property>::const_iterator propIter = m_properties.begin();
    std::vector<Property>::const_iterator propEnd  = m_properties.end();
    for ( ; propIter != propEnd; ++propIter )
        names.push_back( (*propIter).Name );  
    // return stringlist
    return names;
}

// builds a default rule string based on m_defaultCompareType
// used if user supplied an explicit rule string
template<typename FilterChecker>
inline void FilterEngine<FilterChecker>::buildDefaultRuleString(void) {
  
    // set up temp string stream 
    std::stringstream ruleStream("");
  
    // get first filterName
    FilterMap::const_iterator mapIter = m_filters.begin();
    ruleStream << (*mapIter).first;
    
    // if there are more filters present
    // iterate over remaining filters, appending compare operator and filter name
    if ( m_filters.size() > 1 ) {        
        for ( ++mapIter ; mapIter != m_filters.end(); ++mapIter )
            ruleStream << ( (m_defaultCompareType == FilterCompareType::AND) ? " & " : " | " ) 
                       << (*mapIter).first;
    }

    // set m_ruleString from temp stream
    m_ruleString = ruleStream.str();
}

// build expression queue based on ruleString
template<typename FilterChecker>
inline void FilterEngine<FilterChecker>::buildRuleQueue(void) {
  
    // skip if no filters present
    if ( m_filters.empty() ) return;
  
    // clear out any prior expression queue data
    while ( !m_ruleQueue.empty() )
        m_ruleQueue.pop();
  
    // create a rule string, if not provided
    if ( m_ruleString.empty() ) 
        buildDefaultRuleString();
    
    // initialize RuleParser, run, and retrieve results
    RuleParser ruleParser(m_ruleString);
    ruleParser.parse();
    m_ruleQueue = ruleParser.results();
    
    // set flag if rule queue contains any values
    m_isRuleQueueGenerated = (!m_ruleQueue.empty());    
}

// returns whether query value passes filter engine rules
template<class FilterChecker> template<typename T>
bool FilterEngine<FilterChecker>::check(const T& query) {
  
    // return result of querying against filter rules
    return evaluateFilterRules(query);
}

// returns list of property names that are 'enabled' ( only those touched by setProperty() )
template<typename FilterChecker>
inline const std::vector<std::string> FilterEngine<FilterChecker>::enabledPropertyNames(void) {
    // initialize stringlist
    std::vector<std::string> names;
    names.reserve(m_properties.size());
    // iterate over all properties, appending if enabled
    std::vector<Property>::const_iterator propIter = m_properties.begin();
    std::vector<Property>::const_iterator propEnd  = m_properties.end();
    for ( ; propIter != propEnd; ++propIter )
        if ( (*propIter).IsEnabled ) 
            names.push_back( (*propIter).Name );    
    // return stringlist
    return names;
}

// evaluates postfix rule queue - with each filter as an operand, AND|OR|NOT as operators
template<class FilterChecker> template<typename T>
bool FilterEngine<FilterChecker>::evaluateFilterRules(const T& query) {
  
    // build ruleQueue if not done before
    if ( !m_isRuleQueueGenerated ) 
        buildRuleQueue();
    
    std::stack<bool> resultStack;
    FilterMap::const_iterator filterIter;
    std::queue<std::string> ruleQueueCopy = m_ruleQueue;
    while ( !ruleQueueCopy.empty() ) {
        const std::string& token = ruleQueueCopy.front();
        
        // token is NOT_OPERATOR
        if ( token == FilterEngine<FilterChecker>::NOT_OPERATOR ) {
            BAMTOOLS_ASSERT_MESSAGE( !resultStack.empty(), "Empty result stack - cannot apply operator: !" );
            resultStack.top() = !resultStack.top();
        }
        
        // token is AND_OPERATOR
        else if ( token == FilterEngine<FilterChecker>::AND_OPERATOR ) {
            BAMTOOLS_ASSERT_MESSAGE( resultStack.size() >= 2 , "Not enough operands - cannot apply operator: &" );
            bool topResult = resultStack.top();
            resultStack.pop();
            resultStack.top() &= topResult;
        }
        
        // token is OR_OPERATOR
        else if ( token == FilterEngine<FilterChecker>::OR_OPERATOR ) {
            BAMTOOLS_ASSERT_MESSAGE( resultStack.size() >= 2 , "Not enough operands - cannot apply operator: |" );
            bool topResult = resultStack.top();
            resultStack.pop();
            resultStack.top() |= topResult;
        }
        
        // token is an operand 
        else {
            // look up PropertyFilter that matches this token 
            filterIter = m_filters.find(token);
            BAMTOOLS_ASSERT_MESSAGE( (filterIter != m_filters.end() ), "Filter mentioned in rule, not found in FilterEngine" );
            const PropertyFilter& filter = (*filterIter).second;
            bool result = m_checker.check(filter, query);
            resultStack.push( result );
        }
        
        // pop token from ruleQueue
        ruleQueueCopy.pop();
    }
    
    // return last result
    BAMTOOLS_ASSERT_MESSAGE( resultStack.size() == 1, "Result stack should only have one value remaining - cannot return result" );
    return resultStack.top();
}

// return list of current filter names
template<typename FilterChecker>
inline const std::vector<std::string> FilterEngine<FilterChecker>::filterNames(void) {
    // initialize stringlist
    std::vector<std::string> names;
    names.reserve(m_filters.size());
    // iterate over all filters, appending filter name
    FilterMap::const_iterator mapIter = m_filters.begin();
    FilterMap::const_iterator mapEnd  = m_filters.end();
    for ( ; mapIter != mapEnd; ++mapIter )
        names.push_back( (*mapIter).first ); 
    // return stringlist
    return names;
}

// parse a filterValue token string that may contain comparison qualifiers (">50", "*SRR", etc.)
template<class FilterChecker> template<typename T>
bool FilterEngine<FilterChecker>::parseToken(const std::string& token, T& value, PropertyFilterValue::ValueCompareType& type) {
    
    // skip if token is empty
    if ( token.empty() ) return false;
    
    // will store token after special chars are removed
    std::string strippedToken;
    
    // if only single character
    if ( token.length() == 1 ) {
        strippedToken = token;
        type = PropertyFilterValue::EXACT;
    } 
    
    // more than one character, check for special chars
    else {
        const int firstChar = (int)token.at(0);
        switch ( firstChar ) {
          
            case ( FilterEngine<FilterChecker>::NOT_CHAR ) :
                strippedToken = token.substr(1);       
                type = PropertyFilterValue::NOT;
                break;
                
            case ( FilterEngine<FilterChecker>::GREATER_THAN_CHAR ) :
                
                // check for '>=' case
                if ( token.at(1) == FilterEngine<FilterChecker>::EQUAL_CHAR ) {
                    if ( token.length() == 2 ) return false;
                    strippedToken = token.substr(2);
                    type = PropertyFilterValue::GREATER_THAN_EQUAL;
                } 
                
                // otherwise only '>'
                else {
                    strippedToken = token.substr(1);
                    type = PropertyFilterValue::GREATER_THAN;
                }
                
                break;
                
            case ( FilterEngine<FilterChecker>::LESS_THAN_CHAR ) : 
         
                // check for '<=' case
                if ( token.at(1) == FilterEngine<FilterChecker>::EQUAL_CHAR ) {
                    if ( token.length() == 2 ) return false;
                    strippedToken = token.substr(2);
                    type = PropertyFilterValue::LESS_THAN_EQUAL;
                } 
                
                // otherwise only '<'
                else {
                    strippedToken = token.substr(1);
                    type = PropertyFilterValue::LESS_THAN;
                }
                
                break;
                
            case ( FilterEngine<FilterChecker>::WILDCARD_CHAR ) : 
              
                // check for *str* case (CONTAINS)
                if ( token.at( token.length() - 1 ) == FilterEngine<FilterChecker>::WILDCARD_CHAR ) {
                    if ( token.length() == 2 ) return false;
                    strippedToken = token.substr(1, token.length() - 2);
                    type = PropertyFilterValue::CONTAINS;
                }
                
                // otherwise *str case (ENDS_WITH)
                else {
                    strippedToken = token.substr(1);
                    type = PropertyFilterValue::ENDS_WITH;
                }
                
                break;
               
            default :
                // check for str* case (STARTS_WITH)
                if ( token.at( token.length() - 1 ) == FilterEngine<FilterChecker>::WILDCARD_CHAR ) {
                    if ( token.length() == 2 ) return false;
                    strippedToken = token.substr(0, token.length() - 1);
                    type = PropertyFilterValue::STARTS_WITH;
                }
                
                // otherwise EXACT
                else {
                    strippedToken = token;
                    type = PropertyFilterValue::EXACT;
                }
                
                break;
        }
    }
    
    // convert stripped token to value
    std::stringstream stream(strippedToken);
    if ( strippedToken == "true" || strippedToken == "false" )
        stream >> std::boolalpha >> value;
    else 
        stream >> value;
    
    // check for valid CompareType on type T
    Variant variantCheck = value;
    
    // if T is not string AND CompareType is for string values, return false
    if ( !variantCheck.is_type<std::string>() ) {
        if ( type == PropertyFilterValue::CONTAINS || 
             type == PropertyFilterValue::ENDS_WITH || 
             type == PropertyFilterValue::STARTS_WITH )          
            
          return false;
    }
    
    // return success
    return true;
}

// sets comparison operator between filters if no rule string given
// default is to do an OR on each filter
template<typename FilterChecker>
inline void FilterEngine<FilterChecker>::setDefaultCompareType(const FilterCompareType::Type& type) {
    // check for supported compare type
    if ( type == FilterCompareType::AND || type == FilterCompareType::OR ) {
        // if not the current compare type
        if ( m_defaultCompareType != type ) {
            m_defaultCompareType = type;
            buildRuleQueue();
        }
    }
}

// sets property filter (value, type) for propertyName, on a particular filter set 
// setProperty("filter1", "mapQuality", 50, GREATER_THAN_EQUAL)
template<class FilterChecker> template<typename T>
bool FilterEngine<FilterChecker>::setProperty(const std::string& filterName, 
                                              const std::string& propertyName, 
                                              const T& value,
                                              const PropertyFilterValue::ValueCompareType& type)
{
    // lookup filter by name, return false if not found
    FilterMap::iterator filterIter = m_filters.find(filterName);
    if ( filterIter == m_filters.end() ) return false;
      
    // lookup property for filter, add new PropertyFilterValue if not found, modify if already exists
    PropertyFilter& filter = (*filterIter).second;
    PropertyMap::iterator propertyIter = filter.Properties.find(propertyName);
    
    bool success;
    
    // property not found for this filter, create new entry
    if ( propertyIter == filter.Properties.end() )
        success = (filter.Properties.insert(std::make_pair(propertyName, PropertyFilterValue(value, type)))).second;
    
    // property already exists, modify
    else {
        PropertyFilterValue& filterValue = (*propertyIter).second;
        filterValue.Value = value;
        filterValue.Type  = type;
        success = true;
    }
    
    // if error so far, return false
    if ( !success ) return false;
    
    // --------------------------------------------
    // otherwise, set Property.IsEnabled to true
    
    // lookup property
    std::vector<Property>::iterator knownPropertyIter = std::find( m_properties.begin(), m_properties.end(), propertyName);
    
    // if not found, create a new (enabled) entry (& re-sort list)
    if ( knownPropertyIter == m_properties.end() ) {
        m_properties.push_back( Property(propertyName, true) );
        std::sort( m_properties.begin(), m_properties.end() );
    } 
    
    // property already known, set as enabled
    else (*knownPropertyIter).IsEnabled = true;

    // return success
    return true;
}

// sets user-specified rule string & signals update of rule-expression queue
template<typename FilterChecker>
inline void FilterEngine<FilterChecker>::setRule(const std::string& ruleString) {
    if ( m_ruleString != ruleString) {
        m_ruleString = ruleString;
        buildRuleQueue();
    }
}

} // namespace BamTools

#endif // BAMTOOLS_FILTER_ENGINE_H
