// ***************************************************************************
// bamtools_filter_properties.h (c) 2010 Derek Barnett, Erik Garrison
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 10 October 2011
// ---------------------------------------------------------------------------
// Provides support data structures & methods for FilterEngine
//
// The FilterEngine consists, most importantly, of :
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

#ifndef BAMTOOLS_FILTER_PROPERTIES_H
#define BAMTOOLS_FILTER_PROPERTIES_H

#include "utils/utils_global.h"
#include "utils/bamtools_utilities.h"
#include "utils/bamtools_variant.h"
#include <iostream>
#include <map>
#include <string>

namespace BamTools {

// ----------------------------------------------------------
// PropertyFilterValue
  
struct UTILS_EXPORT PropertyFilterValue {
  
    // define valid ValueCompareTypes
    enum ValueCompareType { CONTAINS = 0
                          , ENDS_WITH
                          , EXACT
                          , GREATER_THAN
                          , GREATER_THAN_EQUAL
                          , LESS_THAN
                          , LESS_THAN_EQUAL
                          , NOT
                          , STARTS_WITH
                          };
                   
    // ctor
    PropertyFilterValue(const Variant& value = Variant(),
                        const ValueCompareType& type = PropertyFilterValue::EXACT)
        : Value(value)
        , Type(type)
    { }
          
    // filter check methods      
    template<typename T>
    bool check(const T& query) const;
    bool check(const std::string& query) const;
             
    // data members
    Variant Value;
    ValueCompareType Type;
};

// checks a query against a filter (value, compare type)
template<typename T>
bool PropertyFilterValue::check(const T& query) const {
  
    // ensure filter value & query are same type
    if ( !Value.is_type<T>() ) { 
        std::cerr << "Cannot compare different types!" << std::endl;
        return false;
    }
    
    // string matching
    if ( Value.is_type<std::string>() ) {
        std::cerr << "Cannot compare different types - query is a string!" << std::endl;
        return false;
    } 
    
    // numeric matching based on our filter type
    switch ( Type ) {
        case ( PropertyFilterValue::EXACT)              : return ( query == Value.get<T>() );
        case ( PropertyFilterValue::GREATER_THAN)       : return ( query >  Value.get<T>() ); 
        case ( PropertyFilterValue::GREATER_THAN_EQUAL) : return ( query >= Value.get<T>() ); 
        case ( PropertyFilterValue::LESS_THAN)          : return ( query <  Value.get<T>() );
        case ( PropertyFilterValue::LESS_THAN_EQUAL)    : return ( query <= Value.get<T>() );
        case ( PropertyFilterValue::NOT)                : return ( query != Value.get<T>() );
        default : BAMTOOLS_ASSERT_UNREACHABLE;
    }
    return false;
}

// checks a string query against filter (value, compare type)
inline
bool PropertyFilterValue::check(const std::string& query) const {
  
    // ensure filter value & query are same type
    if ( !Value.is_type<std::string>() ) {
        std::cerr << "Cannot compare different types!" << std::endl;
        return false;
    }
  
    // localize string version of our filter value
    const std::string& valueString = Value.get<std::string>();
    
    // string matching based on our filter type
    switch ( Type ) {
        case ( PropertyFilterValue::CONTAINS)           : return ( query.find(valueString) != std::string::npos );
        case ( PropertyFilterValue::ENDS_WITH)          : return ( query.find(valueString) == (query.length() - valueString.length()) ); 
        case ( PropertyFilterValue::EXACT)              : return ( query == valueString );
        case ( PropertyFilterValue::GREATER_THAN)       : return ( query >  valueString ); 
        case ( PropertyFilterValue::GREATER_THAN_EQUAL) : return ( query >= valueString ); 
        case ( PropertyFilterValue::LESS_THAN)          : return ( query <  valueString );
        case ( PropertyFilterValue::LESS_THAN_EQUAL)    : return ( query <= valueString );
        case ( PropertyFilterValue::NOT)                : return ( query != valueString );
        case ( PropertyFilterValue::STARTS_WITH)        : return ( query.find(valueString) == 0 );
        default : BAMTOOLS_ASSERT_UNREACHABLE;
    }
    return false;
}

inline
const std::string toString(const PropertyFilterValue::ValueCompareType& type) {
  
    switch ( type ) {
        case ( PropertyFilterValue::CONTAINS )           : return std::string( "CONTAINS");
        case ( PropertyFilterValue::ENDS_WITH )          : return std::string( "ENDS_WITH");
        case ( PropertyFilterValue::EXACT )              : return std::string( "EXACT");
        case ( PropertyFilterValue::GREATER_THAN )       : return std::string( "GREATER_THAN");
        case ( PropertyFilterValue::GREATER_THAN_EQUAL ) : return std::string( "GREATER_THAN_EQUAL");
        case ( PropertyFilterValue::LESS_THAN )          : return std::string( "LESS_THAN");
        case ( PropertyFilterValue::LESS_THAN_EQUAL )    : return std::string( "LESS_THAN_EQUAL");
        case ( PropertyFilterValue::NOT )                : return std::string( "NOT");
        case ( PropertyFilterValue::STARTS_WITH )        : return std::string( "STARTS_WITH");
        default : BAMTOOLS_ASSERT_UNREACHABLE;
    }
    return std::string();
}

// property name => property filter value 
// ('name' => ('SSR', STARTS_WITH), 'mapQuality' => (50, GREATER_THAN_EQUAL), etc...)
typedef std::map<std::string, PropertyFilterValue> PropertyMap;

// ----------------------------------------------------------
// PropertyFilter

struct UTILS_EXPORT PropertyFilter {
    // data members
    PropertyMap Properties;
};

// filter name => properties  
// ('filter1' => properties1, 'filter2' => properties2, etc...)
typedef std::map<std::string, PropertyFilter> FilterMap;
  
// ----------------------------------------------------------
// Property
  
// used to store properties known to engine & keep track of enabled state
struct UTILS_EXPORT Property {
    std::string Name;
    bool IsEnabled;
    Property(const std::string& name, bool isEnabled = false) 
        : Name(name)
        , IsEnabled(isEnabled) 
    { }
};

inline bool operator<  (const Property& lhs, const Property& rhs) { return lhs.Name <  rhs.Name; }
inline bool operator== (const Property& lhs, const Property& rhs) { return lhs.Name == rhs.Name; }

} // namespace BamTools

#endif // BAMTOOLS_FILTER_PROPERTIES_H
