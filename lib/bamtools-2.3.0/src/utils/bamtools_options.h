// ***************************************************************************
// bamtools_options.h (c) 2010 Derek Barnett, Erik Garrison
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 10 October 2011
// ---------------------------------------------------------------------------
// Parses command line arguments and creates a help menu
// ---------------------------------------------------------------------------
// Modified from:
// The Mosaik suite's command line parser class: COptions
// (c) 2006 - 2009 Michael Str�mberg
// Marth Lab, Department of Biology, Boston College
// Re-licensed under MIT License with author's permission.
//
// * Modified slightly to fit BamTools, otherwise code is same. 
// *  (BamTools namespace, added stdin/stdout) (DB)
// ***************************************************************************

#ifndef BAMTOOLS_OPTIONS_H
#define BAMTOOLS_OPTIONS_H

#include "utils/bamtools_variant.h"
#include "utils/utils_global.h"

#include <map>
#include <string>
#include <vector>

#ifndef WIN32
    #include <stdint.h>
#endif

namespace BamTools {

#define ARGUMENT_LENGTH       35
#define DESC_LENGTH_FIRST_ROW 30
#define DESC_LENGTH           42
#define MAX_LINE_LENGTH       78

#ifdef WIN32
  #define snprintf _snprintf
  typedef __int64          int64_t;
  typedef unsigned __int64 uint64_t;
  #define strtoui64 _strtoui64
#else
  #define strtoui64 strtoull
#endif

struct UTILS_EXPORT Option {
  
    // data members
    std::string Argument;
    std::string ValueDescription;
    std::string Description;
    bool StoreValue;
    bool HasDefaultValue;
    Variant DefaultValue;

    // constructor
    Option(void)
        : StoreValue(true)
        , HasDefaultValue(false)
    { }
};

struct UTILS_EXPORT OptionValue {
  
    // data members
    bool* pFoundArgument;
    void* pValue;
    std::string ValueTypeDescription;
    bool UseVector;
    bool StoreValue;
    bool IsRequired;
    Variant VariantValue;

    // constructor
    OptionValue(void)
        : pFoundArgument(NULL)
        , pValue(NULL)
        , UseVector(false)
        , StoreValue(true)
        , IsRequired(false)
    { } 
};

struct UTILS_EXPORT OptionGroup {
    std::string Name;
    std::vector<Option> Options;
};

class UTILS_EXPORT Options {
  
    // add option/argument rules
    public:
        // adds a simple option to the parser
        static void AddOption(const std::string& argument, 
                       const std::string& optionDescription, 
                       bool& foundArgument, 
                       OptionGroup* group);
                       
        // adds a value option to the parser
        template<typename T>
        static void AddValueOption(const std::string& argument, 
                            const std::string& valueDescription, 
                            const std::string& optionDescription, 
                            const std::string& valueTypeDescription, 
                            bool& foundArgument, 
                            T& val, 
                            OptionGroup* group);
                            
        // adds a value option to the parser (with a default value)
        template<typename T, typename D>
        static void AddValueOption(const std::string& argument, 
                            const std::string& valueDescription, 
                            const std::string& optionDescription, 
                            const std::string& valueTypeDescription, 
                            bool& foundArgument, 
                            T& val, 
                            OptionGroup* group, 
                            D& defaultValue);
       
    // other API methods
    public:
        // creates an option group
        static OptionGroup* CreateOptionGroup(const std::string& groupName);    
        // displays the help menu
        static void DisplayHelp(void);
        // parses the command line
        static void Parse(int argc, char* argv[], int offset = 0);
        // sets the program info
        static void SetProgramInfo(const std::string& programName,
                                   const std::string& description,
                                   const std::string& arguments);
        // returns string representation of stdin
        static const std::string& StandardIn(void);
        // returns string representation of stdout
        static const std::string& StandardOut(void);
        
    // static data members
    private:
        // the program name
        static std::string m_programName;
        // the main description
        static std::string m_description;
        // the example arguments
        static std::string m_exampleArguments;
        // stores the option groups
        static std::vector<OptionGroup> m_optionGroups;
        // stores the options in a map
        static std::map<std::string, OptionValue> m_optionsMap;
        // string representation of stdin
        static const std::string m_stdin;
        // string representation of stdout
        static const std::string m_stdout;
};

// adds a value option to the parser
template<typename T>
void Options::AddValueOption(const std::string& argument, 
                             const std::string& valueDescription, 
                             const std::string& optionDescription, 
                             const std::string& valueTypeDescription, 
                             bool& foundArgument, 
                             T& val, 
                             OptionGroup* group) 
{
        Option o;
        o.Argument         = argument;
        o.ValueDescription = valueDescription;
        o.Description      = optionDescription;
        group->Options.push_back(o);

        OptionValue ov;
        ov.pFoundArgument       = &foundArgument;
        ov.pValue               = (void*)&val;
        ov.VariantValue         = val;
        ov.IsRequired           = (valueTypeDescription.empty() ? false : true);
        ov.ValueTypeDescription = valueTypeDescription;
        m_optionsMap[argument] = ov;
}

// adds a value option to the parser (with a default value)
template<typename T, typename D>
void Options::AddValueOption(const std::string& argument, 
                             const std::string& valueDescription, 
                             const std::string& optionDescription, 
                             const std::string& valueTypeDescription, 
                             bool& foundArgument, 
                             T& val, 
                             OptionGroup* group, 
                             D& defaultValue) 
{
        Option o;
        o.Argument         = argument;
        o.ValueDescription = valueDescription;
        o.Description      = optionDescription;
        o.DefaultValue     = defaultValue;
        o.HasDefaultValue  = true;
        group->Options.push_back(o);

        OptionValue ov;
        ov.pFoundArgument       = &foundArgument;
        ov.pValue               = (void*)&val;
        ov.VariantValue         = val;
        ov.IsRequired           = (valueTypeDescription.empty() ? false : true);
        ov.ValueTypeDescription = valueTypeDescription;
        m_optionsMap[argument] = ov;
}

} // namespace BamTools

#endif // BAMTOOLS_OPTIONS_H
