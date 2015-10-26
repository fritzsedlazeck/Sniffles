// ***************************************************************************
// bamtools_options.cpp (c) 2010 Derek Barnett, Erik Garrison
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

#include "utils/bamtools_options.h"
using namespace BamTools;

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <sstream>
using namespace std;

string Options::m_programName;                   // the program name
string Options::m_description;                   // the main description
string Options::m_exampleArguments;              // the example arguments
vector<OptionGroup> Options::m_optionGroups;     // stores the option groups
map<string, OptionValue> Options::m_optionsMap;  // stores the options in a map
const string Options::m_stdin  = "stdin";        // string representation of stdin
const string Options::m_stdout = "stdout";       // string representation of stdout

// adds a simple option to the parser
void Options::AddOption(const string& argument,
                        const string& optionDescription,
                        bool& foundArgument,
                        OptionGroup* group)
{
    Option o;
    o.Argument    = argument;
    o.Description = optionDescription;
    o.StoreValue  = false;
    group->Options.push_back(o);

    OptionValue ov;
    ov.pFoundArgument = &foundArgument;
    ov.StoreValue     = false;

    m_optionsMap[argument] = ov;
}

// creates an option group
OptionGroup* Options::CreateOptionGroup(const string& groupName) {
    OptionGroup og;
    og.Name = groupName;
    m_optionGroups.push_back(og);
    return &m_optionGroups[m_optionGroups.size() - 1];
}

// displays the help menu
void Options::DisplayHelp(void) {

    // initialize
    char argumentBuffer[ARGUMENT_LENGTH + 1];
    ostringstream sb;

    char indentBuffer[MAX_LINE_LENGTH - DESC_LENGTH + 1];
    memset(indentBuffer, ' ', MAX_LINE_LENGTH - DESC_LENGTH);
    indentBuffer[MAX_LINE_LENGTH - DESC_LENGTH] = 0;

    // display the menu
    printf("Description: %s.\n\n", m_description.c_str());
    printf("Usage: ");
    printf("%s", m_programName.c_str());
    printf(" %s\n\n", m_exampleArguments.c_str());

    vector<Option>::const_iterator      optionIter;
    vector<OptionGroup>::const_iterator groupIter;
    for (groupIter = m_optionGroups.begin(); groupIter != m_optionGroups.end(); ++groupIter) {
        
        printf("%s:\n", groupIter->Name.c_str());

        for (optionIter = groupIter->Options.begin(); optionIter != groupIter->Options.end(); ++optionIter) {

            if (optionIter->StoreValue) 
                snprintf(argumentBuffer, ARGUMENT_LENGTH + 1, "  %s <%s>", optionIter->Argument.c_str(), optionIter->ValueDescription.c_str());
            else 
                snprintf(argumentBuffer, ARGUMENT_LENGTH + 1, "  %s", optionIter->Argument.c_str());
            printf("%-35s ", argumentBuffer);

            string description = optionIter->Description;

            // handle default values
            if (optionIter->HasDefaultValue) {
                
                sb.str("");
                sb << description << " [";

                if (optionIter->DefaultValue.is_type<unsigned int>()) {
                    sb << (unsigned int)optionIter->DefaultValue;
                } else if (optionIter->DefaultValue.is_type<unsigned char>()) {
                    sb << (unsigned short)(unsigned char)optionIter->DefaultValue;
                } else if (optionIter->DefaultValue.is_type<float>()) {
                    sb << std::fixed << std::setprecision(2) << (float)optionIter->DefaultValue;
                } else if (optionIter->DefaultValue.is_type<double>()) {
                    sb << std::fixed << std::setprecision(4) << (double)optionIter->DefaultValue;
                } else if (optionIter->DefaultValue.is_type<std::string>()) {
                    const std::string stringValue = optionIter->DefaultValue;
                    sb << stringValue;
                } else {
                    printf("ERROR: Found an unsupported data type for argument %s when casting the default value.\n",
                           optionIter->Argument.c_str());
                    exit(1);
                }

                sb << "]";
                description = sb.str(); 
            }

            if ( description.size() <= DESC_LENGTH_FIRST_ROW ) {
                printf("%s\n", description.c_str());
            } else {

                // handle the first row
                const char* pDescription = description.data();
                unsigned int cutIndex = DESC_LENGTH_FIRST_ROW;
                while(pDescription[cutIndex] != ' ') 
                    cutIndex--;
                printf("%s\n", description.substr(0, cutIndex).c_str());
                description = description.substr(cutIndex + 1);

                // handle subsequent rows
                while(description.size() > DESC_LENGTH) {
                    pDescription = description.data();
                    cutIndex = DESC_LENGTH;
                    while(pDescription[cutIndex] != ' ') 
                        cutIndex--;
                    printf("%s%s\n", indentBuffer, description.substr(0, cutIndex).c_str());
                    description = description.substr(cutIndex + 1);
                }

                // handle last row
                printf("%s%s\n", indentBuffer, description.c_str());
            }                       
        }

        printf("\n");
    }

    printf("Help:\n"); 
    printf("  --help, -h                        shows this help text\n");
    exit(1);
}

// parses the command line
void Options::Parse(int argc, char* argv[], int offset) {

    // initialize
    map<string, OptionValue>::const_iterator ovMapIter;
    map<string, OptionValue>::const_iterator checkMapIter;
    const int LAST_INDEX = argc - 1;
    ostringstream errorBuilder;
    bool foundError = false;
    char* end_ptr = NULL;
    const string ERROR_SPACER(7, ' ');

    // check if we should show the help menu
    bool showHelpMenu = false;
    if (argc > 1) {
        for (int i = 1; i < argc; i++) {
            const std::string argument = argv[i];
            if ( (argument == "-h") || (argument == "--help") || (argument == "help") ) 
                showHelpMenu = true;
        }
    } else showHelpMenu = true;

    if (showHelpMenu) 
        DisplayHelp();

    // check each argument
    for (int i = offset+1; i < argc; i++) {
      
        const string argument = argv[i];
        ovMapIter = m_optionsMap.find(argument);

        if (ovMapIter == m_optionsMap.end()) {
            errorBuilder << ERROR_SPACER << "An unrecognized argument was found: " << argument << std::endl;
            foundError = true;
        } else {

            *ovMapIter->second.pFoundArgument = true;

            // grab the value
            if (ovMapIter->second.StoreValue) {

                if (i < LAST_INDEX) {

                    // check if the next argument is really a command line option
                    const string val = argv[i + 1]; 
                    checkMapIter = m_optionsMap.find(val);

                    if (checkMapIter == m_optionsMap.end()) {
                        
                        ++i;
                        
                        if (ovMapIter->second.VariantValue.is_type<unsigned int>()) {
                            const unsigned int uint32 = (unsigned int)strtoul(val.c_str(), &end_ptr, 10);
                            unsigned int* varValue = (unsigned int*)ovMapIter->second.pValue;
                            *varValue = uint32;
                        } else if (ovMapIter->second.VariantValue.is_type<unsigned char>()) {
                            const unsigned char uint8 = (unsigned char)strtoul(val.c_str(), &end_ptr, 10);
                            unsigned char* varValue = (unsigned char*)ovMapIter->second.pValue;
                            *varValue = uint8;
                        } else if (ovMapIter->second.VariantValue.is_type<uint64_t>()) {
                            const uint64_t uint64 = strtoui64(val.c_str(), &end_ptr, 10);
                            uint64_t* varValue = (uint64_t*)ovMapIter->second.pValue;
                            *varValue = uint64;
                        } else if (ovMapIter->second.VariantValue.is_type<double>()) {
                            const double d = strtod(val.c_str(), &end_ptr);
                            double* varValue = (double*)ovMapIter->second.pValue;
                            *varValue = d;
                        } else if (ovMapIter->second.VariantValue.is_type<float>()) {
                            const float f = (float)strtod(val.c_str(), &end_ptr);
                            float* varValue = (float*)ovMapIter->second.pValue;
                            *varValue = f;
                        } else if (ovMapIter->second.VariantValue.is_type<string>()) {
                            string* pStringValue = (string*)ovMapIter->second.pValue;
                            *pStringValue = val;
                        } else if (ovMapIter->second.VariantValue.is_type<vector<string> >()) {
                            vector<string>* pVectorValue = (vector<string>*)ovMapIter->second.pValue;
                            pVectorValue->push_back(val);
                        } else {
                            printf("ERROR: Found an unsupported data type for argument %s when parsing the arguments.\n",
                                   argument.c_str());
                            exit(1);
                        }
                    } else {
                        errorBuilder << ERROR_SPACER << "The argument (" << argument
                                     << ") expects a value, but none was found." << endl;
                        foundError = true;
                    }
                } else {
                    errorBuilder << ERROR_SPACER << "The argument (" << argument
                                 << ") expects a value, but none was found." << endl;
                    foundError = true;
                }
            }
        }
    }

    // check if we missed any required parameters
    for (ovMapIter = m_optionsMap.begin(); ovMapIter != m_optionsMap.end(); ++ovMapIter) {
        if (ovMapIter->second.IsRequired && !*ovMapIter->second.pFoundArgument) {
            errorBuilder << ERROR_SPACER << ovMapIter->second.ValueTypeDescription
                         << " was not specified. Please use the " << ovMapIter->first << " parameter." << endl;
            foundError = true;
        }
    }

    // print the errors if any were found
    if (foundError) {
        printf("ERROR: Some problems were encountered when parsing the command line options:\n");
        printf("%s\n", errorBuilder.str().c_str());
        printf("For a complete list of command line options, type \"%s help %s\"\n", argv[0], argv[1]);
        exit(1);
    }
}

// sets the program info
void Options::SetProgramInfo(const string& programName,
                             const string& description,
                             const string& arguments)
{
    m_programName      = programName;
    m_description      = description;
    m_exampleArguments = arguments;
}

// return string representations of stdin
const string& Options::StandardIn(void) { return m_stdin; }

// return string representations of stdout
const string& Options::StandardOut(void) { return m_stdout; }
