// ***************************************************************************
// bamtools_filter_ruleparser.h (c) 2010 Derek Barnett, Erik Garrison
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 10 October 2011
// ---------------------------------------------------------------------------
// Provides a compound rule parser for FilterEngine.
// ***************************************************************************

#ifndef BAMTOOLS_FILTER_RULEPARSER_H
#define BAMTOOLS_FILTER_RULEPARSER_H

#include "utils/bamtools_utilities.h"
#include <queue>
#include <stack>
#include <string>

namespace BamTools {

// -------------------------------------------
// char constants  
  
const char LEFT_PARENTHESIS_CHAR  = '(';
const char RIGHT_PARENTHESIS_CHAR = ')';
const char AND_OPERATOR_CHAR      = '&';
const char OR_OPERATOR_CHAR       = '|';
const char NOT_OPERATOR_CHAR      = '!';
const char SPACE_CHAR             = ' ';
  
// -------------------------------------------
// RuleToken implementation
  
struct RuleToken {
  
    // enums
    enum RuleTokenType { OPERAND = 0
                       , AND_OPERATOR
                       , OR_OPERATOR
                       , NOT_OPERATOR
                       , LEFT_PARENTHESIS
                       , RIGHT_PARENTHESIS
                       };
    
    // data members
    RuleTokenType Type;
    std::string Value;
};

inline int priority(const RuleToken& token) {
    switch ( token.Type ) {
        case ( RuleToken::NOT_OPERATOR )      : return 3;
        case ( RuleToken::AND_OPERATOR )      : return 2;
        case ( RuleToken::OR_OPERATOR  )      : return 1;
        case ( RuleToken::LEFT_PARENTHESIS )  : return 0;
        case ( RuleToken::RIGHT_PARENTHESIS ) : return 0;
        default:
            BAMTOOLS_ASSERT_UNREACHABLE;
            return -1;
    } 
}

inline bool isRightAssociative(const RuleToken& token) {
    return (token.Type == RuleToken::NOT_OPERATOR || 
            token.Type == RuleToken::LEFT_PARENTHESIS);
}

inline bool isLeftAssociative(const RuleToken& token) {
    return !isRightAssociative(token);
}

inline bool isLeftParenthesis(const RuleToken& token) {
    return ( token.Type == RuleToken::LEFT_PARENTHESIS );
}

inline bool isRightParenthesis(const RuleToken& token) {
    return ( token.Type == RuleToken::RIGHT_PARENTHESIS );
}

inline bool isOperand(const RuleToken& token) {
    return ( token.Type == RuleToken::OPERAND );
}

inline bool isOperator(const RuleToken& token) {
    return ( token.Type == RuleToken::AND_OPERATOR ||
             token.Type == RuleToken::OR_OPERATOR  ||
             token.Type == RuleToken::NOT_OPERATOR);
}
  
// -------------------------------------------
// RuleParser implementation  
  
class RuleParser {

    // ctor & dtor
    public:
        RuleParser(const std::string& ruleString)
            : m_ruleString(ruleString)
        { 
            // initialize char markers
            m_begin = (char*)m_ruleString.c_str();
            m_end   = m_begin + m_ruleString.length();
            ignoreQuotes();
        }
        
        ~RuleParser(void) { }
  
    // public interface
    public:
        void parse(void);
        std::queue<std::string> results(void) const { return m_ruleQueue; }

    // internal methods
    private:
        char getNextChar(void);
        void ignoreQuotes(void);
        bool readToken(RuleToken& token);
        void skipSpaces(void);
      
    // data members
    private:
        std::string m_ruleString;
        char* m_begin;
        char* m_current;
        char* m_end;
        
        std::queue<std::string> m_ruleQueue;
        std::stack<RuleToken> m_operatorStack;
};

inline
char RuleParser::getNextChar(void) {
   if ( m_current == m_end ) return 0;
   return *m_current++;
}

inline
void RuleParser::ignoreQuotes(void) {
    if ( *m_begin == '\"' ) ++m_begin;
    if ( *m_end   == '\"' ) --m_end;
}

inline
void RuleParser::parse(void) {
  
    // clear out any prior data
    while ( !m_ruleQueue.empty() ) 
        m_ruleQueue.pop();
    
    // skip if no rule to parse
    if ( m_ruleString.empty() ) return;
  
    // start at beginning of ruleString
    m_current = m_begin;
    
    // iterate through tokens in rule string
    RuleToken token;
    while ( readToken(token) ) {
      
        if ( token.Value.empty() ) break;
      
        // if token is an operand
        if ( isOperand(token) )
            m_ruleQueue.push(token.Value);

        // if token is an operator 
        else if ( isOperator(token) ) {

            // pop any operators at top of stack with higher priority
            while ( !m_operatorStack.empty() ) {
                const RuleToken& opToken = m_operatorStack.top();
                if ( (isLeftAssociative(token) && (priority(token) <= priority(opToken))) ||
                     (isRightAssociative(token) && (priority(token) < priority(opToken))) 
                    )
                {
                    m_ruleQueue.push(opToken.Value);
                    m_operatorStack.pop();
                }
                else break;
            }
            
            // push current operator token onto stack
            m_operatorStack.push(token);
        }
        
        // if token is left parenthesis
        else if ( isLeftParenthesis(token) )
            m_operatorStack.push(token);
        
        // if token is right parenthesis
        else if ( isRightParenthesis(token) ) {
          
            bool foundLeftParenthesis = false;
          
            // push operators into rule queue until left parenthesis found
            while ( !m_operatorStack.empty() && !foundLeftParenthesis ) {
                const RuleToken& opToken = m_operatorStack.top();
                if ( !isLeftParenthesis(opToken) )
                    m_ruleQueue.push(opToken.Value);
                else 
                    foundLeftParenthesis = true;
                m_operatorStack.pop();
            }
          
            // no left parenthesis found, error
            BAMTOOLS_ASSERT_MESSAGE( foundLeftParenthesis, "ERROR: Mismatched parenthesis in rule string.1");
        }
        
        // error: unknown operand
        else BAMTOOLS_ASSERT_UNREACHABLE;
    }    
    
    // while there are still operators on stack
    while ( !m_operatorStack.empty() ) {
        const RuleToken& token = m_operatorStack.top();
        BAMTOOLS_ASSERT_MESSAGE( (!isLeftParenthesis(token) && !isRightParenthesis(token)), "ERROR: Mismatched parenthesis in rule string.2");
        m_ruleQueue.push(token.Value);
        m_operatorStack.pop();
    }
}

inline
bool RuleParser::readToken(RuleToken& token) {
  
    // skip any preceding whitespace
    skipSpaces();
    if ( m_current == m_end ) return false;

    // clear out prior token value
    token.Value.clear();
    
    // read chars while still in token
    char c = 1;
    bool keepReading = true;
    bool inOperandString = false;
    while ( keepReading && (c != 0) ) {
      
      // get next char
      c = getNextChar();
      switch (c) {
        
          // current char is '('
          case ( LEFT_PARENTHESIS_CHAR ) :
              token.Type = RuleToken::LEFT_PARENTHESIS;
              token.Value.append(1, LEFT_PARENTHESIS_CHAR);
              keepReading = false;
              break;
              
          // current char is ')'
          case ( RIGHT_PARENTHESIS_CHAR ) :
              if ( inOperandString )
                  --m_current;
              else {
                  token.Type = RuleToken::RIGHT_PARENTHESIS;
                  token.Value.append(1, RIGHT_PARENTHESIS_CHAR);
              }
              keepReading = false;
              break;
        
          // current char is '&'
          case ( AND_OPERATOR_CHAR ) :
              if ( inOperandString ) 
                  --m_current;
              else {
                  token.Type = RuleToken::AND_OPERATOR;
                  token.Value.append(1, AND_OPERATOR_CHAR);
              }
              keepReading = false;
              break;
              
          // current char is '|' 
          case ( OR_OPERATOR_CHAR ) :
              if ( inOperandString )
                  --m_current;
              else {  
                  token.Type = RuleToken::OR_OPERATOR;
                  token.Value.append(1, OR_OPERATOR_CHAR);
              }
              keepReading = false;
              break;
              
          // current char is '!'
          case ( NOT_OPERATOR_CHAR ) :
              token.Type = RuleToken::NOT_OPERATOR;
              token.Value.append(1, NOT_OPERATOR_CHAR);
              keepReading = false;
              break;
              
          // current char is ' '
          case ( SPACE_CHAR ) : 
              keepReading = false;
              break;
            
          // current char is a true value token
          default:
              if ( c != 0 ) {
                  token.Type = RuleToken::OPERAND;
                  token.Value.append(1, c);
                  inOperandString = true;
                  keepReading = true;
              }
        }         
    }
      
    return true;
}

inline
void RuleParser::skipSpaces(void) {
    while ( m_current != m_end ) {
        const char c = *m_current;
        if ( c == ' ' || c == '\t' || c == '\r' || c == '\n')
            ++m_current;
        else break;
    }
}

} // namespace BamTools

#endif // BAMTOOLS_FILTER_RULEPARSER_H
