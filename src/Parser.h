/*
 * Parser.h
 *
 *  Created on: May 29, 2012
 *      Author: fritz
 */

#ifndef PARSER_H_
#define PARSER_H_

#include "Alignment.h"

class Parser {

public:
	virtual ~Parser(){};
	virtual Alignment * parseRead(uint16_t mappingQv) = 0;
	virtual void parseReadFast(uint16_t mappingQv,Alignment *& align) = 0;
};

#endif /* PARSER_H_ */
