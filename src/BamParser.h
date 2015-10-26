/*
 * parser.h
 *
 *  Created on: Apr 17, 2012
 *      Author: fritz
 */

#ifndef BAMPARSER_H_
#define BAMPARSER_H_

#include "api/BamMultiReader.h"
#include "api/BamWriter.h"
#include "Alignment.h"
#include "Parser.h"

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>


using namespace BamTools;
using namespace std;

class BamParser: public Parser {
private:
	BamMultiReader reader;

public:
	BamParser(string file);
	~BamParser(){
		reader.Close();
	}
	Alignment * parseRead(uint16_t mappingQv);
	void parseReadFast(uint16_t mappingQv,Alignment*& aln);
	string get_header();
	RefVector get_refInfo();
};

#endif /* PARSER_H_ */
