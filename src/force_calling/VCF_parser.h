/*
 * VCF_parser.h
 *
 *  Created on: Aug 24, 2017
 *      Author: sedlazec
 */

#ifndef FORCE_CALLING_VCF_PARSER_H_
#define FORCE_CALLING_VCF_PARSER_H_

#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <vector>
#include <map>
#include <stdlib.h>
#include <string.h>
#include <iosfwd>
#include <fstream>

struct strcoordinate{
	int pos;
	std::string chr;
};

struct strvcfentry{
	std::string header;
	strcoordinate start;
	strcoordinate stop;
	short type; //0=DEL,1=DUP,2=INV,3=TRA
	std::map<std::string,std::string> calls;
	int sup_lumpy;
	int caller_id;
	std::vector<int> caller_supports;
	std::pair<bool,bool> strands;
	std::pair<int,int> num_reads; //ref alt
	std::string genotype;
	int sv_len;
	//int num_reads;
};


std::vector<strvcfentry> parse_vcf(std::string filename, int min_svs);


#endif /* FORCE_CALLING_VCF_PARSER_H_ */
