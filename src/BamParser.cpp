/*
 * parser.cpp
 *
 *  Created on: Apr 17, 2012
 *      Author: fritz
 */

#include "BamParser.h"

BamParser::BamParser(string file){
	vector<string > tmps;
	tmps.push_back(file);

	if(!reader.Open(tmps)){
		cerr<<"BAM Parser: could not open file: "<<file<<endl;
		exit(EXIT_FAILURE);
	}

}

Alignment* BamParser::parseRead(uint16_t mappingQv){

	Alignment *align = new Alignment();
	BamAlignment* al = new BamAlignment();
	while(reader.GetNextAlignmentCore(al[0])){
		if( al->IsMapped() && al->MapQuality > mappingQv){
			al->BuildCharData();
			align->setAlignment(al);
			return align;
		}
	}
	return align;

}
void BamParser::parseReadFast(uint16_t mappingQv,Alignment*& align){

//	Alignment *align = new Alignment();
	BamAlignment* al = align->getAlignment();
//	getSequence().first
//	align->initSequence();
	align->getQueryBases().clear();
	align->clear_QueryBases();
	while(reader.GetNextAlignmentCore(al[0])){

		if( al->IsMapped() && al->MapQuality > mappingQv){
			al->BuildCharData();
			align->setAlignment(al);
			return;
		}
	}
}
RefVector BamParser::get_refInfo(){
	return reader.GetReferenceData();
}

string BamParser::get_header(){
	return reader.GetHeaderText();
}
