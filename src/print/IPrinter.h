/*
 * IPrinter.h
 *
 *  Created on: Aug 24, 2015
 *      Author: fsedlaze
 */

#ifndef PRINT_IPRINTER_H_
#define PRINT_IPRINTER_H_
#include <vector>
#include <iostream>
#include "api/BamReader.h"
#include "../sub/Breakpoint.h"

class IPrinter{
protected:
	RefVector ref;
	BamParser *mapped_file;
	virtual void print_header()=0;
	virtual void print_body(std::vector<Breakpoint *> &SV, RefVector ref)=0;
	long calc_pos(long pos, RefVector ref,std::string &chr);
	std::string get_chr(long pos, RefVector ref);
	std::string get_type(char type);
public:
	IPrinter(){
		//we just need the ref information:
	}
	virtual ~IPrinter(){
		delete mapped_file;
	};
	void printSV(std::vector<Breakpoint *> &SV){
		print_header();
		print_body(SV,ref);
	}
	void init(){
		BamParser *mapped_file = new BamParser(Parameter::Instance()->bam_files[0]);
		this->ref = mapped_file->get_refInfo();
	}
};


#endif /* PRINT_IPRINTER_H_ */
