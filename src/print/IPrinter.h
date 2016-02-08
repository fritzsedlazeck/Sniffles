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

#include "../tree/Intervall_bed.h"
#include "api/BamReader.h"
#include "../Ignore_Regions.h"
#include "../sub/Breakpoint.h"

class IPrinter {
protected:
	FILE *file;
	uint id;
	RefVector ref;
	BamParser *mapped_file;
	IntervallTree_bed bed_tree;
	Leaf *root;

	virtual void print_header()=0;
	virtual void print_body(Breakpoint * &SV, RefVector ref)=0;
	long calc_pos(long pos, RefVector ref, std::string &chr);
	std::string get_chr(long pos, RefVector ref);
	std::string get_type(char type);

public:
	IPrinter() {
		id = 0;
		root = NULL;
		//we just need the ref information:
	}
	virtual ~IPrinter() {
		delete mapped_file;
		fclose(this->file);
	}
	;
	void printSV(Breakpoint * SV) {
		print_body(SV, ref);
	}
	void init() {
		if(!Parameter::Instance()->output_vcf.empty()){
			file = fopen(Parameter::Instance()->output_vcf.c_str(), "w");
		}else if(!Parameter::Instance()->output_maria.empty()){
			file = fopen(Parameter::Instance()->output_maria.c_str(), "w");
		}
		print_header();
		BamParser *mapped_file = new BamParser(Parameter::Instance()->bam_files[0]);
		this->ref = mapped_file->get_refInfo();

		if (!Parameter::Instance()->ignore_regions_bed.empty()) {
			std::cout << "Cross checking..." << std::endl;
			initialize_bed(bed_tree, root, ref);
		}

	}
};

#endif /* PRINT_IPRINTER_H_ */
