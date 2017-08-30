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
#include <time.h>
#include "../tree/Intervall_bed.h"
#include "api/BamReader.h"
#include "../Ignore_Regions.h"
#include "../sub/Breakpoint.h"
#include "../cluster/Cluster_SVs.h"
#include <math.h>

double const uniform_variance = 0.2886751; //sqrt(1/12) see variance of uniform distribution -> std

class IPrinter {
protected:
	FILE *file;
	FILE *distances;
	FILE *tmp_file;
	uint id;
	RefVector ref;
	BamParser *mapped_file;
	IntervallTree_bed bed_tree;
	Leaf *root;

	virtual void print_header()=0;
	virtual void print_body(Breakpoint * &SV, RefVector ref)=0;
	virtual void print_body_recall(Breakpoint * &SV, RefVector ref)=0;

	long calc_pos(long pos, RefVector ref, std::string &chr);
	std::string get_chr(long pos, RefVector ref);
	std::string get_type(char type);
	void sort_insert(int pos, std::vector<int> & positons);
	bool is_huge_ins(Breakpoint * &SV);
public:

	IPrinter() {
		id = 0;
		root = NULL;
		//we just need the ref information:
	}
	virtual ~IPrinter() {
		delete mapped_file;

	}

	void printSV(Breakpoint * SV) {
		if(Parameter::Instance()->input_vcf.empty()){
			print_body(SV, ref);
		}else{
			print_body_recall(SV,ref);
		}
	}
	void init() {
		try {
			if (!Parameter::Instance()->output_vcf.empty()) {
				file = fopen(Parameter::Instance()->output_vcf.c_str(), "w");
			} else if (!Parameter::Instance()->output_bedpe.empty()) {
				file = fopen(Parameter::Instance()->output_bedpe.c_str(), "w");
			}
		} catch (...) {
			std::cerr << "Output file could not be created. Please check if path exists and if you have write permissions." << std::endl;
			exit(0);
		}
		if (file == NULL) {
			std::cerr << "Output file could not be created. Please check if path exists and if you have write permissions." << std::endl;
			exit(EXIT_FAILURE);
		}


		BamParser *mapped_file = new BamParser(Parameter::Instance()->bam_files[0]);
		this->ref = mapped_file->get_refInfo();
		print_header();
		if (!Parameter::Instance()->ignore_regions_bed.empty()) {
			std::cout << "Cross checking..." << std::endl;
			initialize_bed(bed_tree, root, ref);
		}
		string tmp_name_file = Parameter::Instance()->tmp_file;
		if (Parameter::Instance()->phase) {
			tmp_name_file += "Names";
			tmp_file = fopen(tmp_name_file.c_str(), "wb");
		}
	}
	bool to_print(Breakpoint * &SV, pair<double, double> &std, pair<double, double> & kurtosis, double & std_length);
	void store_readnames(std::vector<long> names, int id);
	void close_file() {
		fclose(this->file);
	}
	void comp_std(Breakpoint * &SV, double & std_start, double & std_stop);
	void comp_std_med(Breakpoint * &SV, double & std_start, double & std_stop);
	pair<double, double> comp_std_quantile(Breakpoint * &SV, pair<double, double>& std, double & std_lenght);
	const std::string currentDateTime();
};

#endif /* PRINT_IPRINTER_H_ */
