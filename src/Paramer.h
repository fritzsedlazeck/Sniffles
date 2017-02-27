/*
 * Paramer.h
 *
 *  Created on: Aug 20, 2015
 *      Author: fsedlaze
 */

#ifndef PARAMER_H_
#define PARAMER_H_

#include <string.h>
#include <string>
#include <vector>
#include <stdlib.h>
#include <iostream>
#include <ctime>
struct region_str {
	std::string chr;
	int start;
	int stop;
};

class Parameter {
private:
	Parameter() {
		window_thresh=10;//TODO check!
		version="1.0.3";
		huge_ins = 2000;//TODO check??
	}
	~Parameter() {

	}
	static Parameter* m_pInstance;
	std::vector<region_str> regions;
public:
	std::string output_vcf;
	std::string output_bedpe;
	std::string ref_seq;
	std::string read_name;
	std::string ignore_regions_bed;
	std::string tmp_file;
	std::string version;

	std::vector<std::string> bam_files;
	short min_mq;
	short report_n_reads;
	short corridor;

	double error_rate;
	double score_treshold;
	double min_allelel_frequency;
	//double min_num_mismatches;

	int window_thresh;
	int min_support;
	int max_splits;
	int max_dist;
	int min_length;
	int min_reads_phase;
	int num_threads;
	int max_readlength;
	int min_grouping_support; //min num reads supporting the overlap of two SVs
	int huge_ins;
	int max_dist_alns;

//	bool splitthreader_output;
	bool debug;
	bool genotype;
	bool phase;

	void set_regions(std::string reg) {
		size_t i = 0;
		region_str tmp;
		short sep;
		while (i < reg.size()) {
			tmp.chr = reg.substr(i, reg.find_first_of(':'));
			i += tmp.chr.size() + 1;
			tmp.start = atoi(&reg[i]);
			i += reg.find_first_of('-') + 1;
			tmp.stop = atoi(&reg[i]);
			i += reg.find_first_of(';') + 1;
			regions.push_back(tmp);
		}
		std::cout << "found regions: " << regions.size() << std::endl;
	}
	bool overlaps(std::string chr, int start, int stop) {
		for (size_t i = 0; i < regions.size(); i++) {
			if (strcmp(chr.c_str(), regions[i].chr.c_str()) == 0 && (abs(start - regions[i].start) < max_dist && abs(stop - regions[i].stop) < max_dist)) {
				return true;
			}
		}
		return false;
	}

	static Parameter* Instance() {
		if (!m_pInstance) {
			m_pInstance = new Parameter;
		}
		return m_pInstance;
	}

	double meassure_time(clock_t begin ,std::string msg){
		clock_t end = clock();
		double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
		std::cout << msg<<" " << elapsed_secs<<std::endl;
		return elapsed_secs;
	}

};

#endif /* PARAMER_H_ */
