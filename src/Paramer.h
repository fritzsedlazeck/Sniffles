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
struct region_str {
	std::string chr;
	int start;
	int stop;
};

class Parameter {
private:
	Parameter() {
		score_treshold = 0;
		min_num_mismatches = 0.3;
	}
	~Parameter() {

	}
	static Parameter* m_pInstance;
	std::vector<region_str> regions;
public:
	std::string output_vcf;
	std::string output_bede;
	std::string ref_seq;
	std::string read_name;
	std::string output_maria;
	std::string ignore_regions_bed;

	std::vector<std::string> bam_files;
	short min_mq;
	short min_cigar_event;
	short report_n_reads;
	short corridor;

	double error_rate;
	double score_treshold;
	double min_num_mismatches;

	int min_support;
	int max_splits;
	int max_dist;
	int min_length;
	int min_reads_phase;
	int num_threads;

	bool realign;
	bool useMD_CIGAR;

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
};

#endif /* PARAMER_H_ */
