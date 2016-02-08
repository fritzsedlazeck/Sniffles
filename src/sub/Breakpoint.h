/*
 * Breakpoint.h
 *
 *  Created on: Jun 23, 2015
 *      Author: fsedlaze
 */

#ifndef SUB_BREAKPOINT_H_
#define SUB_BREAKPOINT_H_

#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <math.h>
#include "../Paramer.h"
#include "../BamParser.h"
#include "../tree/BinTree.h"

const unsigned char DEL = 0x01; // hex for 0000 0001
const unsigned char DUP = 0x02; // hex for 0000 0010
const unsigned char INS = 0x04; // hex for 0000 0100
const unsigned char INV = 0x08; // hex for 0000 1000
const unsigned char TRA = 0x10; // hex for 0001 0000

struct region_ref_str{ //not very nice!
	std::string read_seq;
	long read_aln_pos;
	bool direction;
	long start;
	long stop;
	std::string ref;
};

struct svs_breakpoint_str{
	long min_pos;
	long max_pos;
	long most_support;
};
struct read_str {
	//to identify
	std::string name;
	region_ref_str aln; //maybe we can use this!
	short type; //split reads, cigar or md string
	//for later assessment:
	pair<bool, bool> strand;
	pair<long,long> coordinates; // I could use the bin tree for that!
	char SV; // bit vector
};
struct position_str {
	svs_breakpoint_str start;
	svs_breakpoint_str stop;
	//int pos; //the chromosomes are encoded over the positions.
	std::map<std::string,read_str> support;
	//std::vector<read_str> support; // map?? -> no duplicated reads, easy to catch up which read is included.
	int coverage;
	int lowmq_cov;
	short read_start;
	short read_stop;
};


//TODO define region object  and inherit from that. Plus define avoid region objects for mappability problems.
class Breakpoint {
private:
	position_str positions;
	std::vector<std::string> strand;
	std::string supporting_types;
	char sv_type;
	std::string sv_debug;
	std::string ref_seq;
	std::vector<short> support;
	double cov;
	short type_support;
	//for phasing:
	BinTree grouped;
	tree_node * grouped_node;
	long length;

	void summarize_support(short type);
	//void summarize_strand(pair<bool, bool> strand, std::vector<short>& array);
	void summarize_type(char SV, std::vector<short>& array);
	//std::string translate_strand(short id);
	char eval_type(std::vector<short> SV);
	std::string rev_complement(std::string seq);
	bool is_in(short id);
	std::string translate_strand(pair<bool, bool> strand);
public:
	Breakpoint(position_str sv, int cov,long len) {
		sv_type=' ';
		type_support=-1;
		this->positions = sv;
		//std::cout << "Break1: " << positions.start << " " << positions.stop << std::endl;

		//std::cout << "Break2: " << positions.start << " " << positions.stop<< std::endl;
		this->cov = cov;
		this->grouped_node=NULL;
		this->set_length(len);
	}
	~Breakpoint() {

	}

	int get_support();
	long overlap(Breakpoint * tmp);
	position_str get_coordinates() {
		return this->positions;
	}
	void predict_SV();
	std::string to_string(RefVector ref);

	void add_read(Breakpoint * point);

	double get_cov() {
		return cov;
	}
	std::string get_chr(long pos, RefVector ref);
	long calc_pos(long pos, RefVector ref);
	char get_SVtype();
	std::string get_strand(int num_best);
	std::string get_ref_seq() {
		return this->ref_seq;
	}
	void set_ref_seq(std::string seq) {
		this->ref_seq = seq;
	}
	long get_length(){
		return length;
	}

	void set_length(long len){
		this->length=len;
	}
	std::string get_supporting_types(){
		return this->supporting_types;
	}

	void add_grouped(int id){
		this->grouped.insert(this->grouped_node, id);
	}
	vector<int> get_groupted(){
		vector<int> tmp;
		this->grouped.get_nodes(this->grouped_node,tmp);
		return tmp;
	}
	void calc_support();
};

#endif /* SUB_BREAKPOINT_H_ */
