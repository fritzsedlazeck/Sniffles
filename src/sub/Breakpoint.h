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
#include <vector>
#include <algorithm>
#include "../Paramer.h"
#include "../BamParser.h"
#include "../tree/BinTree.h"


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
//	std::string name;
	long id;
	region_ref_str aln; //maybe we can use this!
	short type; //split reads, cigar or md string
	//for later assessment:
	pair<bool, bool> strand;
	pair<bool, bool> read_strand;
	pair<long,long> coordinates; // I could use the bin tree for that!
	char SV; // bit vector
	int length;
	std::string sequence; //just for indels!
};
struct position_str {
	svs_breakpoint_str start;
	svs_breakpoint_str stop;
	//int pos; //the chromosomes are encoded over the positions.
	std::map<std::string,read_str> support;
	//std::vector<read_str> support; // map?? -> no duplicated reads, easy to catch up which read is included.
	int coverage;
	int lowmq_cov;
	int read_start;
	int read_stop;
};


struct str_types{
	bool is_SR;
	bool is_ALN;
	bool is_Noise;
};

//TODO define region object  and inherit from that. Plus define avoid region objects for mappability problems.
class Breakpoint {
private:
	str_types type;
	position_str positions;
	std::vector<std::string> strand;
	std::string supporting_types;
	char sv_type;
	std::string sv_debug;
	std::string ref_seq;
	//std::vector<short> support;
	short type_support;
	//for phasing:
	BinTree grouped;
	tree_node * grouped_node;
	long length;
	std::string indel_sequence;
	bool should_be_stored;
	int ref_allele;

	void summarize_support(short type);
	//void summarize_strand(pair<bool, bool> strand, std::vector<short>& array);
	void summarize_type(char SV, ushort * array);
	//std::string translate_strand(short id);
	char eval_type(ushort *SV);
	std::string rev_complement(std::string seq);
	bool is_in(short id);
	std::string translate_strand(pair<bool, bool> strand);
	bool is_same_strand(Breakpoint * tmp);
	bool check_SVtype(Breakpoint * break1, Breakpoint * break2);
	bool is_NEST(Breakpoint * next, Breakpoint * curr){
		return (( (*next->get_coordinates().support.begin()).second.SV& NEST) ||( (*curr->get_coordinates().support.begin()).second.SV& NEST) );
	}
public:
	Breakpoint(position_str sv,long len) {
		ref_allele=0;
		should_be_stored=false;
		sv_type |= NA;
		type.is_ALN=((*sv.support.begin()).second.type==0);
		type.is_SR=((*sv.support.begin()).second.type==1);
		type.is_Noise=((*sv.support.begin()).second.type==2);
		type_support=-1;
		this->positions = sv;
		this->grouped_node=NULL;
		this->length=len;
	}
	Breakpoint(position_str sv,long len,char sv_type) {
		ref_allele=0;
		should_be_stored=false;
		this->sv_type =sv_type;
		type.is_ALN=((*sv.support.begin()).second.type==0);
		type.is_SR=((*sv.support.begin()).second.type==1);
		type.is_Noise=((*sv.support.begin()).second.type==2);
		type_support=-1;
		this->positions = sv;
		this->grouped_node=NULL;
		this->length=len;
	}
	~Breakpoint() {

	}

	int get_support();
	long overlap(Breakpoint * tmp);
	long overlap_breakpoint(long start,long stop);
	void set_coordinates(int start, int stop){
		this->positions.start.min_pos=start;
		this->positions.stop.max_pos=stop;
	}
	position_str get_coordinates() {
		return this->positions;
	}
	void predict_SV();
	std::string to_string(RefVector ref);

	void add_read(Breakpoint * point);

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
	str_types get_types(){
		return this->type;
	}
	std::string get_read_names();
	std::vector<long> get_read_ids();
	std::string to_string();
	std::string get_sequence(){
		return this->indel_sequence;
	}
	void set_valid(bool valid){
		this-> should_be_stored=valid;
	}
	bool get_valid(){
		return this->should_be_stored;
	}
	int get_refcount(){
		return this->ref_allele;
	}
	void set_refcount(int value){
		this->ref_allele+=value;
	}
};

#endif /* SUB_BREAKPOINT_H_ */
