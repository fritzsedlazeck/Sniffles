/*
 * Genotyper.h
 *
 *  Created on: Mar 28, 2016
 *      Author: fsedlaze
 */

#ifndef GENOTYPER_H_
#define GENOTYPER_H_
#include "../Paramer.h"
#include "../print/IPrinter.h"
#include "../tree/Breakpoint_Tree.h"
//#include "../sub/Detect_Breakpoints.h"
struct variant_str{
	std::string chr;
	std::string chr2;
	int pos;
	int pos2;
	int len;
};
struct str_breakpoint_slim{
	int pos;
//	int cov;
	std::map<std::string,bool> rnames;
	//bool is_start;
	std::string chr;
};

class Genotyper{
private:
	Breakpoint_Tree tree;
	breakpoint_node * node;
	std::vector<std::string>  read_SVs(Breakpoint_Tree & tree,breakpoint_node *& node );
	void compute_cov(Breakpoint_Tree & tree,breakpoint_node *& node,std::vector<std::string>  ref_dict);
	void update_file(Breakpoint_Tree & tree,breakpoint_node *& node);
	variant_str get_breakpoint_vcf(string buffer);
	variant_str get_breakpoint_bedpe(string buffer);
	std::string mod_breakpoint_vcf(string buffer, int ref_plus, int ref_minus);
	std::string mod_breakpoint_bedpe(string buffer, int ref_plus, int ref_minus);
	void parse_pos(char * buffer, int & pos, std::string & chr);
	void get_breakpoint_vcf(string buffer,str_breakpoint_slim & start,str_breakpoint_slim & stop);
	void read_SVs(std::map<std::string, std::vector<str_breakpoint_slim> > & entries);
	void update_svs_output(std::map<std::string, std::vector<str_breakpoint_slim> > entries);
	variant_str get_breakpoint_bedpe(string buffer, str_breakpoint_slim & start, str_breakpoint_slim & stop);

public:
	Genotyper(){
		node=NULL;
	}
	~Genotyper(){

	}
	void update_SVs();
	void update_SVs2();
	void update_SVs(std::vector<Breakpoint *> & points,long ref_space);
	std::string assess_genotype(int ref, int support);
};
#endif /* GENOTYPER_H_ */
