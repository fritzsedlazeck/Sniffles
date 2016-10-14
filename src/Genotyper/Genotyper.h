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
struct variant_str{
	std::string chr;
	std::string chr2;
	int pos;
	int pos2;
};
class Genotyper{
private:
	Breakpoint_Tree tree;
	breakpoint_node * node;
	void read_SVs(Breakpoint_Tree & tree,breakpoint_node *& node );
	void compute_cov(Breakpoint_Tree & tree,breakpoint_node *& node);
	void update_file(Breakpoint_Tree & tree,breakpoint_node *& node);
	variant_str get_breakpoint_vcf(char *buffer);
	variant_str get_breakpoint_bedpe(char *buffer);
	std::string mod_breakpoint_vcf(char *buffer, int ref);
	std::string mod_breakpoint_bedpe(char *buffer, int ref);

public:
	Genotyper(){
		node=NULL;
	}
	~Genotyper(){

	}
	void update_SVs();
};
#endif /* GENOTYPER_H_ */
