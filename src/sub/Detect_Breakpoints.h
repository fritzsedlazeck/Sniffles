/*
 * Detect_Breakpoints.h
 *
 *  Created on: Jun 19, 2015
 *      Author: fsedlaze
 */

#ifndef SUB_DETECT_BREAKPOINTS_H_
#define SUB_DETECT_BREAKPOINTS_H_

#include "../BamParser.h"
#include "../Parser.h"
#include "../Alignment.h"
#include "../plane-sweep/Plane-sweep.h"
#include "../tree/IntervallTree.h"
#include "../tree/TNode.h"
#include "../tree/IntervallContainer.h"
#include "../tree/IntervallList.h"
#include "../Paramer.h"
#include "../print/IPrinter.h"

#include <iostream>
#include <omp.h>

struct hist_str{
	long position;
	int hits;
	std::vector<std::string> names;
};
void clarify(std::vector<Breakpoint *> & points);
void detect_breakpoints(std::string filename, IPrinter *& printer);
//void screen_for_events(Node * list,IntervallTree & bst ,TNode *&root, int cov, int lowMQ_cov,RefVector ref);
bool screen_for_events(Alignment * tmp, IntervallTree & bst, TNode *&root, RefVector ref, int cov);
void add_events(Alignment *& tmp, std::vector<str_event> events, short type, long ref_space, IntervallTree & bst, TNode *&root,long read_id,bool add);
void add_splits(Alignment *& tmp, std::vector<aln_str> events, short type, RefVector ref, IntervallTree & bst, TNode *&root,long read_id,bool add);
void estimate_parameters(std::string read_filename);
bool overlaps(aln_str prev,aln_str curr);
void detect_merged_svs(Breakpoint * point);

long get_ref_lengths(int id, RefVector ref);
std::string TRANS_type(char type);


#endif /* SUB_DETECT_BREAKPOINTS_H_ */
