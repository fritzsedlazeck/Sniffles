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
#include "../Paramer.h"
#include <iostream>
#include <omp.h>

void clarify(std::vector<Breakpoint *> & points);
std::vector<Breakpoint *> detect_breakpoints(std::string filename);
//void screen_for_events(Node * list,IntervallTree & bst ,TNode *&root, int cov, int lowMQ_cov,RefVector ref);
bool screen_for_events(Alignment * tmp, IntervallTree & bst, TNode *&root, RefVector ref, int cov);
void add_events(Alignment * tmp, std::vector<str_event> events, short type, long ref_space, IntervallTree & bst, TNode *&root, int cov,std::string read_seq);
void add_splits(Alignment * tmp, std::vector<aln_str> events, short type, RefVector ref, IntervallTree & bst, TNode *&root, int cov,std::string read_seq);
void estimate_parameters(std::string read_filename);
#endif /* SUB_DETECT_BREAKPOINTS_H_ */
