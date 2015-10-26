/*
 * IntervallTree.h
 *
 *  Created on: Jun 23, 2015
 *      Author: fsedlaze
 */

#ifndef TREE_INTERVALLTREE_H_
#define TREE_INTERVALLTREE_H_

#include <vector>

#include "TNode.h"
class IntervallTree {
private:
	int max(int, int);
	TNode * srl(TNode *&);
	TNode * drl(TNode *&);
	TNode * srr(TNode *&);
	TNode * drr(TNode *&);
public:
	void insert(Breakpoint * point, TNode *&);
	void del(Breakpoint * point, TNode *&);
	int deletemin(TNode *&);
	void find(Breakpoint * point, TNode *&);
	TNode * findmin(TNode*);
	TNode * findmax(TNode*);
	void makeempty(TNode *&);
	void copy(TNode * &, TNode *&);
	TNode * nodecopy(TNode *&);
	void preorder(TNode*);
	void inorder(TNode*,TNode * root);
	void postorder(TNode*);
	int bsheight(TNode*);
	void get_breakpoints(TNode *p,std::vector<Breakpoint *> & points);
	int nonodes(TNode*);
	void collapse_intervalls(TNode *&p);
};

#endif /* TREE_INTERVALLTREE_H_ */
