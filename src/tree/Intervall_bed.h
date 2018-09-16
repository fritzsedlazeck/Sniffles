/*
 * Intervall_bed.h
 *
 *  Created on: Feb 4, 2016
 *      Author: fsedlaze
 */

#ifndef TREE_INTERVALL_BED_H_
#define TREE_INTERVALL_BED_H_

#include "Leaf.h"
#include <vector>
#include "../Paramer.h"

class IntervallTree_bed {
private:
	int max(int, int);
	Leaf * srl(Leaf *&);
	Leaf * drl(Leaf *&);
	Leaf * srr(Leaf *&);
	Leaf * drr(Leaf *&);
public:
	void insert(long start, long stop, Leaf *&);
	int deletemin(Leaf *&);
	bool is_in(long pos, Leaf *&); //true if found
	Leaf * findmin(Leaf*);
	Leaf * findmax(Leaf*);
	void makeempty(Leaf *&);
	void copy(Leaf * &, Leaf *&);
	Leaf * nodecopy(Leaf *&);
	void preorder(Leaf*);
	void inorder(Leaf*, Leaf * root);
	void postorder(Leaf*);
	int bsheight(Leaf*);
	int nonodes(Leaf*);
	void del(long start, long stop, Leaf * &p);
};

#endif /* TREE_INTERVALL_BED_H_ */
