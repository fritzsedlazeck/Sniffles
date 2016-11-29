/*
 * TNode.h
 *
 *  Created on: Jun 23, 2015
 *      Author: fsedlaze
 */

#ifndef TREE_TNODE_H_
#define TREE_TNODE_H_

#include <iostream>
#include <cstdlib>
#include <cmath>
#include "../sub/Breakpoint.h"
//#include "TNode.h"
class TNode {
private:

	Breakpoint * data;
	//int value;
	int height;
	int MAX_DIST;
	void init() {
		this->parent = NULL;
		this->left = NULL;
		this->right = NULL;
		MAX_DIST=500;
	}
public:
	TNode * parent;
	TNode * left;
	TNode * right;
	TNode() {
		height=0;
		init();
		this->data=NULL;
	}
	TNode(Breakpoint * point) {
		init();
		this->data=point;
		height=0;
	}

	~TNode() {

	}

	Breakpoint * get_data() {
		return data;
	}
	int get_height() {
		return height;
	}
	void set_height(int val) {
		this->height = val;
	}
};

#endif /* TREE_TNODE_H_ */
