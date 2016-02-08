/*
 * Leaf.h
 *
 *  Created on: Feb 4, 2016
 *      Author: fsedlaze
 */

#ifndef TREE_LEAF_H_
#define TREE_LEAF_H_
#include "../Paramer.h"
#include <iostream>
#include <cstdlib>
#include <cmath>
class Leaf {
private:
	long start;
	long stop;
	int height;
	void init() {
		height = 0;
		this->parent = NULL;
		this->left = NULL;
		this->right = NULL;
	}

public:
	Leaf * parent;
	Leaf * left;
	Leaf * right;

	Leaf(long start, long stop) {
		this->start = start;
		this->stop = stop;
		init();
	}
	int overlap(long position) {
		if (abs(position - get_start()) < Parameter::Instance()->max_dist && abs(position - get_stop()) < Parameter::Instance()->max_dist) {
			return 0;
		}
		return (position - start); //((start < position) && (stop > position));
	}
	int overlap(long start, long stop) {
		if (abs(start - get_start()) < Parameter::Instance()->max_dist && abs(stop - get_stop()) < Parameter::Instance()->max_dist) {
			return 0;
		}
		//as abstraction lets try the start+stop coordinate!
		return (start - get_start()); // + (stop-get_stop());
	}

	long get_start() {
		return start;
	}
	long get_stop() {
		return stop;
	}
	int get_height() {
		return height;
	}
	void set_height(int val) {
		this->height = val;
	}
};

#endif /* TREE_LEAF_H_ */
