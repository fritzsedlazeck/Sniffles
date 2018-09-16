/*
 * List.h
 *
 *  Created on: Nov 2, 2016
 *      Author: fsedlaze
 */

#ifndef TREE_INTERVALLLIST_H_
#define TREE_INTERVALLLIST_H_
#include <vector>

#include "TNode.h"
#include "IntervallContainer.h"

class IntervallList:public IntervallContainer {
private:
	std::vector<TNode * > breakpoints;
public:
	IntervallList(){

	}
	~IntervallList(){

	}
	void insert(Breakpoint * point, TNode *&);
	void get_breakpoints(TNode *p,std::vector<Breakpoint *> & points);
	void clear(TNode*&);
	void print(TNode *p);
};


#endif /* TREE_INTERVALLLIST_H_ */
