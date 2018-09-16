/*
 * IntervallContainer.h
 *
 *  Created on: Nov 2, 2016
 *      Author: fsedlaze
 */

#ifndef TREE_INTERVALLCONTAINER_H_
#define TREE_INTERVALLCONTAINER_H_

class IntervallContainer {
protected:

public:
	IntervallContainer() {

	}
	virtual ~IntervallContainer() {

	}
	virtual void insert(Breakpoint * point, TNode *&)=0;
	virtual void get_breakpoints(TNode *p, std::vector<Breakpoint *> & points)=0;
	virtual void clear(TNode*&)=0;
	virtual void print(TNode *p)=0;
};

#endif /* TREE_INTERVALLCONTAINER_H_ */
