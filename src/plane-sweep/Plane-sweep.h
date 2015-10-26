/*
 * Plane-sweep.h
 *
 *  Created on: Apr 14, 2015
 *      Author: fsedlaze
 */

#ifndef PLANESWEEP_H_
#define PLANESWEEP_H_

#include <iostream>
#include <string.h>
#include "IContainer.h"
#include "../Alignment.h"

#include "MyList.h"
#include "MyHeap.h"
//#include "MyHeap.h"

class PlaneSweep {
private:
	MyList* current_reads; // will be used as heap.
//	void release_pos(int new_start);
	int RefID; //the current chr/contig;
	int curr_pos;
public:
	PlaneSweep(){
		curr_pos=0;
		RefID=-1;
		current_reads=new MyList();
		//current_reads=new MyHeap();
	}
	~PlaneSweep(){
		delete current_reads;
	}
	void release_pos(int new_start);
	void add_read(Alignment* read);
	void finalyze();
	Node * get_current_reads(){
		return current_reads->get_start();
	}
	int get_num_reads(){
		return current_reads->get_normal_reads();
	}
	int get_num_SVreads(){
		return current_reads->get_split_reads();
	}
	int get_num_lowMQ_reads(){
		return current_reads->get_lowMqcov();
	}
	int get_RefID(){
		return RefID;
	}
};



#endif /* PLANESWEEP_H_ */
