/*
 * MyList.h
 *
 *  Created on: Apr 14, 2015
 *      Author: fsedlaze
 */

#ifndef MYLIST_H_
#define MYLIST_H_
#include <iostream>
#include "Node.h"
#include "IContainer.h"
class MyList: public IContainer{
private:
	Node * start;
	Node * last;
	int num_nodes;
	int split_read;
	int num_nodes_lowMQ;
public:
	MyList(){
		split_read=0;
		num_nodes=0;
		num_nodes_lowMQ=0;
		start=NULL;
		last=NULL;
	}
	~MyList(){
		delete [] start;
		num_nodes=0;
		split_read=0;
	}
	Node * pop();
	void push(Alignment * read);
	int size();
	Node * get_entries(int new_start);
	Node * get_start(){
		return start;
	}
	int get_normal_reads(){
		return this->num_nodes; //strange value!
	}
	int get_split_reads(){
		return this->split_read;
	}
	int get_lowMqcov(){
		return this->num_nodes_lowMQ;
	}
};

#endif /* MYLIST_H_ */
