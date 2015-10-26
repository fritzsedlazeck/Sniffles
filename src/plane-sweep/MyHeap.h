/*
 * MyHeap.h
 *
 *  Created on: Apr 14, 2015
 *      Author: fsedlaze
 */

#ifndef MYHEAP_H_
#define MYHEAP_H_
#include <vector>
#include <iostream>
#include "IContainer.h"
#include "Node.h"
class MyHeap: public IContainer{
	private:
    std::vector<Node*> heap;
    int left(int parent);
    int right(int parent);
    int parent(int child);
    void heapifyup(int index);
    void heapifydown(int index);
public:
	MyHeap(){};
    ~MyHeap(){};
    void push(Alignment * read);
    Node * get_entries(int new_start);
    int size() { return this->heap.size(); }
    void print();
};

#endif /* MYHEAP_H_ */
