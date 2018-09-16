/*
 * MyHeap.cpp
 *
 *  Created on: Apr 15, 2015
 *      Author: fsedlaze
 */
#include "MyHeap.h"


void MyHeap::push(Alignment * read){
	//std::cout<<"push"<<std::endl;
	Node * tmp=new Node(read);
    heap.push_back(tmp);
    heapifyup(heap.size() - 1);
    //std::cout<<"push end"<<std::endl;
}

Node * MyHeap::get_entries(int new_start)
{

    if(!heap.empty() && heap[0]->get_stop() <= new_start){
    	Node * min = heap[0];
		heap[0] = heap.at(heap.size() - 1);
		heap.pop_back();
		heapifydown(0);
	    std::cout<<"get entries end "<<this->size()<<std::endl;
	    print();
		return min;
    }
    return NULL;
}

void MyHeap::print()
{
   std::vector<Node*>::iterator pos = heap.begin();
    std::cout << "Heap = ";
    while ( pos != heap.end() ) {
        std::cout << (*pos)->get_stop() << " ";
        ++pos;
    }
    std::cout << std::endl;
}

void MyHeap::heapifyup(int index)
{
    while ( ( index > 0 ) && ( parent(index) >= 0 ) &&
            ( heap[parent(index)]->get_stop() > heap[index]->get_stop() ) )
    {
        Node * tmp = heap[parent(index)];
        heap[parent(index)] = heap[index];
        heap[index] = tmp;
        index = parent(index);
    }
}

void MyHeap::heapifydown(int index)
{
    int child = left(index);
    if ( ( child > 0 ) && ( right(index) > 0 ) &&
         ( heap[child]->get_stop() > heap[right(index)]->get_stop() ) )
    {
        child = right(index);
    }
    if ( child > 0 )
    {
        Node* tmp = heap[index];
        heap[index] = heap[child];
        heap[child] = tmp;
        heapifydown(child);
    }
}

int MyHeap::left(int parent)
{
    int i = ( parent << 1 ) + 1; // 2 * parent + 1
    return ( i < heap.size() ) ? i : -1;
}

int MyHeap::right(int parent)
{
    int i = ( parent << 1 ) + 2; // 2 * parent + 2
    return ( i < heap.size() ) ? i : -1;
}

int MyHeap::parent(int child)
{
    if (child != 0)
    {
        int i = (child - 1) >> 1;
        return i;
    }
    return -1;
}
