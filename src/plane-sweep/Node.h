/*
 * Node.h
 *
 *  Created on: Apr 14, 2015
 *      Author: fsedlaze
 */

#ifndef NODE_H_
#define NODE_H_
#include "../Alignment.h"


class Node{
private:
	//Alignment * read;
	int start;
	int stop;
	short mq;
	Node * next;
	short times; //indicate if two reads stop at the same location
	bool processed;
	bool support; //to del
public:
	Node(){
		processed=false;
		next=NULL;
		times=0;
		start=0;
		stop=0;
	//	read=NULL;
		mq=0;
	}
	Node(Alignment * read){
		processed=false;
		next=NULL;
		times=1;
		mq=read->getMappingQual();
		//this->read=read;
		this->start=read->getPosition();
		this->stop=this->start+read->getRefLength();
		support=read->supports_SV();
	}
	~Node(){

	}
	bool supports_SV(){
		return support;
	}
	Node * get_next(){
		return this->next;
	}
	int get_stop(){
		return this->stop;
	}
	void set_next(Node *next){
		this->next=next;
	}
	/*Alignment * get_read(){
		return this->read;
	}
	void set_read(Alignment * read){
		this->start=read->getPosition();
		this->stop=read->getRefLength();
		this->read=read;
	}*/
	void set_processed (bool value){
		this->processed=value;
	}

	bool was_processed(){
		return this->processed;
	}
	uint16_t getmq(){
		return mq;
	}
};


#endif /* NODE_H_ */
