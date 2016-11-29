/*
 * List.cpp
 *
 *  Created on: Nov 2, 2016
 *      Author: fsedlaze
 */

#include "IntervallList.h"

void IntervallList::insert(Breakpoint * point, TNode *& note) {

	for (size_t i = 0; i < this->breakpoints.size(); i++) {
		long score = this->breakpoints[i]->get_data()->overlap(point);
		if (score == 0) {
		//	std::cout<<"overlap: "<<this->breakpoints[i]->get_data()->get_coordinates().support.size()<<std::endl;
			this->breakpoints[i]->get_data()->add_read(point);
			delete point;
			return;
		}/* else if (score < 0) {
			TNode * p = new TNode(point);
			this->breakpoints.insert((this->breakpoints.begin()+i),p);
			return;
		}*/
	}
	TNode * p = new TNode(point);
	this->breakpoints.push_back(p);
}
void IntervallList::get_breakpoints(TNode *p, std::vector<Breakpoint *> & points) {
	for (size_t i = 0; i < this->breakpoints.size(); i++) {
		points.push_back(this->breakpoints[i]->get_data());
	}
}

void IntervallList::clear(TNode*&){
	this->breakpoints.clear();
}
void IntervallList::print(TNode *p){
	std::cout<<"Print:"<<std::endl;
	for (size_t i = 0; i < this->breakpoints.size(); i++) {
		std::cout<<"( "<<this->breakpoints[i]->get_data()->get_coordinates().start.min_pos<<"-"<<this->breakpoints[i]->get_data()->get_coordinates().stop.max_pos<<" )\t";
	}
	std::cout<<std::endl;
}
