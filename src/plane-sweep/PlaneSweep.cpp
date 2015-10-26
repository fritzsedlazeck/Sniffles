/*
 * PlaneSweep.cpp
 *
 *  Created on: Apr 14, 2015
 *      Author: fsedlaze
 */

#include "Plane-sweep.h"

void PlaneSweep::release_pos(int new_start) {

	Node * curr = this->current_reads->get_entries(new_start);
	while (curr != NULL) {
		delete curr;
		curr = NULL;
		curr = this->current_reads->get_entries(new_start);
	}

}

void PlaneSweep::add_read(Alignment* read) {
	//std::cout<<"\tnew: "<<read->getPosition()<<" "<<read->getPosition()+read->getRefLength()<<std::endl;
	if (this->current_reads->size() > 0 && read->getRefID() != this->RefID) {
		finalyze();
		this->RefID = read->getRefID();
	} else if (this->current_reads->size() == 0) {
		this->RefID = read->getRefID();
	}
	//first check if we can release already some positions:
	release_pos(read->getPosition());

	//insert read to our list if it does not support an SV:
	if(!read->supports_SV() && read->getMappingQual()>20){
		current_reads->push(read);
	}
}

void PlaneSweep::finalyze() {
	//report the remaining positions/reads;
	std::cout << "Finalize" << std::endl;
	release_pos(-1); //flag for releasing all;
}

