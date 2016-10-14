/*
 * PlaneSweep_slim.cpp
 *
 *  Created on: Mar 9, 2016
 *      Author: fsedlaze
 */

#include "PlaneSweep_slim.h"

pair_str PlaneSweep_slim::add_mut(int pos,int new_cov, int min_cov) {
	//check if we need to release reads:
	std::vector<pair_str>::iterator j = entries.begin();
	for (size_t i = 0; i < this->entries.size() && pos > entries[i].position; i++) {
		//no need to record ending events. We just search for starting events!
		this->cov-=entries[i].coverage;
		j++;
	}

	//erase old events:
	entries.erase(entries.begin(), j);


	//add current mut:
	//insert the stop coordinate!
	pair_str tmp;
	tmp.position=pos+boundary;
	tmp.coverage=new_cov;
	this->entries.push_back(tmp);
	this->cov+=new_cov;

	//record if we met the threshold:

	tmp.position=-1; //flag for not meeting the threshold
	if(this->cov>min_cov){
		tmp.coverage = this->cov;
		tmp.position = pos;
	}
	return tmp;
}
void PlaneSweep_slim::finalyze() {
	entries.clear();
	this->cov=0;
}
