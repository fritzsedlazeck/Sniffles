/*
 * PlaneSweep_slim.h
 *
 *  Created on: Mar 9, 2016
 *      Author: fsedlaze
 */

#ifndef PLANE_SWEEP_PLANESWEEP_SLIM_H_
#define PLANE_SWEEP_PLANESWEEP_SLIM_H_
#include <vector>
#include <iostream>
#include <string.h>
#include <list>
#include "../Paramer.h" //for testing/debug
struct pair_str{
	int position;
	int coverage;
};


class PlaneSweep_slim {
private:
	int boundary;
	int cov;
	std::vector<pair_str> entries;
public:
	PlaneSweep_slim() {
		boundary=100;
		cov=0;
	}
	~PlaneSweep_slim(){
		entries.clear();
	}
	void release_pos(int new_start);
	pair_str add_mut(int pos,int cov, int min_cov);
	void finalyze();

};

#endif /* PLANE_SWEEP_PLANESWEEP_SLIM_H_ */
