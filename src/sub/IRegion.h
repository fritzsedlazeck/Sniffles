/*
 * IRegion.h
 *
 *  Created on: Aug 27, 2015
 *      Author: fsedlaze
 */

#ifndef SUB_IREGION_H_
#define SUB_IREGION_H_
#include "../Paramer.h"
#include "../BamParser.h"
struct read_str {
	//to identify
	std::string name;
	short type;
	//for later assessment:
	pair<bool, bool> strand;
	char SV; // bits vector
};
struct position_str {
	long start;
	long stop;
	//int pos; //the chromosomes are encoded over the positions.
	std::vector<read_str> support;
	int coverage;
	int lowmq_cov;
	short read_start;
	short read_stop;
};

class IRegion {
protected:
	position_str start;

public:
	IRegion(position_str reg) {
		this->start = reg;
		//std::cout << "Break1: " << start.start << " " << start.stop << std::endl;
		if (reg.start > reg.stop) {
			this->start.start = reg.stop;
			this->start.stop = reg.start;
		}
	}
	virtual ~IRegion() {

	}

	virtual std::string to_string(RefVector ref)=0;
	virtual long overlap(IRegion * tmp) =0;
	position_str get_coordinates() {
		return this->start;
	}
	virtual int support()=0;
};

#endif /* SUB_IREGION_H_ */
