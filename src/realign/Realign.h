/*
 * Realign.h
 *
 *  Created on: Aug 24, 2015
 *      Author: fsedlaze
 */

#ifndef REALIGN_REALIGN_H_
#define REALIGN_REALIGN_H_
#include "../sub/Breakpoint.h"
#include <iostream>
#include <fstream>
#include <vector>
#include "IAlignment.h"
#include "SWCPU.h"
struct ref_str {
	long length;
	streampos file_pos;
};

class Realigner {
private:
	std::string read_new_part(long start, long stop);
	vector<ref_str> meta_info;
	size_t buffer_size;
	char*buffer;
	ifstream myfile;
	void init();
	std::string read_chr(short id);
public:
	Realigner() {
		init();
	}
	~Realigner() {
		delete[] buffer;
	}
	void align(std::vector<Breakpoint *> SV);
};

#endif /* REALIGN_REALIGN_H_ */
