/*
 * IPrinter.cpp
 *
 *  Created on: Aug 24, 2015
 *      Author: fsedlaze
 */

#include "IPrinter.h"

std::string IPrinter::get_chr(long pos, RefVector ref) {
	//	std::cout << "pos: " << pos << std::endl;
	size_t id = 0;
	while (id < ref.size() && pos >= 0) {
		pos -= ((long) ref[id].RefLength + (long) Parameter::Instance()->max_dist);
		//	std::cout << id << std::endl;
		id++;
	}

	return ref[id - 1].RefName;
}

long IPrinter::calc_pos(long pos, RefVector ref, std::string &chr) {
	size_t i = 0;
	pos -= (ref[i].RefLength + Parameter::Instance()->max_dist);

	while (i < ref.size() && pos >= 0) {
		i++;
		pos -= ((long) ref[i].RefLength + (long) Parameter::Instance()->max_dist);
	}
	chr = ref[i].RefName;
	return pos + ref[i].RefLength + (long) Parameter::Instance()->max_dist;
}

std::string IPrinter::get_type(char type) {
	if (type & DEL) {
		return "DEL";
	}
	if (type & INV) {
		return "INV";
	}
	if (type & DUP) {
		return "DUP";
	}
	if (type & INS) {
		return "INS";
	}
	if (type & TRA) {
		return "TRA";
	}
	return "WTF"; // should not occur!
}
