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

void IPrinter::store_readnames(std::vector<int> names, int id) {
	name_str tmp;
	tmp.svs_id = id; //stays the same
	for (size_t i = 0; i < names.size(); i++) {
		tmp.read_name = names[i];
		fwrite(&tmp, sizeof(struct name_str), 1, this->tmp_file);
	}
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
	string tmp;
	if (type & DEL) {
		tmp += "DEL";
	}
	if (type & INV) {
		if (!tmp.empty()) {
			tmp += '/';
		}
		tmp += "INV";
	}
	if (type & DUP) {
		if (!tmp.empty()) {
			tmp += '/';
		}
		tmp += "DUP";
	}
	if (type & INS) {
		if (!tmp.empty()) {
			tmp += '/';
		}
		tmp += "INS";
	}
	if (type & TRA) {
		if (!tmp.empty()) {
			tmp += '/';
		}
		//tmp += "BND";
		tmp += "TRA";
	}

	return tmp; // should not occur!
}
// Get current date/time, format is YYYY-MM-DD.HH:mm:ss
const std::string IPrinter::currentDateTime() {
	time_t now = time(0);
	struct tm tstruct;
	char buf[80];
	tstruct = *localtime(&now);
	// Visit http://en.cppreference.com/w/cpp/chrono/c/strftime
	// for more information about date/time format
	strftime(buf, sizeof(buf), "%Y%m%d", &tstruct);
	return buf;
}
