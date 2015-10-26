/*
 * Breakpoint.cpp
 *
 *  Created on: Sep 1, 2015
 *      Author: fsedlaze
 */

/*
 * Breakpoint.h
 *
 *  Created on: Jun 23, 2015
 *      Author: fsedlaze
 */
#include "Breakpoint.h"

//TODO define region object  and inherit from that. Plus define avoid region objects for mappability problems.

std::string Breakpoint::translate_strand(pair<bool, bool> strand_pair) {
	if (strand_pair.first && strand_pair.second) {
		return "++";
	} else if (strand_pair.first && !strand_pair.second) {
		return "+-";
	} else if (!strand_pair.first && strand_pair.second) {
		return "-+";
	} else if (!strand_pair.first && !strand_pair.second) {
		return "--";
	}
	return " ";
}

void Breakpoint::summarize_type(char SV, std::vector<short>& array) {
	//std::string ss;
	if (SV & DEL) {
		//	ss += "DEL; ";
		array[0]++;
	}

	if (SV & DUP) {
		//	ss += "DUP; ";
		array[1]++;
	}

	if (SV & INS) {
		//	ss += "INS; ";
		array[2]++;
	}

	if (SV & INV) {
		//	ss += "INV; ";
		array[3]++;
	}

	if (SV & TRA) {
		//	ss += "TRA; ";
		array[4]++;
	}
	//return ss;
}
/*std::string Breakpoint::translate_strand(short id) {
 switch (id) {
 case 0:
 return "++";
 break;
 case 1:
 return "+-";
 break;
 case 2:
 return "-+";
 break;
 case 3:
 return "--";
 break;
 }
 return "";
 }*/
char Breakpoint::get_SVtype() {
	if (sv_type == ' ') {
		//	std::cout<<"was not set"<<std::endl;
		calc_support();
		predict_SV();
	}
	return this->sv_type;
}

void Breakpoint::calc_support() {
	std::vector<short> SV;
	for (size_t i = 0; i < 5; i++) {
		SV.push_back(0);
	}
	for (std::map<std::string, read_str>::iterator i = positions.support.begin(); i != positions.support.end(); i++) {
		summarize_type((*i).second.SV, SV);
	}
	this->sv_type = eval_type(SV);
}
char Breakpoint::eval_type(std::vector<short> SV) {
	std::stringstream ss;

	if (SV[0] != 0) {
		ss << " DEL(";
		ss << SV[0];
		ss << ")";
	}
	if (SV[1] != 0) {
		ss << " DUP(";
		ss << SV[1];
		ss << ")";
	}
	if (SV[2] != 0) {
		ss << " INS(";
		ss << SV[2];
		ss << ")";
	}
	if (SV[3] != 0) {
		ss << " INV(";
		ss << SV[3];
		ss << ")";
	}
	if (SV[4] != 0) {
		ss << " TRA(";
		ss << SV[4];
		ss << ")";
	}
	this->sv_debug = ss.str(); //only for debug!
	//std::cout << sv_debug << std::endl;

	int maxim = 0;
	int id = 0;
	for (size_t i = 0; i < SV.size(); i++) {
		if (maxim < SV[i]) {
			maxim = SV[i];
			id = i;
		} else if (maxim == SV[i]) {
			id = -1;
		}
	}
	this->type_support = maxim;
	switch (id) {
	case 0:
		return DEL;
		break;
	case 1:
		return DUP;
		break;
	case 2:
		return INS;
		break;
	case 3:
		return INV;
		break;
	case 4:
		return TRA;
		break;
	}
	return 'n';
}

long Breakpoint::overlap(Breakpoint * tmp) {
	if (abs(tmp->get_coordinates().start.min_pos - positions.start.min_pos) < Parameter::Instance()->max_dist && abs(tmp->get_coordinates().stop.max_pos - positions.stop.max_pos) < Parameter::Instance()->max_dist) {
		return 0;
	}
//as abstraction lets try the start+stop coordinate!
	return (tmp->get_coordinates().start.min_pos - positions.start.min_pos) + (tmp->get_coordinates().stop.max_pos - positions.stop.max_pos);
}

void Breakpoint::predict_SV() {
	bool md = false;
	bool cigar = false;
	bool split = false;
	int num = 0;
	std::map<long, int> starts;
	std::map<long, int> stops;
	std::map<std::string, int> strands;

	for (std::map<std::string, read_str>::iterator i = positions.support.begin(); i != positions.support.end(); i++) {
		if ((*i).second.SV & this->sv_type) { ///check type
			if (starts.find((*i).second.coordinates.first) == starts.end()) {
				starts[(*i).second.coordinates.first] = 1;
			} else {
				starts[(*i).second.coordinates.first]++;
			}

			if (stops.find((*i).second.coordinates.second) == stops.end()) {
				stops[(*i).second.coordinates.second] = 1;
			} else {
				stops[(*i).second.coordinates.second]++;
			}

			std::string tmp = translate_strand((*i).second.strand); //std::string tmp=
			//std::cout << tmp << std::endl;
			if (strands.find(tmp) == strands.end()) {
				strands[tmp] = 1;
			} else {
				strands[tmp]++;
			}

			if ((*i).second.type == 0) {
				cigar = true;
			} else if ((*i).second.type == 1) {
				md = true;
			} else if ((*i).second.type == 2) {
				split = true;
			} else {
				std::cerr << "Type " << (*i).second.type << std::endl;
			}
			num++;
		}
	}

	int maxim = 0;
	long coord = 0;
	for (map<long, int>::iterator i = starts.begin(); i != starts.end(); i++) {
		if ((*i).second > maxim) {
			coord = (*i).first;
		}
	}
	this->positions.start.most_support = coord;

	maxim = 0;
	coord = 0;
	for (map<long, int>::iterator i = stops.begin(); i != stops.end(); i++) {
		if ((*i).second > maxim) {
			coord = (*i).first;
		}
	}
	this->positions.stop.most_support = coord;
	starts.clear();
	stops.clear();

	for (size_t i = 0; i < strands.size(); i++) {
		maxim = 0;
		std::string id;
		for (std::map<std::string, int>::iterator j = strands.begin(); j != strands.end(); j++) {
			if (maxim < (*j).second) {
				maxim = (*j).second;
				id = (*j).first;
				//std::cout << '\t' << id << std::endl;
			}
		}
		if (maxim > 0) {
			this->strand.push_back(id);
			strands[id] = 0;
		}
	}
	strands.clear();

	this->supporting_types = "";
	if (md) {
		this->supporting_types += "MD";
	}
	if (cigar) {
		if (!supporting_types.empty()) {
			this->supporting_types += ",";
		}
		this->supporting_types += "CI";
	}
	if (split) {
		if (!supporting_types.empty()) {
			this->supporting_types += ",";
		}
		this->supporting_types += "SR";
	}
	//this->strand = eval_strand(strand);
}

std::string Breakpoint::to_string(RefVector ref) {

	std::stringstream ss;
	ss << "(";
	ss << get_chr(get_coordinates().start.min_pos, ref);
	ss << ":";
	ss << calc_pos(get_coordinates().start.min_pos, ref);
	ss << "-";
	ss << get_chr(get_coordinates().stop.max_pos, ref);
	ss << ":";
	ss << calc_pos(get_coordinates().stop.max_pos, ref);
	ss << " ";
	ss << positions.support.size();
	ss << " ";
	ss << this->sv_debug;
	ss << " ";
	ss << this->get_strand(2);
	ss << "\n";
	int num = 0;
	for (std::map<std::string, read_str>::iterator i = positions.support.begin(); i != positions.support.end() && num < Parameter::Instance()->report_n_reads; i++) {
		ss << "\t";
		ss << (*i).first;
		ss << " ";
		ss << (*i).second.type;
		if ((*i).second.strand.first) {
			ss << "+";
		} else {
			ss << "-";
		}
		if ((*i).second.strand.second) {
			ss << "+";
		} else {
			ss << "-";
		}
		num++;
		ss << "\n";
	}
	ss << " ";
	return ss.str();
}

void Breakpoint::add_read(Breakpoint * point) {
	if (point != NULL) {
		this->positions.start.min_pos = min(this->positions.start.min_pos, point->get_coordinates().start.min_pos);
		this->positions.start.max_pos = max(this->positions.start.max_pos, point->get_coordinates().start.max_pos);

		this->positions.stop.min_pos = min(this->positions.stop.min_pos, point->get_coordinates().stop.min_pos);
		this->positions.stop.max_pos = max(this->positions.stop.max_pos, point->get_coordinates().stop.max_pos);

		//this->start.start = (double)(this->start.start+point->get_coordinates().start)/2;
		//this->start.stop = (double)(this->start.stop+point->get_coordinates().stop)/2;

		std::map<std::string, read_str> support = point->get_coordinates().support;
		for (std::map<std::string, read_str>::iterator i = support.begin(); i != support.end(); i++) {
			this->positions.support[(*i).first] = (*i).second;
		}
	}
}

std::string Breakpoint::get_chr(long pos, RefVector ref) {
//	std::cout << "pos: " << pos << std::endl;
	size_t id = 0;
	while (id < ref.size() && pos >= 0) {
		pos -= (long) ref[id].RefLength;
		//	std::cout << id << std::endl;
		id++;
	}

	return ref[id - 1].RefName;
}

long Breakpoint::calc_pos(long pos, RefVector ref) {
	size_t i = 0;
	pos -= ref[i].RefLength;
	while (i < ref.size() && pos >= 0) {
		i++;
		pos -= ref[i].RefLength;
	}
	return pos + ref[i].RefLength;
}

int Breakpoint::get_support() {
	return type_support;
}
char complement(char nuc) {
	switch (nuc) {
	case 'A':
		return 'T';
		break;
	case 'C':
		return 'G';
		break;
	case 'G':
		return 'C';
		break;
	case 'T':
		return 'A';
		break;
	default:
		return nuc;
		break;
	}
}
std::string Breakpoint::rev_complement(std::string seq) {
	std::string tmp;
	for (std::string::reverse_iterator i = seq.rbegin(); i != seq.rend(); i++) {
		tmp += complement((*i));
	}
	return tmp;
}

std::string Breakpoint::get_strand(int num_best) {
	//if(this->strand.empty()){
	//	predict_SV();
	//}
	if (sv_type == ' ') {
		//	std::cout<<"was not set"<<std::endl;
		calc_support();
		predict_SV();
	}
	std::string tmp = this->strand[0];
	for (int i = 1; i < num_best; i++) {
		tmp += '\t';
		if (i < (int) this->strand.size()) {
			tmp += this->strand[i];
		} else {
			tmp += ' ';
		}
	}
	return tmp;
}

