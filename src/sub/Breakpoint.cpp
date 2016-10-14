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
#include "../print/IPrinter.h"
#include "Breakpoint.h"

///////////////////////////////// MERGING////////////////////////////////////////////
bool Breakpoint::check_SVtype(Breakpoint * break1, Breakpoint * break2) { //todo check that!
	char SV1 = (*break1->get_coordinates().support.begin()).second.SV;
	char SV2 = (*break2->get_coordinates().support.begin()).second.SV;
	//we have to check it that way,because we can have multiple types!
	if ((SV1 & INS) && (SV2 & INS)) {
		return true;
	} else if ((SV1 & INV) && (SV2 & INV)) {
		return true;
	} else if ((SV1 & TRA) && (SV2 & TRA)) {
		return true;
	} else if ((SV1 & DUP) && (SV2 & DUP)) {
		return true;
	} else if ((SV1 & DEL) && (SV2 & DEL)) {
		return true;
	} /*else if (((SV1 & DUP) && (SV2 & INS))  || ((SV1 & INS) && (SV2 & DUP))) {
		return true;
	}*/
	return false;
}
bool Breakpoint::is_same_strand(Breakpoint * tmp) {
	//todo check for same SVtype? except for INV+DEL
	if (check_SVtype(tmp, this)) {
		if ((*tmp->get_coordinates().support.begin()).second.SV & TRA) { //only for tra since we get otherwise a problem with the cigar events
			return ((*tmp->get_coordinates().support.begin()).second.strand.first == (*this->get_coordinates().support.begin()).second.strand.first && (*tmp->get_coordinates().support.begin()).second.strand.second == (*this->get_coordinates().support.begin()).second.strand.second);
		}
		return true;
	}
	return false;
}
int get_dist(Breakpoint * tmp) {
	position_str pos = tmp->get_coordinates();
	//return Parameter::Instance()->max_dist;
	if ((*tmp->get_coordinates().support.begin()).second.SV & TRA || (*tmp->get_coordinates().support.begin()).second.SV & INS) {
		return Parameter::Instance()->max_dist; //TODO: change!
	} else {
		int max_val = Parameter::Instance()->max_dist;// * 0.5; //TOOD maybe 0.5?
		return std::max(std::min((int) (pos.stop.max_pos - pos.start.min_pos) / 80, Parameter::Instance()->max_dist), max_val); //20% of the lenght of the SV
	}
}

long Breakpoint::overlap(Breakpoint * tmp) {

	int max_dist =  get_dist(tmp);
	//std::cout<<"\tOverlap: "<<max_dist<< " start: "<<abs(tmp->get_coordinates().start.min_pos - positions.stop.min_pos) << " stop :" <<abs(tmp->get_coordinates().stop.max_pos - positions.start.max_pos);

	//check type. ALN could be part of the SPlit read event! Not merge two split reads!
	if (is_same_strand(tmp) && (abs(tmp->get_coordinates().start.min_pos - positions.start.min_pos) < max_dist && abs(tmp->get_coordinates().stop.max_pos - positions.stop.max_pos) < max_dist)) {
//		std::cout << "\tfound hit!" << std::endl;
		return 0;
	}

	//TODO:: if one of the two SVs is based on ALN it might be good to join them: to merge noisy flanking regions
	if (((tmp->get_types().is_ALN || this->get_types().is_ALN) && !(tmp->get_types().is_ALN && this->get_types().is_ALN)) && (abs(tmp->get_coordinates().start.min_pos - positions.stop.min_pos) < max_dist/2 || abs(tmp->get_coordinates().stop.max_pos - positions.start.max_pos) < max_dist/2)) { //TODO maybe add SV type check!
//		std::cout << "\t hit!" << std::endl;
		return 0;
	}

	//std::cout << "no hit? "<< std::endl;
//as abstraction lets try the start+stop coordinate!
	long diff=(tmp->get_coordinates().start.min_pos - positions.start.min_pos);
	if(abs(diff)<max_dist){
		return (tmp->get_coordinates().stop.max_pos - positions.stop.max_pos);
	}
	return diff;// + (tmp->get_coordinates().stop.max_pos - positions.stop.max_pos);
}

void Breakpoint::add_read(Breakpoint * point) {
	if (point != NULL) {
		if ((*point->get_coordinates().support.begin()).second.type == 0) {
			this->type.is_ALN = true;
		}
		if ((*point->get_coordinates().support.begin()).second.type == 1) {
			this->type.is_SR = true;
		}

		if ((point->get_types().is_ALN || this->get_types().is_ALN) && !(point->get_types().is_ALN && this->get_types().is_ALN)) { //
			this->type.is_ALN = false; //TODO trick ;)
			//THis is to merger noisy flanking regions:
			if (abs(point->get_coordinates().start.min_pos - positions.stop.min_pos) < Parameter::Instance()->max_dist / 2) { //left of the current position
			//start stays; stop changes
				this->positions.stop.min_pos = positions.stop.min_pos;
				this->positions.stop.max_pos = positions.stop.max_pos;
				//Now we make the stop positions invalid:
				std::map<std::string, read_str> support = this->positions.support;
				for (std::map<std::string, read_str>::iterator i = support.begin(); i != support.end(); i++) {
					(*i).second.coordinates.second = -1;
				}
				support = point->get_coordinates().support;
				for (std::map<std::string, read_str>::iterator i = support.begin(); i != support.end(); i++) {
					(*i).second.coordinates.second = -1;
				}
			} else if (abs(point->get_coordinates().stop.max_pos - positions.start.max_pos) < Parameter::Instance()->max_dist / 2) { //right of the current position
			//	cout << "right " << endl;
			//stop stays; start changes:
				this->positions.start.min_pos = positions.start.min_pos;
				this->positions.start.max_pos = positions.start.max_pos;

				//Now we make the start positions invalid:
				std::map<std::string, read_str> support = this->positions.support;
				for (std::map<std::string, read_str>::iterator i = support.begin(); i != support.end(); i++) {
					(*i).second.coordinates.first = -1;
				}
				support = point->get_coordinates().support;
				for (std::map<std::string, read_str>::iterator i = support.begin(); i != support.end(); i++) {
					(*i).second.coordinates.first = -1;
				}
			}
		} else {
			this->positions.start.min_pos = min(this->positions.start.min_pos, point->get_coordinates().start.min_pos);
			this->positions.start.max_pos = max(this->positions.start.max_pos, point->get_coordinates().start.max_pos);

			if (point->get_coordinates().support.begin()->second.SV & INS) {
				this->positions.stop.min_pos = this->positions.start.max_pos;
				this->positions.stop.max_pos = this->positions.start.max_pos;
			} else {
				this->positions.stop.min_pos = min(this->positions.stop.min_pos, point->get_coordinates().stop.min_pos);
				this->positions.stop.max_pos = max(this->positions.stop.max_pos, point->get_coordinates().stop.max_pos);
			}
		}
		bool flag = false;
		std::map<std::string, read_str> support = point->get_coordinates().support;
		for (std::map<std::string, read_str>::iterator i = support.begin(); i != support.end(); i++) {
			if (this->positions.support[(*i).first].SV & INS) {
				flag = true;
			}
			this->positions.support[(*i).first] = (*i).second;
		}
		if (flag) {
			this->length = (this->length + point->get_length()) / 2;
		}
	}
}

///////////////////////////////// MERGING////////////////////////////////////////////
std::vector<std::string> Breakpoint::get_read_names(int maxim) {
	std: vector<std::string> read_names;
	std::map<std::string, read_str> support = this->positions.support;
	int num = 0;
	for (std::map<std::string, read_str>::iterator i = support.begin(); (num < maxim || maxim == -1) && i != support.end(); i++) {
		read_names.push_back((*i).first);
		num++;
	}
	return read_names;
}

std::vector<int> Breakpoint::get_read_ids() {
	std: vector<int> read_names;
	std::map<std::string, read_str> support = this->positions.support;
	int num = 0;
	for (std::map<std::string, read_str>::iterator i = support.begin(); i != support.end(); i++) {
		read_names.push_back((*i).second.id);
		num++;
	}
	return read_names;
}

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

char Breakpoint::get_SVtype() {

	if (sv_type == ' ') {
		std::cerr << "was not set" << std::endl;
		calc_support();
		predict_SV();
	}
	return this->sv_type;
}

void Breakpoint::calc_support() {
	//if (positions.support.size() > Parameter::Instance()->min_support) {
	std::vector<short> SV;
	for (size_t i = 0; i < 5; i++) {
		SV.push_back(0);
	}
	//run over all supports and check the majority type:
	for (std::map<std::string, read_str>::iterator i = positions.support.begin(); i != positions.support.end(); i++) {
		summarize_type((*i).second.SV, SV);
	}
	//given the majority type get the stats:
	this->sv_type = eval_type(SV);
	//}
}

char Breakpoint::eval_type(std::vector<short> SV) {
	/*	std::stringstream ss;
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
	 std::cout << sv_debug << std::endl;
	 */
	int maxim = 0;
	int id = 0;
	for (size_t i = 0; i < SV.size(); i++) {
		if (maxim < SV[i]) {
			maxim = SV[i];
		}
	}
	this->type_support = maxim;
	char max_SV = 0;
	if (maxim == SV[0]) {
		max_SV |= DEL;
	}
	if (maxim == SV[1]) {
		max_SV |= DUP;
	}
	if (maxim == SV[2]) {
		max_SV |= INS;
	}
	if (maxim == SV[3]) {
		max_SV |= INV;
	}
	if (maxim == SV[4]) {
		max_SV |= TRA;
	}

	return max_SV;
}
/*void trans_strand(pair<bool, bool> tmp) {
 if (tmp.first) {
 std::cout << "+";
 } else {
 std::cout << "-";
 }

 if (tmp.second) {
 std::cout << "+";
 } else {
 std::cout << "-";
 }
 std::cout << std::endl;
 }*/

/*
 char direction(bool dir) {
 if (dir) {
 return '+';
 }
 return '-';
 }*/

void Breakpoint::predict_SV() {
	bool aln = false;
	bool split = false;
	int num = 0;
	std::map<long, int> starts;
	std::map<long, int> stops;
	std::map<std::string, int> strands;

	for (std::map<std::string, read_str>::iterator i = positions.support.begin(); i != positions.support.end(); i++) {
		if ((*i).second.SV & this->sv_type) { ///check type
			//cout << "Hit" << endl;
			if ((*i).second.coordinates.first != -1) {
				if (starts.find((*i).second.coordinates.first) == starts.end()) {
					starts[(*i).second.coordinates.first] = 1;
				} else {
					starts[(*i).second.coordinates.first]++;
				}
			}
			if ((*i).second.coordinates.second != -1) { //TODO test
				if (stops.find((*i).second.coordinates.second) == stops.end()) {
					stops[(*i).second.coordinates.second] = 1;
				} else {
					stops[(*i).second.coordinates.second]++;
				}
			}
			if (!((*i).second.type == 0 && ((*i).second.SV & INV))) {

				std::string tmp = translate_strand((*i).second.strand);

				//std::cout << tmp << std::endl;
				if (strands.find(tmp) == strands.end()) {
					strands[tmp] = 1;
				} else {
					strands[tmp]++;
				}
			}
			if ((*i).second.type == 0) {
				aln = true;
			} else if ((*i).second.type == 1) {
				split = true;
			} else {
				std::cerr << "Type " << (*i).second.type << std::endl;
			}
			num++;
		}
	}

	long mean = 0;
	long counts = 0;
	int maxim = 0;
	long coord = 0;
	for (map<long, int>::iterator i = starts.begin(); i != starts.end(); i++) {
		//cout<<"start:\t"<<(*i).first<<" "<<(*i).second<<endl;
		if ((*i).second > maxim) {
			coord = (*i).first;
			maxim = (*i).second;
		}
		mean += ((*i).first * (*i).second);
		counts += (*i).second;
	}
	if (maxim < 5) {
		this->positions.start.most_support = mean / counts;
	} else {
		this->positions.start.most_support = coord;
	}

	maxim = 0;
	coord = 0;
	mean = 0;
	counts = 0;
	for (map<long, int>::iterator i = stops.begin(); i != stops.end(); i++) {
		if ((*i).second > maxim) {
			coord = (*i).first;
			maxim = (*i).second;
		}
		mean += ((*i).first * (*i).second);
		counts += (*i).second;
	}
	if (maxim < 5) {
		this->positions.stop.most_support = mean / counts;
	} else {
		this->positions.stop.most_support = coord;
	}
	if (!(this->get_SVtype() & INS)) {
		this->length = this->positions.stop.most_support - this->positions.start.most_support;
	}
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
	if (aln) {
		this->supporting_types += "AL";
	}
	if (split) {
		if (!supporting_types.empty()) {
			this->supporting_types += ",";
		}
		this->supporting_types += "SR";
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
	if (this->strand.empty()) {
		return "UNDEF";
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
#include "Detect_Breakpoints.h"
std::string Breakpoint::to_string(){
	std::stringstream ss;
	ss << "\t\tTREE: ";
	ss << TRANS_type(this->get_SVtype());
	ss << " ";
	ss << get_coordinates().start.min_pos;
	ss << ":";
	ss << get_coordinates().stop.max_pos;
	ss << " ";
	ss << positions.support.size();
	ss << " ";
	ss << get_strand(2);
	return ss.str();
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

