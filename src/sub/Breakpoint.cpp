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
std::string print_type(char SV1) {
	std::string type = "";
	if ((SV1 & INS)) {
		type += "INS";
	}
	if ((SV1 & TRA)) {
		type += "TRA";
	}
	if ((SV1 & INV)) {
		type += "INV";
	}
	if ((SV1 & DEL)) {
		type += "DEL";
	}
	if ((SV1 & DUP)) {
		type += "DUP";
	}
	if ((SV1 & NEST)) {
		type += "INVDUP";
	}
	if ((SV1 & NEST)) {
		type += "INS_open";
	}
	return type;
}
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
	} else if (((SV1 & NEST) && (SV2 & NEST)) || (((SV1 & NEST) && (SV2 & INV)) || ((SV1 & INV) && (SV2 & NEST)))) {
		return true;
	} else if (((SV1 & DUP) && (SV2 & INS)) || ((SV1 & INS) && (SV2 & DUP))) { //DUP and ins have often the same signal for alignments.
		return true;
	}
	//std::cout<<"S1: "<<print_type(SV1)<<" S2 "<<print_type(SV2)<<" "<<(*break2->get_coordinates().support.begin()).second.type<<std::endl;
	return (SV1 == SV2);
}
bool Breakpoint::is_same_strand(Breakpoint * tmp) {
	//todo check for same SVtype? except for INV+DEL
	if (check_SVtype(tmp, this)) {

		//if ((*tmp->get_coordinates().support.begin()).second.SV & TRA) { //only for tra since we get otherwise a problem with the cigar events
		//	//std::string readname= (*tmp->get_coordinates().support.begin()).first;
		//	return ((*tmp->get_coordinates().support.begin()).second.strand.first == (*this->get_coordinates().support.begin()).second.strand.first && (*tmp->get_coordinates().support.begin()).second.strand.second == (*this->get_coordinates().support.begin()).second.strand.second);
		//}
		return true;
	}
	return false;
}
int get_dist(Breakpoint * tmp) {
	position_str pos = tmp->get_coordinates();

	//return Parameter::Instance()->max_dist;

	if ((*tmp->get_coordinates().support.begin()).second.SV & TRA) {
		return Parameter::Instance()->max_dist; //TODO: change!
	}
	long dist = (pos.stop.max_pos - pos.start.min_pos);

	if ((*pos.support.begin()).second.length != 1) {
		dist = (*pos.support.begin()).second.length;
	}

	/*	if(dist<10){
	 std::cout<<"DIST <10 ! "<<dist <<std::endl;
	 }*/
	dist = std::max((long) Parameter::Instance()->min_length * 2, dist);

	/*if(dist <10){//TODO
	 dist=(long)Parameter::Instance()->max_dist;
	 std::cout<<"LEN SMALLER! "<<(*pos.support.begin()).first<<" "<<Parameter::Instance()->max_dist<<" "<<(pos.stop.max_pos - pos.start.min_pos)<< std::endl;
	 }*/
	return std::min((int) (dist * 4), Parameter::Instance()->max_dist);
}

long Breakpoint::overlap(Breakpoint * tmp) {
	bool flag = false;
//flag = ((*tmp->get_coordinates().support.begin()).second.SV & DEL);
	int max_dist = std::min(get_dist(tmp), get_dist(this)); // Parameter::Instance()->max_dist
	if (flag) {
		std::cout << "\t Overlap: " << max_dist << " start: " << tmp->get_coordinates().start.min_pos << " " << positions.start.min_pos << " stop :" << tmp->get_coordinates().stop.max_pos << " " << positions.stop.max_pos;
		if ((*positions.support.begin()).second.SV & DEL) {
			std::cout << " Is DEL";
		} else if ((*positions.support.begin()).second.SV & INS) {
			std::cout << " Is Ins ";
		} else if ((*positions.support.begin()).second.SV & DUP) {
			std::cout << " Is Dup ";
		} else if ((*positions.support.begin()).second.SV & INV) {
			std::cout << " Is Inv ";
		}
		std::cout << " Support: " << positions.support.size();
		std::cout << std::endl;

	}
//merging two robust calls:
	/*if (is_same_strand(tmp) && (abs(tmp->get_coordinates().start.min_pos - positions.start.min_pos) < max_dist || abs(tmp->get_coordinates().stop.max_pos - positions.stop.max_pos) < max_dist)) {
	 if (tmp->get_coordinates().stop.max_pos - tmp->get_coordinates().start.min_pos == Parameter::Instance()->huge_ins || positions.stop.max_pos - positions.start.min_pos == Parameter::Instance()->huge_ins) {
	 if (flag) {
	 cout << "\tHIT" << endl;
	 }
	 return 0;
	 }
	 }*/

	//Standard merging.
	if (is_same_strand(tmp) && (abs(tmp->get_coordinates().start.min_pos - positions.start.min_pos) < max_dist && abs(tmp->get_coordinates().stop.max_pos - positions.stop.max_pos) < max_dist)) {
		if (flag) {
			cout << "\tHIT" << endl;
		}
		return 0;
	}

	//merging "huge ins" and observed ins:
	if (is_same_strand(tmp) && (abs(tmp->get_coordinates().start.min_pos - tmp->get_coordinates().stop.max_pos) == Parameter::Instance()->huge_ins || abs(positions.start.min_pos - positions.stop.max_pos) == Parameter::Instance()->huge_ins) && (abs(tmp->get_coordinates().start.min_pos - positions.start.min_pos) < max_dist || abs(tmp->get_coordinates().stop.max_pos - positions.stop.max_pos) < max_dist)) {
		return 0;
	}

	//If there is an INVDUP:
	if ((is_NEST(tmp, this) && is_same_strand(tmp)) && (abs(tmp->get_coordinates().start.min_pos - positions.start.min_pos) < Parameter::Instance()->max_dist || abs(tmp->get_coordinates().start.min_pos - positions.stop.min_pos) < Parameter::Instance()->max_dist)) {
		return 0;
	}

//extend Split read by noisy region: //not longer needed??
	/*	if (((tmp->get_types().is_Noise || this->get_types().is_Noise) && !(tmp->get_types().is_Noise && this->get_types().is_Noise))
	 && (abs(tmp->get_coordinates().start.min_pos - positions.stop.min_pos) < max_dist / 2 || abs(tmp->get_coordinates().stop.max_pos - positions.start.max_pos) < max_dist / 2)) { //TODO maybe add SV type check!
	 if (flag) {
	 cout << "\tHIT Noise" << endl;
	 }
	 return 0;
	 }*/

//as abstraction lets try the start+stop coordinate!
	long diff = (tmp->get_coordinates().start.min_pos - positions.start.min_pos);
//if (abs(diff) < max_dist) {
//return (tmp->get_coordinates().stop.max_pos - positions.stop.max_pos);
//}
	if (diff == 0) {
		return 1;
	}
	return diff; // + (tmp->get_coordinates().stop.max_pos - positions.stop.max_pos);
}

long Breakpoint::overlap_breakpoint(long start, long stop) {
	int max_dist = get_dist(this); // Parameter::Instance()->max_dist
	if ((start < positions.start.min_pos && positions.start.min_pos < stop) || (start < positions.stop.max_pos && positions.stop.max_pos < stop)) {
		return 0;
	}
	long diff = (start - positions.start.min_pos);
	//if (abs(diff) < max_dist) {
	//return (tmp->get_coordinates().stop.max_pos - positions.stop.max_pos);
	//}
	if (diff == 0) {
		return 1;
	}
	return diff;

}
void Breakpoint::add_read(Breakpoint * point) { //point = one read support!

	if (point != NULL) {
		//merge the support:
		std::map<std::string, read_str> support = point->get_coordinates().support;
		this->set_refcount(max(point->get_refcount(), this->get_refcount())); //set ref count!
		for (std::map<std::string, read_str>::iterator i = support.begin(); i != support.end(); i++) {
			this->positions.support[(*i).first] = (*i).second;
		}
	}
}

///////////////////////////////// MERGING////////////////////////////////////////////
std::string Breakpoint::get_read_names() {
	std::string read_names;
	int num = Parameter::Instance()->report_n_reads;
	if (num == -1) {
		num = this->positions.support.size();
	}
	for (std::map<std::string, read_str>::iterator i = this->positions.support.begin(); num != 0 && i != this->positions.support.end(); i++) {
		read_names += ",";
		read_names += (*i).first;
		num--;
	}
	return read_names.substr(1);
}

std::vector<long> Breakpoint::get_read_ids() {
	std: vector<long> read_names;
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

void Breakpoint::summarize_type(char SV, ushort * array) {
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
	if (SV & NEST) {
		//	ss += "NEST; ";
		array[5]++;
	}
//return ss;
}

char Breakpoint::get_SVtype() {

	if (sv_type & NA) {
//		std::cerr << "was not set" << std::endl;
		calc_support();
		predict_SV();
	}
	return this->sv_type;
}

void Breakpoint::calc_support() {

	ushort sv[6] = { 0, 0, 0, 0, 0, 0 };

//run over all supports and check the majority type:
	for (std::map<std::string, read_str>::iterator i = positions.support.begin(); i != positions.support.end(); i++) {
		summarize_type((*i).second.SV, sv);
	}
//given the majority type get the stats:
	this->sv_type = eval_type(sv);

	if (get_SVtype() & TRA) { // we cannot make assumptions abut support yet.
		set_valid((bool) (get_support() >= 1)); // this is needed as we take each chr independently and just look at the primary alignment
	} else if (get_support() >= Parameter::Instance()->min_support) {
		predict_SV();
		set_valid((bool) (get_length() > Parameter::Instance()->min_length));
	}
}

char Breakpoint::eval_type(ushort* SV) {
	int maxim = 0;
	int id = 0;
	for (size_t i = 0; i < 6; i++) {
		if (maxim < SV[i]) {
			maxim = SV[i];
		}
	}
	this->type_support = maxim;
	if (!Parameter::Instance()->input_vcf.empty()) {
		this->type_support--; // this is needed since we introduce a pseudo element
	}
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
	if (maxim == SV[5]) {
		max_SV |= NEST;
	}
	return max_SV;
}

long get_median(std::vector<long> corrds) {

	sort(corrds.begin(), corrds.end());
	if (corrds.size() % 2 == 0) {
		return (corrds[corrds.size() / 2 - 1] + corrds[corrds.size() / 2]) / 2;
	}
	if (corrds.size() == 1) {
		return corrds[0];
	}
	return corrds[(int) corrds.size() / 2];
}


void Breakpoint::predict_SV() {
	bool aln = false;
	bool split = false;
	bool noise = false;
	int num = 0;
	std::map<long, int> starts;
	std::map<long, int> stops;
	std::map<long, int> lengths;			//ins!
	std::map<std::string, int> strands;
	std::vector<long> start2;
	std::vector<long> stops2;
	std::vector<long> lengths2;

	bool scan_reads = true;
	if (positions.support.find("input") != positions.support.end() && !Parameter::Instance()->change_coords) {
		scan_reads = false;
	}

	for (std::map<std::string, read_str>::iterator i = positions.support.begin(); i != positions.support.end() && scan_reads; i++) {

		if (((*i).second.SV & this->sv_type) && strncmp((*i).first.c_str(), "input", 5) != 0) {			// && !((*i).second.SV & INS && (*i).second.length==Parameter::Instance()->huge_ins)) { ///check type
			if ((*i).second.coordinates.first != -1) {
				if ((*i).second.length != Parameter::Instance()->huge_ins) {
					if (starts.find((*i).second.coordinates.first) == starts.end()) {
						starts[(*i).second.coordinates.first] = 1;
					} else {
						starts[(*i).second.coordinates.first]++;
					}
				}
				start2.push_back((*i).second.coordinates.first);
			}

			if ((*i).second.coordinates.second != -1) { //TODO test
				if ((*i).second.length != Parameter::Instance()->huge_ins) {
					if (stops.find((*i).second.coordinates.second) == stops.end()) {
						stops[(*i).second.coordinates.second] = 1;
					} else {
						stops[(*i).second.coordinates.second]++;
					}
				}
				stops2.push_back((*i).second.coordinates.second);
			}

			if (((*i).second.SV & INS)) { //check lenght for ins only!
				if ((*i).second.length != Parameter::Instance()->huge_ins) {
					if (lengths.find((*i).second.length) == lengths.end()) {
						lengths[(*i).second.length] = 1;
					} else {

						lengths[(*i).second.length]++;
					}
				}
				lengths2.push_back((*i).second.length);
			}

			if (!((*i).second.type == 0 && ((*i).second.SV & INV))) {
				std::string tmp = translate_strand((*i).second.strand);
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
			} else if ((*i).second.type == 2) {
				noise = true;
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
	long min_start = 0;
	long max_stop = 0;
	if (num > 0) {
		min_start = (*starts.begin()).first;
		//find the start coordinate:
		for (map<long, int>::iterator i = starts.begin(); i != starts.end(); i++) {
			if ((*i).second > maxim) {
				coord = (*i).first;
				maxim = (*i).second;
			}
			if ((*i).first < min_start) {
				min_start = (*i).first;
			}
		}
		if (maxim < 5) {
			this->positions.start.most_support = get_median(start2);	//check if "input"!
		} else {
			this->positions.start.most_support = coord;
		}
		this->indel_sequence = "";

		//find the stop coordinate:
		maxim = 0;
		coord = 0;
		mean = 0;
		counts = 0;
		for (map<long, int>::iterator i = stops.begin(); i != stops.end(); i++) {
			if ((*i).second > maxim) {
				coord = (*i).first;
				maxim = (*i).second;
			}
			if ((*i).first > max_stop) {
				max_stop = (*i).first;
			}
		}
		if (maxim < 5) {
			this->positions.stop.most_support = get_median(stops2);	// mean / counts;
		} else {
			this->positions.stop.most_support = coord;
		}

		//find the length:
		if (!(this->get_SVtype() & INS)) { //all types but Insertions:
			this->length = this->positions.stop.most_support - this->positions.start.most_support;

		} else { //compute most supported length for insertions:
			maxim = 0;
			coord = 0;
			mean = 0;
			counts = 0;
			for (map<long, int>::iterator i = lengths.begin(); i != lengths.end(); i++) {
				if ((*i).second > maxim) {
					coord = (*i).first;
					maxim = (*i).second;
				}
			}
			if (maxim == 0) { //for forcecalling
				this->length = this->positions.stop.most_support - this->positions.start.most_support;
			} else if (maxim < 3) {
				this->length = get_median(lengths2);
			} else {
				this->length = coord;
			}
			//	cout << "third len: " << this->length << endl; // problem!
		}
		starts.clear();
		stops.clear();

		//strand information:
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
	}

	if (!scan_reads) {
		//	std::cout<<"Test strand"<<std::endl;
		//if forcecalling we need to define the strands:
		std::map<std::string, read_str>::iterator it = positions.support.find("input");

		//no output!??
		/*if ((*it).second.strand.first) {
		 std::cout << "s1: true" << endl;
		 }
		 if ((*it).second.strand.second) {
		 std::cout << "s2: true" << endl;
		 }*/
		if (it != positions.support.end()) {
			this->strand.push_back(translate_strand((*it).second.strand));
		}
	}
	std::map<std::string, read_str>::iterator tmp = positions.support.begin();
	int start_prev_dist = 1000;
	int stop_prev_dist = 1000;

	while (tmp != positions.support.end()) {
		int start_dist = abs((*tmp).second.coordinates.first - this->positions.start.most_support);
		int stop_dist = abs((*tmp).second.coordinates.second - this->positions.stop.most_support);
		if (((*tmp).second.SV & this->sv_type) && ((start_dist < start_prev_dist && stop_dist < stop_prev_dist) && (strcmp((*tmp).second.sequence.c_str(), "NA") != 0 && !(*tmp).second.sequence.empty()))) {
			this->indel_sequence = (*tmp).second.sequence;

			stop_prev_dist = stop_dist; //new to test!
			start_prev_dist = start_dist;
		}
		tmp++;
	}
	//}
	this->supporting_types = "";
	if (positions.support.find("input") != positions.support.end()) {
		if (num == 0) {
			//cout << "Force called: "<< this->positions.support["input"].length << endl;
			this->positions.start.most_support = this->positions.support["input"].coordinates.first;
			this->positions.stop.most_support = this->positions.support["input"].coordinates.second;
			//this->length = this->positions.support["input"].length;
		}
		this->supporting_types += "input";
	}

	if (aln) {
		this->type.is_ALN = true;
		this->supporting_types += "AL";
	}
	if (split) {
		this->type.is_SR = true;
		if (!supporting_types.empty()) {
			this->supporting_types += ",";
		}
		this->supporting_types += "SR";
	}
	if (noise) {
		this->type.is_Noise = true;
		if (!supporting_types.empty()) {
			this->supporting_types += ",";
		}
		this->supporting_types += "NR";
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
	if (sv_type & NA) {
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
std::string Breakpoint::to_string() {
	std::stringstream ss;
	if (positions.support.size() >= 1) {
		ss << "\t\tTREE: ";
		ss << TRANS_type(this->get_SVtype());
		ss << " ";
		ss << get_coordinates().start.min_pos;
		ss << ":";
		ss << get_coordinates().stop.max_pos;
		ss << " ";
		ss << this->length;
		ss << " PE=";
		ss << positions.support.size();
		ss << " ";
		ss << get_strand(1);
	}
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
	for (std::map<std::string, read_str>::iterator i = positions.support.begin(); i != positions.support.end(); i++) {
		if (num < Parameter::Instance()->report_n_reads && Parameter::Instance()->report_n_reads != -1) {
			break;
		}
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

