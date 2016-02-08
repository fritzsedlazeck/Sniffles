/*
 * Alignments.cpp
 *
 *  Created on: May 25, 2012
 *      Author: fritz
 */

#include "Alignment.h"

void Alignment::setRef(string sequence) {
	alignment.second = sequence;
}
void Alignment::initAlignment() {
	al = new BamAlignment();
}
void Alignment::setAlignment(BamAlignment * align) {
	al = align;
	alignment.first.clear();
	alignment.second.clear();
	is_computed = false;

	orig_length = al->QueryBases.size();
	for (size_t i = 0; i < al->QueryBases.size(); i++) {
		alignment.first += toupper(al->QueryBases[i]);
	}
	stop = this->getPosition() + this->getRefLength();
}
void Alignment::computeAlignment() {
	int pos = 0;
	for (size_t i = 0; i < al->CigarData.size(); i++) {
		if (al->CigarData[i].Type == 'I') {
			for (uint32_t t = 0; t < al->CigarData[i].Length; t++) {
				alignment.second.insert(pos, "-");
				alignment.second.erase(alignment.second.size() - 1, 1);
				pos++;
			}
		} else if (al->CigarData[i].Type == 'D') {
			for (uint32_t t = 0; t < al->CigarData[i].Length; t++) {
				alignment.first.insert(pos, "-");
				pos++;
			}
		} else if (al->CigarData[i].Type == 'S') {

			if (pos == 0) { //front side
				alignment.second.erase(((int) alignment.second.size()) - al->CigarData[i].Length, al->CigarData[i].Length);
			} else { //backside
				alignment.second.erase(pos, al->CigarData[i].Length);
			}
			alignment.first.erase(pos, al->CigarData[i].Length);

		} else if (al->CigarData[i].Type == 'M') {
			pos += al->CigarData[i].Length;
		} else if (al->CigarData[i].Type == 'H') {

		} else if (al->CigarData[i].Type == 'N') {
			alignment.second.erase(pos, al->CigarData[i].Length);
		}
	}
	for (size_t i = 0; i < alignment.first.size(); i++) {
		if (alignment.first[i] == '=') {
			alignment.first[i] = alignment.second[i];
		}
	}

	is_computed = true;

	if (alignment.first.size() != alignment.second.size()) {
		cerr << "Error alignment has different length" << endl;
		cerr << " ignoring alignment " << al->Name << endl;
		cerr << al->Position << endl;

		cerr << endl;
		cerr << "read: " << alignment.first << endl;
		cerr << endl;
		cerr << " ref: " << alignment.second << endl;
		cerr << endl;
		cerr << orig_length << endl;
		vector<CigarOp> cig = getCigar();

		for (size_t i = 0; i < cig.size(); i++) {
			cerr << cig[i].Length << cig[i].Type << " ";
		}
		cerr << endl;
		exit(0);
		return;
	}
}
int32_t Alignment::getPosition() {
	return al->Position;
}
int32_t Alignment::getRefID() {
	return al->RefID;
}
bool Alignment::getStrand() {
	return !al->IsReverseStrand();
}
vector<CigarOp> Alignment::getCigar() {
	return al->CigarData;
}
string Alignment::getQualitValues() {
	return al->Qualities;
}
size_t Alignment::get_length(std::vector<CigarOp> CigarData) {
	size_t len = 0; //orig_length;
	for (size_t i = 0; i < CigarData.size(); i++) {
		if (CigarData[i].Type == 'D' || CigarData[i].Type == 'M' || CigarData[i].Type == 'N') {

			len += CigarData[i].Length;
		}
	}
	return len;
}
size_t Alignment::getRefLength() {
	return get_length(this->al->CigarData);
}
size_t Alignment::getOrigLen() {
	return orig_length;
}
pair<string, string> Alignment::getSequence() {
	return alignment;
}
BamAlignment * Alignment::getAlignment() {
	return al;
}
string Alignment::getName() {
	return al->Name;
}
uint16_t Alignment::getMappingQual() {
	return al->MapQuality;
}
float Alignment::getIdentity() {
	if (is_computed) {
		float match = 0;
		for (size_t i = 0; i < alignment.first.size(); i++) {
			if (alignment.first[i] == alignment.second[i]) {
				match++;
			}
		}
		return match / (float) alignment.first.size();
	}
	return -1;
}
int Alignment::getAlignmentFlag() {
	return al->AlignmentFlag;
}
string Alignment::getQueryBases() {
	return al->QueryBases;
}
string Alignment::getQualities() {
	return al->Qualities;
}
string convertInt(int number) {
	stringstream ss; //create a stringstream
	ss << number; //add number to the stream
	return ss.str(); //return a string with the contents of the stream
}
string Alignment::getTagData() {
	vector<string> tags;

	uint32_t i = 0;
	if (al->GetTag("AS", i)) {
		string tmp = "AS:i:";
		tmp += convertInt(i);
		tags.push_back(tmp);

	}
	i = 0;
	if (al->GetTag("NM", i)) {
		string tmp = "NM:i:";
		tmp += convertInt(i);
		tags.push_back(tmp);
	}

	string md;
	if (al->GetTag("MD", md)) {
		string tmp = "MD:Z:";
		tmp += md;
		tags.push_back(tmp);
	}

	i = 0;
	if (al->GetTag("UQ", i)) {
		string tmp = "UQ:i:";
		tmp += convertInt(i);
		tags.push_back(tmp);
	}
	string sa;
	if (al->GetTag("SA", sa)) {
		string tmp = "SA:Z:";
		tmp += sa;
		tags.push_back(tmp);
	}

	string res;
	for (size_t i = 0; i < tags.size(); i++) {
		res += tags[i];
		if (i + 1 < tags.size()) {
			res += '\t';
		}
	}
	return res;
}
void Alignment::initSequence() {
	this->alignment.first.clear();
	this->alignment.second.clear();
}

int Alignment::get_id(RefVector ref, std::string chr) {
	for (size_t i = 0; i < ref.size(); i++) {
		if (strcmp(ref[i].RefName.c_str(), chr.c_str()) == 0) {
			return i;
		}
	}
	return -1; //should not happen!
}
void Alignment::get_coords(aln_str tmp, int & start, int &stop) {

	if (tmp.cigar[0].Type == 'S' || tmp.cigar[0].Type == 'H') {
		start = tmp.cigar[0].Length;
	} else {
		start = 0;
	}

	size_t index = tmp.cigar.size() - 1;
	if (tmp.cigar[index].Type == 'S' || tmp.cigar[index].Type == 'H') {
		stop = tmp.cigar[index].Length;
	} else {
		stop = 0;
	}

	if (!tmp.strand) {
		int h = start;
		start = stop;
		stop = h;
	}

	/*start = 0;
	 stop = 0;
	 for (size_t i = 0; i < cigar.size(); i++) {
	 if (start == 0 && (cigar[i].Type == 'H' || cigar[i].Type == 'S')) {
	 start += cigar[i].Length;
	 stop += cigar[i].Length;
	 }
	 if (cigar[i].Type == 'I' || cigar[i].Type == 'M') {
	 stop += cigar[i].Length;
	 }
	 }*/
}
void sort_insert(aln_str tmp, vector<aln_str> &entries) {

	for (vector<aln_str>::iterator i = entries.begin(); i != entries.end(); i++) {
		if (tmp.read_pos_start < (*i).read_pos_start) {
			entries.insert(i, tmp);
			return;
		}
	}
	entries.push_back(tmp);
}
vector<aln_str> Alignment::getSA(RefVector ref) {
	string sa;

	vector<aln_str> entries;
	if (al->GetTag("SA", sa) && !sa.empty()) {

		//store the main aln:
		aln_str tmp;
		tmp.RefID = this->getRefID();
		tmp.cigar = this->getCigar();
		tmp.length = this->getRefLength();
		tmp.mq = this->getMappingQual();
		tmp.pos = this->getPosition(); //+get_ref_lengths(tmp.RefID, ref);
		tmp.strand = getStrand();
		bool flag = strcmp(getName().c_str(), Parameter::Instance()->read_name.c_str()) == 0;

		get_coords(tmp, tmp.read_pos_start, tmp.read_pos_stop);
		entries.push_back(tmp);
		if (flag) {
			std::cout << "Main Read: " << tmp.read_pos_start << " REF: " << tmp.pos << " " << tmp.RefID << std::endl;
		}

		//parse the rest:
		size_t i = 0;
		int count = 0;

		std::string cigar;
		std::string chr;
		while (i < sa.size()) {
			if (count == 0 && sa[i] != ',') {
				chr += sa[i];
			}
			if (count == 1 && sa[i - 1] == ',') {
				tmp.pos = atoi(&sa[i]);
			}
			if (count == 2 && sa[i - 1] == ',') {
				tmp.strand = (bool) (sa[i] == '+');
			}
			if (count == 3 && sa[i] != ',') {
				cigar += sa[i];
			}
			if (count == 4 && sa[i - 1] == ',') {
				tmp.mq = atoi(&sa[i]);
			}
			if (count == 5 && sa[i] != ';') {
				tmp.nm = atoi(&sa[i]);
			}

			if (sa[i] == ',') {
				count++;
			}
			if (sa[i] == ';') {
				if (tmp.mq > Parameter::Instance()->min_mq && entries.size() <= Parameter::Instance()->max_splits) {
					//check this!
					tmp.cigar = translate_cigar(cigar); //translates the cigar (string) to a type vector
					get_coords(tmp, tmp.read_pos_start, tmp.read_pos_stop); //get the coords on the read.
					tmp.length = get_length(tmp.cigar); //gives the length on the reference.
					tmp.RefID = get_id(ref, chr); //translates back the chr to the id of the chr;
					//TODO: should we do something about the MD string?
					if (flag) {
						std::cout << "Read: " << tmp.read_pos_start << " REF: " << tmp.pos << " " << tmp.RefID;
						if (tmp.strand) {
							std::cout << "+" << std::endl;
						} else {
							std::cout << "+" << std::endl;
						}
					}
					//tmp.pos+=get_ref_lengths(tmp.RefID, ref);
					//insert sorted:
					includes_SV = true;
					sort_insert(tmp, entries);
				}
				chr.clear();
				cigar.clear();
				tmp.cigar.clear();
				count = 0;
				tmp.mq = 0;
			}
			i++;
		}
	}
	return entries;
}

//returns -1 if flags are not set!
double Alignment::get_scrore_ratio() {
	uint score = 0;
	uint subscore = 0;
	if (al->GetTag("AS", score) && al->GetTag("XS", subscore)) {
		if(subscore==0){
			return 40;// -1;
		}
		if (strcmp(getName().c_str(), Parameter::Instance()->read_name.c_str()) == 0) {
			std::cout << getName()<<" score: "<<(double) score / (double) subscore << std::endl;
		}
		return (double) score / (double) subscore;
	}
	return 100; //TODO: -1
}
bool Alignment::get_is_save() {
	string sa;

	double score = get_scrore_ratio();

	/*if((al->GetTag("XA", sa) && !sa.empty())){
		std::cout<<this->getName()<<"XA"<<std::endl;
	}
	if( (al->GetTag("XT", sa) && !sa.empty()) ){
		std::cout<<this->getName()<<"XT"<<std::endl;
	}*/

	return !((al->GetTag("XA", sa) && !sa.empty()) || (al->GetTag("XT", sa) && !sa.empty()) || (score == -1 || score< Parameter::Instance()->score_treshold)); //TODO: 7.5
}

std::vector<CigarOp> Alignment::translate_cigar(std::string cigar) {

	std::vector<CigarOp> new_cigar;

	size_t i = 0;
	bool first = true;
	CigarOp tmp;
	tmp.Length = -1;
	while (i < cigar.size()) {
		if (tmp.Length == -1) {
			tmp.Length = atoi(&cigar[i]);
		} else if (tmp.Length != -1 && atoi(&cigar[i]) == 0 && cigar[i] != '0') {
			tmp.Type = cigar[i];
			new_cigar.push_back(tmp);

			tmp.Length = -1;
			first = false;
		}
		i++;
	}
	return new_cigar;
}

double Alignment::get_avg_indel_length_Cigar() {
	double len=0;
	double num=0;
	for (size_t i = 0; i < al->CigarData.size(); i++) {
		if ((al->CigarData[i].Type == 'I'||al->CigarData[i].Type == 'D') && al->CigarData[i].Length > 1) {
			len+= al->CigarData[i].Length;
			num++;
		}
	}

	return len/num;
}

vector<str_event> Alignment::get_events_CIGAR() {

	size_t read_pos = 0;
	size_t pos = this->getPosition(); //orig_length;
	vector<str_event> events;
	for (size_t i = 0; i < al->CigarData.size(); i++) {
		if (al->CigarData[i].Type == 'H' || (al->CigarData[i].Type == 'S' || al->CigarData[i].Type == 'M')) {
			read_pos += al->CigarData[i].Length;
		}
		if (al->CigarData[i].Type == 'D' && al->CigarData[i].Length > Parameter::Instance()->min_cigar_event) {
			str_event ev;
			ev.read_pos = read_pos;
			ev.length = al->CigarData[i].Length; //deletion
			ev.pos = pos;
			includes_SV = true;
			events.push_back(ev);
		}
		if (al->CigarData[i].Type == 'I' && al->CigarData[i].Length > Parameter::Instance()->min_cigar_event) {
		//	std::cout<<"CIGAR: "<<al->CigarData[i].Length<<" "<<this->getName()<<std::endl;
			str_event ev;
			ev.length = al->CigarData[i].Length * -1; //insertion;
			ev.pos = pos;
			ev.read_pos = read_pos;
			includes_SV = true;
			events.push_back(ev);
			read_pos += al->CigarData[i].Length;
		}
		if (al->CigarData[i].Type == 'D' || al->CigarData[i].Type == 'M' || al->CigarData[i].Type == 'N') {
			pos += al->CigarData[i].Length;
		}
	}

	return events;
}

double Alignment::get_num_mismatches(std::string md) {
	bool deletion = false;
	bool match = false;
	vector<int> helper;
	double mis = 0;
	double len = 0;
	double maxim = 0;
	for (size_t i = 0; i < md.size(); i += 20) {
		mis = 0;
		len = 0;
		for (size_t j = 0; len < 100 && j + i < md.size(); j++) {
			if (match && atoi(&md[i + j]) == 0 && md[i + j] != '0') { //is not a number:
				if (md[i] == '^') {
					deletion = true;
				} else {
					len++;
				}
				if (!deletion) {
					//mistmatch!!
					mis++;
					match = false;
				}
			} else {
				len += atoi(&md[i + j]);

				match = true;
				deletion = false;
			}
		}

		if (strcmp(getName().c_str(), Parameter::Instance()->read_name.c_str()) == 0) {
			std::cout << (mis / len) << std::endl;
		}
		if ((mis / len) > maxim) {
			maxim = (mis / len);
		}
	}
	return maxim; // 0.03);
}
std::string Alignment::get_md() {
	std::string md;
	if (al->GetTag("MD", md)) {
		return md;
	}
	return md;
}
vector<str_event> Alignment::get_events_MD(int min_mis) {
	vector<str_event> events;
	std::string md;
	if (al->GetTag("MD", md)) {
		//TODO: remove:
		bool flag = strcmp(getName().c_str(), Parameter::Instance()->read_name.c_str()) == 0;

		if (flag) {
			std::cout << "found!" << std::endl;
		}
		//TODO think of a good threshold!
		if (get_num_mismatches(md) > Parameter::Instance()->min_num_mismatches) {
			if (flag) {
				std::cout << "is_strange!" << std::endl;
			}
			//generate a vector that holds the positions of the read
			std::vector<int> aln;
			int pos = getPosition();

			for (size_t i = 0; i < al->CigarData.size(); i++) {
				if (al->CigarData[i].Type == 'I') { //TODO check
				}
				if (al->CigarData[i].Type == 'D') {
					pos += al->CigarData[i].Length;
				}
				if (al->CigarData[i].Type == 'M') {
					for (size_t j = 0; j < al->CigarData[i].Length; j++) {
						aln.push_back(pos);
						pos++;
						//aln += "=";
					}
				}
			}
			//fill in the mismatches:
			bool deletion = false;
			bool match = false;
			double mis = 0;
			double len = 0;
			//MD:Z:106G48^C33C41T0G5^G15G0C52^T15^C5^GC16G0A48^T3^C9^G3^G17^G23^G3^C2^G40C0G6
			for (size_t i = 0; i < md.size(); i++) {
				if ((atoi(&md[i]) == 0 && md[i] != '0')) { //is not a number:
					if (md[i] == '^') {
						deletion = true;
					}
					if (!deletion) {
						//mistmatch!!
						mis++;
						aln[len] = aln[len] * -1;
						len++;
					}
					match = false;
				} else if (!match) {
					len += atoi(&md[i]);
					match = true;
					deletion = false;
				}
			}

		/*	if (flag) {
				for (size_t i = 0; i < aln.size(); i++) {
					std::cout << aln[i] << " ";
				}
				std::cout << endl;
			}
*/
			int runlength = 100;
			str_event ev;
			ev.pos = -1;
			ev.length = -1;
			ev.read_pos = 0;
			int start = 0;
			int last = 0;
			for (size_t i = 0; i < aln.size(); i += 50) { //+=runlength/2 ??
				//std::cout<<aln[i]<<";";
				int mis = 0;
				int first = 0;

				for (size_t j = 0; (j + i) < aln.size() && j < runlength; j++) {
					if (aln[(i + j)] < 0) {
						if (first == 0) {
							first = abs(aln[(i + j)]);
						}
						last = abs(aln[(i + j)]);
						mis++;
					}
				}
				if (mis > min_mis) { //TOOD ratio?
					if (ev.pos == -1) {
						start = i;
						ev.pos = first;
						ev.read_pos = ev.pos - getPosition();
					}
				} else {
					if ((start > 20 && abs((int) (i + runlength) - (int) aln.size()) > 20) && ev.pos != -1) {
						if (flag) {
							std::cout << i << " " << (i + runlength) << " " << aln.size() << std::endl;
							std::cout << ev.pos << " " << last << " " << std::endl;
						}
						includes_SV = true;
						ev.length = last - ev.pos;
						if (flag) {
							std::cout << ev.pos << " " << ev.length << std::endl;
						}
						if (ev.length > runlength) {
							events.push_back(ev);
						}
						last = 0;
						ev.pos = -1;
					} else {
						ev.pos = -1;
					}
				}
			}
		}

	}
	return events;
}
