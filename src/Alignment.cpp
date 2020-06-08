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
}

void update_aln(std::string & alignment, int & i, int pos_to_modify) {
	int ref_pos = 0;
	while (i < alignment.size() && ref_pos != pos_to_modify) {
		if (alignment[i] != '-') {
			ref_pos++;
		}
		i++;
	}
	alignment[i] = 'Y';
}

void add_event(int pos, list<differences_str>::iterator & i, list<differences_str> & events) {
	//insert sorted into vector:
	while (i != events.end() && pos > (*i).position) {
		++i;
	}
	differences_str ev;
	ev.position = pos;
	ev.type = 0; //mismatch
	events.insert(i, ev);
}

void add_event(int pos, size_t & i, vector<differences_str> & events) {
	//insert sorted into vector:
	while (i < events.size() && pos > events[i].position) {
		i++;
	}
	differences_str ev;
	ev.position = pos;
	ev.type = 0; //mismatch
	ev.readposition = -1;
	events.insert(events.begin() + i, ev);
}

/*
 * :6-ata:10+gtc:4*at:3,
 * -ata: del
 * +gtc: ins
 * *at a->t
 * :10 match bp
 * cs:Z::5+g:7+t:4+c:7-a:6*tg:15*tg:2+g:6+c:1+c:1+c:4+a:15+g:4+c:12*gt*ta:2+c:7+t:8+t:17+t:8+g:5+g:10-g:4+g*tc:5+t*ag:10+a:4-ttt
 */

vector<differences_str> Alignment::summarize_csstring(std::vector<indel_str> &dels) {

	string cs = this->get_cs();
	int pos = this->getPosition();
	int corr = 0;
	int ref_pos = 0;
	size_t pos_events = 0;
	int max_size = (this->getRefLength() * 0.9) + getPosition();
	//	comp_aln = clock();

	indel_str del;
	del.sequence = "";
	del.pos = -1;

	vector<differences_str> events;
	differences_str ev;
	std::cout << "CS: " << cs << std::endl;
	for (size_t i = 0; i < cs.size() && pos < max_size; i++) {
		ev.position = pos;
		if (cs[i] == ':') { //match (can be numbers or bp!)
			i++;
			int num = atoi(&cs[i]); //is 0 if letter!
			bool is_long = ((num == 0 && cs[i] != '0'));
			while ((i < cs.size() && cs[i] != '+') && (cs[i] != '-' && cs[i] != '+')) {
				i++;
				if (is_long) {
					num++;
				}
			}
			cout << "\tMatch: " << num << std::endl;
			pos += num;
		} else if (cs[i] == '*') { //mismatch (ref,alt) pairs
			//only every second char counts!
			add_event(pos, pos_events, events);
			pos++;
			i += 2;
			ev.type = 0; //mismatch
			cout << "\tMiss: " << 1 << std::endl;
		} else if (cs[i] == '-') { //del
			//collet del seq in dels!
			indel_str del;
			del.sequence = "";
			del.pos = pos;
			while ((i < cs.size() && cs[i] != '+') && (cs[i] != '-' && cs[i] != '+')) {
				del.sequence += cs[i];
			}
			dels.push_back(del);
			pos += del.sequence.size();
			ev.type = del.sequence.size();
			cout << "\tDEL: " << del.sequence.size() << std::endl;

		} else if (cs[i] == '+') { //ins
			int num = 0;
			while ((i < cs.size() && cs[i] != '+') && (cs[i] != '-' && cs[i] != '+')) {
				i++;
				num++;
			}
			ev.type = num * -1;
			cout << "\tINS: " << num << std::endl;
		}
		events.push_back(ev);

	}
	std::cout << "end CS: " << cs << std::endl;
	return events;
}

vector<differences_str> Alignment::summarizeAlignment(std::vector<indel_str> &dels) {
//	clock_t comp_aln = clock();
	vector<differences_str> events;
	int pos = this->getPosition();
	differences_str ev;
	bool flag = (strcmp(this->getName().c_str(), Parameter::Instance()->read_name.c_str()) == 0);
	int read_pos = 0;
	if (al->CigarData[0].Type == 'S') {
		read_pos += al->CigarData[0].Length;
	}

	int sum_mis = 0;
	int sum_events = 0;
	int sum_single = 0;
	for (size_t i = 0; i < al->CigarData.size(); i++) {
		if (al->CigarData[i].Type == 'M' || (al->CigarData[i].Type == '=' || al->CigarData[i].Type == 'X')) {
			pos += al->CigarData[i].Length;
			read_pos += al->CigarData[i].Length;
		} else if (al->CigarData[i].Type == 'D') {
			ev.position = pos;
			ev.type = al->CigarData[i].Length; //deletion
			ev.readposition = read_pos;
			ev.resolved = true;
			if (al->CigarData[i].Length > 2) {
				sum_events++;
			} else {
				sum_single++;
			}
			events.push_back(ev);
			pos += al->CigarData[i].Length;
		} else if (al->CigarData[i].Type == 'I') {
			ev.position = pos;
			ev.resolved = true;
			ev.readposition = read_pos;
			ev.type = al->CigarData[i].Length * -1; //insertion
			if (al->CigarData[i].Length > 2) {
				sum_events++;
			} else {
				sum_single++;
			}

			events.push_back(ev);
			read_pos += al->CigarData[i].Length;
		} else if (al->CigarData[i].Type == 'N') {
			pos += al->CigarData[i].Length;
			ev.resolved = true;
			read_pos += al->CigarData[i].Length;
		} else if ((al->CigarData[i].Type == 'S' || al->CigarData[i].Type == 'H') && al->CigarData[i].Length > Parameter::Instance()->huge_ins) { /// Used for reads ranging into an inser
			string sa;
			al->GetTag("SA", sa);
			uint32_t sv;
			if ((al->GetTag("SV", sv) && sa.empty()) && (!(sv & Ns_CLIPPED) && !(sv & FULLY_EXPLAINED))) { // TODO remove last )
				ev.position = pos; // - Parameter::Instance()->huge_ins;
				if (i == 0) {
					ev.readposition = 0;
				} else {
					ev.readposition = read_pos;
				}
				ev.resolved = false;
				ev.type = Parameter::Instance()->huge_ins * -1; //insertion: WE have to fix the length since we cannot estimate it!]
				events.push_back(ev);
			} else if (!al->GetTag("SV", sv) && sa.empty()) {
				if (flag) {
					cout << "HIT ALN" << endl;
				}
				ev.position = pos; // - Parameter::Instance()->huge_ins;
				if (i == 0) {
					ev.readposition = 0;
				} else {
					ev.readposition = read_pos;
				}
				ev.resolved = false;
				ev.type = Parameter::Instance()->huge_ins * -1; //insertion: WE have to fix the length since we cannot estimate it!]
				events.push_back(ev);
			}
		}
	}

	//exit(0);
	if (flag) {
		/*	std::cout << "FIRST:" << std::endl;
		 for (size_t i = 0; i < events.size(); i++) {
		 // if (abs(events[i].type) > 200) {
		 cout << events[i].position << " " << events[i].type << endl;
		 // }
		 }
		 cout << endl;*/
	}

//set ref length requ. later on:
	this->ref_len = pos - getPosition(); //TODO compare to get_length!
			//Parameter::Instance()->meassure_time(comp_aln, "\t\tCigar: ");

	string md = this->get_md();
	pos = this->getPosition();
	int corr = 0;
	bool match = false;
	bool gap = false;
	int ref_pos = 0;
	size_t pos_events = 0;
	int max_size = (this->getRefLength() * 0.9) + getPosition();
//	comp_aln = clock();
	indel_str del;
	del.sequence = "";
	del.pos = -1;
	for (size_t i = 0; i < md.size() && pos < max_size; i++) {
		if (md[i] == '^') {
			gap = true;
		}
		if ((atoi(&md[i]) == 0 && md[i] != '0')) { //is not a number
			if (!gap) { // only mismatches are stored. We should have the rest from CIGAR
				//correct for shift in position with respect to the ref:
				while (ref_pos < events.size() && pos > events[ref_pos].position) {
					if (events[ref_pos].type > 0) {
						pos += events[ref_pos].type;
					}
					ref_pos++;
				}
				sum_mis++;
				//store in sorted order:
				add_event(pos, pos_events, events);
				pos++;			//just the pos on ref!
			} else if (Parameter::Instance()->print_seq) { //can only be a deletion:
				if (del.pos == -1) {
					del.pos = pos;
				} else { //avoid first string position;
					del.sequence += md[i];
				}
			}
			match = false;
		} else if (!match) {
			match = true;
			pos += atoi(&md[i]);
			gap = false;
			if (Parameter::Instance()->print_seq && del.sequence.size() > Parameter::Instance()->min_length) {
				dels.push_back(del);
			}
			del.sequence = "";
			del.pos = -1;
		}
	}

	//if (flag) {
	//	std::cout << this->getName() << " " << (double) sum_mis << " " << (double) sum_events<<" "<<sum_single <<" "<< (double) sum_mis / (double) (sum_single+sum_mis)<< endl;
	//}
	//if (Parameter::Instance()->ccs_reads && (double) sum_mis / (double) (sum_events+sum_mis) > 0.95) {
	//	events.clear();
	//}
//	Parameter::Instance()->meassure_time(comp_aln, "\t\tMD string: ");

	return events;
}
void Alignment::computeAlignment() {
	cout << "COMP ALN!" << endl;

	clock_t comp_aln = clock();
	int to_del = 0;
	int pos = 0;

	for (size_t i = 0; i < al->CigarData.size(); i++) {
		if (al->CigarData[i].Type == 'I') {
			to_del += al->CigarData[i].Length;
			alignment.second.insert(pos, al->CigarData[i].Length, '-');
			pos += al->CigarData[i].Length;
		} else if (al->CigarData[i].Type == 'D') {

			alignment.first.insert(pos, al->CigarData[i].Length, '-');
			alignment.second.insert(pos, al->CigarData[i].Length, 'X');
			pos += al->CigarData[i].Length;
			/*for (uint32_t t = 0; t < al->CigarData[i].Length; t++) {
			 alignment.first.insert(pos, "-");
			 alignment.second.insert(pos, "X");
			 pos++;
			 }*/
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
			//nothing todo
		} else if (al->CigarData[i].Type == 'N') {
			alignment.second.erase(pos, al->CigarData[i].Length);
		}
	}
	if (to_del > 0) {
		alignment.second = alignment.second.substr(0, alignment.second.size() - to_del);
//alignment.second.erase(alignment.second.size() - to_del, to_del);
	}
	Parameter::Instance()->meassure_time(comp_aln, "\t\tCIGAR opterations ");
	comp_aln = clock();
//Apply MD string:
	string md = this->get_md();
	pos = 0;
	int corr = 0;
	bool match = false;
	int last_pos_string = 0;
	int last_pos_ref = 0;

	for (size_t i = 0; i < md.size(); i++) {
		if (atoi(&md[i]) == 0 && md[i] != '0') { //is not a number!
			if (md[i] != '^') {
				update_aln(alignment.second, last_pos_string, pos - last_pos_ref);
				last_pos_ref = pos;
				pos++;
			}
			match = false;
		} else if (!match) {
			match = true;
			pos += atoi(&md[i]);
		}
	}
	Parameter::Instance()->meassure_time(comp_aln, "\t\tMD opterations ");

	if (alignment.first.size() != alignment.second.size()) { // || strcmp(this->getName().c_str(),"IIIIII_10892000")==0) {
			//if(al->CigarData[0].Length!=100){
		cout << "Error alignment has different length" << endl;
		cout << " ignoring alignment " << al->Name << endl;
		cout << al->Position << endl;

		cout << endl;
		cout << "read: " << alignment.first << endl;
		cout << " ref: " << alignment.second << endl;
		cout << endl;
		cout << orig_length << endl;
		vector<CigarOp> cig = getCigar();
		for (size_t i = 0; i < cig.size(); i++) {
			cout << cig[i].Length << cig[i].Type << " ";
		}
		cout << endl;

		cout << this->get_md() << endl;

//	exit(0);
//	return;
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
	if (this->ref_len < 0) {
		return this->ref_len;
	}
//	return get_length(this->getCigar());

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

/*float Alignment::getIdentity() {
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
 }*/
int Alignment::getAlignmentFlag() {
	return al->AlignmentFlag;
}
string Alignment::getQueryBases() {
	if (al != NULL) {
		return al->QueryBases;
	} else {
		return "";
	}
}
void Alignment::clear_QueryBases() {
	al->QueryBases.clear();
	al->QueryBases = "";
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

int get_readlen(std::vector<CigarOp> cigar) {
	int pos = 0;
	for (size_t i = 0; i < cigar.size(); i++) {
		if (cigar[i].Type == 'I') {
			pos += cigar[i].Length;
		} else if (cigar[i].Type == 'D') {
			//pos += cigar[i].Length;
		} else if (cigar[i].Type == 'M') {
			pos += cigar[i].Length;
		}
	}
	return pos;
}
void Alignment::get_coords(aln_str tmp, int & start, int &stop) {

	size_t index = 0;
	if (!tmp.strand) {
		index = tmp.cigar.size() - 1;
	}
//	cout<<"Cigar: "<<this->getName()<<" "<<tmp.cigar.size()<<" "<<index<<endl;
	if (tmp.cigar[index].Type == 'S' || tmp.cigar[index].Type == 'H') {
		start = tmp.cigar[index].Length;
	} else {
		start = 0;
	}
	stop = get_readlen(tmp.cigar) + start;
}
void Alignment::check_entries(vector<aln_str> &entries) {

	bool flag = (strcmp(this->getName().c_str(), Parameter::Instance()->read_name.c_str()) == 0);

	if (flag) {
		std::cout << "Nested? " << std::endl;
		for (size_t i = 0; i < entries.size(); i++) {
			std::cout << entries[i].pos << "-" << entries[i].pos + entries[i].length << "(" << entries[i].read_pos_start << "-" << entries[i].read_pos_stop << ")";
			if (entries[i].strand) {
				std::cout << "+ ";
			} else {
				std::cout << "- ";
			}
			//sort_insert_ref(entries[i], new_entries);
		}
		std::cout << std::endl;

	}

	int chr = entries[0].RefID;
	bool strand = entries[0].strand;
	int strands = 1;
	int valid = 1;
	double read_gaps = 0;
	double ref_gaps = 0;

	int ref_size = 0;
	int read_size = 0;

	for (size_t i = 1; i < entries.size(); i++) {
		if (entries[i].read_pos_stop - entries[i].read_pos_start > 200) { //only consider segments > 200bp.
			ref_size = min((int) abs((entries[i - 1].pos + entries[i - 1].length) - entries[i].pos), (int) abs(entries[i - 1].pos - (entries[i].pos + entries[i].length)));

			read_size = abs(entries[i - 1].read_pos_stop - entries[i].read_pos_start);
			if (abs(ref_size - read_size) > Parameter::Instance()->min_length) {
				valid++;
			}
			if (flag) {
				cout << "Read: " << read_size << " Ref: " << ref_size << " " << this->getName() << std::endl;
			}

			if (chr != entries[i].RefID) {
				return;
			}
			if (strand != entries[i].strand) {
				strands++;
				strand = entries[i].strand;
			}
		}
	}

	if (flag) {
		std::cout << "summary: " << strands << " " << valid << " " << std::endl;
	}
	if (strands < 3 || valid < 2) { //check!
		if (flag) {
			std::cout << "Return" << std::endl;
		}
		return;
	}

	for (size_t i = 1; i < entries.size(); i++) {
		int ref_dist = 0;
		int read_dist = 0;
		if (entries[i - 1].strand) {
			ref_dist = abs((entries[i - 1].pos + entries[i - 1].length) - entries[i].pos);
			read_dist = abs(entries[i - 1].read_pos_stop - entries[i].read_pos_start);
		} else {
			ref_dist = abs((entries[i - 1].pos) - (entries[i].pos + entries[i].length));
			read_dist = abs(entries[i - 1].read_pos_stop - entries[i].read_pos_start);
		}

		if (flag) {
			std::cout << "REF DIST: " << ref_dist << " READ DIST: " << read_dist << std::endl;
		}
		if (abs(entries[i - 1].pos - entries[i].pos) < 100) { //inv dup:
			aln_str tmp;
			tmp.RefID = entries[i].RefID;
			tmp.strand = !entries[i].strand;
			tmp.mq = 60;
			tmp.length = 1;
			tmp.pos = entries[i].pos + entries[i].length;
			tmp.read_pos_start = entries[i].read_pos_stop; //fake...

			if (entries[0].strand) {
				tmp.pos = entries[i - 1].pos + entries[i - 1].length;
				tmp.read_pos_start = entries[i - 1].read_pos_stop; //fake...
				tmp.strand = !tmp.strand;
			} else {
				tmp.pos = entries[i].pos + entries[i].length;
				tmp.read_pos_start = entries[i].read_pos_stop; //fake...
			}
			tmp.read_pos_stop = tmp.read_pos_start + 1;
			entries.insert(entries.begin() + (i), tmp);
			break;
		}

		if (abs(ref_dist - read_dist) > Parameter::Instance()->min_length) { //distances between the inversion and the other split reads!
			aln_str tmp;
			tmp.RefID = entries[i].RefID;
			tmp.strand = !entries[i].strand;
			tmp.length = 1;
			tmp.mq = 60;

			//before the current element:

			tmp.pos = entries[i].pos - 1;
			tmp.read_pos_start = entries[i].read_pos_start - 1;
			tmp.read_pos_stop = tmp.read_pos_start + 1;

			//sort_insert(tmp, new_entries); //read_pos_start
			aln_str tmp2;
			tmp2 = tmp;
			//after the current element:
			tmp2.pos = entries[i].pos + entries[i].length;
			tmp2.read_pos_start = entries[i].read_pos_stop; //fake...
			tmp2.read_pos_stop = tmp2.read_pos_start + 1;
			//sort_insert(tmp, new_entries);
			if (entries[i - 1].strand) {
				entries.insert(entries.begin() + (i + 1), tmp2);
				entries.insert(entries.begin() + (i), tmp);
			} else {
				int start = tmp.read_pos_start;
				tmp.read_pos_start = tmp2.read_pos_start;
				tmp2.read_pos_start = start;
				tmp2.read_pos_stop = tmp2.read_pos_start + 1;
				tmp.read_pos_stop = tmp.read_pos_start + 1;
				entries.insert(entries.begin() + (i + 1), tmp);
				entries.insert(entries.begin() + (i), tmp2);
			}
			break;
		}

	}
	if (flag) {
		for (size_t i = 0; i < entries.size(); i++) {
			std::cout << entries[i].pos << "-" << entries[i].pos + entries[i].length << "(" << entries[i].read_pos_start << "-" << entries[i].read_pos_stop << ")";
			if (entries[i].strand) {
				std::cout << "+ ";
			} else {
				std::cout << "- ";
			}
		}
		std::cout << std::endl;
	}
}

void Alignment::sort_insert_ref(aln_str tmp, vector<aln_str> &entries) {

	for (vector<aln_str>::iterator i = entries.begin(); i != entries.end(); i++) {
		if ((tmp.pos < (*i).pos)) { //insert before
			entries.insert(i, tmp);
			return;
		}
	}
	entries.push_back(tmp);
}

void Alignment::sort_insert(aln_str tmp, vector<aln_str> &entries) {

	for (vector<aln_str>::iterator i = entries.begin(); i != entries.end(); i++) {
		if ((tmp.read_pos_start < (*i).read_pos_start)) { //insert before
			entries.insert(i, tmp);
			return;
		}
	}
	entries.push_back(tmp);
}

bool Alignment::overlapping_segments(vector<aln_str> entries) {
	bool flag = (strcmp(this->getName().c_str(), Parameter::Instance()->read_name.c_str()) == 0);
	if (flag) {
		std::cout << "HO: " << entries.size() << std::endl;
		for (size_t i = 0; i < entries.size(); i++) {
			std::cout << "Seg: " << i << " " << entries[i].pos << " " << entries[i].length << std::endl;
		}
	}
	return (entries.size() == 2 && abs(entries[0].pos - entries[1].pos) < 100);
}
vector<aln_str> Alignment::getSA(RefVector ref) {

	string sa;
	vector<aln_str> entries;
	if (al->GetTag("SA", sa) && !sa.empty()) {
//store the main aln:
		aln_str tmp;
		tmp.RefID = this->getRefID();
		tmp.cigar = this->getCigar();
		tmp.length = (long) get_length(tmp.cigar);
		tmp.mq = this->getMappingQual();
		tmp.pos = (long) this->getPosition(); //+get_ref_lengths(tmp.RefID, ref);
		tmp.strand = getStrand();
		uint32_t sv;
		al->GetTag("SV", sv);
		tmp.cross_N = ((sv & Ns_CLIPPED));
		bool flag = false; // strcmp(getName().c_str(), Parameter::Instance()->read_name.c_str()) == 0;

		get_coords(tmp, tmp.read_pos_start, tmp.read_pos_stop);
		if (flag) {
			cout << "\t read " << tmp.read_pos_start << " stop " << tmp.read_pos_stop << endl;
		}
		entries.push_back(tmp);
		if (flag) {
			std::cout << "Main Read: read start:" << tmp.read_pos_start << " REF: " << tmp.pos << " RefID: " << tmp.RefID << std::endl;
		}
		size_t i = 0;
		int count = 0;

		std::string cigar;
		std::string chr;
		bool nested = true;
		while (i < sa.size()) {
			if (count == 0 && sa[i] != ',') {
				chr += sa[i];
			}
			if (count == 1 && sa[i - 1] == ',') {
				tmp.pos = (long) atoi(&sa[i]);
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
			if (sa[i] == ';' && !cigar.empty()) {
				//TODO: maybe check how often this happens per read!
				if ((tmp.mq > Parameter::Instance()->min_mq || sv & FULLY_EXPLAINED) && entries.size() <= Parameter::Instance()->max_splits) {
					//TODO: check this!
					tmp.cigar = translate_cigar(cigar); //translates the cigar (string) to a type vector
					get_coords(tmp, tmp.read_pos_start, tmp.read_pos_stop); //get the coords on the read.
					if (flag) {
						cout << "\t read " << tmp.read_pos_start << " stop " << tmp.read_pos_stop << endl;
					}
					tmp.length = (long) get_length(tmp.cigar); //gives the length on the reference.
					tmp.RefID = get_id(ref, chr); //translates back the chr to the id of the chr;
					//TODO: should we do something about the MD string?
					if (flag) {
						std::cout << "Read: " << tmp.read_pos_start << " " << tmp.read_pos_stop << " REF: " << tmp.pos << " " << tmp.RefID;
						if (tmp.strand) {
							std::cout << "+" << std::endl;
						} else {
							std::cout << "-" << std::endl;
						}
					}
					//tmp.pos+=get_ref_lengths(tmp.RefID, ref);
					//insert sorted:
					includes_SV = true;
					sort_insert(tmp, entries);

					//al->GetTag("SV", sv);   <-get that involved

				} else if (tmp.mq < Parameter::Instance()->min_mq) {
					nested = false;
				} else {					//Ignore read due to too many splits
					entries.clear();
					return entries;
				}
				chr.clear();
				cigar.clear();
				tmp.cigar.clear();
				count = 0;
				tmp.mq = 0;
			}
			i++;
		}
		if (nested && (entries.size() > 2 || overlapping_segments(entries))) {
			check_entries(entries);
		}
		if (flag) {
			for (size_t i = 0; i < entries.size(); i++) {
				cout << "ENT: " << entries[i].pos << " " << entries[i].pos + entries[i].length << " Read: " << entries[i].read_pos_start << " " << entries[i].read_pos_stop << " ";
				if (entries[i].strand) {
					cout << "+" << endl;
				} else {
					cout << "-" << endl;
				}
			}
		}
	}
	return entries;
}

//returns -1 if flags are not set!
double Alignment::get_scrore_ratio() {
	uint score = -1;
	uint subscore = -1;
	if (al->GetTag("AS", score) && al->GetTag("XS", subscore)) {
		if (subscore == 0) {
			subscore = 1;
		}
		return (double) score / (double) subscore;
	}
	return -1;
}
bool Alignment::get_is_save() {
	string sa;

	double score = get_scrore_ratio(); //TODO should I use this again for bwa?
//	cout<<score<<endl;

	return !((al->GetTag("XA", sa) && !sa.empty()) || (al->GetTag("XT", sa) && !sa.empty())) && (score == -1 || score > Parameter::Instance()->score_treshold); //|| //TODO: 7.5
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
	double len = 0;
	double num = 0;
	for (size_t i = 0; i < al->CigarData.size(); i++) {
		if ((al->CigarData[i].Type == 'I' || al->CigarData[i].Type == 'D') && al->CigarData[i].Length > 1) {
			len += al->CigarData[i].Length;
			num++;
		}
	}

	return len / num;
}

vector<str_event> Alignment::get_events_CIGAR() {

	size_t read_pos = 0;
	size_t pos = this->getPosition(); //orig_length;
	vector<str_event> events;
	for (size_t i = 0; i < al->CigarData.size(); i++) {
		if (al->CigarData[i].Type == 'H' || (al->CigarData[i].Type == 'S' || al->CigarData[i].Type == 'M')) {
			read_pos += al->CigarData[i].Length;
		}
		if (al->CigarData[i].Type == 'D' && al->CigarData[i].Length > Parameter::Instance()->min_length) {
			str_event ev;
			ev.read_pos = read_pos;
			ev.length = al->CigarData[i].Length; //deletion
			ev.pos = pos;
			includes_SV = true;
			events.push_back(ev);
		}
		if (al->CigarData[i].Type == 'I' && al->CigarData[i].Length > Parameter::Instance()->min_length) {
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
	} else {
		std::cerr << "No MD string detected! Check bam file! Otherwise generate using e.g. samtools." << std::endl;
		cout << "MD: TEST" << this->getName() << endl;
		exit(EXIT_FAILURE);
	}
	return md;
}

std::string Alignment::get_cs() {
	std::string cs;
	if (al->GetTag("cs", cs)) {
		return cs;
	} else {
		std::cerr << "No CS string detected! Check bam file!" << std::endl;
		exit(EXIT_FAILURE);
	}
	return cs;
}

vector<str_event> Alignment::get_events_MD(int min_mis) {
	vector<str_event> events;
	/*std::string md;
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

	 }*/
	return events;
}

vector<int> Alignment::get_avg_diff(double & dist, double & avg_del, double & avg_ins) {

//computeAlignment();
//cout<<alignment.first<<endl;
//cout<<alignment.second<<endl;
	avg_del = 0;
	avg_ins = 0;
	vector<int> mis_per_window;
	std::vector<indel_str> dels;
	vector<differences_str> event_aln = summarizeAlignment(dels);
	if (event_aln.empty()) {
		dist = 0;
		return mis_per_window;
	}

	PlaneSweep_slim * plane = new PlaneSweep_slim();
	int min_tresh = 5; //reflects a 10% error rate.
//compute the profile of differences:
	double del = 0;
	double ins = 0;
	double mis = 0;
	if (event_aln.size() > 1) {
		double length = event_aln[event_aln.size() - 1].position - event_aln[0].position;
		for (size_t i = 0; i < event_aln.size(); i++) {
			if (i != 0) {
				dist += event_aln[i].position - event_aln[i - 1].position;
			}

			pair_str tmp;
			tmp.position = -1;
			if (event_aln[i].type == 0) {
				tmp = plane->add_mut(event_aln[i].position, 1, min_tresh);
			} else {
				tmp = plane->add_mut(event_aln[i].position, abs(event_aln[i].type), min_tresh);
			}
			if (tmp.position != -1) { //check that its not the prev event!
				mis_per_window.push_back(tmp.coverage); //store #mismatch per window each time it exceeds. (which might be every event position!)
			}
			if (event_aln[i].type > 0) {
				avg_del += event_aln[i].type;
			} else if (event_aln[i].type < 0) {
				avg_ins += event_aln[i].type * -1;
			}
		}
		//cout << "len: " << length << endl;
		avg_ins = avg_ins / length;
		avg_del = avg_del / length;

		dist = dist / (double) event_aln.size();
	}
	plane->finalyze();

	return mis_per_window;	//total_num /num;
}

vector<str_event> Alignment::get_events_Aln() {

	bool flag = (strcmp(this->getName().c_str(), Parameter::Instance()->read_name.c_str()) == 0);

//clock_t comp_aln = clock();
	std::vector<indel_str> dels;
	vector<differences_str> event_aln;
//	if (Parameter::Instance()->cs_string) {/
//		cout << "run cs check " << std::endl;
//		event_aln = summarize_csstring(dels);
//	} else {
	event_aln = summarizeAlignment(dels);
	if (flag) {
		cout << "\tALN events " << event_aln.size() << endl;
	}
//	}
//double time2 = Parameter::Instance()->meassure_time(comp_aln, "\tcompAln Events: ");

	vector<str_event> events;
	PlaneSweep_slim * plane = new PlaneSweep_slim();
	vector<pair_str> profile;
//	comp_aln = clock();

	bool is_N_region = false;
	uint32_t sv;
	if (al->GetTag("SV", sv) && (!(sv & Ns_CLIPPED) && !(sv & FULLY_EXPLAINED))) {
		is_N_region = true;
	}

	int noise_events = 0;
//compute the profile of differences:
	for (size_t i = 0; i < event_aln.size(); i++) {
		pair_str tmp;
		tmp.position = -1;
		if (event_aln[i].type == 0) { //substitutions.
			tmp = plane->add_mut(event_aln[i].position, 1, Parameter::Instance()->window_thresh);
		} else {
			tmp = plane->add_mut(event_aln[i].position, 1, Parameter::Instance()->window_thresh);	// abs(event_aln[i].type)
		}
		if (tmp.position != -1 && (profile.empty() || (tmp.position - profile[profile.size() - 1].position) > 100)) {	//for noisy events;
			profile.push_back(tmp);
		} else if (abs(event_aln[i].type) > Parameter::Instance()->min_length) {	//for single events like NGM-LR would produce them.
			tmp.position = event_aln[i].position;
			profile.push_back(tmp);
		}
	}

//comp_aln = clock();
	int stop = 0;
	size_t start = 0;
	for (size_t i = 0; i < profile.size() && stop < event_aln.size(); i++) {
		if (profile[i].position >= event_aln[stop].position) {
			//find the postion:
			size_t pos = 0;
			while (pos < event_aln.size() && event_aln[pos].position != profile[i].position) {
				pos++;
			}
			//run back to find the start:
			start = pos;
			int prev = event_aln[pos].position;
			start = pos;
			int prev_type = 1;
			//todo it is actually pos + type and not *type
			while (start > 0 && (prev - event_aln[start].position) < (Parameter::Instance()->max_dist_alns)) {	//13		//} * abs(event_aln[start].type) + 1)) { //TODO I  dont like 13!??
				prev = event_aln[start].position;
				prev_type = abs(event_aln[start].type);
				start--;

				if (prev_type == 0) {
					prev_type = 1;
				}
				prev += prev_type;
			}

			if (start + 1 < event_aln.size()) { //TODO do some testing!
				start++; //we are running one too far!
			}
			//run forward to identify the stop:
			prev = event_aln[pos].position;
			stop = pos;
			prev_type = 1;
			while (stop < event_aln.size() && (event_aln[stop].position - prev) < (Parameter::Instance()->max_dist_alns)) {		// * abs(event_aln[stop].type) + 1)) {
				prev = event_aln[stop].position;

				prev_type = abs(event_aln[stop].type);
				stop++;
				if (prev_type == 0) {
					prev_type = 1;
				}
				prev += prev_type;
			}
			if (stop > 0) {
				stop--;
			}

			//	cout<<start<<" events: "<<event_aln[start].type <<" pos "<<event_aln[start].readposition<<endl;
			int insert_max_pos = 0;
			int insert_max = 0;

			if (event_aln[start].type < 0) {
				insert_max_pos = event_aln[start].position;
				insert_max = abs(event_aln[start].type);
			}

			int del_max = 0;
			int del_max_pos = 0;

			if (event_aln[start].type > 0) {
				//	cout<<"HIT"<<endl;
				del_max_pos = event_aln[start].position;
				del_max = event_aln[start].type;

			}

			double insert = 0;
			double del = 0;
			double mismatch = 0;

			for (size_t k = start; k <= stop; k++) {
				if (event_aln[k].type == 0) {
					mismatch++;
				} else if (event_aln[k].type > 0) {
					del += abs(event_aln[k].type);
					if (del_max < abs(event_aln[k].type)) {
						del_max = abs(event_aln[k].type);
						del_max_pos = event_aln[k].position;
					}
				} else if (event_aln[k].type < 0) {
					insert += abs(event_aln[k].type);
					if (insert_max < abs(event_aln[k].type)) {
						insert_max = abs(event_aln[k].type);
						insert_max_pos = event_aln[k].position;
					}
				}
			}

			//	cout << "DELMAX: " << del_max << " " << Parameter::Instance()->avg_del << endl;
			str_event tmp;
			tmp.pos = event_aln[start].position;

			tmp.length = event_aln[stop].position;
			if (event_aln[stop].type > 1) {		//because of the way we summarize mutations to one location
				tmp.length += event_aln[stop].type;
			}
			tmp.length = (tmp.length - event_aln[start].position);

			tmp.type = 0;
			if (insert_max > Parameter::Instance()->min_length && insert > (del + del)) { //we have an insertion! //todo check || vs. &&
				if (is_N_region && insert_max - (insert_max * Parameter::Instance()->avg_ins) < Parameter::Instance()->min_length) {
					tmp.type = 0;
				} else {

					if (flag) {
						cout << "Is INS" << endl;
					}

					tmp.length = insert_max; //TODO not sure!
					while (start < stop && event_aln[start].readposition == -1) {
						start++;
					}
					if (flag) {
						cout << event_aln[start].readposition << " " << event_aln[start].type << endl;
					}
					tmp.read_pos = event_aln[start].readposition;
					if (Parameter::Instance()->print_seq) {
						//	if (flag) {
						//		std::cout << "Seq+:" << this->getAlignment()->QueryBases.substr(tmp.read_pos, tmp.length) << std::endl;

						//	}
						if (this->getAlignment()->QueryBases.size() < tmp.read_pos) {
							cerr << "Read sequence is shorter than expected. Please check your bam file if the read sequence is reported!" << endl;
							exit(-1);
						}
						tmp.sequence = this->getAlignment()->QueryBases.substr(tmp.read_pos, tmp.length);
					} else {
						tmp.sequence = "NA";
					}
					tmp.pos = insert_max_pos;
					tmp.type |= INS;
					tmp.is_noise = false;
				}
			} else if (del_max > Parameter::Instance()->min_length && (insert + insert) < del) { //deletion
				if (is_N_region && del_max - (del_max * Parameter::Instance()->avg_del) < Parameter::Instance()->min_length) {
					tmp.type = 0;
				} else {
					if (Parameter::Instance()->print_seq) {
						for (size_t del_pos = 0; del_pos < dels.size(); del_pos++) {
							if (abs(dels[del_pos].pos - tmp.pos) < 10) {
								tmp.sequence = dels[del_pos].sequence;
							}
						}
					} else {
						tmp.sequence = "NA";
					}
					tmp.length = del_max;
					tmp.type |= DEL;
					tmp.is_noise = false;
				}
			} else if ((mismatch + del + insert) / 2 > Parameter::Instance()->min_length) { //TODO
				if (is_N_region || ((del_max > Parameter::Instance()->min_length && insert_max > Parameter::Instance()->min_length) && (del_max / insert_max) < Parameter::Instance()->min_length)) {
					tmp.type = 0;
				} else {
					noise_events++;
					tmp.type |= DEL;
					tmp.type |= INV;
					tmp.sequence = "NA";
					tmp.is_noise = true;
				}
			}

			if (flag) {
				cout << "Read: " << " " << (double) this->getRefLength() << " events: " << event_aln.size() << " " << this->al->Name << std::endl;
				cout << "INS max " << insert_max << " del_max " << del_max << std::endl;
				cout << "INS:" << insert << " DEL: " << del << " MIS: " << mismatch << endl;
				cout << event_aln[start].position << " " << event_aln[stop].position << endl;
				cout << "store: " << tmp.pos << " " << tmp.pos + abs(tmp.length) << " " << tmp.length << endl;
				cout << tmp.sequence << endl;
				cout << tmp.type << endl;
				cout << endl;
			}

			if (tmp.type != 0) {
				if (flag) {
					cout << "ADDED" << endl;
				}
				events.push_back(tmp);
			} else {
				if (flag) {
					cout << "NOT ADDED" << endl;
				}
			}
		}
	}
//	Parameter::Instance()->meassure_time(comp_aln, "\tcompPosition: ");
	if (noise_events > 4) {
		if (flag) {
			cout << "!dumped" << endl;
		}
		events.clear();
	}

	if (flag) {
		cout << "events" << events.size() << endl;
	}
	return events;
}

