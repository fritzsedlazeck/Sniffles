/*
 * Realign.cpp
 *
 *  Created on: Aug 24, 2015
 *      Author: fsedlaze
 */

#include "Realign.h"

void Realigner::init() {
	//run through ref sequence and store a file * at the begining of each chr;
	myfile.open(Parameter::Instance()->ref_seq.c_str(), ifstream::in);
	if (!myfile.good()) {
		cout << "Fastq Parser: could not open file: " << Parameter::Instance()->ref_seq.c_str() << endl;
		exit(0);
	}

	buffer_size = 20000;
	buffer = new char[buffer_size];

	myfile.getline(buffer, buffer_size);

	long len = 0;
	while (!myfile.eof()) {
		if (buffer[0] == '>') {
			ref_str tmp;
			tmp.length = len;
			tmp.file_pos = myfile.tellg();
			meta_info.push_back(tmp);
		} else {
			for (size_t i = 0; i < buffer_size && buffer[i] != '\0'; i++) {
				len++;
			}
		}
		myfile.getline(buffer, buffer_size);
	}
	myfile.close();
}
std::string Realigner::read_new_part(long start, long stop) {
	long pos = start;
	int i = 0;
	for (; i < (int) meta_info.size() && start - meta_info[i].length > 0; i++) {
	}
	i--; //one step back
	start -= meta_info[i].length;
	stop -= meta_info[i].length;
	myfile.open(Parameter::Instance()->ref_seq.c_str(), ifstream::in);
	myfile.seekg(meta_info[i].file_pos);
	myfile.getline(buffer, buffer_size);
	string seq;
	pos = 0;
	while (!myfile.eof() && buffer[0] != '>') {
		for (size_t i = 0; i < buffer_size && buffer[i] != '\0' && buffer[i] != '\n' && buffer[0] != '>'; i++) {
			if (pos >= start && pos <= stop) {
				seq += toupper(buffer[i]);
			}
			pos++;
		}
		myfile.getline(buffer, buffer_size);
	}
	if (buffer[0] != '>') {
		for (size_t i = 0; i < buffer_size && buffer[i] != '\0' && buffer[i] != '\n' && buffer[0] != '>'; i++) {
			if (pos >= start && pos <= stop) {
				seq += toupper(buffer[i]);
			}
			pos++;
		}
	}

	myfile.close();
	return seq;
}

std::string Realigner::read_chr(short id) {

	myfile.open(Parameter::Instance()->ref_seq.c_str(), ifstream::in);
	myfile.seekg(meta_info[id].file_pos);
	myfile.getline(buffer, buffer_size);
	string seq;
	while (!myfile.eof() && buffer[0] != '>') {
		for (size_t i = 0; i < buffer_size && buffer[i] != '\0' && buffer[i] != '\n'; i++) {
			seq += toupper(buffer[i]);
		}
		myfile.getline(buffer, buffer_size);
	}
	myfile.close();
	return seq;
}

void get_coords_DEL(Breakpoint *& sv) {
	//not much todo unless size is large -> tra???
	region_ref_str tmp;

	if(sv->get_length()<5000){
		//take normal ref region;
		//tmp.start=sv->get_coordinates().start;
	}else{
		//chop:
	}
}
void get_coords_DUP(Breakpoint *& sv) {
	//duplicate the dup region next to each other?? -> define 2 overlapping regions
}
void get_coords_TRA(Breakpoint *& sv) {
	//define 2 regions on the chr
}
void get_coords_INS(Breakpoint *& sv) {
	//nothing todo
}
void get_coords_INV(Breakpoint *& sv) {
	// Estimate 2 breakpoints and set directions!
}
void Realigner::align(std::vector<Breakpoint *> sv) {
	long len = 0;
	std::string seq;
	//1 collect regions from ref
	for (size_t i = 0; i < sv.size(); i++) {//parallel

		//check if split read-> concatenate sequences;
		//else just a standard alignment.

		//split to ease job:
		if (sv[i]->get_SVtype() & DEL) {
			get_coords_DEL(sv[i]);
		} else if (sv[i]->get_SVtype() & DUP) {
			get_coords_DUP(sv[i]);
		} else if (sv[i]->get_SVtype() & TRA) {
			get_coords_TRA(sv[i]);
		} else if (sv[i]->get_SVtype() & INV) {
			get_coords_INV(sv[i]);
		} else if (sv[i]->get_SVtype() & INS) {
			get_coords_INS(sv[i]);
		}
	}
	//2: collect regions from ref:
	//for each chr run through all SV??
	for (size_t i = 0; i < this->meta_info.size(); i++) {
		long curr_length=meta_info[i].length;

		std::string ref = this->read_chr(i);
		long next_length=curr_length+ref.size();

		for(size_t j = 0; j < sv.size(); j++) { //parallel
			/*for(size_t k=0;k<sv[j]->get_ref_coord().size();k++){
				if((sv[j]->get_ref_coord()[k].start-curr_length) >0 && (sv[j]->get_ref_coord()[k].stop-next_length)<0){
					//extract region:
					int start=sv[j]->get_ref_coord()[k].start-curr_length;
					int stop=sv[j]->get_ref_coord()[k].stop-curr_length;
					sv[j]->set_ref_seq(k,ref.substr(start,stop-start));
				}
			}*/
		}
	}
	for (size_t i = 0; i < sv.size(); i++) {//parallel
		//2 send SV+ regions to alignment using OpenMP
		IAlignment * aligner = new SWCPUCor(0);
		Align align;
		align.pBuffer1 = new char[400];
		align.pBuffer2 = new char[400];
		char * refSeq = new char[400];
		char * readSeq = new char[400];
		int mode = 0;
		int cigarLength = aligner->SingleAlign(mode, Parameter::Instance()->corridor, refSeq, readSeq, align, 0);

		//3 Backtrack information

		//4 filter out SV if necessary

	}


}
