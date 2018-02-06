/*
 * Alignments.h
 *
 *  Created on: May 25, 2012
 *      Author: fritz
 */

#ifndef ALIGNMENTS_H_
#define ALIGNMENTS_H_
#include <ctime>
#include <string.h>
#include <vector>
#include <math.h>
#include "api/BamAux.h"
#include "api/BamMultiReader.h"
#include "api/BamWriter.h"
#include "Paramer.h"
#include "plane-sweep/PlaneSweep_slim.h"

const unsigned char DEL = 0x01; // hex for 0000 0001
const unsigned char DUP = 0x02; // hex for 0000 0010
const unsigned char INS = 0x04; // hex for 0000 0100
const unsigned char INV = 0x08; // hex for 0000 1000
const unsigned char TRA = 0x10; // hex for 0001 0000
const unsigned char NEST =0x20; // hex for 0010 0000
const unsigned char NA = 0x80;  // hex for 1000 0000


//NGM: choped alns:
const unsigned int NOT_CLIPPED = 0x0;
const unsigned int Ns_CLIPPED = 0x1;
const unsigned int FULLY_EXPLAINED = 0x2;

using namespace BamTools;
using namespace std;
typedef unsigned short ushort;
typedef unsigned int uint;

struct differences_str{
	int position;
	int readposition;
	short type;
	bool resolved;
};
struct indel_str{
	int pos;
	std::string sequence;
};

struct str_event{
	short length;
	int pos;
	int read_pos;
	char type;
	bool is_noise;
	std::string sequence; //just for indels;
};
struct aln_str{
	int RefID;
	long pos;
	bool strand;
	std::vector<CigarOp> cigar;
	ushort mq;
	ushort nm;
	long length;
	int read_pos_start;
	int read_pos_stop;
	bool cross_N;
};

class Alignment {

private:
	int ref_len;
	BamAlignment * al;
	bool includes_SV;
	pair<string,string> alignment;
	bool is_computed;
	int32_t orig_length;
	int stop;
	 std::vector<CigarOp> translate_cigar(std::string cigar);
	 size_t get_length(std::vector<CigarOp> CigarData);
	 int get_id(RefVector ref, std::string chr);
	 vector<differences_str> summarizeAlignment(std::vector<indel_str> &dels);
	 void sort_insert(aln_str tmp, vector<aln_str> &entries);

	 void sort_insert_ref(aln_str tmp, vector<aln_str> &entries);
	 void check_entries(vector<aln_str> &entries);
	 bool overlapping_segments(vector<aln_str> entries);
public:
	Alignment(){
		al=NULL;
		ref_len=0;
	    stop=0;
		orig_length=0;
		al=NULL;
		is_computed=false;
		includes_SV=false;
	}
	~Alignment(){
		alignment.first.clear();
		alignment.second.clear();
		delete al;
	}
	void setAlignment(BamAlignment * al);
	void setRef(string sequence);
	void computeAlignment();
	void clear_QueryBases();
	pair<string,string> getSequence();
	int32_t getPosition();
	int32_t getRefID();
	bool getStrand();
	uint16_t getMappingQual();
	string getName();
	vector<CigarOp> getCigar();
	string getQualitValues();
	size_t getRefLength();
	size_t getOrigLen();
	BamAlignment * getAlignment();
	//float getIdentity();
	void initAlignment();
	int getAlignmentFlag();
	string getQueryBases();
	string getQualities();
	string getTagData();
	vector<aln_str> getSA(RefVector ref);
	void initSequence();
	vector<str_event> get_events_Aln();
	int get_stop(){
		return stop;
	}
	vector<str_event> get_events_CIGAR();
	vector<str_event> get_events_MD(int min_length);
	 void get_coords(aln_str tmp,int & start,int &stop);
	 bool supports_SV(){
		 return this->includes_SV;
	 }
	 void set_supports_SV(bool flag){
		 this->includes_SV=flag;
	 }
	 bool get_is_save();
	 double get_num_mismatches(std::string md);
	 double get_scrore_ratio();
	 std::string get_md();
	 double get_avg_indel_length_Cigar();
	 vector<int> get_avg_diff(double & dist,double & avg_del, double & avg_len);

};


#endif /* ALIGNMENTS_H_ */
