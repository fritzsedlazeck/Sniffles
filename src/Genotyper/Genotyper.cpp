/*
 / * Genotyper.cpp
 *
 *  Created on: Mar 28, 2016
 *      Author: fsedlaze
 */

#include "Genotyper.h"

long get_ref_lengths3(int id, RefVector ref) {
	long length = 0;

	for (size_t i = 0; i < (size_t) id && i < ref.size(); i++) {
		length += (long) ref[i].RefLength + (long) Parameter::Instance()->max_dist;
	}
	return length;
}

void update_entries(std::vector<str_breakpoint_slim> &entries, int read_start, int read_stop, size_t & current_pos, int wobble, std::string rname, bool strand) { //TODO room for optimization!

	bool flag = false;//(strcmp(rname.c_str(), "dba0f5fd-bd4b-417f-8a57-b45de4ca017a") == 0);
	if (flag) {
		cout << "FOUND" << endl;
	}

	if (entries.empty() || read_stop + wobble < entries[0].pos) {
		if (flag) {
			cout << "Quit Return" << endl;
		}
		return;
	}

	for (size_t i = current_pos; i < entries.size(); i++) { //runs over the SV:
		if (entries[i].pos < read_start - wobble) {
			current_pos = i;
		}
		//that if is the issue.
		//if the read is rangin inside on the left:

		if ((read_start - wobble < entries[i].pos && read_stop + wobble > entries[i].pos)) {		// && (abs(entries[i].pos - read_start) > wobble/2 && abs(entries[i].pos - read_stop) > wobble/2)) {	//TODO not sure if I cannot combine these two.
			entries[i].rnames[rname] = strand; //TOOD maybe just a normal vector!
			//if (entries[i].pos == 2) {
			//	cout << "\tHIT: " << entries[i].pos << " " << entries[i].rnames.size() << endl;
			//}
		}
		if (entries[i].pos > read_stop + wobble) {
			if (flag) {
				cout << "Return" << endl;
			}
			break;
		}
	}
}

void update_coverage(std::map<std::string, std::vector<str_breakpoint_slim> > & entries) {

	//run across BAM file to parse intervals and compute coverage
	//tree.overlaps <- strand 1 = true , else = false??
	//do I need to think about filtering etc?

	cout << "\tReopening Bam file for parsing coverage " << endl;
	RefVector ref;
	BamParser * mapped_file = 0;
	if (Parameter::Instance()->bam_files[0].find("bam") != string::npos) {
		mapped_file = new BamParser(Parameter::Instance()->bam_files[0]);
		ref = mapped_file->get_refInfo();
	} else {
		cerr << "File Format not recognized. File must be a sorted .bam file!" << endl;
		exit(EXIT_FAILURE);
	}

	Alignment * tmp_aln = mapped_file->parseRead((uint16_t) Parameter::Instance()->min_mq);
	long num_reads = 0;

	int current_RefID = tmp_aln->getRefID();
	//std::vector<str_breakpoint_slim> tmp = entries[ref[tmp_aln->getRefID()].RefName];

	size_t current_pos = 0;	// index of entries vector;

	while (!tmp_aln->getQueryBases().empty()) {
		if (current_RefID != tmp_aln->getRefID()) {
			current_pos = 0;
			//	entries[ref[current_RefID].RefName] = tmp;
			std::cout << "\t\tSwitch Chr " << ref[tmp_aln->getRefID()].RefName << std::endl;
			current_RefID = tmp_aln->getRefID();
			//	tmp = entries[ref[current_RefID].RefName];
		}

		int start = (int) tmp_aln->getPosition();
		int stop = (int) start + tmp_aln->getRefLength();
		//cout<<"RNAME: "<<tmp_aln->getName() <<endl;
		update_entries(entries[ref[current_RefID].RefName], start, stop, current_pos, 5, tmp_aln->getName(), tmp_aln->getStrand());

		mapped_file->parseReadFast((uint16_t) Parameter::Instance()->min_mq, tmp_aln);
	}

	std::cout << "\tFinalizing  .." << std::endl;
	delete mapped_file;
}

str_breakpoint_slim init_breakpoint() {
	str_breakpoint_slim tmp;
	tmp.chr = "";
//	tmp.cov = 0;
//	tmp.is_start = false;
	tmp.pos = 0;
	return tmp;
}

void Genotyper::get_breakpoint_vcf(string buffer, str_breakpoint_slim & start, str_breakpoint_slim & stop) {
	size_t i = 0;
	int count = 0;

	start = init_breakpoint();
	stop = init_breakpoint();

	while (buffer[i] != '\0' && buffer[i] != '\n') {
		if (count == 0 && buffer[i] != '\t') {
			start.chr += buffer[i];
		}
		if (count == 1 && buffer[i - 1] == '\t') {
			start.pos = atoi(&buffer[i]);

		}
		if (stop.pos == -1 && (count == 4 && (buffer[i - 1] == '[' || buffer[i - 1] == ']'))) { //FOR TRA
			parse_pos(&buffer[i - 1], stop.pos, stop.chr);
		}

		if (count > 6 && strncmp(";CHR2=", &buffer[i], 6) == 0) {
			i += 6;
			while (buffer[i] != ';') {
				stop.chr += buffer[i];
				i++;
			}
		}
		if (count > 6 && strncmp(";END=", &buffer[i], 5) == 0) {
			stop.pos = atoi(&buffer[i + 5]); //stores right most breakpoint
			break;
		}

		if (buffer[i] == '\t') {
			count++;
		}
		i++;
	}
}

variant_str Genotyper::get_breakpoint_bedpe(string buffer, str_breakpoint_slim & start, str_breakpoint_slim & stop) {
	size_t i = 0;
	int count = 0;
	std::string chr;
	variant_str tmp;
	start = init_breakpoint();
	stop = init_breakpoint();

	while (buffer[i] != '\0' && buffer[i] != '\n') {
		if (count == 12 && buffer[i] != '\t') {
			start.chr += buffer[i];
		}

		if (count == 13 && buffer[i - 1] == '\t') {
			start.pos = atoi(&buffer[i]);
		}
		if (count == 14 && buffer[i] != '\t') {
			stop.chr += buffer[i];
		}
		if (count == 15 && buffer[i - 1] == '\t') {
			stop.pos = atoi(&buffer[i]);
			break;
		}
		if (buffer[i] == '\t') {
			count++;
		}
		i++;
	}
	return tmp;
}

void insert_sorted_entry(std::map<std::string, std::vector<str_breakpoint_slim> > & entries, str_breakpoint_slim pos) {
	size_t i = 0;
	while (i < entries[pos.chr].size() && entries[pos.chr][i].pos < pos.pos) {
		i++;
	}
	if (entries.size() == i) {
		entries[pos.chr].push_back(pos);
	} else {
		std::vector<str_breakpoint_slim>::iterator it = entries[pos.chr].begin();
		entries[pos.chr].insert(it + i, pos);
	}
}
void Genotyper::read_SVs(std::map<std::string, std::vector<str_breakpoint_slim> > & entries) {

	//cout<<"READ SV "<<endl;
	std::ifstream myfile;
	bool is_vcf = !Parameter::Instance()->output_vcf.empty();

	if (!Parameter::Instance()->output_vcf.empty()) {
		myfile.open(Parameter::Instance()->output_vcf.c_str(), std::ifstream::in);
	} else if (!Parameter::Instance()->output_bedpe.empty()) {
		myfile.open(Parameter::Instance()->output_bedpe.c_str(), std::ifstream::in);
	}

	if (!myfile.good()) {
		std::cout << "SVParse: could not open file: " << std::endl;
		exit(EXIT_FAILURE);
	}
	string buffer;

	getline(myfile, buffer);

	int num_sv = 0;
	while (!myfile.eof()) {
		if (buffer[0] != '#') {

			str_breakpoint_slim start;
			str_breakpoint_slim stop;
			if (is_vcf) {
				get_breakpoint_vcf(buffer, start, stop);
			} else {
				//get_breakpoint_bedpe(buffer); //TODO
			}

			//we want a sorted inclusion!
			insert_sorted_entry(entries, start);
			insert_sorted_entry(entries, stop);

			num_sv++;
			if (num_sv % 5000 == 0) {
				cout << "\t\tRead in SV: " << num_sv << endl;
			}
		}
		getline(myfile, buffer);
	}
	myfile.close();
}

void Genotyper::update_svs_output(std::map<std::string, std::vector<str_breakpoint_slim> > entries) {
	std::ifstream myfile;
	bool is_vcf = !Parameter::Instance()->output_vcf.empty();

	string file_name;
	if (!Parameter::Instance()->output_vcf.empty()) {
		file_name = Parameter::Instance()->output_vcf;
		myfile.open(Parameter::Instance()->output_vcf.c_str(), std::ifstream::in);
	} else if (!Parameter::Instance()->output_bedpe.empty()) {
		file_name = Parameter::Instance()->output_bedpe;
		myfile.open(Parameter::Instance()->output_bedpe.c_str(), std::ifstream::in);
	}

	FILE*file = fopen(Parameter::Instance()->tmp_file.c_str(), "w"); //

	if (!myfile.good()) {
		std::cout << "SVParse: could not open file: " << std::endl;
		exit(EXIT_FAILURE);
	}

	string buffer;
	getline(myfile, buffer);
//parse SVs breakpoints in file

	while (!myfile.eof()) {
		if (buffer[0] != '#') {

			std::string to_print;
			//GO PER ENTRY IN VCF FILE:
			variant_str tmp;
			str_breakpoint_slim start;
			str_breakpoint_slim stop;
			if (is_vcf) {
				get_breakpoint_vcf(buffer, start, stop);
			} else {
				//get_breakpoint_bedpe(buffer); //TODO
			}
			map<std::string, bool> tmp_names;
			pair<int, int> strands;
			strands.first = 0; //+
			strands.second = 0; //-
			bool flag = false;//(start.pos == 78257723);
			bool control = false;
			int counts_plus = 0;
			int counts_minus = 0;

			for (size_t i = 0; i < entries[start.chr].size(); i++) {
				if (start.pos == entries[start.chr][i].pos) { //found match with next in line VCF entry
					control = true;
					counts_plus = entries[start.chr][i].rnames.size();
					for (map<std::string, bool>::iterator t = entries[start.chr][i].rnames.begin(); t != entries[start.chr][i].rnames.end(); t++) {
						tmp_names[(*t).first] = (*t).second;
					}
					break;
				}
			}

			for (size_t i = 0; i < entries[stop.chr].size(); i++) {
				if (stop.pos == entries[stop.chr][i].pos) {
					control = true;
					counts_minus = entries[stop.chr][i].rnames.size();
					//final_ref.second = entries[stop.chr][i].cov;
					for (map<std::string, bool>::iterator t = entries[stop.chr][i].rnames.begin(); t != entries[stop.chr][i].rnames.end(); t++) {
						tmp_names[(*t).first] = (*t).second;
					}
					break;
				}
			}

			if (!control) {
				std::cerr << "Error in GT: Tree node not found. Exiting." << std::endl;
				exit(EXIT_FAILURE);
			}
			for (map<std::string, bool>::iterator t = tmp_names.begin(); t != tmp_names.end(); t++) {
				if (flag) {
					cout << (*t).first << " ,";
				}
				if ((*t).second) {
					strands.first++;
				} else {
					strands.second++;
				}
			}

			if (flag) {
				cout << endl;
				if (tmp_names["9cff7725-df41-441b-8bec-a7e378fd4864"]) {
					cout << " + " << endl;
				} else if (!tmp_names["9cff7725-df41-441b-8bec-a7e378fd4864"]) {
					cout << " + " << endl;
				} else {
					cout << " NA " << endl;
				}
				cout << "FOUND coverage: " << start.pos << " total coverage +: " << strands.first << " -: " << strands.second << " tot+: " << counts_minus + counts_plus << endl;
			}
			if (is_vcf) {
				to_print = mod_breakpoint_vcf(buffer, strands.first, strands.second); //(int) tmp_names.size());
			} else {
				to_print = mod_breakpoint_bedpe(buffer, strands.first, strands.second); //(int) tmp_names.size());
			}
			if (!to_print.empty()) {
				fprintf(file, "%s", to_print.c_str());
				fprintf(file, "%c", '\n');
			}
		} else {
			fprintf(file, "%s", buffer.c_str());
			fprintf(file, "%c", '\n');
		}

		getline(myfile, buffer);
	}
	myfile.close();
	fclose(file);

	string move = "mv ";
	move += Parameter::Instance()->tmp_file;
	move += " ";
	move += file_name;
	system(move.c_str());
}

void Genotyper::update_SVs2() {
//parse SVs not just breakpoints and update with the coverage info
	//per CHR:
	//	Reconstruct the tree: similar to the IVCF parser?
	//	Compute the coverage per chr (<can be in parallele to 1)
	//	Intersect result from 1 + 2 -> print out.
	//	Done;
	RefVector ref;
	BamParser * mapped_file = 0;
	if (Parameter::Instance()->bam_files[0].find("bam") != string::npos) {
		mapped_file = new BamParser(Parameter::Instance()->bam_files[0]);
		ref = mapped_file->get_refInfo();
	} else {
		cerr << "File Format not recognized. File must be a sorted .bam file!" << endl;
		exit(EXIT_FAILURE);
	}

	std::map<std::string, std::vector<str_breakpoint_slim> > entries;

	//initialize for chrs
	std::vector<str_breakpoint_slim> tmp;
	for (size_t i = 0; i < ref.size(); i++) {
		entries[ref[i].RefName] = tmp;
	}

	//parse + include sorted
	read_SVs(entries);

	//obtain coverage
	update_coverage(entries);

	//Update VCF file for these entries and put it in a tmp file.
	update_svs_output(entries);

	string del = "rm ";
	del += Parameter::Instance()->tmp_genotyp;
	system(del.c_str());
}

//================================================================================================
//================================================================================================
//================================================================================================
//================================================================================================

std::string Genotyper::assess_genotype(int total_coverage, int SV_support) {
	double allele = 0;

	if (total_coverage < SV_support) {
		if (double((double)SV_support / (double)total_coverage) < 2.0) {
			cerr << "Warning @Genotyper:refrence coverage. Please report this! " << endl;
		}
		total_coverage = SV_support;
	}
//	if (!Parameter::Instance()->testing) {
//		ref += support; // just for now!
//	}

	//cout << "REF " << total_coverage << endl;
	if (total_coverage > 0) {
		allele = (double) SV_support / (double) total_coverage;	//(support + ref);
	}
	if (allele < Parameter::Instance()->min_allelel_frequency) {
		return "";
	}
	//cout << "ALLELE " << allele << endl;
	std::stringstream ss;
	ss << ";AF=";
	ss << allele;
	ss << "\tGT:DR:DV\t";
	if (total_coverage < 2) {
		ss << "./."; //we cannot define it.
	} else if (allele > Parameter::Instance()->homfreq) {
		ss << "1/1";
	} else if (allele > Parameter::Instance()->hetfreq) {
		ss << "0/1";
	} else {
		ss << "0/0";
	}
	ss << ":";
	ss << total_coverage - SV_support;
	ss << ":";
	ss << SV_support;
	return ss.str();
}

//// ============== Fisher exact test for strandness ===========
void initLogFacs(double* logFacs, int n) {
	logFacs[0] = 0;
	for (int i = 1; i < n + 1; ++i) {
		logFacs[i] = logFacs[i - 1] + log((double) i); // only n times of log() calls
	}
}

double logHypergeometricProb(double* logFacs, int a, int b, int c, int d) {
	return logFacs[a + b] + logFacs[c + d] + logFacs[a + c] + logFacs[b + d] - logFacs[a] - logFacs[b] - logFacs[c] - logFacs[d] - logFacs[a + b + c + d];
}

double logFac(int n) {
	double ret;
	for (ret = 0.; n > 0; --n) {
		ret += log((double) n);
	}
	return ret;
}
double logHypergeometricProb(int a, int b, int c, int d) {
	return logFac(a + b) + logFac(c + d) + logFac(a + c) + logFac(b + d) - logFac(a) - logFac(b) - logFac(c) - logFac(d) - logFac(a + b + c + d);
}

double fisher_exact(int sv_plus, int sv_minus, int ref_plus, int ref_minus) {
	int n = sv_plus + sv_minus + ref_plus + ref_minus;
	double* logFacs = new double[n + 1]; // *** dynamically allocate memory logFacs[0..n] ***
	initLogFacs(logFacs, n); // *** initialize logFacs array ***
	double logpCutoff = logHypergeometricProb(logFacs, sv_plus, sv_minus, ref_plus, ref_minus); // *** logFacs added
	double pFraction = 0;
	for (int x = 0; x <= n; ++x) {
		if (sv_plus + sv_minus - x >= 0 && sv_plus + ref_plus - x >= 0 && ref_minus - sv_plus + x >= 0) {
			double l = logHypergeometricProb(logFacs, x, sv_plus + sv_minus - x, sv_plus + ref_plus - x, ref_minus - sv_plus + x);
			if (l <= logpCutoff)
				pFraction += exp(l - logpCutoff);
		}
	}
	double logpValue = logpCutoff + log(pFraction);
//	std::cout << "Two-sided log10-p-value is " << logpValue / log(10.) << std::endl;
//	std::cout << "Two-sided p-value is " << exp(logpValue) << std::endl;
	delete[] logFacs;
	return exp(logpValue);
}

std::string Genotyper::mod_breakpoint_vcf(string buffer, int ref_plus, int ref_minus) {
//find last of\t
//parse #reads supporting
//print #ref

	string entry;
	size_t pos = 0;
	pair<int, int> read_strands;
	pos = buffer.find("STRANDS2=");
	if (pos != string::npos) {
		read_strands.second = 0;
		while (pos < buffer.size() && buffer[pos] != '\t') {
			if (buffer[pos - 1] == '=') {
				read_strands.first = atoi(&buffer[pos]);
			}
			if (buffer[pos - 1] == ',') {
				read_strands.second = atoi(&buffer[pos]);
				break;
			}
			pos++;
		}
	}
	//cout<<"start2 "<<read_strands.first<<" "<< read_strands.second <<endl;
	double pval = fisher_exact(read_strands.first, read_strands.second, ref_plus, ref_minus);
	//cout<<"next"<<endl;
	pos = buffer.find_last_of("GT");
//tab
	entry = buffer.substr(0, pos - 2);
	std::stringstream ss;
	ss << ";REF_strand=";
	buffer = buffer.substr(pos + 1);		// the right part is only needed:
	pos = buffer.find_last_of(':');
	int support = atoi(buffer.substr(pos + 1).c_str()); //parse RE in GT field.

	int total_coverage = (ref_plus + ref_minus);		// TODO not nice but just to make sure.
	//cout << "\t support: " << support << " ref_plus " << ref_plus << " ref_minus " << ref_minus << endl;

	ss << ref_plus << "," << ref_minus;
	ss << ";Strandbias_pval=" << pval;
	entry += ss.str();

	if (read_strands.first + read_strands.second > 5 && pval < 0.01) {
		pos = entry.find("PASS");
		if (pos != string::npos) {
			entry.erase(pos, 4);
			entry.insert(pos, "STRANDBIAS");
		}

		/*	pos = 0;
		 int count = 0;
		 for (size_t i = 0; i < entry.size(); i++) {

		 if (entry[i] == '.') {
		 pos = i + 2; //for avoiding . and \t
		 }
		 if (entry[i] == '\t' && pos != 0) {
		 count++;
		 if (count == 2) {
		 entry.erase(pos, i - pos);
		 entry.insert(pos, "STRANDBIAS");
		 }
		 }
		 }*/

	}

	string msg = assess_genotype(total_coverage, support);
	if (msg.empty()) {
		return "";
	}
	entry += msg;
	//cout<<"done"<<endl;
	return entry;

}

std::string Genotyper::mod_breakpoint_bedpe(string buffer, int ref_plus, int ref_minus) {

	std::string tmp = buffer;
	std::string entry = tmp;
	entry += '\t';
//int ref = max(tree.get_ref(node,var.chr,var.pos),tree.get_ref(node,var.chr2,var.pos2));

	int pos = tmp.find_last_of('\t');		//TODO!!
	int support = atoi(tmp.substr(pos + 1).c_str());
	double allele = (double) support / (double) (ref_plus + ref_minus);

	if (allele < Parameter::Instance()->min_allelel_frequency) {
		return "";
	}

	std::stringstream ss;
	ss << ref_plus << "," << ref_minus;
	ss << "\t";
	ss << support;
	entry += ss.str();
	return entry;
}

void Genotyper::parse_pos(char * buffer, int & pos, std::string & chr) {
	chr = "";
	pos = -1;
	size_t i = 0;
	int count = 0;
	while (buffer[i] != '\t') {
		if (count == 1 && ((buffer[i] != '[' || buffer[i] != ']') && buffer[i] != ':')) {
			chr += buffer[i];
		}
		if (count == 2 && buffer[i - 1] == ':') {
			pos = atoi(&buffer[i]);
		}
		if ((buffer[i] == ']' || buffer[i] == '[') || buffer[i] == ':') {
			count++;
		}
		i++;
	}
}

variant_str Genotyper::get_breakpoint_vcf(string buffer) {
//TODO extend for TRA!
	size_t i = 0;
	int count = 0;

	variant_str tmp;

	while (buffer[i] != '\0' && buffer[i] != '\n') {
		if (count == 0 && buffer[i] != '\t') {
			tmp.chr += buffer[i];
		}
		if (count == 1 && buffer[i - 1] == '\t') {
			tmp.pos = atoi(&buffer[i]);
		}
		if (tmp.pos2 == -1 && (count == 4 && (buffer[i - 1] == '[' || buffer[i - 1] == ']'))) {
			parse_pos(&buffer[i - 1], tmp.pos2, tmp.chr2);
		}

		if (count > 6 && strncmp(";CHR2=", &buffer[i], 6) == 0) {
			i += 6;
			while (buffer[i] != ';') {
				tmp.chr2 += buffer[i];
				i++;
			}
		}
		if (count > 6 && strncmp(";END=", &buffer[i], 5) == 0) {
			tmp.pos2 = atoi(&buffer[i + 5]); //stores right most breakpoint
			break;
		}

		if (buffer[i] == '\t') {
			count++;
		}
		i++;
	}
	return tmp;
}
variant_str Genotyper::get_breakpoint_bedpe(string buffer) {
	size_t i = 0;
	int count = 0;
	std::string chr;
	variant_str tmp;

	while (buffer[i] != '\0' && buffer[i] != '\n') {
		if (count == 12 && buffer[i] != '\t') {
			tmp.chr += buffer[i];
		}

		if (count == 13 && buffer[i - 1] == '\t') {
			tmp.pos = atoi(&buffer[i]);
		}
		if (count == 14 && buffer[i] != '\t') {
			tmp.chr2 += buffer[i];
		}
		if (count == 15 && buffer[i - 1] == '\t') {
			tmp.pos2 = atoi(&buffer[i]);
			break;
		}
		if (buffer[i] == '\t') {
			count++;
		}
		i++;
	}
	return tmp;
}

void Genotyper::update_file(Breakpoint_Tree & tree, breakpoint_node *& node) {
	std::ifstream myfile;
	bool is_vcf = !Parameter::Instance()->output_vcf.empty();

	string file_name;
	if (!Parameter::Instance()->output_vcf.empty()) {
		file_name = Parameter::Instance()->output_vcf;
		myfile.open(Parameter::Instance()->output_vcf.c_str(), std::ifstream::in);
	} else if (!Parameter::Instance()->output_bedpe.empty()) {
		file_name = Parameter::Instance()->output_bedpe;
		myfile.open(Parameter::Instance()->output_bedpe.c_str(), std::ifstream::in);
	}

	FILE*file = fopen(Parameter::Instance()->tmp_file.c_str(), "w"); //

	if (!myfile.good()) {
		std::cout << "SVParse: could not open file: " << std::endl;
		exit(EXIT_FAILURE);
	}

	string buffer;
	getline(myfile, buffer);
//parse SVs breakpoints in file

	while (!myfile.eof()) { // TODO:if first -> we need to define AF!
		if (buffer[0] != '#') {
			std::string to_print;

			//Go variant per variant in the VCF obtain positions:
			variant_str tmp;
			if (is_vcf) {
				tmp = get_breakpoint_vcf(buffer);
			} else {
				tmp = get_breakpoint_bedpe(buffer);
			}

			std::pair<int, int> first_node = tree.get_ref(node, tmp.chr, tmp.pos);
			std::pair<int, int> second_node = tree.get_ref(node, tmp.chr2, tmp.pos2);

			std::pair<int, int> final_ref;
			if (first_node.first + first_node.second > second_node.first + second_node.second) {
				final_ref = first_node;
			} else {
				final_ref = second_node;
			}

			if (final_ref.first == -1) {
				std::cerr << "Error in GT: Tree node not found. Exiting." << std::endl;
				exit(EXIT_FAILURE);
			}

			//int ref = final_ref.first + final_ref.second;

			if (is_vcf) {
				to_print = mod_breakpoint_vcf(buffer, final_ref.first, final_ref.second);
			} else {
				to_print = mod_breakpoint_bedpe(buffer, final_ref.first, final_ref.second);
			}
			if (!to_print.empty()) {
				fprintf(file, "%s", to_print.c_str());
				fprintf(file, "%c", '\n');
			}
		} else {
			fprintf(file, "%s", buffer.c_str());
			fprintf(file, "%c", '\n');
		}

		getline(myfile, buffer);
	}
	myfile.close();
	fclose(file);

	string move = "mv ";
	move += Parameter::Instance()->tmp_file;
	move += " ";
	move += file_name;
	system(move.c_str());
}

std::vector<std::string> Genotyper::read_SVs(Breakpoint_Tree & tree, breakpoint_node * &node) {
	std::vector<std::string> ref_dict;
	std::ifstream myfile;
	bool is_vcf = !Parameter::Instance()->output_vcf.empty();

	if (!Parameter::Instance()->output_vcf.empty()) {
		myfile.open(Parameter::Instance()->output_vcf.c_str(), std::ifstream::in);
	} else if (!Parameter::Instance()->output_bedpe.empty()) {
		myfile.open(Parameter::Instance()->output_bedpe.c_str(), std::ifstream::in);
	}

	if (!myfile.good()) {
		std::cout << "SVParse: could not open file: " << std::endl;
		exit(EXIT_FAILURE);
	}
//size_t buffer_size = 250000000;
	string buffer;

	getline(myfile, buffer);
//char* buffer = new char[buffer_size];
//myfile.getline(buffer, buffer_size);
//parse SVs breakpoints in file

	int num_sv = 0;
	int prev_pos1 = 0;
	int prev_pos2 = 0;
	while (!myfile.eof()) {
		//cout << buffer << endl;
		if (buffer[0] != '#') {
			// create binary tree to hold breakpoints!

			variant_str tmp;
			if (is_vcf) {
				tmp = get_breakpoint_vcf(buffer);
			} else {
				tmp = get_breakpoint_bedpe(buffer);
			}
			//	std::cout << "SV: " << tmp.pos << " " << tmp.pos2 << std::endl;
			tree.insert(node, tmp.chr, tmp.pos, true); //true: start;
			tree.insert(node, tmp.chr2, tmp.pos2, false); //false: stop;//
			num_sv++;
			if (num_sv % 5000 == 0) {
				cout << "\t\tRead in SV: " << num_sv << endl;
			}
		} else if (buffer[2] == 'c' && buffer[3] == 'o') { //##contig=<ID=chr1,length=699930>
		//fill the refdict.
			std::string id = "";
			for (size_t i = 13; i < buffer.size() && buffer[i] != ','; i++) {
				id += buffer[i];
			}
			ref_dict.push_back(id);
		}
		getline(myfile, buffer);
		//myfile.getline(buffer, buffer_size);
	}
	myfile.close();
	return ref_dict;
//tree.inorder(node);
}
void Genotyper::compute_cov(Breakpoint_Tree & tree, breakpoint_node *& node, std::vector<std::string> ref_dict) {

	/*FILE * ref_allel_reads = fopen(Parameter::Instance()->tmp_genotyp.c_str(), "r");
	 if (ref_allel_reads == NULL) {
	 std::cerr << "CovParse: could not open file: " << Parameter::Instance()->tmp_genotyp << std::endl;
	 }
	 //check if we want to compute the full coverage!

	 size_t nbytes = fread(&tmp, sizeof(struct str_read), 1, ref_allel_reads);*/
	std::string buffer;
	std::ifstream myfile;

	str_read tmp;
//cout<<"File: "<< Parameter::Instance()->tmp_genotyp.c_str()<<endl;
	myfile.open(Parameter::Instance()->tmp_genotyp.c_str(), std::ifstream::in);
	if (!myfile.good()) {
		std::cout << "Genotype Parser: could not open file: " << Parameter::Instance()->tmp_genotyp << std::endl;
		exit(0);
	}
	int prev_id = -1;
	int num_reads = 0;

	getline(myfile, buffer);
//cout<<buffer<<endl;
	while (!myfile.eof()) {
		int count = 0;

		tmp.chr_id = atoi(&buffer[0]);
		for (size_t i = 0; i < buffer.size() && buffer[i] != '\n' && buffer[i] != '\0'; i++) {
			if (count == 1 && buffer[i - 1] == '\t') {
				tmp.start = atoi(&buffer[i]);
			}
			if (count == 2 && buffer[i - 1] == '\t') {
				tmp.length = atoi(&buffer[i]);
			}
			if (count == 3 && buffer[i - 1] == '\t') {
				tmp.strand = atoi(&buffer[i]);
				break;
			}
			if (buffer[i] == '\t') {
				count++;
			}
		}

		//while (nbytes != 0) {
		//	std::cout << "Read: " << " " << tmp.chr_id << ":" << ref_dict[tmp.chr_id] << " " << tmp.start << " " << tmp.length << std::endl;
		if (prev_id != tmp.chr_id) {
			cout << "\t\tScanning CHR " << ref_dict[tmp.chr_id] << endl;
			prev_id = tmp.chr_id;
		}
		if (tmp.strand == 1) {		//strand of read
			tree.overalps(tmp.start, tmp.start + tmp.length, ref_dict[tmp.chr_id], node, true);
		} else {
			tree.overalps(tmp.start, tmp.start + tmp.length, ref_dict[tmp.chr_id], node, false);
		}
		//nbytes = fread(&tmp, sizeof(struct str_read), 1, ref_allel_reads);
		num_reads++;
		getline(myfile, buffer);
	}

//cout << "Num: " << num_reads << endl;
	myfile.close();
//fclose (ref_allel_reads);
//	tree.inorder(node);
}

void Genotyper::update_SVs() {
//parse SVs not just breakpoints and update with the coverage info
	cout << "\tConstruct tree" << endl;
	std::vector<std::string> ref_dict = read_SVs(this->tree, this->node);
	cout << "\tUpdate reference alleles" << endl;
	compute_cov(this->tree, this->node, ref_dict);
	cout << "\tWriting SV calls" << endl;
	update_file(this->tree, this->node);
	cout << "\tCleaning tmp files" << endl;
	string del = "rm ";
	del += Parameter::Instance()->tmp_genotyp;
	system(del.c_str());
}

void Genotyper::update_SVs(std::vector<Breakpoint *> & svs, long ref_space) { //refspace for the ref reads!!
	FILE * ref_allel_reads = fopen(Parameter::Instance()->tmp_genotyp.c_str(), "r");
	if (ref_allel_reads == NULL) {
		std::cerr << "Genotype Parser: could not open file: " << Parameter::Instance()->tmp_genotyp << std::endl;
	}
	str_read tmp;
	size_t nbytes = fread(&tmp, sizeof(struct str_read), 1, ref_allel_reads);
	int num_reads = 0;

//Do the tree per chr  to reduce the complexity?
//
	while (nbytes != 0) {
		for (size_t i = 0; i < svs.size(); i++) {
			if (svs[i]->get_valid()) {
				long start = tmp.start + ref_space;
				long stop = start + (long) tmp.length;
				//start - 100 orig!
				if ((svs[i]->get_coordinates().start.min_pos - 100 > start && svs[i]->get_coordinates().start.min_pos + 100 < stop)) { //found
					svs[i]->set_refcount(svs[i]->get_refcount() + 1);
				}
				//stop coordinate
				if ((svs[i]->get_coordinates().stop.max_pos - 100 > start + 100 && svs[i]->get_coordinates().stop.max_pos + 100 < stop - 100)) { //found
					svs[i]->set_refcount(svs[i]->get_refcount() + 1);
				}
			}
		}
		//if reads should be included-> Planesweep for +- breakpoint (Maybe hit -> extra function for that region around the breakpoint!
		num_reads++;
		if (num_reads % 1000 == 0) {
			cout << "\tProcessed " << num_reads << endl;
		}
		nbytes = fread(&tmp, sizeof(struct str_read), 1, ref_allel_reads);
	}
	fclose(ref_allel_reads);
}

