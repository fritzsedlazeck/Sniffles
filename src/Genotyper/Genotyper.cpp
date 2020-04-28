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

void update_entries(std::vector<str_breakpoint_slim> &entries, int start, int stop, size_t & current_pos, int wobble, std::string rname) { //TODO room for optimization!

	if (entries.empty() || stop + wobble < entries[0].pos) {
		return;
	}
//	cout<<"PS: "<<current_pos<<endl;
//	cout << "cov: " << start << " " << stop << endl;
	for (size_t i = current_pos; i < entries.size(); i++) {
		if (entries[i].pos < start - wobble) {
			current_pos = i;
		}
		if ((start - wobble < entries[i].pos && stop + wobble > entries[i].pos) && (abs(entries[i].pos - start) > wobble && abs(entries[i].pos - stop) > wobble)) {	//TODO not sure if I cannot combine these two.
			//entries[i].cov++;
			entries[i].rnames[rname] = true;
//			cout << "\tHIT: " << entries[i].rnames.size() << endl;
		}
		if (entries[i].pos > stop + wobble) {
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

	Alignment * tmp_aln = mapped_file->parseRead(Parameter::Instance()->min_mq);
	long num_reads = 0;

	int current_RefID = tmp_aln->getRefID();
	std::vector<str_breakpoint_slim> tmp = entries[ref[tmp_aln->getRefID()].RefName];

	size_t current_pos = 0;	// index of entries vector;

	while (!tmp_aln->getQueryBases().empty()) {

		if (current_RefID != tmp_aln->getRefID()) {
			current_pos = 0;
			entries[ref[current_RefID].RefName] = tmp;
			std::cout << "\t\tSwitch Chr " << ref[tmp_aln->getRefID()].RefName << std::endl;
			current_RefID = tmp_aln->getRefID();
			tmp = entries[ref[current_RefID].RefName];
		}

		int start = (int) tmp_aln->getPosition();
		int stop = (int) start + tmp_aln->getRefLength();
		//	cout << "Ref: " << ref[current_RefID].RefName << endl;
		update_entries(tmp, start, stop, current_pos, 5, tmp_aln->getName());

		mapped_file->parseReadFast(Parameter::Instance()->min_mq, tmp_aln);
	}

	entries[ref[current_RefID].RefName] = tmp;
	/*	cout << "Check:" << endl;
	 for (std::map<std::string, std::vector<str_breakpoint_slim> >::iterator i = entries.begin(); i != entries.end(); i++) {
	 for (size_t j = 0; j < (*i).second.size(); j++) {
	 cout << (*i).second[j].chr << " " << (*i).second[j].pos<<" "<< (*i).second[j].cov<< endl;
	 }
	 }*/

	std::cout << "\tFinalizing  .." << std::endl;
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

			//	if (start.pos == 10441961) {
			//		cout << "Found: " << start.chr << " " << start.pos << " " << stop.chr << " " << stop.pos << endl;
			//	}

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

	/*cout << "Check:" << endl;
	 for (std::map<std::string, std::vector<str_breakpoint_slim> >::iterator i = entries.begin(); i != entries.end(); i++) {
	 for (size_t j = 0; j < (*i).second.size(); j++) {
	 //if ((*i).second[j].pos == 10441961) {
	 cout << (*i).second[j].chr << " " << (*i).second[j].pos << endl;
	 //}
	 }
	 }*/
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
			// create binary tree to hold breakpoints!
			variant_str tmp;
			str_breakpoint_slim start;
			str_breakpoint_slim stop;
			if (is_vcf) {
				get_breakpoint_vcf(buffer, start, stop);
			} else {
				//get_breakpoint_bedpe(buffer); //TODO
			}
			map<std::string, bool> tmp_names;
			for (size_t i = 0; i < entries[start.chr].size(); i++) {
				if (start.pos == entries[start.chr][i].pos) {
					//	final_ref.first = entries[start.chr][i].cov;
					for (map<std::string, bool>::iterator t = entries[start.chr][i].rnames.begin(); t != entries[start.chr][i].rnames.end(); t++) {
						tmp_names[(*t).first] = true;
					}
					break;
				}
			}

			for (size_t i = 0; i < entries[stop.chr].size(); i++) {
				if (stop.pos == entries[stop.chr][i].pos) {
					//final_ref.second = entries[stop.chr][i].cov;
					for (map<std::string, bool>::iterator t = entries[stop.chr][i].rnames.begin(); t != entries[stop.chr][i].rnames.end(); t++) {
						tmp_names[(*t).first] = true;
					}

					break;
				}
			}

			if (tmp_names.empty() == -1) {
				std::cerr << "Error in GT: Tree node not found. Exiting." << std::endl;
				exit(EXIT_FAILURE);
			}
			if (is_vcf) {
				to_print = mod_breakpoint_vcf(buffer, (int) tmp_names.size());
			} else {
				to_print = mod_breakpoint_bedpe(buffer, (int) tmp_names.size());
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

std::string Genotyper::assess_genotype(int ref, int support) {
	double allele = 0;
	if (!Parameter::Instance()->testing) {
		ref += support; // just for now!
	}
	if ((ref) > 0) {
		allele = (double) support / (double) ref;	//(support + ref);
	}
	if (allele < Parameter::Instance()->min_allelel_frequency) {
		return "";
	}
	std::stringstream ss;
	ss << ";AF=";
	ss << allele;
	ss << "\tGT:DR:DV\t";
	if (ref < 2) {
		ss << "./."; //we cannot define it.
	} else if (allele > Parameter::Instance()->homfreq) {
		ss << "1/1";
	} else if (allele > Parameter::Instance()->hetfreq) {
		ss << "0/1";
	} else {
		ss << "0/0";
	}
	ss << ":";
	if (ref < support) {
		cerr << "Warning @Genotyper:refrence coverage. Please report this! " << endl;
		ss << ref - support;
	} else {
		ss << ref - support;
	}
	ss << ":";
	ss << support;
	return ss.str();
}

std::string Genotyper::mod_breakpoint_vcf(string buffer, int ref_strand) {
//find last of\t
//parse #reads supporting
//print #ref
	string entry;
	int pos = 0;

	pos = buffer.find_last_of("GT");
//tab
	entry = buffer.substr(0, pos - 2);
	std::stringstream ss;
	ss << ";REF_strand=";
	buffer = buffer.substr(pos + 1);		// the right part is only needed:
	pos = buffer.find_last_of(':');
	int support = atoi(buffer.substr(pos + 1).c_str());
	ref_strand = max(ref_strand, support);		// TODO not nice but just to make sure.
	ss << max(ref_strand, support);
	entry += ss.str();

	string msg = assess_genotype(ref_strand, support);
	if (msg.empty()) {
		return "";
	}
	entry += msg;
	return entry;

}

std::string Genotyper::mod_breakpoint_bedpe(string buffer, int ref) {

	std::string tmp = buffer;
	std::string entry = tmp;
	entry += '\t';
//int ref = max(tree.get_ref(node,var.chr,var.pos),tree.get_ref(node,var.chr2,var.pos2));

	int pos = tmp.find_last_of('\t');		//TODO!!
	int support = atoi(tmp.substr(pos + 1).c_str());
	double allele = (double) support / (double) (ref);

	if (allele < Parameter::Instance()->min_allelel_frequency) {
		return "";
	}

	std::stringstream ss;
	ss << ref;
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
			// create binary tree to hold breakpoints!
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

			int ref = final_ref.first + final_ref.second;

			if (is_vcf) {
				to_print = mod_breakpoint_vcf(buffer, ref);
			} else {
				to_print = mod_breakpoint_bedpe(buffer, ref);
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

