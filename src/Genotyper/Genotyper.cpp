/*
 / * Genotyper.cpp
 *
 *  Created on: Mar 28, 2016
 *      Author: fsedlaze
 */

#include "Genotyper.h"

std::string Genotyper::assess_genotype(int ref, int support) {
	double allele = (double) support / (double) (support + ref);

	if (allele < Parameter::Instance()->min_allelel_frequency) {
		return "";
	}
	if ((support + ref) == 0) {
		allele = 0;
	}

	std::stringstream ss;
	ss << ";AF=";
	ss << allele;
	ss << "\tGT:DR:DV\t";
	if (ref == 0 && support == 0) {
		ss << "./."; //we cannot define it.
	} else if (allele > Parameter::Instance()->homfreq) {
		ss << "1/1";
	} else if (allele > Parameter::Instance()->hetfreq) {
		ss << "0/1";
	} else {
		ss << "0/0";
	}
	ss << ":";
	ss << ref;
	ss << ":";
	ss << support;
	return ss.str();
}

std::string Genotyper::mod_breakpoint_vcf(string buffer, std::pair<int, int> ref_strand) {
	int ref = ref_strand.first + ref_strand.second;
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
	ss << ref_strand.first;
	ss << ",";
	ss << ref_strand.second;
	entry += ss.str();

	buffer = buffer.substr(pos + 1);		// the right part is only needed:
	pos = buffer.find_last_of(':');
	int support = atoi(buffer.substr(pos + 1).c_str());
	string msg = assess_genotype(ref, support);
	if (msg.empty()) {
		return "";
	}
	entry += msg;
	return entry;

}

std::string Genotyper::mod_breakpoint_bedpe(string buffer, std::pair<int, int> ref_strand) {

	int ref = ref_strand.first + ref_strand.second;
	std::string tmp = buffer;
	std::string entry = tmp;
	entry += '\t';
	//int ref = max(tree.get_ref(node,var.chr,var.pos),tree.get_ref(node,var.chr2,var.pos2));

	int pos = tmp.find_last_of('\t');		//TODO!!
	int support = atoi(tmp.substr(pos + 1).c_str());
	double allele = (double) support / (double) (support + ref);

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
			if (is_vcf) {
				to_print = mod_breakpoint_vcf(buffer, final_ref);
			} else {
				to_print = mod_breakpoint_bedpe(buffer, final_ref);
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

	/*fprintf(ref_allel_reads, "%i",tmp_aln->getRefID());
	 fprintf(ref_allel_reads, "%i",tmp_aln->getPosition());
	 fprintf(ref_allel_reads, "%i",tmp_aln->getRefLength());
	 if (tmp_aln->getStrand()) {
	 fprintf(ref_allel_reads, "%c",'1');
	 } else {
	 fprintf(ref_allel_reads, "%c",'2');
	 }*/
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

