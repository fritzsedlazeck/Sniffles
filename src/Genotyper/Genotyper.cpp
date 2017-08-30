/*
 / * Genotyper.cpp
 *
 *  Created on: Mar 28, 2016
 *      Author: fsedlaze
 */

#include "Genotyper.h"

std::string Genotyper::mod_breakpoint_vcf(char *buffer, int ref) {
	//find last of\t
	//parse #reads supporting
	//print #ref
	string entry;
	string tmp = string(buffer);
	int pos = 0;

	pos = tmp.find_last_of("GT");
	//tab
	entry = tmp.substr(0, pos - 2);

	tmp = tmp.substr(pos + 1);	// the right part is only needed:
	pos = tmp.find_last_of(':');
	int support = atoi(tmp.substr(pos + 1).c_str());
	double allele = (double) support / (double) (support + ref);

	if (allele < Parameter::Instance()->min_allelel_frequency) {
		return "";
	}
	std::stringstream ss;
	ss << ";AF=";
	ss << allele;
	ss << "\tGT:DR:DV\t";

	if (allele > 0.8) {
		ss << "1/1:";
	} else if (allele > 0.3) {
		ss << "0/1:";
	} else {
		ss << "0/0:";
	}

	ss << ref;
	ss << ":";
	ss << support;

	entry += ss.str();
	return entry;

}

std::string Genotyper::mod_breakpoint_bedpe(char *buffer, int ref) {

	std::string tmp = string(buffer);
	std::string entry = tmp;
	entry += '\t';
	//int ref = max(tree.get_ref(node,var.chr,var.pos),tree.get_ref(node,var.chr2,var.pos2));

	int pos = tmp.find_last_of('\t'); //TODO!!
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

variant_str Genotyper::get_breakpoint_vcf(char *buffer) {
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
			tmp.pos2 = atoi(&buffer[i+5]); //stores right most breakpoint
			break;
		}

		if (buffer[i] == '\t') {
			count++;
		}
		i++;
	}
	return tmp;
}
variant_str Genotyper::get_breakpoint_bedpe(char *buffer) {
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

	FILE*file = fopen(Parameter::Instance()->tmp_file.c_str(), "w");

	if (!myfile.good()) {
		std::cout << "SVParse: could not open file: " << std::endl;
		exit(0);
	}
	size_t buffer_size = 25000;
	char* buffer = new char[buffer_size];
	myfile.getline(buffer, buffer_size);
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
			int ref = max(tree.get_ref(node, tmp.chr, tmp.pos), tree.get_ref(node, tmp.chr2, tmp.pos2));
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
			fprintf(file, "%s", buffer);
			fprintf(file, "%c", '\n');
		}

		myfile.getline(buffer, buffer_size);
	}
	myfile.close();
	fclose(file);

	string move = "mv ";
	move += Parameter::Instance()->tmp_file;
	move += " ";
	move += file_name;
	system(move.c_str());
}

void Genotyper::read_SVs(Breakpoint_Tree & tree, breakpoint_node *& node) {

	std::ifstream myfile;
	bool is_vcf = !Parameter::Instance()->output_vcf.empty();

	if (!Parameter::Instance()->output_vcf.empty()) {
		myfile.open(Parameter::Instance()->output_vcf.c_str(), std::ifstream::in);
	} else if (!Parameter::Instance()->output_bedpe.empty()) {
		myfile.open(Parameter::Instance()->output_bedpe.c_str(), std::ifstream::in);
	}

	if (!myfile.good()) {
		std::cout << "SVParse: could not open file: " << std::endl;
		exit(0);
	}
	size_t buffer_size = 25000;
	char* buffer = new char[buffer_size];
	myfile.getline(buffer, buffer_size);
	//parse SVs breakpoints in file

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
			//std::cout<<"SV: "<<tmp.pos<<" "<<tmp.pos2<<std::endl;
			tree.insert(node, tmp.chr, tmp.pos,true); //true: start;
			tree.insert(node, tmp.chr2, tmp.pos2,false);//false: stop;//
		}
		myfile.getline(buffer, buffer_size);
	}
	myfile.close();
	//tree.inorder(node);
}
void Genotyper::compute_cov(Breakpoint_Tree & tree, breakpoint_node *& node) {
	std::string output = Parameter::Instance()->tmp_file.c_str();
	output += "ref_allele";

	FILE * ref_allel_reads = fopen(output.c_str(), "r");
	if (ref_allel_reads == NULL) {
		std::cerr << "CovParse: could not open file: " << output.c_str() << std::endl;
	}

	//check if we want to compute the full coverage!
	str_read tmp;
	size_t nbytes = fread(&tmp, sizeof(struct str_read), 1, ref_allel_reads);
	while (nbytes != 0) {
		if (!tmp.SV_support){
			//std::cout<<"Read: "<<tmp.start<<" "<<tmp.length<<std::endl;
			//if reads should be included-> Planesweep for +- breakpoint (Maybe hit -> extra function for that region around the breakpoint!
			tree.overalps(tmp.start, tmp.start + tmp.length, tmp.chr, node, tmp.SV_support);
		}
		nbytes = fread(&tmp, sizeof(struct str_read), 1, ref_allel_reads);
	}
	fclose(ref_allel_reads);
//	tree.inorder(node);
}

void Genotyper::update_SVs() {
	//parse SVs not just breakpoints and update with the coverage info
	read_SVs(this->tree, this->node);
	compute_cov(this->tree, this->node);
	update_file(this->tree, this->node);
	cout << "Cleaning tmp files" << endl;
	string del = "rm ";
	del += Parameter::Instance()->tmp_file;
	del += "ref_allele";
	system(del.c_str());
}

