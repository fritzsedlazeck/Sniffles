/*
 * Cluster_SVs.cpp
 *
 *  Created on: Apr 28, 2016
 *      Author: fsedlaze
 */

#include "Cluster_SVs.h"

std::map<long, std::vector<int> > Cluster_SVS::parse_names_ids(int & max_ID) {
	FILE * alt_allel_reads = fopen(Parameter::Instance()->tmp_phasing.c_str(), "r");
	if (alt_allel_reads == NULL) {
		std::cerr << "ClusterParse: could not open tmp file: " << Parameter::Instance()->tmp_phasing << std::endl;
	}

	std::map<long, std::vector<int> > names;
	name_str tmp;
	size_t nbytes = fread(&tmp, sizeof(struct name_str), 1, alt_allel_reads);
	while (nbytes != 0) {
		max_ID = std::max(max_ID, tmp.svs_id); //needs to be a long as we need to know the size prior to storing!
		names[tmp.read_name].push_back(tmp.svs_id);
		nbytes = fread(&tmp, sizeof(struct name_str), 1, alt_allel_reads);
	}
	fclose(alt_allel_reads);
	return names;
}

void Cluster_SVS::update_SVs(std::vector<combine_str> & ids) {
	std::ifstream myfile;
	bool is_vcf = !Parameter::Instance()->output_vcf.empty();
	std::string filename;
	int col;
	if (is_vcf) {
		col = 2;
		filename = Parameter::Instance()->output_vcf;
	} else {
		col = 6;
		filename = Parameter::Instance()->output_bedpe;
	}
	myfile.open(filename.c_str(), std::ifstream::in);
	if (!myfile.good()) {
		std::cout << "Cluster Parse: could not open file: " << std::endl;
		exit(EXIT_FAILURE);
	}

	std::string tmp_name_file = filename;
	tmp_name_file += ".tmp";
	FILE*file = fopen(tmp_name_file.c_str(), "w");

	size_t buffer_size = 2500000;
	char* buffer = new char[buffer_size];
	myfile.getline(buffer, buffer_size);
	//parse SVs breakpoints in file
	while (!myfile.eof()) {
		if (buffer[0] != '#') {
			int count = 0;
			for (size_t i = 0; i < buffer_size && (buffer[i] != '\0' && buffer[i] != '\n'); i++) {
				if (count == col) { //if colum of id:
					if (buffer[i - 1] == '\t') {
						int id = atoi(&buffer[i]);
						fprintf(file, "%s", find_id(id, ids).c_str());
						fprintf(file, "%c", '\t');
					}
				} else {
					fprintf(file, "%c", buffer[i]);
				}
				if (buffer[i] == '\t') {
					count++;
				}
			}
		} else {
			fprintf(file, "%s", buffer);
		}
		fprintf(file, "%c", '\n');
		myfile.getline(buffer, buffer_size);
	}
	myfile.close();
	fclose(file);

	std::string move = "mv ";
	move += tmp_name_file;
	move += " ";
	move += filename;
	system(move.c_str());
}
void Cluster_SVS::add_id(int curr_id, int new_id, std::vector<combine_str> & ids, int subkey) {

	for (size_t i = 0; i < ids.size(); i++) { //check if already in the new array
		if (ids[i].curr_id == curr_id && ids[i].alt_id == new_id) {
			ids[i].support++;
			ids[i].hit = subkey;
			return;
		}
	}

	//make new entry:
	combine_str tmp;
	tmp.curr_id = curr_id;
	tmp.alt_id = new_id; //smallest ID of SVs
	tmp.support = 1;
	tmp.hit = subkey;
	ids.push_back(tmp);
}
std::string Cluster_SVS::find_id(int curr_id, std::vector<combine_str> & ids) {
	std::stringstream ss;
	for (size_t i = 0; i < ids.size(); i++) {
		if (ids[i].support > Parameter::Instance()->min_grouping_support) {
			if (ids[i].curr_id == curr_id) {

				ss << ids[i].alt_id;
				ss << '_';
				ss << ids[i].hit;
				return ss.str();
			}
		}
	}
	ss << curr_id;
	return ss.str();
}
void Cluster_SVS::update_SVs() {
//1: read in names + IDs -> store in map!
	int max_ID = 0;

	//TODO: restructure!
	//id=svs_id;

	std::map<long, std::vector<int> > names = parse_names_ids(max_ID); //key = read_id values: SVs id's

//2: make array with ID as entry and value is the smalles ID in the colum of all storred readnames.
	std::vector<combine_str> ids;
	for (std::map<long, std::vector<int> >::iterator i = names.begin(); i != names.end(); i++) {
		if ((*i).second.size() > 1) {
			int min_id = max_ID + 1;
			for (size_t j = 0; j < (*i).second.size(); j++) { // get the smallest ID of SVs associated with the read id.
				min_id = std::min(min_id, (*i).second[j]);
			}
			//min_id is now the smallest SVs id that this read is associated with.
			int subkey = 0;
			for (size_t j = 0; j < (*i).second.size(); j++) { // update the other SVs IDs
				//if ((*i).second[j] != min_id) {
				add_id((*i).second[j], min_id, ids, subkey);
				subkey++;
				//}
			}
		}
	}
	names.clear();

//3: Update the IDS in the VCF/Bedpe files.
	update_SVs(ids);
	std::cout << "\tCleaning tmp files" << std::endl;
	std::string del = "rm ";
	del += Parameter::Instance()->tmp_phasing;
	system(del.c_str());

}
