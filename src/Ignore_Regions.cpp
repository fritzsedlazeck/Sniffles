/*
 * Ignore_Regions.cpp
 *
 *  Created on: Feb 4, 2016
 *      Author: fsedlaze
 */

#include "Ignore_Regions.h"
#include "sub/Detect_Breakpoints.h"

long get_ref_coords(std::string chr, RefVector ref) {
	long length = 0;
	for (size_t i = 0; i < ref.size(); i++) {
		if (strcmp(ref[i].RefName.c_str(), chr.c_str()) == 0) {
			return length;
		}
		length += ref[i].RefLength + Parameter::Instance()->max_dist;
	}
	return -1; //should not happen

}


long get_ref_lengths2(int id, RefVector ref) {
	long length = 0;

	for (size_t i = 0; i < (size_t) id && i < ref.size(); i++) {
		length += ref[i].RefLength + Parameter::Instance()->max_dist;
	}
	return length;
}

int get_id2(RefVector ref, std::string chr) {
	for (size_t i = 0; i < ref.size(); i++) {
		if (strcmp(ref[i].RefName.c_str(), chr.c_str()) == 0) {
			return i;
		}
	}
	return -1; //should not happen!
}


void initialize_bed(IntervallTree_bed &bed_tree, Leaf *&root,RefVector ref) {
	//bst.insert(point, root);
	size_t buffer_size = 2000000;
	char*buffer = new char[buffer_size];
	std::ifstream myfile;

	myfile.open(Parameter::Instance()->ignore_regions_bed.c_str(), std::ifstream::in);
	if (!myfile.good()) {
		std::cout << "SAM Parser: could not open file: " << Parameter::Instance()->ignore_regions_bed.c_str() << std::endl;
		exit(0);
	}

	myfile.getline(buffer, buffer_size);

	while (!myfile.eof()) {
		int count = 0;
		string chr;
		int p1;
		int p2;
		for (size_t i = 0; i < buffer_size && buffer[i] != '\0' && buffer[i] != '\n'; i++) {
			if (count == 0 && buffer[i] != '\t') {
				chr += buffer[i];
			}
			if (count == 1 && buffer[i - 1] == '\t') {
				p1 = atoi(&buffer[i]);
			}
			if (count == 2 && buffer[i - 1] == '\t') {
				p2 = atoi(&buffer[i]);
				break;
			}
			if (buffer[i] == '\t') {
				count++;
			}
		}
		//transfer coordinates:
		long ref_dist = get_ref_coords(chr, ref);
		std::cout << (long) p1 + ref_dist << " " << (long) p2 + ref_dist << std::endl;
		bed_tree.insert((long) p1 + ref_dist, (long) p2 + ref_dist, root);

		myfile.getline(buffer, buffer_size);
	}

}

void ignore_regions(std::vector<Breakpoint *> & final_SV,RefVector ref) {

	IntervallTree_bed bed_tree;
	Leaf *root = NULL;
	initialize_bed(bed_tree, root,ref);

	bed_tree.postorder(root);
	size_t i = 0;
	while (i < final_SV.size()) {
		if (final_SV[i]->get_SVtype() & DUP) {
			std::cout << final_SV[i]->get_coordinates().start.most_support << " ";
		}
		if (bed_tree.is_in(final_SV[i]->get_coordinates().start.most_support, root) || bed_tree.is_in(final_SV[i]->get_coordinates().stop.most_support, root)) {
			if (final_SV[i]->get_SVtype() & DUP) {
				std::cout << "erase" << endl;
			}
			final_SV.erase(final_SV.begin() + i);
		} else {
			if (final_SV[i]->get_SVtype() & DUP) {
				std::cout << "keep" << endl;
			}
			i++;
		}
	}
}

