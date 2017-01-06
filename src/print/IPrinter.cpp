/*
 * IPrinter.cpp
 *
 *  Created on: Aug 24, 2015
 *      Author: fsedlaze
 */

#include "IPrinter.h"

bool IPrinter::to_print(Breakpoint * &SV, double & std_start, double & std_stop) {

	std_start = 0;
	std_stop = 0;
	//comp_std(SV, std_start, std_stop);
	comp_std_quantile(SV, std_start, std_stop);
	bool to_print = true;

	if(((SV->get_SVtype() & INS) && SV->get_coordinates().stop.most_support- SV->get_coordinates().start.most_support == Parameter::Instance()->huge_ins) && (SV->get_types().is_ALN)){
		return (!SV->get_types().is_SR && (std_start<5 && std_stop<5));
	}

	if ((SV->get_SVtype() & INS) || (SV->get_SVtype() & DEL)) { //for insertions  + deletions:
		double dist = (double) (SV->get_coordinates().stop.most_support - SV->get_coordinates().start.most_support);
		//dist = (double)std::min(((int)dist * 4), Parameter::Instance()->max_dist);
		dist = dist * 4.0 * (uniform_variance / 2); //because we test against corrected value!
		//std::cout<<"DIST: "<<(SV->get_coordinates().stop.most_support - SV->get_coordinates().start.most_support)<<" "<<dist<<" STD: "<<std_start<<" "<<std_stop<<std::endl;
		return ((std_start < dist && std_stop < dist)); //0.2886751
	}

	//TODO test!
	if (SV->get_SVtype() & NEST) {
		return true;
	}
	double max_allowed=4*Parameter::Instance()->max_dist*(uniform_variance / 2);

	return (std_start < max_allowed && std_stop < max_allowed);

}

std::string IPrinter::get_chr(long pos, RefVector ref) {
	//	std::cout << "pos: " << pos << std::endl;
	size_t id = 0;
	while (id < ref.size() && pos >= 0) {
		pos -= ((long) ref[id].RefLength + (long) Parameter::Instance()->max_dist);
		//	std::cout << id << std::endl;
		id++;
	}

	return ref[id - 1].RefName;
}

void IPrinter::store_readnames(std::vector<long> names, int id) {
	name_str tmp;
	tmp.svs_id = id; //stays the same
	for (size_t i = 0; i < names.size(); i++) {
		tmp.read_name = names[i];
		fwrite(&tmp, sizeof(struct name_str), 1, this->tmp_file);
	}
}

long IPrinter::calc_pos(long pos, RefVector ref, std::string &chr) {
	size_t i = 0;
	pos -= (ref[i].RefLength + Parameter::Instance()->max_dist);

	while (i < ref.size() && pos >= 0) {
		i++;
		pos -= ((long) ref[i].RefLength + (long) Parameter::Instance()->max_dist);
	}
	chr = ref[i].RefName;
	return pos + ref[i].RefLength + (long) Parameter::Instance()->max_dist;
}

std::string IPrinter::get_type(char type) {
	string tmp;
	if (type & DEL) {
		tmp += "DEL";
	}
	if (type & INV) {
		if (!tmp.empty()) {
			tmp += '/';
		}
		tmp += "INV";
	}
	if (type & DUP) {
		if (!tmp.empty()) {
			tmp += '/';
		}
		tmp += "DUP";
	}
	if (type & INS) {
		if (!tmp.empty()) {
			tmp += '/';
		}
		tmp += "INS";
	}
	if (type & TRA) {
		if (!tmp.empty()) {
			tmp += '/';
		}
		//tmp += "BND";
		tmp += "TRA";
	}
	if (type & NEST) {
		if (!tmp.empty()) {
			tmp += '/';
		}
		tmp += "INVDUP";
	}
	return tmp; // should not occur!
}
// Get current date/time, format is YYYY-MM-DD.HH:mm:ss
const std::string IPrinter::currentDateTime() {
	time_t now = time(0);
	struct tm tstruct;
	char buf[80];
	tstruct = *localtime(&now);
	// Visit http://en.cppreference.com/w/cpp/chrono/c/strftime
	// for more information about date/time format
	strftime(buf, sizeof(buf), "%Y%m%d", &tstruct);
	return buf;
}

void IPrinter::sort_insert(int pos, std::vector<int> & positions) {

	size_t i = 0;
	while (i < positions.size() && positions[i] < pos) {
		i++;
	}
	positions.insert(positions.begin() + i, pos);

}
void IPrinter::comp_std_med(Breakpoint * &SV, double & std_start, double & std_stop) {
	std::vector<int> std_start_dists;
	std::vector<int> std_stop_dists;

	std::map<std::string, read_str> support = SV->get_coordinates().support;
	for (std::map<std::string, read_str>::iterator i = support.begin(); i != support.end(); i++) {
		if ((*i).second.SV & SV->get_SVtype()) {
			if ((*i).second.coordinates.first != -1) {
				long diff = (SV->get_coordinates().start.most_support - (*i).second.coordinates.first);
				//	std::cout << "DIFF Start: " << diff << std::endl;
				sort_insert(std::pow((double) diff, 2.0), std_start_dists);
				//std_start += std::pow((double) diff, 2.0);
			}
			if ((*i).second.coordinates.second != -1) {
				long diff = (SV->get_coordinates().stop.most_support - (*i).second.coordinates.second);
				//	std::cout << "DIFF Stop: " << diff << std::endl;
				sort_insert(std::pow((double) diff, 2.0), std_stop_dists);
				//std_stop += std::pow((double) diff, 2.0);
			}
		}
	}
	int median = std_stop_dists.size() / 2;
	std_start = std::sqrt(std_start_dists[median]);
	std_stop = std::sqrt(std_stop_dists[median]);
}

void IPrinter::comp_std_quantile(Breakpoint * &SV, double & std_start, double & std_stop) {
	double count = 0;
	std::vector<int> std_start_dists;
	std::vector<int> std_stop_dists;

	std::stringstream ss;

	std::map<std::string, read_str> support = SV->get_coordinates().support;

	fprintf(distances, "%s", IPrinter::get_type(SV->get_SVtype()).c_str());
	fprintf(distances, "%c", '\t');
	fprintf(distances, "%i", SV->get_coordinates().stop.most_support-SV->get_coordinates().start.most_support);
	fprintf(distances, "%c", '\t');
	fprintf(distances, "%i", support.size());
	fprintf(distances, "%c", '\t');
	for (std::map<std::string, read_str>::iterator i = support.begin(); i != support.end(); i++) {
		if ((*i).second.SV & SV->get_SVtype()) {
			if ((*i).second.coordinates.first != -1) {
				long diff = (SV->get_coordinates().start.most_support - (*i).second.coordinates.first);
				ss << '\t';
				ss << diff;
				sort_insert(std::pow((double) diff, 2.0), std_start_dists);
			}
			if ((*i).second.coordinates.second != -1) {
				long diff = (SV->get_coordinates().stop.most_support - (*i).second.coordinates.second);
				ss << '\t';
				ss << diff;

				sort_insert(std::pow((double) diff, 2.0), std_stop_dists);
			}
		}
	}

	count = 0;
	for (int i = 0; i < std::max((int)std_stop_dists.size() / 2,8); i++) {
		std_start += std_start_dists[i];
		std_stop += std_stop_dists[i];
		count++;
	}

	std_start = std::sqrt(std_start / count);
	std_stop = std::sqrt(std_stop / count);
	fprintf(distances, "%f", std_start);
	fprintf(distances, "%s",ss.str().c_str());
	fprintf(distances, "%c", '\n');


}
void IPrinter::comp_std(Breakpoint * &SV, double & std_start, double & std_stop) {
	double count = 0;
	std_start = 0;
	std_stop = 0;
	std::map<std::string, read_str> support = SV->get_coordinates().support;
	for (std::map<std::string, read_str>::iterator i = support.begin(); i != support.end(); i++) {
		if ((*i).second.SV & SV->get_SVtype()) {
			count++;
			if ((*i).second.coordinates.first != -1) {
				long diff = (SV->get_coordinates().start.most_support - (*i).second.coordinates.first);
				//	std::cout << "DIFF Start: " << diff << std::endl;
				std_start += std::pow((double) diff, 2.0);
			}
			if ((*i).second.coordinates.second != -1) {
				long diff = (SV->get_coordinates().stop.most_support - (*i).second.coordinates.second);
				//	std::cout << "DIFF Stop: " << diff << std::endl;
				std_stop += std::pow((double) diff, 2.0);
			}

		}
	}
	std_start = std::sqrt(std_start / count);
	std_stop = std::sqrt(std_stop / count);
}
