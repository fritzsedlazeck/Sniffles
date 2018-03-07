/*
 * BedePrinter.cpp
 *
 *  Created on: Aug 24, 2015
 *      Author: fsedlaze
 */

#include "BedpePrinter.h"

void BedpePrinter::print_header() {
	fprintf(file, "%s", "#Chrom\tstart\tstop\tchrom2\tstart2\tstop2\tvariant_name/ID\tscore (smaller is better)\tstrand1\tstrand2\ttype\tnumber_of_split_reads\tbest_chr1\tbest_start\tbest_chr2\tbest_stop\tpredicted_length\tFILTER\n");
}
void BedpePrinter::print_body(Breakpoint * &SV, RefVector ref) {
	if (!this->bed_tree.is_in(SV->get_coordinates().start.most_support, this->root) && !this->bed_tree.is_in(SV->get_coordinates().stop.most_support, this->root)) {
		//temp. store read names supporting this SVs to later group the SVs together.
		if (Parameter::Instance()->phase) {
			store_readnames(SV->get_read_ids(), id);
		}
		double std_quant_start = 0;
		double std_quant_stop = 0;

		pair<double, double> kurtosis;
		pair<double, double> std_quant;
		double std_length = 0;
		int zmws = 0;
		bool ok_to_print = (to_print(SV, std_quant, kurtosis, std_length, zmws) || Parameter::Instance()->ignore_std);
		if (ok_to_print && (zmws == 0 || zmws >= Parameter::Instance()->min_zmw)) {

			std::string chr;
			std::string strands = SV->get_strand(2);
			int pos = IPrinter::calc_pos(SV->get_coordinates().start.min_pos, ref, chr);

			//start coordinates:
			fprintf(file, "%s", chr.c_str());
			fprintf(file, "%c", '\t');
			fprintf(file, "%i", pos);
			fprintf(file, "%c", '\t');
			fprintf(file, "%i", IPrinter::calc_pos(SV->get_coordinates().start.max_pos, ref, chr));
			fprintf(file, "%c", '\t');

			//stop coordinates
			string chr_start;
			int start = IPrinter::calc_pos(SV->get_coordinates().start.most_support, ref, chr_start);

			long end_coord = SV->get_coordinates().stop.min_pos;
			if (((SV->get_SVtype() & INS) && SV->get_length() == Parameter::Instance()->huge_ins) && SV->get_types().is_ALN) {
				end_coord = std::max((SV->get_coordinates().stop.min_pos -(long)SV->get_length()), (long)start);
			}

			pos = IPrinter::calc_pos(end_coord, ref, chr);

			fprintf(file, "%s", chr.c_str());
			fprintf(file, "%c", '\t');
			fprintf(file, "%i", pos);
			fprintf(file, "%c", '\t');
			end_coord = SV->get_coordinates().stop.max_pos;
			if (((SV->get_SVtype() & INS) && SV->get_length() == Parameter::Instance()->huge_ins) && SV->get_types().is_ALN) {
				end_coord = std::max((SV->get_coordinates().stop.max_pos - (long) SV->get_length()), (long)start);
			}
			fprintf(file, "%i", IPrinter::calc_pos(end_coord, ref, chr));
			fprintf(file, "%c", '\t');
			fprintf(file, "%i", id);
			id++;
			fprintf(file, "%c", '\t');
			fprintf(file, "%i", -1); //TODO: score
			fprintf(file, "%c", '\t');
			fprintf(file, "%c", strands[0]);
			fprintf(file, "%c", '\t');
			fprintf(file, "%c", strands[1]);
			fprintf(file, "%c", '\t');
			fprintf(file, "%s", IPrinter::get_type(SV->get_SVtype()).c_str());
			fprintf(file, "%c", '\t');
			fprintf(file, "%i", SV->get_support());
			fprintf(file, "%c", '\t');
			fprintf(file, "%s", chr_start.c_str());
			fprintf(file, "%c", '\t');
			fprintf(file, "%i", start);
			fprintf(file, "%c", '\t');

			end_coord = SV->get_coordinates().stop.most_support;
			if (((SV->get_SVtype() & INS) && SV->get_length() == Parameter::Instance()->huge_ins) && SV->get_types().is_ALN) {
				end_coord =  std::max((SV->get_coordinates().stop.most_support - (long) SV->get_length()), (long)start);
			}

			pos = IPrinter::calc_pos(end_coord, ref, chr);
			fprintf(file, "%s", chr.c_str());
			fprintf(file, "%c", '\t');
			fprintf(file, "%i", pos);
			fprintf(file, "%c", '\t');


			if (((SV->get_SVtype() & INS) && SV->get_length() == Parameter::Instance()->huge_ins) && SV->get_types().is_ALN) {//!
				fprintf(file, "%s", "NA");
			} else {
				fprintf(file, "%i", SV->get_length());
			}

			if (((SV->get_SVtype() & INS) && SV->get_length() == Parameter::Instance()->huge_ins) && SV->get_types().is_ALN) {
				fprintf(file, "%s", "\tUNRESOLVED\t");
			} else 	if (std_quant.first < 10 && std_quant.second < 10) {
				fprintf(file, "%s", "PRECISE");
			} else {
				fprintf(file, "%s", "IMPRECISE");
			}

			//fprintf(file, "%c", '\t');
			//fprintf(file, "%i", SV->get_support());
			fprintf(file, "%c", '\n');
		}
	}
}

void BedpePrinter::print_body_recall(Breakpoint * &SV, RefVector ref) {
	std::string chr;
	std::string strands = SV->get_strand(2);
	int pos = IPrinter::calc_pos(SV->get_coordinates().start.min_pos, ref, chr);
	fprintf(file, "%s", chr.c_str());
	fprintf(file, "%c", '\t');
	fprintf(file, "%i", pos);
	fprintf(file, "%c", '\t');
	fprintf(file, "%i", IPrinter::calc_pos(SV->get_coordinates().start.max_pos, ref, chr));
	fprintf(file, "%c", '\t');
	pos = IPrinter::calc_pos(SV->get_coordinates().stop.min_pos, ref, chr);
	fprintf(file, "%s", chr.c_str());
	fprintf(file, "%c", '\t');
	fprintf(file, "%i", pos);
	fprintf(file, "%c", '\t');
	fprintf(file, "%i", IPrinter::calc_pos(SV->get_coordinates().stop.max_pos, ref, chr));
	fprintf(file, "%c", '\t');
	fprintf(file, "%i", id);
	id++;
	fprintf(file, "%c", '\t');
	fprintf(file, "%i", -1); //TODO: score
	fprintf(file, "%c", '\t');
	fprintf(file, "%c", strands[0]);
	fprintf(file, "%c", '\t');
	fprintf(file, "%c", strands[1]);
	fprintf(file, "%c", '\t');
	fprintf(file, "%s", IPrinter::get_type(SV->get_SVtype()).c_str());
	fprintf(file, "%c", '\t');
	fprintf(file, "%i", SV->get_support());
	fprintf(file, "%c", '\t');
	pos = IPrinter::calc_pos(SV->get_coordinates().start.most_support, ref, chr);
	fprintf(file, "%s", chr.c_str());
	fprintf(file, "%c", '\t');
	fprintf(file, "%i", pos);
	fprintf(file, "%c", '\t');
	pos = IPrinter::calc_pos(SV->get_coordinates().stop.most_support, ref, chr);
	fprintf(file, "%s", chr.c_str());
	fprintf(file, "%c", '\t');
	fprintf(file, "%i", pos);
	fprintf(file, "%c", '\t');
	if (((SV->get_SVtype() & INS) && SV->get_length() == Parameter::Instance()->huge_ins) && !SV->get_types().is_SR) {
		fprintf(file, "%s", "NA");
	} else {
		fprintf(file, "%i", SV->get_length());
	}
	//fprintf(file, "%c", '\t');
	//fprintf(file, "%i", SV->get_support());
	fprintf(file, "%c", '\n');
}

