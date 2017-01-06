/*
 * BedePrinter.cpp
 *
 *  Created on: Aug 24, 2015
 *      Author: fsedlaze
 */

#include "BedpePrinter.h"

void BedpePrinter::print_header() {
	fprintf(file, "%s", "#Chrom\tstart\tstop\tchrom2\tstart2\tstop2\tvariant_name/ID\tscore (smaller is better)\tstrand1\tstrand2\ttype\tnumber_of_split_reads\tbest_chr1\tbest_start\tbest_chr2\tbest_stop\tpredicted_length\n");
}
void BedpePrinter::print_body(Breakpoint * &SV, RefVector ref) {
	if (!this->bed_tree.is_in(SV->get_coordinates().start.most_support, this->root) && !this->bed_tree.is_in(SV->get_coordinates().stop.most_support, this->root)) {
		//temp. store read names supporting this SVs to later group the SVs together.
		if (Parameter::Instance()->phase) {
			store_readnames(SV->get_read_ids(), id);
		}

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
		fprintf(file, "%i", SV->get_length());
		fprintf(file, "%c", '\t');
		fprintf(file, "%i", SV->get_support());
		fprintf(file, "%c", '\n');
	}
}
