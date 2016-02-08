/*
 * MariaPrinter.cpp
 *
 *  Created on: Sep 4, 2015
 *      Author: fsedlaze
 */

#include "MariaPrinter.h"

void MariaPrinter::print_header() {
	fprintf(file, "%s", "Chrom\tstart\tstop\tchrom2\tstart2\tstop2\tvariant_name/ID\tscore (smaller is better)\tstrand1\tstrand2\ttype\tnumber_of_split_reads\tbest_chrom\tbest_pos\tbest_chrom2\tbest_pos2\n");
}
void MariaPrinter::print_body(Breakpoint *& SV, RefVector ref) {
	//"Chrom\tstart\tstop\tchrom2\tstart2\tstop2\tvariant_name/ID\tscore (smaller is better)\tstrand1\tstrand2\ttype\tnumber_of_split_reads\n"
	if (!this->bed_tree.is_in(SV->get_coordinates().start.most_support, this->root) && !this->bed_tree.is_in(SV->get_coordinates().stop.most_support, this->root)) {
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
		fprintf(file, "%c", '\n');
	}
}
