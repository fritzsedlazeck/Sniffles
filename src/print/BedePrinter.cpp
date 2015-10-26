/*
 * BedePrinter.cpp
 *
 *  Created on: Aug 24, 2015
 *      Author: fsedlaze
 */

#include "BedePrinter.h"

void BedePrinter::print_header() {
//nothing to be done!
}
void BedePrinter::print_body(std::vector<Breakpoint *> &SV, RefVector ref) {
	FILE *file;
	file = fopen(Parameter::Instance()->output_bede.c_str(), "w");
	//chr1    934247  934273  chr1    934690  934692  1       1.36749e-107    +       -       TYPE:DELETION   IDS:2,36        STRANDS:+-,36   MAX:chr1:934248;chr1:934692     95:chr1:934248-934252;chr1:934692-934692
	std::string chr;
	for (size_t i = 0; i < SV.size(); i++) {
		int pos = IPrinter::calc_pos(SV[i]->get_coordinates().start.max_pos, ref, chr);

		fprintf(file, "%s", chr.c_str());
		fprintf(file, "%c", '\t');

		fprintf(file, "%i", pos);
		fprintf(file, "%c", '\t');
		fprintf(file, "%i", pos + 1);
		fprintf(file, "%c", '\t');
		pos = IPrinter::calc_pos(SV[i]->get_coordinates().stop.max_pos, ref, chr);
		fprintf(file, "%s", chr.c_str());
		fprintf(file, "%c", '\t');
		fprintf(file, "%i", pos);
		fprintf(file, "%c", '\t');
		fprintf(file, "%i", pos + 1);
		fprintf(file, "%c", '\t');
		fprintf(file, "%i", SV[i]->get_id());
		fprintf(file, "%c", '\t');
		fprintf(file, "%c", '0'); //TODO: think about eval!
		fprintf(file, "%c", '\t');
		fprintf(file, "%s", SV[i]->get_strand(1).c_str());
		fprintf(file, "%c", '\t');
		fprintf(file, "%s", IPrinter::get_type(SV[i]->get_SVtype()).c_str()); //TODO
		fprintf(file, "%c", '\t');
		fprintf(file, "%s", "IDS:??");
		fprintf(file, "%c", '\t');
		fprintf(file, "%s", "STRANDS:??,666");
		fprintf(file, "%c", '\t');
		fprintf(file, "%s", "MAX:");
		pos = IPrinter::calc_pos(SV[i]->get_coordinates().start.most_support, ref, chr);
		fprintf(file, "%s", chr.c_str());
		fprintf(file, "%c", ':');
		fprintf(file, "%i", pos);
		fprintf(file, "%c", ',');
		pos = IPrinter::calc_pos(SV[i]->get_coordinates().stop.most_support, ref, chr);
		fprintf(file, "%s", chr.c_str());
		fprintf(file, "%c", ':');
		fprintf(file, "%i", pos);
		fprintf(file, "%c", '\n');
	}
	fclose(file);
}
