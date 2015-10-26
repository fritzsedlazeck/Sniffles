/*
 * MariaPrinter.cpp
 *
 *  Created on: Sep 4, 2015
 *      Author: fsedlaze
 */

#include "MariaPrinter.h"

void MariaPrinter::print_header() {
	FILE *file;
	file = fopen(Parameter::Instance()->output_maria.c_str(), "w");
	fprintf(file, "%s", "Chrom\tstart\tstop\tchrom2\tstart2\tstop2\tvariant_name/ID\tscore (smaller is better)\tstrand1\tstrand2\ttype\tnumber_of_split_reads\tbest_chrom\tbest_pos\tbest_chrom2\tbest_pos2\n");
	fclose(file);
}
void MariaPrinter::print_body(std::vector<Breakpoint *>& SV, RefVector ref) {
	FILE *file;
	file = fopen(Parameter::Instance()->output_maria.c_str(), "a");
	//"Chrom\tstart\tstop\tchrom2\tstart2\tstop2\tvariant_name/ID\tscore (smaller is better)\tstrand1\tstrand2\ttype\tnumber_of_split_reads\n"
	std::string chr;
	for (size_t i = 0; i < SV.size(); i++) {
		std::string strands = SV[i]->get_strand(2);
		if (SV[i]->get_support() > Parameter::Instance()->min_support) {

			int pos = IPrinter::calc_pos(SV[i]->get_coordinates().start.min_pos, ref, chr);
			fprintf(file, "%s", chr.c_str());
			fprintf(file, "%c", '\t');
			fprintf(file, "%i", pos);
			fprintf(file, "%c", '\t');
			fprintf(file, "%i", IPrinter::calc_pos(SV[i]->get_coordinates().start.max_pos, ref, chr));
			fprintf(file, "%c", '\t');
			pos = IPrinter::calc_pos(SV[i]->get_coordinates().stop.min_pos, ref, chr);
			fprintf(file, "%s", chr.c_str());
			fprintf(file, "%c", '\t');
			fprintf(file, "%i", pos);
			fprintf(file, "%c", '\t');
			fprintf(file, "%i", IPrinter::calc_pos(SV[i]->get_coordinates().stop.max_pos, ref, chr));
			fprintf(file, "%c", '\t');
			fprintf(file, "%i", SV[i]->get_id());
			fprintf(file, "%c", '\t');
			fprintf(file, "%i", -1); //TODO: score
			fprintf(file, "%c", '\t');
			//std::cout<<IPrinter::get_type(SV[i]->get_SVtype()).c_str()<<std::endl;
			fprintf(file, "%c", strands[0]);
			fprintf(file, "%c", '\t');
			fprintf(file, "%c", strands[1]);
//			fprintf(file, "%c", '\t');
//			fprintf(file, "%c", strands[3]);
//			fprintf(file, "%c", '\t');
//			fprintf(file, "%c", strands[4]);
			fprintf(file, "%c", '\t');
			fprintf(file, "%s", IPrinter::get_type(SV[i]->get_SVtype()).c_str());
			fprintf(file, "%c", '\t');
			fprintf(file, "%i", SV[i]->get_support());
			fprintf(file, "%c", '\t');
			pos = IPrinter::calc_pos(SV[i]->get_coordinates().start.most_support, ref, chr);
			fprintf(file, "%s", chr.c_str());
			fprintf(file, "%c", '\t');
			fprintf(file, "%i", pos);
			fprintf(file, "%c", '\t');

			pos = IPrinter::calc_pos(SV[i]->get_coordinates().stop.most_support, ref, chr);
			fprintf(file, "%s", chr.c_str());
			fprintf(file, "%c", '\t');
			fprintf(file, "%i", pos);
			fprintf(file, "%c", '\n');
		}
	}
	fclose(file);
}
