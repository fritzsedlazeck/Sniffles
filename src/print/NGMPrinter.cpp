/*
 * NGMPrinter.cpp
 *
 *  Created on: Sep 23, 2015
 *      Author: fsedlaze
 */

/*
 * MariaPrinter.cpp
 *
 *  Created on: Sep 4, 2015
 *      Author: fsedlaze
 */

#include "NGMPrinter.h"

void NGMPrinter::print_header() {

}
void NGMPrinter::print_body(Breakpoint *& SV, RefVector ref) {
	//"Chrom\tstart\tstop\tchrom2\tstart2\tstop2\tvariant_name/ID\tscore (smaller is better)\tstrand1\tstrand2\ttype\tnumber_of_split_reads\n"

	if (!this->bed_tree.is_in(SV->get_coordinates().start.most_support, this->root) && !this->bed_tree.is_in(SV->get_coordinates().stop.most_support, this->root)) {
		std::string chr;
		std::string strands = SV->get_strand(2);
		if ((SV->get_SVtype() & TRA) || SV->get_length() > 1000000) { //1MB??
			int pos = IPrinter::calc_pos(SV->get_coordinates().start.min_pos, ref, chr) - Parameter::Instance()->max_dist;
			fprintf(file, "%s", chr.c_str());
			fprintf(file, "%c", '\t');
			if (pos > 0) {
				fprintf(file, "%i", pos);
			} else {
				fprintf(file, "%i", 0);
			}
			fprintf(file, "%c", '\t');
			pos = IPrinter::calc_pos(SV->get_coordinates().start.max_pos, ref, chr) + Parameter::Instance()->max_dist;
			fprintf(file, "%i", pos);
			fprintf(file, "%c", '\n');

			pos = IPrinter::calc_pos(SV->get_coordinates().stop.min_pos, ref, chr) - Parameter::Instance()->max_dist;
			fprintf(file, "%s", chr.c_str());
			fprintf(file, "%c", '\t');
			if (pos > 0) {
				fprintf(file, "%i", pos);
			} else {
				fprintf(file, "%i", 0);
			}
			fprintf(file, "%c", '\t');
			pos = IPrinter::calc_pos(SV->get_coordinates().stop.max_pos, ref, chr) - Parameter::Instance()->max_dist;
			fprintf(file, "%i", pos);
			fprintf(file, "%c", '\n');
		} else { //smaller SV:

			int pos = IPrinter::calc_pos(SV->get_coordinates().start.min_pos, ref, chr) - Parameter::Instance()->max_dist;
			fprintf(file, "%s", chr.c_str());
			fprintf(file, "%c", '\t');
			if (pos > 0) {
				fprintf(file, "%i", pos);
			} else {
				fprintf(file, "%i", 0);
			}
			fprintf(file, "%c", '\t');
			pos = IPrinter::calc_pos(SV->get_coordinates().stop.max_pos, ref, chr) + Parameter::Instance()->max_dist;
			fprintf(file, "%i", pos);
			fprintf(file, "%c", '\n');

		}

	}
}
