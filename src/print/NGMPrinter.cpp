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
void NGMPrinter::print_body(std::vector<Breakpoint *>& SV, RefVector ref) {
	FILE *file;
	file = fopen(Parameter::Instance()->output_vcf.c_str(), "w");
	//"Chrom\tstart\tstop\tchrom2\tstart2\tstop2\tvariant_name/ID\tscore (smaller is better)\tstrand1\tstrand2\ttype\tnumber_of_split_reads\n"
	std::string chr;
	for (size_t i = 0; i < SV.size(); i++) {
		std::string strands = SV[i]->get_strand(2);
		if (SV[i]->get_support() > Parameter::Instance()->min_support) {
			fprintf(file, "%s", IPrinter::get_type(SV[i]->get_SVtype()).c_str());
			fprintf(file, "%c", '\t');

			int pos = IPrinter::calc_pos(SV[i]->get_coordinates().start.most_support, ref, chr);
			fprintf(file, "%s", chr.c_str());
			fprintf(file, "%c", '\t');
			fprintf(file, "%i", pos);
			fprintf(file, "%c", '\t');

			pos = IPrinter::calc_pos(SV[i]->get_coordinates().stop.most_support, ref, chr);
			fprintf(file, "%s", chr.c_str());
			fprintf(file, "%c", '\t');
			fprintf(file, "%i", pos);

			std::map<std::string, read_str> support = SV[i]->get_coordinates().support;
			int num = 0;
			for (std::map<std::string, read_str>::iterator read = support.begin(); read != support.end() && num < Parameter::Instance()->report_n_reads; read++) {
				if (SV[i]->get_SVtype() & (*read).second.SV) {
					fprintf(file, "%c", '\t');
					fprintf(file, "%s", (*read).first.c_str());
					num++;
				}
			}
			fprintf(file, "%c", '\n');

		}
	}
	fclose(file);
}
