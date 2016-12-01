/*
 * VCFPrinter.cpp
 *
 *  Created on: Aug 24, 2015
 *      Author: fsedlaze
 */

#include "VCFPrinter.h"

void VCFPrinter::print_header() {

	fprintf(file, "%s", "##fileformat=VCFv4.2\n");
	fprintf(file, "%s", "##source=Sniffles\n");
	string time = currentDateTime();
	fprintf(file, "%s", "##fileDate=");
	fprintf(file, "%s", time.c_str());
	fprintf(file, "%s", "\n"); //TODO change!
	fprintf(file, "%s", "##ALT=<ID=DEL,Description=\"Deletion\">\n");
	fprintf(file, "%s", "##ALT=<ID=DUP,Description=\"Duplication\">\n");
	fprintf(file, "%s", "##ALT=<ID=INV,Description=\"Inversion\">\n");
	fprintf(file, "%s", "##ALT=<ID=TRA,Description=\"Translocation\">\n");
	fprintf(file, "%s", "##ALT=<ID=INS,Description=\"Insertion\">\n");
	fprintf(file, "%s", "##INFO=<ID=CHR2,Number=1,Type=String,Description=\"Chromosome for END coordinate in case of a translocation\">\n");
	fprintf(file, "%s", "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the structural variant\">\n");
	fprintf(file, "%s", "##INFO=<ID=MAPQ,Number=1,Type=Integer,Description=\"Median mapping quality of paired-ends\">\n");
	fprintf(file, "%s", "##INFO=<ID=RE,Number=1,Type=Integer,Description=\"read support\">\n");
	fprintf(file, "%s", "##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description=\"Imprecise structural variation\">\n");
	fprintf(file, "%s", "##INFO=<ID=PRECISE,Number=0,Type=Flag,Description=\"Precise structural variation\">\n");
	fprintf(file, "%s", "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Length of the SV\">\n");
	fprintf(file, "%s", "##INFO=<ID=SVMETHOD,Number=1,Type=String,Description=\"Type of approach used to detect SV\">\n");
	fprintf(file, "%s", "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n");
	fprintf(file, "%s", "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
	fprintf(file, "%s", "##FORMAT=<ID=DR,Number=1,Type=Integer,Description=\"# high-quality reference reads\">\n");
	fprintf(file, "%s", "##FORMAT=<ID=DV,Number=1,Type=Integer,Description=\"# high-quality variant reads\">\n");

	fprintf(file, "%s", "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
	for (size_t i = 0; i < Parameter::Instance()->bam_files.size(); i++) {
		fprintf(file, "%c", '\t');
		fprintf(file, "%s", Parameter::Instance()->bam_files[i].c_str());
	}
	fprintf(file, "%c", '\n');
}
void VCFPrinter::print_body(Breakpoint * &SV, RefVector ref) {
	if (!this->bed_tree.is_in(SV->get_coordinates().start.most_support, this->root) && !this->bed_tree.is_in(SV->get_coordinates().stop.most_support, this->root)) {
		//temp. store read names supporting this SVs to later group the SVs together.

		if (Parameter::Instance()->phase) {
			store_readnames(SV->get_read_ids(), id);
		}
		std::string chr;
		int start = IPrinter::calc_pos(SV->get_coordinates().start.most_support, ref, chr);
		fprintf(file, "%s", chr.c_str());
		fprintf(file, "%c", '\t');
		fprintf(file, "%i", start);
		fprintf(file, "%c", '\t');
		fprintf(file, "%i", id);
		id++;

		int end = IPrinter::calc_pos(SV->get_coordinates().stop.most_support, ref, chr);
		std::string strands = SV->get_strand(1);
		fprintf(file, "%s", "\tN\t");
		/*if (SV->get_SVtype() & TRA) {	//TODO check for INV!
			//N[22:36765684[ +-
			//]21:10540232]N -+
			if (strands[0] == '+') {
				fprintf(file, "%s", "N[");
				fprintf(file, "%s", chr.c_str());
				fprintf(file, "%c", ':');
				fprintf(file, "%i", pos);

			} else {
				fprintf(file, "%c", ']');
				fprintf(file, "%s", chr.c_str());
				fprintf(file, "%c", ':');
				fprintf(file, "%i", pos);
			}
			if (strands[1] == '+') {
				fprintf(file, "%c", '[');
			} else {
				fprintf(file, "%c", ']');
			}
			if (strands[0] == '-') {
				fprintf(file, "%c", 'N');
			}
			fprintf(file, "%s", "\t.\tPASS\tIMPRECISE;SVMETHOD=Snifflesv");
			fprintf(file, "%s", Parameter::Instance()->version.c_str());

		} else {*/
			fprintf(file, "%c", '<');
			fprintf(file, "%s", IPrinter::get_type(SV->get_SVtype()).c_str());
			fprintf(file, "%c", '>');
			fprintf(file, "%s", "\t.\tPASS\tIMPRECISE;SVMETHOD=Snifflesv");
			fprintf(file, "%s", Parameter::Instance()->version.c_str());
			fprintf(file, "%s", ";CHR2=");
			fprintf(file, "%s", chr.c_str());
			fprintf(file, "%s", ";END=");
			if (SV->get_SVtype() & INS) {
				fprintf(file, "%i", std::max((int)(end - SV->get_length()) , start));
			} else {
				fprintf(file, "%i", end);
			}

			double std_start=0;
			double std_stop=0;
//			comp_std(SV,std_start,std_stop);
//			fprintf(file, "%s", ";STD_start=");
//			fprintf(file, "%f",std_start);
//			fprintf(file, "%s", ";STD_end=");
//			fprintf(file, "%f",std_stop);
	//	}
		fprintf(file, "%s", ";SVTYPE=");
		fprintf(file, "%s", IPrinter::get_type(SV->get_SVtype()).c_str());
		if (Parameter::Instance()->report_n_reads > 0) {
			std::vector<std::string> names = SV->get_read_names(Parameter::Instance()->report_n_reads);
			fprintf(file, "%s", ";RNAMES=");
			for (size_t i = 0; i < names.size() && i < Parameter::Instance()->report_n_reads; i++) {
				fprintf(file, "%s", names[i].c_str());
				fprintf(file, "%s", ";");
			}
		} else {
			fprintf(file, "%s", ";");
		}
		fprintf(file, "%s", "SUPTYPE=");
		fprintf(file, "%s", SV->get_supporting_types().c_str());
		fprintf(file, "%s", ";SVLEN=");
		fprintf(file, "%i", SV->get_length());
		fprintf(file, "%s", ";STRANDS=");
		fprintf(file, "%s", strands.c_str());
		fprintf(file, "%s", ";RE=");
		fprintf(file, "%i", SV->get_support());
		fprintf(file, "%s", "\tGT:DR:DV\t./.:.:");
		fprintf(file, "%i", SV->get_support());
		fprintf(file, "%c", '\n');
	}
}

