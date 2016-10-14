/*
 * VCFPrinter.cpp
 *
 *  Created on: Aug 24, 2015
 *      Author: fsedlaze
 */

#include "VCFPrinter.h"

void VCFPrinter::print_header() {
	fprintf(file, "%s", "##fileformat=VCFv4.1\n");
	//string time = currentDateTime();
	fprintf(file, "%s", "##fileDate=");
//	fprintf(file, "%s", time.c_str());
	fprintf(file, "%s", "\n"); //TODO change!
	fprintf(file, "%s", "##ALT=<ID=DEL,Description=\"Deletion\">\n");
	fprintf(file, "%s", "##ALT=<ID=DUP,Description=\"Duplication\">\n");
	fprintf(file, "%s", "##ALT=<ID=INV,Description=\"Inversion\">\n");
	fprintf(file, "%s", "##ALT=<ID=TRA,Description=\"Translocation\">\n");
	fprintf(file, "%s", "##ALT=<ID=INS,Description=\"Insertion\">\n");
	//rintf(file, "%s", "##FILTER=<ID=LowQual,Description=\"PE support below 3 or mapping quality below 20.\">\n");
	//printf(file, "%s", "##INFO=<ID=CIEND,Number=2,Type=Integer,Description=\"PE confidence interval around END\">\n");
	//fprintf(file, "%s", "##INFO=<ID=CIPOS,Number=2,Type=Integer,Description=\"PE confidence interval around POS\">\n");
	fprintf(file, "%s", "##INFO=<ID=CHR2,Number=1,Type=String,Description=\"Chromosome for END coordinate in case of a translocation\">\n");
	fprintf(file, "%s", "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the structural variant\">\n");
//	fprintf(file, "%s", "##INFO=<ID=PE,Number=1,Type=Integer,Description=\"Paired-end support of the structural variant\">\n");
	fprintf(file, "%s", "##INFO=<ID=MAPQ,Number=1,Type=Integer,Description=\"Median mapping quality of paired-ends\">\n");
	fprintf(file, "%s", "##INFO=<ID=RE,Number=1,Type=Integer,Description=\"read support\">\n");
//	fprintf(file, "%s", "##INFO=<ID=SRQ,Number=1,Type=Float,Description=\"Split-read consensus alignment quality\">\n");
//	fprintf(file, "%s", "##INFO=<ID=CONSENSUS,Number=1,Type=String,Description=\"Split-read consensus sequence\">\n");
//	fprintf(file, "%s", "##INFO=<ID=CT,Number=1,Type=String,Description=\"Paired-end signature induced connection type\">\n");
	fprintf(file, "%s", "##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description=\"Imprecise structural variation\">\n");
	fprintf(file, "%s", "##INFO=<ID=PRECISE,Number=0,Type=Flag,Description=\"Precise structural variation\">\n");
	fprintf(file, "%s", "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Length of the SV\">\n");
//	fprintf(file, "%s", "##INFO=<ID=SUBTYPE,Number=1,Type=String,Description=\"What supports the SV.\">\n");
//	fprintf(file, "%s", "##INFO=<SUBTYPE=SP,Description=\"SV supported by split reads\">\n");
//	fprintf(file, "%s", "##INFO=<SUBTYPE=MD,Description=\"SV supported by MD string\">\n");
//	fprintf(file, "%s", "##INFO=<SUBTYPE=CI,Description=\"SV supported by Cigar string\">\n");
	//fprintf(file, "%s", "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n");
	fprintf(file, "%s", "##INFO=<ID=SVMETHOD,Number=1,Type=String,Description=\"Type of approach used to detect SV\">\n");
	fprintf(file, "%s", "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
	//fprintf(file, "%s", "##FORMAT=<ID=GL,Number=G,Type=Float,Description=\"Log10-scaled genotype likelihoods for RR,RA,AA genotypes\">\n");
	//fprintf(file, "%s", "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">\n");
	//fprintf(file, "%s", "##FORMAT=<ID=FT,Number=1,Type=String,Description=\"Per-sample genotype filter\">\n");
//	fprintf(file, "%s", "##FORMAT=<ID=RC,Number=1,Type=Integer,Description=\"Normalized high-quality read count for the SV\">\n");
	//fprintf(file, "%s", "##FORMAT=<ID=DR,Number=1,Type=Integer,Description=\"# high-quality reference pairs\">\n");
	fprintf(file, "%s", "##FORMAT=<ID=DR,Number=1,Type=Integer,Description=\"# high-quality reference reads\">\n");
	fprintf(file, "%s", "##FORMAT=<ID=DV,Number=1,Type=Integer,Description=\"# high-quality variant reads\">\n");

	//fprintf(file, "%s", "##FORMAT=<ID=RR,Number=1,Type=Integer,Description=\"# high-quality reference junction reads\">\n");
	//fprintf(file, "%s", "##FORMAT=<ID=RV,Number=1,Type=Integer,Description=\"# high-quality variant junction reads\">\n");
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
			store_readnames(SV->get_read_ids(),id);
		}
		std::string chr;
		int pos = IPrinter::calc_pos(SV->get_coordinates().start.most_support, ref, chr);
		fprintf(file, "%s", chr.c_str());
		fprintf(file, "%c", '\t');
		fprintf(file, "%i", pos);
		fprintf(file, "%c", '\t');
		fprintf(file, "%i", id);
		id++;
		fprintf(file, "%s", "\tN\t<");
		fprintf(file, "%s", IPrinter::get_type(SV->get_SVtype()).c_str());
		fprintf(file, "%s", ">\t.\tPASS\tIMPRECISE;SVMETHOD=Snifflesv0.0.1;CHR2=");
		pos = IPrinter::calc_pos(SV->get_coordinates().stop.most_support, ref, chr);
		fprintf(file, "%s", chr.c_str());
		fprintf(file, "%s", ";END=");
		if (SV->get_SVtype() & INS) {
			fprintf(file, "%i", pos+SV->get_length());
		}else{
			fprintf(file, "%i", pos);
		}
		if (Parameter::Instance()->debug) {
			std::vector< std::string> names=SV->get_read_names(Parameter::Instance()->report_n_reads);
			fprintf(file, "%s", ";RNAMES=");
			for(size_t i=0;i<names.size()&& i<Parameter::Instance()->report_n_reads;i++){
				fprintf(file, "%s", names[i].c_str());
				fprintf(file, "%s", ";");
			}
		}else{
			fprintf(file, "%s", ";");
		}
		fprintf(file, "%s", "SUPTYPE=");
		fprintf(file, "%s", SV->get_supporting_types().c_str());
		fprintf(file, "%s", ";SVLEN=");
		fprintf(file, "%i", SV->get_length());
		fprintf(file, "%s", ";STRANDS=");
		fprintf(file, "%s", SV->get_strand(1).c_str());
		fprintf(file, "%s", ";RE=");
		fprintf(file, "%i", SV->get_support());
		fprintf(file, "%s", "\tGT:DR:DV\t./.:");
		fprintf(file, "%i", SV->get_support());
		fprintf(file, "%c", '\n');
	}
}

