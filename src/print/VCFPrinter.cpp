/*
 * VCFPrinter.cpp
 *
 *  Created on: Aug 24, 2015
 *      Author: fsedlaze
 */

#include "VCFPrinter.h"

void VCFPrinter::print_header() {
	FILE *file;
	file = fopen(Parameter::Instance()->output_vcf.c_str(), "w");
	fprintf(file, "%s", "##fileformat=VCFv4.1\n");
	fprintf(file, "%s", "##fileDate=20150221\n"); //TODO change!
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
	fprintf(file, "%s", "##INFO=<ID=SUPTYPE,Number=1,Type=String,Description=\"What supports the SV.\">\n");
	fprintf(file, "%s", "##INFO=<SUPTYPE=SP,Description=\"SV supported by split reads\">\n");
	fprintf(file, "%s", "##INFO=<SUPTYPE=MD,Description=\"SV supported by MD string\">\n");
	fprintf(file, "%s", "##INFO=<SUPTYPE=CI,Description=\"SV supported by Cigar string\">\n");
	//fprintf(file, "%s", "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n");
	fprintf(file, "%s", "##INFO=<ID=SVMETHOD,Number=1,Type=String,Description=\"Type of approach used to detect SV\">\n");
	fprintf(file, "%s", "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
	//fprintf(file, "%s", "##FORMAT=<ID=GL,Number=G,Type=Float,Description=\"Log10-scaled genotype likelihoods for RR,RA,AA genotypes\">\n");
	//fprintf(file, "%s", "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">\n");
	//fprintf(file, "%s", "##FORMAT=<ID=FT,Number=1,Type=String,Description=\"Per-sample genotype filter\">\n");
//	fprintf(file, "%s", "##FORMAT=<ID=RC,Number=1,Type=Integer,Description=\"Normalized high-quality read count for the SV\">\n");
	//fprintf(file, "%s", "##FORMAT=<ID=DR,Number=1,Type=Integer,Description=\"# high-quality reference pairs\">\n");
	fprintf(file, "%s", "##FORMAT=<ID=DV,Number=1,Type=Integer,Description=\"# high-quality variant pairs\">\n");
	//fprintf(file, "%s", "##FORMAT=<ID=RR,Number=1,Type=Integer,Description=\"# high-quality reference junction reads\">\n");
	//fprintf(file, "%s", "##FORMAT=<ID=RV,Number=1,Type=Integer,Description=\"# high-quality variant junction reads\">\n");
	fprintf(file, "%s", "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
	for (size_t i = 0; i < Parameter::Instance()->bam_files.size(); i++) {
		fprintf(file, "%c", '\t');
		fprintf(file, "%s", Parameter::Instance()->bam_files[i].c_str());
	}
	fprintf(file, "%c", '\n');
	fclose(file);
}
void VCFPrinter::print_body(std::vector<Breakpoint *> &SV, RefVector ref) {
	FILE *file;
	file = fopen(Parameter::Instance()->output_vcf.c_str(), "a");
	//MT      5289    DEL.MT:5289..6499       N       <DEL>   .       LowQual IMPRECISE;CIEND=0,0;CIPOS=0,0;SVTYPE=DEL;SVMETHOD=EMBL.DELLYv0.5.9;CHR2=MT;END=6499;SVLEN=1210;CT=3to5;PE=26;MAPQ=0     GT:GL:GQ:FT:RC:DR:DV:RR:RV      0/0:0,-1710.68,-34129.3:17107:PASS:0:5721:2:314:0
	std::string chr;
	for (size_t i = 0; i < SV.size(); i++) {
		int pos = IPrinter::calc_pos(SV[i]->get_coordinates().start.most_support, ref,chr);
		fprintf(file, "%s", chr.c_str());
		fprintf(file, "%c", '\t');


		fprintf(file, "%i", pos);
		fprintf(file, "%c", '\t');
		fprintf(file, "%i", SV[i]->get_id());
		fprintf(file, "%s", "\tN\t<");
		fprintf(file, "%s", IPrinter::get_type(SV[i]->get_SVtype()).c_str());
		fprintf(file, "%s", ">\t.\tPASS\tIMPRECISE;SVMETHOD=Snifflesv0.0.1;CHR2=");
		pos = IPrinter::calc_pos(SV[i]->get_coordinates().stop.most_support, ref,chr);
		fprintf(file, "%s", chr.c_str());
		fprintf(file, "%s", ";END=");
		fprintf(file, "%i", pos);
		fprintf(file, "%s", ";SUPTYPE=");
		fprintf(file, "%s", SV[i]->get_supporting_types().c_str());
		fprintf(file, "%s", ";SVLEN=");
		fprintf(file, "%i", SV[i]->get_length());
		fprintf(file, "%s", ";STRANDS=");
		fprintf(file, "%s", SV[i]->get_strand(2).c_str());
		fprintf(file, "%s", ";RE=");
		fprintf(file, "%i", SV[i]->get_support());
		fprintf(file, "%s", "\tGT:DV\t./.:");
		fprintf(file, "%i", SV[i]->get_support());
		fprintf(file, "%c", '\n');
	}
	fclose(file);
}

