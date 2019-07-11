/*
 * VCFPrinter.cpp
 *
 *  Created on: Aug 24, 2015
 *      Author: fsedlaze
 */

#include "VCFPrinter.h"

void VCFPrinter::print_header() {
	fprintf(file, "%s", "##fileformat=VCFv4.1\n");
	fprintf(file, "%s", "##source=Sniffles\n");
	string time = currentDateTime();
	fprintf(file, "%s", "##fileDate=");
	fprintf(file, "%s", time.c_str());

	//REport over all chrs:
	for (size_t i = 0; i < this->ref.size(); i++) {
		fprintf(file, "%s", "\n");
		fprintf(file, "%s", "##contig=<ID=");
		fprintf(file, "%s", ref[i].RefName.c_str());
		fprintf(file, "%s", ",length=");
		fprintf(file, "%i", (int) ref[i].RefLength);
		fprintf(file, "%c", '>');
	}

	fprintf(file, "%s", "\n");
	fprintf(file, "%s", "##ALT=<ID=DEL,Description=\"Deletion\">\n");
	fprintf(file, "%s", "##ALT=<ID=DUP,Description=\"Duplication\">\n");
	fprintf(file, "%s", "##ALT=<ID=INV,Description=\"Inversion\">\n");
	fprintf(file, "%s", "##ALT=<ID=INVDUP,Description=\"InvertedDUP with unknown boundaries\">\n");
	fprintf(file, "%s", "##ALT=<ID=TRA,Description=\"Translocation\">\n");
	fprintf(file, "%s", "##ALT=<ID=INS,Description=\"Insertion\">\n");
	fprintf(file, "%s", "##FILTER=<ID=UNRESOLVED,Description=\"An insertion that is longer than the read and thus we cannot predict the full size.\">\n");
	fprintf(file, "%s", "##INFO=<ID=CHR2,Number=1,Type=String,Description=\"Chromosome for END coordinate in case of a translocation\">\n");
	fprintf(file, "%s", "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the structural variant\">\n");
	fprintf(file, "%s", "##INFO=<ID=MAPQ,Number=1,Type=Integer,Description=\"Median mapping quality of paired-ends\">\n");
	fprintf(file, "%s", "##INFO=<ID=RE,Number=1,Type=Integer,Description=\"read support\">\n");
	fprintf(file, "%s", "##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description=\"Imprecise structural variation\">\n");
	fprintf(file, "%s", "##INFO=<ID=PRECISE,Number=0,Type=Flag,Description=\"Precise structural variation\">\n");

	//##FILTER=<ID=LowQual,Description="PE/SR support below 3 or mapping quality below 20.">
	fprintf(file, "%s", "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Length of the SV\">\n");
	fprintf(file, "%s", "##INFO=<ID=REF_strand,Number=2,Type=Integer,Description=\"Length of the SV\">\n");
	fprintf(file, "%s", "##INFO=<ID=SVMETHOD,Number=1,Type=String,Description=\"Type of approach used to detect SV\">\n");
	fprintf(file, "%s", "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n");
	if (Parameter::Instance()->report_n_reads > 0 || Parameter::Instance()->report_n_reads == -1) {
		fprintf(file, "%s", "##INFO=<ID=RNAMES,Number=.,Type=String,Description=\"Names of reads supporting SVs (comma separated)\">\n");
	}
	if (Parameter::Instance()->print_seq) {
		fprintf(file, "%s", "##INFO=<ID=SEQ,Number=1,Type=String,Description=\"Extracted sequence from the best representative read.\">\n");
	}

	if (Parameter::Instance()->read_strand) {
		fprintf(file, "%s", "##INFO=<ID=STRANDS2,Number=4,Type=Integer,Description=\"alt reads first + ,alt reads first -,alt reads second + ,alt reads second -.\">\n");
		fprintf(file, "%s", "##INFO=<ID=REF_strand,Number=2,Type=Integer,Description=\"plus strand ref, minus strand ref.\">\n");
	}

	fprintf(file, "%s", "##INFO=<ID=STD_quant_start,Number=A,Type=Float,Description=\"STD of the start breakpoints across the reads.\">\n");
	fprintf(file, "%s", "##INFO=<ID=STD_quant_stop,Number=A,Type=Float,Description=\"STD of the stop breakpoints across the reads.\">\n");
	fprintf(file, "%s", "##INFO=<ID=Kurtosis_quant_start,Number=A,Type=Float,Description=\"Kurtosis value of the start breakpoints across the reads.\">\n");
	fprintf(file, "%s", "##INFO=<ID=Kurtosis_quant_stop,Number=A,Type=Float,Description=\"Kurtosis value of the stop breakpoints across the reads.\">\n");
	fprintf(file, "%s", "##INFO=<ID=SUPTYPE,Number=A,Type=String,Description=\"Type by which the variant is supported.(SR,ALN,NR)\">\n");
	fprintf(file, "%s", "##INFO=<ID=STRANDS,Number=A,Type=String,Description=\"Strand orientation of the adjacency in BEDPE format (DEL:+-, DUP:-+, INV:++/--)\">\n");
	fprintf(file, "%s", "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency.\">\n");
	fprintf(file, "%s", "##INFO=<ID=ZMW,Number=A,Type=Integer,Description=\"Number of ZMWs (Pacbio) supporting SV.\">\n");
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

map<std::string, vector<int> > init_motives2() {
	map<std::string, vector<int> > motives;
	motives["TGAA"].push_back(0); // = 0;
	motives["ATCT"].push_back(0); //  = 0; //rev comp

	motives["TCTA"].push_back(0); //  = 0;
	motives["GGAA"].push_back(0); //  = 0;
	motives["GGCA"].push_back(0); //  = 0;
	motives["TCTA"].push_back(0); //  = 0;
	motives["AAGG"].push_back(0); //  = 0;
	motives["AGAA"].push_back(0); //  = 0;
	motives["AGAT"].push_back(0); //  = 0;
	motives["CTAT"].push_back(0); //  = 0;
	motives["TCTA"].push_back(0); //  = 0;
	motives["TGAA"].push_back(0); //  = 0;
	motives["GATA"].push_back(0); //  = 0;
	motives["GACA"].push_back(0); //  = 0;
	motives["TAGA"].push_back(0); //  = 0;
	motives["CAGA"].push_back(0); //  = 0;
	motives["TATC"].push_back(0); //  = 0;
	motives["TTTTC"].push_back(0); //  = 0;
	motives["GATA"].push_back(0); //  = 0;
	motives["AGAA"].push_back(0); //  = 0;
	motives["TCCT"].push_back(0); //  = 0;
	motives["TATC"].push_back(0); //  = 0;
	motives["AAAGA"].push_back(0); //  = 0;
	motives["ATT"].push_back(0); //  = 0;

	motives["TAGA"].push_back(0); //  = 0;
	motives["GAAA"].push_back(0); //  = 0;
	motives["TCTG"].push_back(0); //  = 0;
	motives["TAT"].push_back(0); //  = 0;
	motives["AGAT"].push_back(0); //  = 0;
	motives["TGGA"].push_back(0); //  = 0;
	motives["TTTTC"].push_back(0); //  = 0;
	motives["GATA"].push_back(0); //  = 0;
	motives["AGAGAT"].push_back(0); //  = 0;
	motives["AGAT"].push_back(0); //  = 0;
	motives["GAAA"].push_back(0); //  = 0;
	motives["CTT"].push_back(0); //  = 0;
	motives["ATCT"].push_back(0); //  = 0;
	motives["GATA"].push_back(0); //  = 0;
	motives["TTTC"].push_back(0); //  = 0;
	motives["AAAG"].push_back(0); //  = 0;
	motives["CTTTT"].push_back(0); //  = 0;

	return motives;
}

map<std::string, int> init_motives() {
	map<std::string, int> motives;
	motives["TGAA"] = 0;
	motives["ATCT"] = 0; //rev comp

	motives["TCTA"] = 0;
	motives["GGAA"] = 0;
	motives["GGCA"] = 0;
	motives["TCTA"] = 0;
	motives["AAGG"] = 0;
	motives["AGAA"] = 0;
	motives["AGAT"] = 0;
	motives["CTAT"] = 0;
	motives["TCTA"] = 0;
	motives["TGAA"] = 0;
	motives["GATA"] = 0;
	motives["GACA"] = 0;
	motives["TAGA"] = 0;
	motives["CAGA"] = 0;
	motives["TATC"] = 0;
	motives["TTTTC"] = 0;
	motives["GATA"] = 0;
	motives["AGAA"] = 0;
	motives["TCCT"] = 0;
	motives["TATC"] = 0;
	motives["AAAGA"] = 0;
	motives["ATT"] = 0;

	motives["TAGA"]= 0;
	motives["GAAA"]= 0;
	motives["TCTA"]= 0;
	motives["TCTG"]= 0;
	motives["TAT"]= 0;
	motives["AGAT"]= 0;
	motives["TGGA"]= 0;
	motives["AGAGAT"]= 0;
	motives["AGAT"]= 0;
	motives["GAAA"]= 0;
	motives["CTT"]= 0;
	motives["ATCT"]= 0;
	motives["GATA"]= 0;
	motives["TTTC"]= 0;
	motives["AAAG"]= 0;
	motives["CTTTT"]= 0;
	return motives;
}

void VCFPrinter::report_STR(Breakpoint * &SV, RefVector ref) {

	//STR code!:::::
	//=============

	if (Parameter::Instance()->str && ((SV->get_SVtype() & INS) || (SV->get_SVtype() & DEL))) {
		map<std::string, std::vector<int> > motives2;// = init_motives2();
		std::string chr;
		int start = IPrinter::calc_pos(SV->get_coordinates().start.most_support, ref, chr);
		cout << "NEW REGION: " << chr << ":" << start << " " << IPrinter::get_type(SV->get_SVtype()) << endl;
		std::map<std::string, read_str> support = SV->get_coordinates().support;
		for (std::map<std::string, read_str>::iterator i = support.begin(); i != support.end(); i++) {
			if (abs(SV->get_coordinates().start.most_support - (*i).second.coordinates.first) < 10) {
				map<std::string, int> motives;// = init_motives();
				int counts = 0;
				std::string sequence = (*i).second.sequence;
				//cout<<"check seq"<<endl;
				for (size_t p = 0; p < sequence.size(); p++) {
				//	if (motives.find(sequence.substr(p, 3)) != motives.end()) {
						motives[sequence.substr(p, 3)]++;
				//	}
				//	if (motives.find(sequence.substr(p, 4)) != motives.end()) {
						motives[sequence.substr(p, 4)]++;
				//	}
				//	if (motives.find(sequence.substr(p, 5)) != motives.end()) {
						motives[sequence.substr(p, 5)]++;
				//	}
				}
				//cout<<"Summarize"<<endl;
				for (map<std::string, int>::iterator p = motives.begin(); p != motives.end(); p++) {
					if ((*p).second > 0) {
						while ((*p).second + 1 > motives2[(*p).first].size()) {
							motives2[(*p).first].push_back(0);
						}
						motives2[(*p).first][(*p).second]++;
					}
				}
			}
		}
		//cout<<"Print"<<endl;
		std::stringstream ss2;
		int num_entries = 0;
		for (map<std::string, vector<int> >::iterator p = motives2.begin(); p != motives2.end(); p++) {
			//if ((*p).second.size() ) {
			int max = 0;
			std::stringstream ss;
			ss << (*p).first << ":";
			for (size_t i = 0; i < (*p).second.size(); i++) {
				ss << (*p).second[i] << ";";
				if ((*p).second[i] > max) {
					max = (*p).second[i];
				}
			}
			if (max > 10) {
				cout << ss.str() << endl;
			}
		}
	}

}

void VCFPrinter::print_body(Breakpoint * &SV, RefVector ref) {

	if (!this->bed_tree.is_in(SV->get_coordinates().start.most_support, this->root) && !this->bed_tree.is_in(SV->get_coordinates().stop.most_support, this->root)) {
		//temp. store read names supporting this SVs to later group the SVs together.
		double std_quant_start = 0;
		double std_quant_stop = 0;

		pair<double, double> kurtosis;
		pair<double, double> std_quant;
		double std_length = 0;
		int zmws = 0;
		bool ok_to_print = (to_print(SV, std_quant, kurtosis, std_length, zmws) || Parameter::Instance()->ignore_std);
		if (Parameter::Instance()->str) {
			report_STR(SV, ref);
		}
		//std::cout << "Print check: " << std_quant.first << " " << std_quant.second << endl;
		if (ok_to_print && (zmws == 0 || zmws >= Parameter::Instance()->min_zmw)) {
			if (Parameter::Instance()->phase) {
				store_readnames(SV->get_read_ids(), id);
			}
			std::string chr;
			int start = IPrinter::calc_pos(SV->get_coordinates().start.most_support, ref, chr);
			fprintf(file, "%s", chr.c_str());
			fprintf(file, "%c", '\t');
			if (start < 1) {
				start = 1;
			}
			fprintf(file, "%i", start);
			fprintf(file, "%c", '\t');
			fprintf(file, "%i", id);
			id++;

			long end_coord = SV->get_coordinates().stop.most_support;
			if (((SV->get_SVtype() & INS))) { // && SV->get_length() == Parameter::Instance()->huge_ins) && SV->get_types().is_ALN) {
				end_coord = std::max((SV->get_coordinates().stop.most_support - (long) SV->get_length()), (long) start);
			}

			int end = IPrinter::calc_pos(end_coord, ref, chr);
			if (end < 1) {
				end = 1;
			}
			std::string strands = SV->get_strand(1);

			if (Parameter::Instance()->reportBND && (SV->get_SVtype() & TRA)) {
				//N[22:36765684[ +-
				//]21:10540232]N -+
				fprintf(file, "%s", "\tN\t");
				if (strands[0] == '-') { //&&
					fprintf(file, "%s", "]");
					fprintf(file, "%s", chr.c_str());
					fprintf(file, "%c", ':');
					fprintf(file, "%i", end);
					fprintf(file, "%s", "]N");

				} else {
					fprintf(file, "%s", "N[");
					fprintf(file, "%s", chr.c_str());
					fprintf(file, "%c", ':');
					fprintf(file, "%i", end);
					fprintf(file, "%c", '[');
				}
			} else if (!SV->get_sequence().empty() && ((SV->get_SVtype() & INS) || (SV->get_SVtype() & DEL))) {
				fprintf(file, "%c", '\t');
				if ((SV->get_SVtype() & DEL)) {
					fprintf(file, "%s", SV->get_sequence().c_str());
				} else {
					fprintf(file, "%c", 'N');
				}
				fprintf(file, "%c", '\t');
				if ((SV->get_SVtype() & INS)) {
					fprintf(file, "%s", SV->get_sequence().c_str());
				} else {
					fprintf(file, "%c", 'N');
				}

			} else {
				fprintf(file, "%s", "\tN\t");
				fprintf(file, "%c", '<');
				fprintf(file, "%s", IPrinter::get_type(SV->get_SVtype()).c_str());
				fprintf(file, "%c", '>');
			}

			if (((SV->get_SVtype() & INS) && SV->get_length() == Parameter::Instance()->huge_ins) && SV->get_types().is_ALN) {
				fprintf(file, "%s", "\t.\tUNRESOLVED\t");
			} else {
				fprintf(file, "%s", "\t.\tPASS\t");
			}

			if (std_quant.first < 10 && std_quant.second < 10) {
				fprintf(file, "%s", "PRECISE");
			} else {
				fprintf(file, "%s", "IMPRECISE");
			}

			fprintf(file, "%s", ";SVMETHOD=Snifflesv");
			fprintf(file, "%s", Parameter::Instance()->version.c_str());

			if (!(Parameter::Instance()->reportBND && (SV->get_SVtype() & TRA))) {

				fprintf(file, "%s", ";CHR2=");
				fprintf(file, "%s", chr.c_str());
				fprintf(file, "%s", ";END=");

				if (SV->get_SVtype() & INS) {
					fprintf(file, "%i", std::max((int) end, start));
				} else {

					fprintf(file, "%i", end);
				}
			}
			if (zmws != 0) {
				fprintf(file, "%s", ";ZMW=");
				fprintf(file, "%i", zmws);
			}
			fprintf(file, "%s", ";STD_quant_start=");
			fprintf(file, "%f", std_quant.first);
			fprintf(file, "%s", ";STD_quant_stop=");
			fprintf(file, "%f", std_quant.second);
			fprintf(file, "%s", ";Kurtosis_quant_start=");
			fprintf(file, "%f", kurtosis.first);
			fprintf(file, "%s", ";Kurtosis_quant_stop=");
			fprintf(file, "%f", kurtosis.second);

			fprintf(file, "%s", ";SVTYPE=");
			if (Parameter::Instance()->reportBND && (SV->get_SVtype() & TRA)) {
				fprintf(file, "%s", "BND");
			} else {
				fprintf(file, "%s", IPrinter::get_type(SV->get_SVtype()).c_str());
			}

			if (Parameter::Instance()->report_n_reads > 0 || Parameter::Instance()->report_n_reads == -1) {
				fprintf(file, "%s", ";RNAMES=");
				fprintf(file, "%s", SV->get_read_names().c_str());
			}
			fprintf(file, "%s", ";SUPTYPE=");
			fprintf(file, "%s", SV->get_supporting_types().c_str());
			fprintf(file, "%s", ";SVLEN=");

			if (((SV->get_SVtype() & INS) && SV->get_length() == Parameter::Instance()->huge_ins) && SV->get_types().is_ALN) {
				if (SV->get_sequence().size() != 0) { //!
					fprintf(file, "%i", SV->get_sequence().size());
				} else {
					fprintf(file, "%i", 1);
				}
			} else if (SV->get_SVtype() & TRA) {
				fprintf(file, "%i", 1);
			} else if (SV->get_SVtype() & DEL) {
				fprintf(file, "%i", SV->get_length() * -1);
			} else {
				fprintf(file, "%i", SV->get_length());
			}
			//	}
			fprintf(file, "%s", ";STRANDS=");
			fprintf(file, "%s", strands.c_str());
			if (Parameter::Instance()->read_strand) {
				fprintf(file, "%s", ";STRANDS2=");
				std::map<std::string, read_str> support = SV->get_coordinates().support;
				pair<int, int> tmp_start;
				pair<int, int> tmp_stop;
				tmp_start.first = 0;
				tmp_start.second = 0;
				tmp_stop.first = 0;
				tmp_stop.second = 0;
				for (std::map<std::string, read_str>::iterator i = support.begin(); i != support.end(); i++) {
					if ((*i).second.read_strand.first) {
						tmp_start.first++;
					} else {
						tmp_start.second++;
					}
					if ((*i).second.read_strand.second) {
						tmp_stop.first++;
					} else {
						tmp_stop.second++;
					}
				}
				fprintf(file, "%i", tmp_start.first);
				fprintf(file, "%s", ",");
				fprintf(file, "%i", tmp_start.second);
				fprintf(file, "%s", ",");
				fprintf(file, "%i", tmp_stop.first);
				fprintf(file, "%s", ",");
				fprintf(file, "%i", tmp_stop.second);
			}

			//	if (Parameter::Instance()->print_seq && !SV->get_sequence().empty()) {
			//		fprintf(file, "%s", ";SEQ=");
			//		fprintf(file, "%s", SV->get_sequence().c_str());
			//	}
			fprintf(file, "%s", ";RE=");
			fprintf(file, "%i", SV->get_support());
			//if(Parameter::Instance()->genotype){
			fprintf(file, "%s", "\tGT:DR:DV\t./.:.:");
			fprintf(file, "%i", SV->get_support());
			//}else{
			//	fprintf(file, "%s",this->assess_genotype(SV->get_refcount(),SV->get_support()).c_str());
			//}
			fprintf(file, "%c", '\n');
		}
	}
}

void VCFPrinter::print_body_recall(Breakpoint * &SV, RefVector ref) {
	if (Parameter::Instance()->phase) {
		store_readnames(SV->get_read_ids(), id);
	}

	pair<double, double> kurtosis;
	pair<double, double> std_quant;
	double std_length = 0;
	int zmws = 0;
	bool ok_to_print = to_print(SV, std_quant, kurtosis, std_length, zmws);

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
	if (Parameter::Instance()->reportBND && (SV->get_SVtype() & TRA)) {
		//N[22:36765684[ +-
		//]21:10540232]N -+
		if (strands[0] == '-' && strands[0] == '+') {
			fprintf(file, "%s", "]");
			fprintf(file, "%s", chr.c_str());
			fprintf(file, "%c", ':');
			fprintf(file, "%i", end);
			fprintf(file, "%s", "]N");

		} else {
			fprintf(file, "%s", "N[");
			fprintf(file, "%s", chr.c_str());
			fprintf(file, "%c", ':');
			fprintf(file, "%i", end);
			fprintf(file, "%c", '[');
		}
	} else {

		fprintf(file, "%c", '<');
		fprintf(file, "%s", IPrinter::get_type(SV->get_SVtype()).c_str());
		fprintf(file, "%c", '>');
	}

	fprintf(file, "%s", "\t.\tPASS\t");
	fprintf(file, "%s", "IMPRECISE");
	fprintf(file, "%s", ";SVMETHOD=Snifflesv");
	fprintf(file, "%s", Parameter::Instance()->version.c_str());
	if (!(Parameter::Instance()->reportBND && (SV->get_SVtype() & TRA))) {
		fprintf(file, "%s", ";CHR2=");
		fprintf(file, "%s", chr.c_str());
		fprintf(file, "%s", ";END=");

		if (SV->get_SVtype() & INS) {
			fprintf(file, "%i", std::max((int) (end - SV->get_length()), start));
		} else {
			fprintf(file, "%i", end);
		}
	}

	fprintf(file, "%s", ";SVTYPE=");
	if (Parameter::Instance()->reportBND && (SV->get_SVtype() & TRA)) {
		fprintf(file, "%s", "BND");
	} else {
		fprintf(file, "%s", IPrinter::get_type(SV->get_SVtype()).c_str());
	}
	if (Parameter::Instance()->report_n_reads > 0 || Parameter::Instance()->report_n_reads == -1) {
		fprintf(file, "%s", ";RNAMES=");
		fprintf(file, "%s", SV->get_read_names().c_str());
	}
	fprintf(file, "%s", ";SUPTYPE=");
	fprintf(file, "%s", SV->get_supporting_types().c_str());
	fprintf(file, "%s", ";SVLEN=");

	if (((SV->get_SVtype() & INS) && SV->get_length() == Parameter::Instance()->huge_ins) && !SV->get_types().is_SR) {
		fprintf(file, "%s", "NA");
	} else {
		fprintf(file, "%i", SV->get_length());
	}
	//	}
	fprintf(file, "%s", ";STRANDS=");
	fprintf(file, "%s", strands.c_str());
	fprintf(file, "%s", ";SEQ=");
	fprintf(file, "%s", SV->get_sequence().c_str());
	fprintf(file, "%s", ";RE=");
	fprintf(file, "%i", SV->get_support());
	fprintf(file, "%s", "\tGT:DR:DV\t./.:.:");
	fprintf(file, "%i", SV->get_support());
	fprintf(file, "%c", '\n');

}

