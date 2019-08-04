/*
 * Force_calling.cpp
 *
 *  Created on: Aug 24, 2017
 *      Author: sedlazec
 */

#include "Force_calling.h"

char assign_type(short type) {

	switch (type) {
	case 0: //DEL
		return DEL;
	case 1: //DUP
		return DUP;
	case 2: //INV
		return INV;
	case 3: //TRA
		return TRA;
	case 4: //INS
		return INS;
	case 6:
		return NEST;
	}
	return ' '; //TODO check default. Should not happen!
}
void fill_tree(IntervallTree & final, TNode *& root_final, RefVector ref, std::map<std::string, long>& ref_lens) {
	//prepare lookup:

	long length = 0;
	for (size_t i = 0; i < ref.size(); i++) {
		ref_lens[ref[i].RefName.c_str()] = length;
		length += (long) ref[i].RefLength + (long) Parameter::Instance()->max_dist;
	}

	//sometimes the stop coordinates are off especially for smaller chrs!??

	//parse VCF file
	std::vector<strvcfentry> entries = parse_vcf(Parameter::Instance()->input_vcf, 0);
	std::cout << "\t\t" << entries.size() << " SVs found in input." << std::endl;
	int invalid_svs = 0;
	for (size_t i = 0; i < entries.size(); i++) {
		if (entries[i].type != -1) {
			position_str svs;

			if (ref_lens.find(entries[i].start.chr) == ref_lens.end()) { // check why this is not called!
				cerr << "Error undefined CHR in VCF vs. BAM header: " << entries[i].start.chr << endl;
				exit(EXIT_FAILURE);
			}
			if (ref_lens.find(entries[i].stop.chr) == ref_lens.end()) {
				cerr << "Error undefined CHR in VCF vs. BAM header: " << entries[i].stop.chr << endl;
				exit(EXIT_FAILURE);
			}
			svs.start.min_pos = (long) entries[i].start.pos + ref_lens[entries[i].start.chr];
			svs.stop.max_pos = (long) entries[i].stop.pos + ref_lens[entries[i].stop.chr];
			read_str read;

			read.coordinates.first = (long) entries[i].start.pos + ref_lens[entries[i].start.chr];
			read.coordinates.second = (long) entries[i].stop.pos + ref_lens[entries[i].stop.chr];
			if (entries[i].type == 4) { //ins?
				if (entries[i].sv_len == Parameter::Instance()->huge_ins) {
					entries[i].sv_len++; // bad hack!
				}

				svs.stop.max_pos += (long) entries[i].sv_len;
				read.coordinates.second += (long) entries[i].sv_len;
				//	cout << "Parse: " << entries[i].start.pos << " " << entries[i].stop.pos << " " << svs.start.min_pos  <<" "<<svs.stop.max_pos  << endl;
			}

			read.SV = assign_type(entries[i].type);
			read.strand = entries[i].strands;
			read.type = 2; //called
			read.length = entries[i].sv_len; //svs.stop.max_pos-svs.start.min_pos;//try
			svs.support["input"] = read;
			//	cout<<"Submit: "<<entries[i].type <<endl;
			Breakpoint * br = new Breakpoint(svs, (long) entries[i].sv_len, read.SV);
			final.insert(br, root_final);
		} else {
			invalid_svs++;
		}
	}
	cerr << "\t\tInvalid types found skipping " << invalid_svs << " entries." << endl;
	//std::cout << "Print:" << std::endl;
//	final.print(root_final);
	entries.clear();
	//exit(0);
}

void force_calling(std::string bam_file, IPrinter *& printer) {
	cout << "Force calling SVs resetting parameters" << endl;
	Parameter::Instance()->min_mq=0;

	//parse reads
	//only process reads overlapping SV
	estimate_parameters(Parameter::Instance()->bam_files[0]);
	BamParser * mapped_file = 0;
	RefVector ref;
	std::string read_filename = Parameter::Instance()->bam_files[0];
	if (read_filename.find("bam") != string::npos) {
		mapped_file = new BamParser(read_filename);
		ref = mapped_file->get_refInfo();
	} else {
		cerr << "File Format not recognized. File must be a sorted .bam file!" << endl;
		exit(EXIT_FAILURE);
	}
	std::cout <<"\tConstruct Tree..." << std::endl;

	//construct the tree:
	IntervallTree final;
	TNode * root_final = NULL;
	std::map<std::string, long> ref_lens;
	fill_tree(final, root_final, ref, ref_lens);

	int current_RefID = 0;
	std::cout << "Start parsing: Chr " << ref[current_RefID].RefName << std::endl;

	//FILE * alt_allel_reads;
	FILE * ref_allel_reads;
	if (Parameter::Instance()->genotype) {
		ref_allel_reads = fopen(Parameter::Instance()->tmp_genotyp.c_str(), "w");
//		//ref_allel_reads = fopen(Parameter::Instance()->tmp_genotyp.c_str(), "wb");
	}
	Alignment * tmp_aln = mapped_file->parseRead(Parameter::Instance()->min_mq);

	long ref_space = ref_lens[ref[tmp_aln->getRefID()].RefName];
	long num_reads = 0;
	while (!tmp_aln->getQueryBases().empty()) {
		if ((tmp_aln->getAlignment()->IsPrimaryAlignment()) && (!(tmp_aln->getAlignment()->AlignmentFlag & 0x800) && tmp_aln->get_is_save())) {
			//change CHR:
			if (current_RefID != tmp_aln->getRefID()) {
				current_RefID = tmp_aln->getRefID();
				ref_space = ref_lens[ref[tmp_aln->getRefID()].RefName];
				std::cout << "\tSwitch Chr " << ref[tmp_aln->getRefID()].RefName << std::endl;				//" " << ref[tmp_aln->getRefID()].RefLength
			}

			//check if overlap with any breakpoint!!
			long read_start_pos = (long) tmp_aln->getPosition() - (long) Parameter::Instance()->max_dist;
			read_start_pos += ref_space;
			long read_stop_pos = read_start_pos + (long) tmp_aln->getAlignment()->Length + (long) Parameter::Instance()->max_dist;	//getRefLength();//(long) tmp_aln->getPosition();

		//	cout<<"Check overlap: "<<read_start_pos<<" "<<read_stop_pos<<endl;;
			if (final.overlaps(read_start_pos, read_stop_pos, root_final)) {
			//	cout<<" found "<<endl;
				//SCAN read:
				std::vector<str_event> aln_event;
				std::vector<aln_str> split_events;
				if (tmp_aln->getMappingQual() > Parameter::Instance()->min_mq) {
					double score = tmp_aln->get_scrore_ratio();
#pragma omp parallel // starts a new team
					{
#pragma omp sections
						{
							{
								//	clock_t begin = clock();
								if ((score == -1 || score > Parameter::Instance()->score_treshold)) {
									aln_event = tmp_aln->get_events_Aln();
								}
								//	Parameter::Instance()->meassure_time(begin, " Alignment ");
							}
#pragma omp section
							{
								//		clock_t begin_split = clock();
								split_events = tmp_aln->getSA(ref);
								//		Parameter::Instance()->meassure_time(begin_split," Split reads ");
							}
						}
					}
					//tmp_aln->set_supports_SV(aln_event.empty() && split_events.empty());

					//Store reference supporting reads for genotype estimation:
				//	str_read tmp;
					if ((Parameter::Instance()->genotype && (aln_event.empty() && split_events.empty()))){//}&& (score == -1 || score > Parameter::Instance()->score_treshold)))) {
						//write read:
						write_read(tmp_aln, ref_allel_reads);
					}

					//store the potential SVs:
					if (!aln_event.empty()) {
					//	cout<<"\t adding aln: "<<endl;
						add_events(tmp_aln, aln_event, 0, ref_space, final, root_final, num_reads, true);
					}
					if (!split_events.empty()) {
					//	cout<<"\t adding split: "<<endl;
						add_splits(tmp_aln, split_events, 1, ref, final, root_final, num_reads, true);
					}
				}
			}
		//	cout<<" none "<<endl;
		}
		//get next read:
		mapped_file->parseReadFast(Parameter::Instance()->min_mq, tmp_aln);

		num_reads++;

		if (num_reads % 10000 == 0) {
			cout << "\t\t# Processed reads: " << num_reads << endl;
		}
	}

	//std::cout << "Print:" << std::endl;
	//final.print(root_final);

	//filter and copy results:
	std::cout << "Finalizing  .." << std::endl;

	if (Parameter::Instance()->genotype) {
		fclose(ref_allel_reads);
	}
	//	sweep->finalyze();

	std::vector<Breakpoint*> points;
	final.get_breakpoints(root_final, points);

	//std::cout<<"fin up"<<std::endl;
	for (size_t i = 0; i < points.size(); i++) {
		points[i]->calc_support();
		points[i]->predict_SV();
		printer->printSV(points[i]); //redo! Ignore min support + STD etc.
	}

}
