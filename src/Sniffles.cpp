//============================================================================
// Name        : Sniffles.cpp
// Author      : Fritz Sedlazeck
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================
// phil: cd ~/hetero/philipp/pacbio/example-svs/reads
//For mac: cmake -D CMAKE_C_COMPILER=/opt/local/bin/gcc-mp-4.7 -D CMAKE_CXX_COMPILER=/opt/local/bin/g++-mp-4.7 ..
#include <iostream>
#include "Paramer.h"
#include <tclap/CmdLine.h>
#include <omp.h>
#include "realign/Realign.h"
#include "sub/Detect_Breakpoints.h"
#include "print/IPrinter.h"
#include "print/BedePrinter.h"
#include "print/VCFPrinter.h"
#include "print/MariaPrinter.h"
#include "phasing/PhaserSV.h"
#include "print/NGMPrinter.h"

Parameter* Parameter::m_pInstance = NULL;
//cmake -D CMAKE_C_COMPILER=/opt/local/bin/gcc-mp-4.7 -D CMAKE_CXX_COMPILER=/opt/local/bin/g++-mp-4.7 ..

//TODO: Friday: redo the positioning. Min,Max as representation, but collect the position and give the most likely!

// TODO: Finish virtual Region class, Inherit Breakpoint and strange regions; Implement the region into the intervall tree.

// Think of method to filter out strange SV.
// flushing the interval tree after a certain amount of SV/ bp? Would speed up things.
// Allow for Illumina split read data?
// Think again of computing the coverage?
// VCF and bede output
// Think of multiple bam files -> setting genotypes
// Think about overlapping SV, maybe flag to report if they share the same read -> phasing info?
// Regular scan through the SV and move those where the end point lies far behind the current pos or reads. Eg. 1MB?

//TODO: look into the gene fusions.


void read_parameters(int argc, char *argv[]) {
	TCLAP::CmdLine cmd("Sniffles version 0.0.1", ' ', "0.0.1");

	TCLAP::ValueArg<std::string> arg_bamfile("m", "mapped_reads", "Bam File", true, "", "string");
	TCLAP::ValueArg<std::string> arg_vcf("v", "vcf", "VCF output file name", false, "", "string");
	//TCLAP::ValueArg<std::string> arg_bede("b", "bede", "Bede output file name", false, "", "string");
	TCLAP::ValueArg<std::string> arg_maria("", "bede", "Simplified format of bede Format.", false, "", "string");
	//TCLAP::ValueArg<std::string> arg_ref("r", "reference", "Reference fasta sequence for realign step", false, "", "string");
	//TCLAP::ValueArg<std::string> arg_region("", "regions", "List of regions CHR:start-stop; to check", false, "", "string");

	TCLAP::ValueArg<int> arg_support("s", "min_support", "Minimum number of reads that support a SV. Default: 10", false, 10, "int");
	TCLAP::ValueArg<int> arg_splits("", "max_num_splits", "Maximum number of splits per read to be still taken into account. Default: 4", false, 4, "int");
	TCLAP::ValueArg<int> arg_dist("d", "max_distance", "Maximum distance to group SV together. Default: 1kb", false, 1000, "int");
	TCLAP::ValueArg<int> arg_threads("t", "threads", "Number of threads to use. Default: 3", false, 3, "int");
	//TCLAP::ValueArg<int> arg_corridor("", "corridor", "Maximum size of corridor for realignment. Default: 2000", false, 2000, "int");
	TCLAP::ValueArg<int> arg_minlength("l", "min_length", "Minimum length of SV to be reported. Default:0", false, 0, "int");
	TCLAP::ValueArg<int> arg_mq("q", "minmapping_qual", "Minimum Mapping Quality. Default: 20", false, 20, "int");
	TCLAP::ValueArg<int> arg_cigar("c", "min_cigar_event", "Minimum Cigar Event (e.g. Insertion, deletion) to take into account. Default:50 ", false, 50, "int");
	TCLAP::ValueArg<int> arg_numreads("n", "num_reads_report", "Report up to N reads that support the SV. Default: 0", false, 0, "int");
//	TCLAP::ValueArg<int> arg_phase_minreads("", "min_reads_phase", "Minimum reads overlapping two SV to phase them together. Default: 1", false, 1, "int");

	//TCLAP::SwitchArg arg_realign("", "re-align", "Enables the realignment of reads at predicted SV sites. Leads to more accurate breakpoint predictions.", cmd, false);
	TCLAP::SwitchArg arg_MD_cigar("", "use_MD_Cigar", "Enables Sniffles to use the alignment information to screen for suspicious regions.", cmd, false);




	cmd.add(arg_numreads);
	cmd.add(arg_dist);
	//cmd.add(arg_bede);
	cmd.add(arg_threads);
	cmd.add(arg_cigar);
	cmd.add(arg_maria);
	cmd.add(arg_vcf);
	//cmd.add(arg_ref);
	//cmd.add(arg_corridor);
	//cmd.add(arg_region);
	cmd.add(arg_minlength);
	//cmd.add(arg_phase_minreads);
	cmd.add(arg_mq);
	cmd.add(arg_splits);
	cmd.add(arg_support);
	cmd.add(arg_bamfile);
	//parse cmd:
	cmd.parse(argc, argv);

	Parameter::Instance()->read_name = " " ;//"m141129_192924_00118_c100715122550000001823152704301506_s1_p0/69574/0_17658"; //TODO: just for debuging reasons!
	Parameter::Instance()->bam_files.push_back(arg_bamfile.getValue());
	Parameter::Instance()->min_mq = arg_mq.getValue();
	Parameter::Instance()->output_vcf = arg_vcf.getValue();
	Parameter::Instance()->min_cigar_event = arg_cigar.getValue();
	Parameter::Instance()->report_n_reads = arg_numreads.getValue();
	Parameter::Instance()->min_support = arg_support.getValue();
	Parameter::Instance()->max_splits = arg_splits.getValue();
	Parameter::Instance()->max_dist = arg_dist.getValue();
	//Parameter::Instance()->ref_seq = arg_ref.getValue();
	//Parameter::Instance()->realign = arg_realign.getValue();
	//Parameter::Instance()->corridor = arg_corridor.getValue();
//	Parameter::Instance()->output_bede = arg_bede.getValue();
	Parameter::Instance()->min_length = arg_minlength.getValue();
	//Parameter::Instance()->min_reads_phase = arg_phase_minreads.getValue();
	Parameter::Instance()->useMD_CIGAR = arg_MD_cigar.getValue();
	Parameter::Instance()->num_threads=arg_threads.getValue();
	Parameter::Instance()->output_maria=arg_maria.getValue();
}

int main(int argc, char *argv[]) {

	try {
		//init parameter and reads user defined parameter from command line.
		read_parameters(argc,argv);

		//init openmp:
		omp_set_dynamic(0);
		omp_set_num_threads(Parameter::Instance()->num_threads);

		//init printer:
		IPrinter * printer;
		if(Parameter::Instance()->realign){
			printer = new NGMPrinter();
		}else if (!Parameter::Instance()->output_vcf.empty()) {
			printer = new VCFPrinter();
			std::cout<<"VCF parser"<<std::endl;
		} else if (!Parameter::Instance()->output_bede.empty() ) {
			printer = new BedePrinter();
		}else if (!Parameter::Instance()->output_maria.empty()){
			printer =new MariaPrinter();
		} else {
			std::cerr << "Please specify an output file using -v or -b" << std::endl;
			return -1;
		}

		//TODO add support of multiple files!
		std::vector<Breakpoint *> final_SV;
		final_SV = detect_breakpoints(Parameter::Instance()->bam_files[0]);


		//realignment: Using NGM???
	/*if (Parameter::Instance()->realign) {
			if (!Parameter::Instance()->ref_seq.empty()) {
				Realigner * re = new Realigner();
				re->align(final_SV);
			} else {
				std::cerr << "You need to specify the reference sequence that was used." << std::endl;
			}
		}*/

		//grouping SV together:
		//PhaserSV * phaser= new PhaserSV();
		//phaser->phase(final_SV);
		//delete phaser;

		//printing results:
		printer->init();
		printer->printSV(final_SV);

	} catch (TCLAP::ArgException &e)  // catch any exceptions
	{
		std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
	}
	return 0;
}
/*
 * Input:
 * Regions to check
 *
 * Parameters
 *
 *
 */

/*
 * 1. Detect strange regions
 * 		Using MD, Cigar, Split reads
 *
 * 2. Extract reads from region (?) The main read holds the whole sequence for BWA-MEM
 *
 * 3. Realign regions (?)
 */
