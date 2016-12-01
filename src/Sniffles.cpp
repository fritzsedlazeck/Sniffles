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
#include "Genotyper/Genotyper.h"
#include "realign/Realign.h"
#include "sub/Detect_Breakpoints.h"
#include "print/IPrinter.h"
#include "print/VCFPrinter.h"
#include "print/BedpePrinter.h"
#include "print/NGMPrinter.h"
#include "Ignore_Regions.h"
#include "plane-sweep/PlaneSweep_slim.h"
#include "print/BedpePrinter.h"

Parameter* Parameter::m_pInstance = NULL;
//cmake -D CMAKE_C_COMPILER=/opt/local/bin/gcc-mp-4.7 -D CMAKE_CXX_COMPILER=/opt/local/bin/g++-mp-4.7 ..

// Think of method to filter out strange SV.
// Think of multiple bam files -> setting genotypes
// Think about overlapping SV, maybe flag to report if they share the same read -> phasing info?
// Regular scan through the SV and move those where the end point lies far behind the current pos or reads. Eg. 1MB?


void read_parameters(int argc, char *argv[]) {

	std::string lable="Sniffles version ";
	lable+=Parameter::Instance()->version;
	TCLAP::CmdLine cmd("Sniffles version 1.0.0", ' ', Parameter::Instance()->version);
	TCLAP::ValueArg<std::string> arg_bamfile("m", "mapped_reads", "Sorted bam File", true, "", "string");
	TCLAP::ValueArg<std::string> arg_vcf("v", "vcf", "VCF output file name", false, "", "string");
	TCLAP::ValueArg<std::string> arg_bedpe("b", "bedpe", " bedpe output file name", false, "", "string");
	TCLAP::ValueArg<int> arg_support("s", "min_support", "Minimum number of reads that support a SV. Default: 10", false, 10, "int");
	TCLAP::ValueArg<int> arg_splits("", "max_num_splits", "Maximum number of splits per read to be still taken into account. Default: 7", false, 7, "int");
	TCLAP::ValueArg<int> arg_dist("d", "max_distance", "Maximum distance to group SV together. Default: 1kb", false, 1000, "int");
	TCLAP::ValueArg<int> arg_threads("t", "threads", "Number of threads to use. Default: 3", false, 3, "int");
	TCLAP::ValueArg<int> arg_minlength("l", "min_length", "Minimum length of SV to be reported. Default: 30", false, 30, "int");
	TCLAP::ValueArg<int> arg_mq("q", "minmapping_qual", "Minimum Mapping Quality. Default: 20", false, 20, "int");
	TCLAP::ValueArg<int> arg_numreads("n", "num_reads_report", "Report up to N reads that support the SV. Default: 0", false, 0, "int");
	TCLAP::ValueArg<std::string> arg_tmp_file("", "tmp_file", "patht to temporary file otherwise Sniffles will use the current directory.", false, "", "string");
	TCLAP::SwitchArg arg_MD_cigar("", "use_MD_Cigar", "Enables Sniffles to use the alignment information to screen for suspicious regions.", cmd, false);
	TCLAP::SwitchArg arg_genotype("", "genotype", "Enables Sniffles to compute the genotypes.", cmd, false);
	TCLAP::SwitchArg arg_cluster("", "cluster", "Enables Sniffles to phase SVs that occure on the same reads", cmd, false);

	cmd.add(arg_numreads);
	cmd.add(arg_tmp_file);
	cmd.add(arg_dist);
	cmd.add(arg_threads);
	cmd.add(arg_bedpe);
	cmd.add(arg_vcf);
	cmd.add(arg_minlength);
	cmd.add(arg_mq);
	cmd.add(arg_splits);
	cmd.add(arg_support);
	cmd.add(arg_bamfile);
	//parse cmd:
	cmd.parse(argc, argv);

	Parameter::Instance()->debug = true;
	Parameter::Instance()->score_treshold = 10;
	Parameter::Instance()->read_name = " ";//m141225_195223_00118_c100750772550000001823151707081513_s1_p0/62524/6888_19885";//0076373e-d278-4316-9967-9b4c0b74df57_Basecall_Alignment_template";//21_30705246";	;//just for debuging reasons!
	Parameter::Instance()->bam_files.push_back(arg_bamfile.getValue());
	Parameter::Instance()->min_mq = arg_mq.getValue();
	Parameter::Instance()->output_vcf = arg_vcf.getValue();
	Parameter::Instance()->report_n_reads = arg_numreads.getValue();
	Parameter::Instance()->min_support = arg_support.getValue();
	Parameter::Instance()->max_splits = arg_splits.getValue();
	Parameter::Instance()->max_dist = arg_dist.getValue();
	Parameter::Instance()->min_length = arg_minlength.getValue();
	Parameter::Instance()->genotype = arg_genotype.getValue();
	Parameter::Instance()->phase = arg_cluster.getValue();
	Parameter::Instance()->num_threads = arg_threads.getValue();
	Parameter::Instance()->output_bedpe = arg_bedpe.getValue();
	Parameter::Instance()->tmp_file = "test.tmp";	//arg_tmp_file.getValue();
	Parameter::Instance()->min_grouping_support = 1;


	if (Parameter::Instance()->tmp_file.empty()) {
		std::stringstream ss;
		//ss<<"."; //TODO: User does not need to see this!
		ss << rand();
		ss << "_tmp";
		Parameter::Instance()->tmp_file = ss.str();
	}
}

void parse_binary() {
	std::string tmp_name_file = Parameter::Instance()->tmp_file; // this file is created in IPrinter and stores the names and ID of SVS.
	tmp_name_file += "Names";

	FILE * alt_allel_reads = fopen(tmp_name_file.c_str(), "r");
	if (alt_allel_reads == NULL) {
		std::cerr << "ClusterParse: could not open tmp file: " << tmp_name_file.c_str() << std::endl;
	}
	std::cout << "start" << std::endl;
	name_str tmp;
	size_t nbytes = fread(&tmp, sizeof(struct name_str), 1, alt_allel_reads);
	std::cout << tmp.read_name << std::endl;
	while (nbytes != 0) {
		int max_ID = std::max(max_ID, tmp.svs_id);

		if (tmp.svs_id == 34 || tmp.svs_id == 35) {
			std::cout << "Cluster: " << tmp.svs_id << " " << tmp.read_name << std::endl;
		}
		//	std::cout << tmp.read_name << std::endl;
		nbytes = fread(&tmp, sizeof(struct name_str), 1, alt_allel_reads);
	}
	fclose(alt_allel_reads);
}

int main(int argc, char *argv[]) {

	try {
		//init parameter and reads user defined parameter from command line.
		read_parameters(argc, argv);

		//init openmp:
		omp_set_dynamic(0);
		omp_set_num_threads(Parameter::Instance()->num_threads);

		if((!Parameter::Instance()->output_vcf.empty()) && (!Parameter::Instance()->output_bedpe.empty())){
			std::cerr<<"Please select only vcf OR bedpe output format!"<<std::endl;
			exit(1);
		}
		//init printer:
		IPrinter * printer;
		if (!Parameter::Instance()->ref_seq.empty()) {
			printer = new NGMPrinter();
		} else if (!Parameter::Instance()->output_vcf.empty()) {
			printer = new VCFPrinter();
		} else if (!Parameter::Instance()->output_bedpe.empty()) {
			printer = new BedpePrinter();
		} else {
			std::cerr << "Please specify an output file using -v or -b" << std::endl;
			return -1;
		}

		printer->init();

		detect_breakpoints(Parameter::Instance()->bam_files[0], printer); //we could write out all read names for each sVs
		printer->close_file();

		//cluster the SVs together:
		if (Parameter::Instance()->phase) {
			std::cout << "Start phasing: " << std::endl;
			Cluster_SVS *cluster = new Cluster_SVS();
			cluster->update_SVs();
		}

		//determine genotypes:
		if (Parameter::Instance()->genotype) {
			std::cout << "Start genotype calling" << std::endl;
			Genotyper * go = new Genotyper();
			go->update_SVs();
		}

		//realignment: Using NGM???
		if (!Parameter::Instance()->ref_seq.empty()) {
			std::cout << "Realignment step activated:" << std::endl;
			//1. parse in output file(s)
			//2. extract reads from bam files (read groups?) //shellscript: picard/samtools??
			//3. realign
			//4. call Sniffles
		}

		cout << "Cleaning tmp files" << endl;
		string del = "rm ";
		del += Parameter::Instance()->tmp_file;
		//system(del.c_str());

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
