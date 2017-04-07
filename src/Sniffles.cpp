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
#include <unistd.h>
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

//cmake -D CMAKE_C_COMPILER=/opt/local/bin/gcc-mp-4.7 -D CMAKE_CXX_COMPILER=/opt/local/bin/g++-mp-4.7 ..

//TODO:
// Check pseudo event introduction.
// AF tag in genotyping.
// strand bias??

Parameter* Parameter::m_pInstance = NULL;

void read_parameters(int argc, char *argv[]) {

	TCLAP::CmdLine cmd("Sniffles version ", ' ', Parameter::Instance()->version);
	TCLAP::ValueArg<std::string> arg_bamfile("m", "mapped_reads", "Sorted bam File", true, "", "string");
	TCLAP::ValueArg<std::string> arg_vcf("v", "vcf", "VCF output file name", false, "", "string");
	TCLAP::ValueArg<std::string> arg_bedpe("b", "bedpe", " bedpe output file name", false, "", "string");
	TCLAP::ValueArg<int> arg_support("s", "min_support", "Minimum number of reads that support a SV. Default: 10", false, 10, "int");
	TCLAP::ValueArg<int> arg_splits("", "max_num_splits", "Maximum number of splits per read to be still taken into account. Default: 7", false, 7, "int");
	TCLAP::ValueArg<int> arg_dist("d", "max_distance", "Maximum distance to group SV together. Default: 1kb", false, 1000, "int");
	TCLAP::ValueArg<int> arg_threads("t", "threads", "Number of threads to use. Default: 3", false, 3, "int");
	TCLAP::ValueArg<int> arg_minlength("l", "min_length", "Minimum length of SV to be reported. Default: 30", false, 30, "int");
	TCLAP::ValueArg<int> arg_mq("q", "minmapping_qual", "Minimum Mapping Quality. Default: 20", false, 20, "int");
	TCLAP::ValueArg<int> arg_numreads("n", "num_reads_report", "Report up to N reads that support the SV in the vcf file. Default: 0", false, 0, "int");
	TCLAP::ValueArg<int> arg_segsize("r","min_seq_size","Discard read if non of its segment is larger then this. Default: 2kb",false,2000,"int");
	TCLAP::ValueArg<std::string> arg_tmp_file("", "tmp_file", "path to temporary file otherwise Sniffles will use the current directory.", false, "", "string");
	TCLAP::SwitchArg arg_genotype("", "genotype", "Enables Sniffles to compute the genotypes.", cmd, false);
	TCLAP::SwitchArg arg_cluster("", "cluster", "Enables Sniffles to phase SVs that occur on the same reads", cmd, false);
	TCLAP::ValueArg<int> arg_cluster_supp("", "cluster_support", "Minimum number of reads supporting clustering of SV. Default: 1", false, 1, "int");
	TCLAP::ValueArg<float> arg_allelefreq("f", "allelefreq", "Threshold on allele frequency (0-1).", false, 0.0, "float");

	cmd.add(arg_cluster_supp);
	cmd.add(arg_numreads);
	cmd.add(arg_segsize);
	cmd.add(arg_tmp_file);
	cmd.add(arg_dist);
	cmd.add(arg_threads);
	cmd.add(arg_minlength);
	cmd.add(arg_mq);
	cmd.add(arg_splits);
	cmd.add(arg_bedpe);
	cmd.add(arg_vcf);
	cmd.add(arg_allelefreq);
	cmd.add(arg_support);
	cmd.add(arg_bamfile);

	//parse cmd:
	cmd.parse(argc, argv);

	Parameter::Instance()->debug = true;
	Parameter::Instance()->score_treshold = 10;
	Parameter::Instance()->read_name = "m151104_233737_42291_c100924012550000001823194105121673_s1_p0/122462/0_24344";//m151102_123142_42286_c100922632550000001823194205121665_s1_p0/80643/0_20394"; //"22_36746138"; //just for debuging reasons!
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
	Parameter::Instance()->tmp_file = arg_tmp_file.getValue();
	Parameter::Instance()->min_grouping_support = arg_cluster_supp.getValue();
	Parameter::Instance()->min_allelel_frequency = arg_allelefreq.getValue();
	Parameter::Instance()->min_segment_size = arg_segsize.getValue();
	if (Parameter::Instance()->min_allelel_frequency > 0) {
		std::cerr << "Automatically enabling genotype mode" << std::endl;
		Parameter::Instance()->genotype = true;
	}

	if (Parameter::Instance()->tmp_file.empty()) {
		std::stringstream ss;
		srand(time(NULL));
		//ss << rand();
		//sleep(5);
		//ss << rand();
		ss << arg_bamfile.getValue();
		ss << "_tmp";
		Parameter::Instance()->tmp_file = ss.str(); //check if file exists! -> if yes throw the dice again
	}
}

//some toy/test functions:
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

double comp_std(std::vector<int> pos, int start) {
	double count = 0;
	double std_start = 0;

	for (size_t i = 0; i < pos.size(); i++) {
		count++;
		if (pos[i] != -1) {
			long diff = (start - pos[i]);
			//	std::cout << "DIFF Start: " << diff << std::endl;
			std_start += std::pow((double) diff, 2.0);
		}
	}
	return std::sqrt(std_start / count);
}

void test_sort_insert(int pos, std::vector<int> & positions) {

	size_t i = 0;
	while (i < positions.size() && positions[i] < pos) {
		i++;
	}
	positions.insert(positions.begin() + i, pos);

}

double test_comp_std_quantile(std::vector<int> positions, int position) {
	double count = 0;
	std::vector<int> std_start_dists;
	double std_start = 0;

	for (std::vector<int>::iterator i = positions.begin(); i != positions.end(); i++) {

		long diff = (position - (*i));
		//	std::cout << "DIFF Start: " << diff << std::endl;
		test_sort_insert(std::pow((double) diff, 2.0), std_start_dists);
		//std_start += std::pow((double) diff, 2.0);

	}

	count = 0;
	for (size_t i = 0; i < std_start_dists.size() / 2; i++) {
		std_start += std_start_dists[i];
		count++;
	}

	return std::sqrt(std_start / count);

}

void test_std() {
	srand(time(NULL));
	int start = rand() % 100000; /// sqrt(1/12) for ins. Plot TRA std vs. cov/support.
	std::vector<int> positions;
	double avg = 0;
	double num = 0;

	for (int border = 100; border < 9001; border = border * 10) {
		for (int t = 0; t < 10; t++) {
			for (int cov = 2; cov < 5; cov += 1) {

				for (size_t i = 0; i < cov; i++) {
					int pos = (rand() % border) + (start - (border / 2));
					positions.push_back(pos);
				}
				avg += comp_std(positions, start) / test_comp_std_quantile(positions, start);
				std::cout << "Cov: " << cov + 1 << " border: " << border << " STD: " << comp_std(positions, start) << std::endl; // / test_comp_std_quantile(positions, start) << std::endl;
				positions.clear();
				num++;
			}
		}
	}
	std::cout << "AVG: " << avg / num << std::endl;
}

void get_rand(int mean, int num, vector<int> & positions, int interval) {
//std::cout << "sim " << num << std::endl;
	for (size_t i = 0; i < num; i++) {
		int pos = (rand() % interval) + (mean - (interval / 2));
		positions.push_back(pos);
	}
}
#include <stdlib.h>
std::vector<int> sort_distance(std::vector<int> positions, int mean) {
	std::vector<int> distances;
	for (size_t i = 0; i < positions.size(); i++) {
		int dist = std::abs(mean - positions[i]);
		size_t j = 0;
		while (j < distances.size()) {
			if (std::abs(mean - distances[j]) < dist) {
				distances.insert(distances.begin() + j, positions[i]);
				break;
			}
			j++;
		}
		if (j == distances.size()) {
			distances.push_back(positions[i]);
		}
	}
	return distances;
}
void test_slimming() {
	double fract = 0.2;
	srand(time(NULL));
	int mean = rand() % 100000; /// sqrt(1/12) for ins. Plot TRA std vs. cov/support.
	int intervall = 1000;

	std::vector<std::vector<double> > stds;
	int key = 0;
	int cov = 100;
	for (double fract = 0.1; fract < 1; fract += 0.1) {

		//std::cout<<fract<<std::endl;
		std::vector<int> positions;
		get_rand(mean, round(cov * fract), positions, intervall); //random process
		get_rand(mean, round(cov * (1 - fract)), positions, 10); //focused calls
		//	std::cout << "Cov: " << cov << " border: " << intervall << " STD: " << comp_std(positions, mean) << std::endl;
		std::vector<int> dists;
		dists = sort_distance(positions, mean);

		/*		for (size_t i = 0; i < dists.size(); i++) {
		 std::cout << abs(mean - dists[i]) << std::endl;
		 }
		 */
		std::vector<double> std_tmp;
		for (size_t i = 0; i < dists.size(); i++) {
			std::vector<int> tmp;
			tmp.assign(dists.rbegin(), dists.rend() - i);
			double std = comp_std(tmp, mean);
			//std::cout << "Points: " << tmp.size() << " STD: " << std << std::endl;
			std_tmp.push_back(std);
		}
		stds.push_back(std_tmp);
	}

	for (size_t i = 0; i < stds.size(); i++) {
		for (size_t j = 0; j < stds[i].size(); j++) {
			std::cout << stds[i][j] << "\t";
		}
		std::cout << std::endl;
	}
}

int main(int argc, char *argv[]) {

	try {
		//	test_slimming();
		//	test_std();
		//
		//exit(0);
		//init parameter and reads user defined parameter from command line.
		read_parameters(argc, argv);

		//init openmp:
		omp_set_dynamic(0);
		omp_set_num_threads(Parameter::Instance()->num_threads);

		if ((!Parameter::Instance()->output_vcf.empty()) && (!Parameter::Instance()->output_bedpe.empty())) {
			std::cerr << "Please select only vcf OR bedpe output format!" << std::endl;
			exit(1);
		}
		//init printer:
		IPrinter * printer;
		if (!Parameter::Instance()->output_vcf.empty()) {
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


	} catch (TCLAP::ArgException &e)  // catch any exceptions
	{
		std::cerr << "Sniffles error: " << e.error() << " for arg " << e.argId() << std::endl;
	}
	return 0;
}
