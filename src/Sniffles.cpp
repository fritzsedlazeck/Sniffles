//============================================================================
// Name        : Sniffles.cpp
// Author      : Fritz Sedlazeck
// Version     :
// Copyright   : MIT License
// Description : Detection of SVs for long read data.
//============================================================================
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
#include "ArgParseOutput.h"
#include "force_calling/Force_calling.h"

//cmake -D CMAKE_C_COMPILER=/usr/local/bin/gcc-8 -D CMAKE_CXX_COMPILER=/usr/local/bin/g++-8 ..

//TODO:
// Check the calibration again in the beginnig (Iceland study)

//check strand headers.
// strand bias??
// I think you could make your performance on PacBio reads even better with a few modifications:
//b. In pbsv, I use a simply mononucleotide consistency check to determine whether to cluster insertions from different reads as supporting the "same" events.  In addition to looking at the similarity of length and breakpoints,
//you could measure [min(Act)+min(Cct)+min(Gct)+min(Tct) / max(Act)+max(Cct)+max(Gct)+max(Tct)]  Even a lax criterion (>0.25)
//can avoid clustering phantom insertions (where one is say all A and the another is G+T).
//[min(A1,A2)+min(C1,C2)+min(G1,G2)+min(T1,T2)[/[max...]/
Parameter* Parameter::m_pInstance = NULL;

template<typename T>
void printParameter(std::stringstream & usage, TCLAP::ValueArg<T> & arg) {

	usage << "    " << arg.longID() << std::endl;
	usage << "        " << arg.getDescription();
	if (!arg.isRequired()) {
		usage << " [" << arg.getValue() << "]";
	}
	usage << std::endl;
}

void printParameter(std::stringstream & usage, TCLAP::SwitchArg & arg) {

	usage << "    " << arg.longID() << std::endl;
	usage << "        " << arg.getDescription();
	if (!arg.isRequired()) {
		usage << " [" << (arg.getValue() ? "true" : "false") << "]";
	}
	usage << std::endl;
}

void read_parameters(int argc, char *argv[]) {

//	TCLAP::CmdLine cmd("", ' ', "", true);
	TCLAP::CmdLine cmd("Sniffles version ", ' ', Parameter::Instance()->version);

	TCLAP::ValueArg<std::string> arg_readname("", "test_read", "readname", false, "", "string", cmd);// for testing only!

	TCLAP::ValueArg<std::string> arg_bamfile("m", "mapped_reads", "Sorted bam File", true, "", "string", cmd);
	TCLAP::ValueArg<std::string> arg_vcf("v", "vcf", "VCF output file name", false, "", "string", cmd);
	TCLAP::ValueArg<std::string> arg_input_vcf("", "Ivcf", "Input VCF file name. Enable force calling", false, "", "string", cmd);
	TCLAP::ValueArg<std::string> arg_bedpe("b", "bedpe", " bedpe output file name", false, "", "string", cmd);
	TCLAP::ValueArg<std::string> arg_tmp_file("", "tmp_file", "path to temporary file otherwise Sniffles will use the current directory.", false, "", "string", cmd);

	//TCLAP::ValueArg<std::string> arg_chrs("c", "chrs", " comma seperated list of chrs to scan", false, "", "string");
	TCLAP::ValueArg<int> arg_support("s", "min_support", "Minimum number of reads that support a SV.", false, 10, "int", cmd);
	TCLAP::ValueArg<int> arg_splits("", "max_num_splits", "Maximum number of splits per read to be still taken into account.", false, 7, "int", cmd);
	TCLAP::ValueArg<int> arg_dist("d", "max_distance", "Maximum distance to group SV together.", false, 1000, "int", cmd);
	TCLAP::ValueArg<int> arg_threads("t", "threads", "Number of threads to use.", false, 3, "int", cmd);
	TCLAP::ValueArg<int> arg_minlength("l", "min_length", "Minimum length of SV to be reported.", false, 30, "int", cmd);
	TCLAP::ValueArg<int> arg_mq("q", "minmapping_qual", "Minimum Mapping Quality.", false, 20, "int", cmd);
	TCLAP::ValueArg<int> arg_numreads("n", "num_reads_report", "Report up to N reads that support the SV in the vcf file. -1: report all.", false, 0, "int", cmd);
	TCLAP::ValueArg<int> arg_segsize("r", "min_seq_size", "Discard read if non of its segment is larger then this.", false, 2000, "int", cmd);
	TCLAP::ValueArg<int> arg_zmw("z", "min_zmw", "Discard SV that are not supported by at least x zmws. This applies only for PacBio recognizable reads.", false, 0, "int", cmd);
	TCLAP::ValueArg<int> arg_cluster_supp("", "cluster_support", "Minimum number of reads supporting clustering of SV.", false, 1, "int", cmd);
	TCLAP::ValueArg<int> arg_parameter_maxdist("", "max_dist_aln_events", "Maximum distance between alignment (indel) events.", false, 4, "int", cmd);
	TCLAP::ValueArg<int> arg_parameter_maxdiff("", "max_diff_per_window", "Maximum differences per 100bp.", false, 50, "int", cmd);

	TCLAP::SwitchArg arg_genotype("", "genotype", "Inactivated: Automatically true.", cmd, true);
	TCLAP::SwitchArg arg_cluster("", "cluster", "Enables Sniffles to phase SVs that occur on the same reads", cmd, false);
	TCLAP::SwitchArg arg_std("", "ignore_sd", "Ignores the sd based filtering. ", cmd, false);
	TCLAP::SwitchArg arg_bnd("", "report_BND", "Dont report BND instead use Tra in vcf output. ", cmd, true);
	TCLAP::SwitchArg arg_seq_old("", "report-seq", "Inactivated (see not_report_seq).", cmd, false);
	TCLAP::SwitchArg arg_seq("", "not_report_seq", "Dont report sequences for indels in vcf output. (Beta version!) ", cmd, false);
	TCLAP::SwitchArg arg_coords("", "change_coords", "Adopt coordinates for force calling if finding evidence. ", cmd, false);
	TCLAP::SwitchArg arg_parameter("", "skip_parameter_estimation", "Enables the scan if only very few reads are present. ", cmd, false);
	TCLAP::SwitchArg arg_cs_string("", "cs_string", "Enables the scan of CS string instead of Cigar and MD. ", cmd, false);
	//TCLAP::SwitchArg arg_read_strand("", "report_read_strands", "Enables the report of the strand categories per read. (Beta) ", cmd, false);
	TCLAP::SwitchArg arg_str("", "report_str", "Enables the report of str. (alpha testing) ", cmd, false);

	TCLAP::SwitchArg arg_ccs("", "ccs_reads", "Preset CCS Pacbio setting. (Beta) ", cmd, false);

	TCLAP::ValueArg<float> arg_allelefreq("f", "allelefreq", "Threshold on allele frequency (0-1). ", false, 0.0, "float", cmd);
	TCLAP::ValueArg<float> arg_hetfreq("", "min_het_af", "Threshold on allele frequency (0-1). ", false, 0.3, "float", cmd);
	TCLAP::ValueArg<float> arg_homofreq("", "min_homo_af", "Threshold on allele frequency (0-1). ", false, 0.8, "float", cmd);
	TCLAP::ValueArg<float> arg_delratio("", "del_ratio", "Estimated ration of deletions per read (0-1). ", false, 0.0458369, "float", cmd);
	TCLAP::ValueArg<float> arg_insratio("", "ins_ratio", "Estimated ratio of insertions per read (0-1). ", false, 0.049379, "float", cmd);

	std::stringstream usage;
	usage << "" << std::endl;
	usage << "Usage: sniffles [options] -m <sorted.bam> -v <output.vcf> " << std::endl;
	usage << "Version: " << Parameter::Instance()->version << std::endl;
	usage << "Contact: fritz.sedlazeck@gmail.com" << std::endl;
	usage << std::endl;
	usage << "Input/Output:" << std::endl;

	printParameter<std::string>(usage, arg_bamfile);
	printParameter<std::string>(usage, arg_vcf);
	printParameter<std::string>(usage, arg_bedpe);
	printParameter<std::string>(usage, arg_input_vcf);
	printParameter<std::string>(usage, arg_tmp_file);

	usage << "" << std::endl;
	usage << "General:" << std::endl;
	printParameter<int>(usage, arg_support);
	printParameter<int>(usage, arg_splits);
	printParameter<int>(usage, arg_dist);
	printParameter<int>(usage, arg_threads);
	printParameter<int>(usage, arg_minlength);
	printParameter<int>(usage, arg_mq);
	printParameter<int>(usage, arg_numreads);
	printParameter<int>(usage, arg_segsize);
	printParameter<int>(usage, arg_zmw);
	printParameter(usage, arg_cs_string);

	usage << "" << std::endl;
	usage << "Clustering/phasing and genotyping:" << std::endl;
	printParameter(usage, arg_genotype);
	printParameter(usage, arg_cluster);
	printParameter<int>(usage, arg_cluster_supp);
	printParameter<float>(usage, arg_allelefreq);
	printParameter<float>(usage, arg_homofreq);
	printParameter<float>(usage, arg_hetfreq);

	usage << "" << std::endl;
	usage << "Advanced:" << std::endl;
	printParameter(usage, arg_bnd);
	printParameter(usage, arg_seq);
	printParameter(usage,arg_seq_old);
	printParameter(usage, arg_std);
	//printParameter(usage, arg_read_strand);
	printParameter(usage, arg_ccs);
	printParameter(usage, arg_str);

	usage << "" << std::endl;
	usage << "Parameter estimation:" << std::endl;
	printParameter(usage, arg_parameter);
	printParameter<float>(usage, arg_delratio);
	printParameter<float>(usage, arg_insratio);
	printParameter<int>(usage, arg_parameter_maxdiff);
	printParameter<int>(usage, arg_parameter_maxdist);

	cmd.setOutput(new ArgParseOutput(usage.str(), ""));

	/*	cmd.add(arg_homofreq);
	 cmd.add(arg_hetfreq);
	 cmd.add(arg_input_vcf);
	 cmd.add(arg_cluster_supp);
	 cmd.add(arg_numreads);
	 cmd.add(arg_zmw);
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
	 //	cmd.add(arg_chrs);*/
	//parse cmd:
	cmd.parse(argc, argv);

	Parameter::Instance()->change_coords = arg_coords.getValue();
	//Parameter::Instance()->debug = true;
	Parameter::Instance()->score_treshold = 10;
	Parameter::Instance()->read_name =arg_readname.getValue(); //just for debuging reasons!
	Parameter::Instance()->bam_files.push_back(arg_bamfile.getValue());
	Parameter::Instance()->min_mq = arg_mq.getValue();
	Parameter::Instance()->output_vcf = arg_vcf.getValue();
	Parameter::Instance()->report_n_reads = arg_numreads.getValue();
	Parameter::Instance()->min_support = arg_support.getValue();
	Parameter::Instance()->max_splits = arg_splits.getValue();
	Parameter::Instance()->max_dist = arg_dist.getValue();
	Parameter::Instance()->min_length = arg_minlength.getValue();
	Parameter::Instance()->genotype = true;//arg_genotype.getValue();
	Parameter::Instance()->phase = arg_cluster.getValue();
	Parameter::Instance()->num_threads = arg_threads.getValue();
	Parameter::Instance()->output_bedpe = arg_bedpe.getValue();
	Parameter::Instance()->tmp_file = arg_tmp_file.getValue();
	Parameter::Instance()->min_grouping_support = arg_cluster_supp.getValue();
	Parameter::Instance()->min_allelel_frequency = arg_allelefreq.getValue();
	Parameter::Instance()->min_segment_size = arg_segsize.getValue();
	Parameter::Instance()->reportBND = arg_bnd.getValue();
	Parameter::Instance()->input_vcf = arg_input_vcf.getValue();
	Parameter::Instance()->print_seq = !arg_seq.getValue();
	Parameter::Instance()->ignore_std = arg_std.getValue();
	Parameter::Instance()->min_zmw = arg_zmw.getValue();
	Parameter::Instance()->homfreq = arg_homofreq.getValue();
	Parameter::Instance()->hetfreq = arg_hetfreq.getValue();
	Parameter::Instance()->skip_parameter_estimation = arg_parameter.getValue();
	Parameter::Instance()->cs_string = arg_cs_string.getValue();
	Parameter::Instance()->read_strand = true;// arg_read_strand.getValue();
	Parameter::Instance()->ccs_reads = arg_ccs.getValue();
	Parameter::Instance()->str = arg_str.getValue();

	if (Parameter::Instance()->ccs_reads) {
		Parameter::Instance()->skip_parameter_estimation = true;
		Parameter::Instance()->ignore_std = false;

	}
	if (Parameter::Instance()->skip_parameter_estimation) {
		cout << "\tSkip parameter estimation." << endl;
		Parameter::Instance()->score_treshold = 2;
		Parameter::Instance()->window_thresh = arg_parameter_maxdiff.getValue();
		Parameter::Instance()->max_dist_alns = arg_parameter_maxdist.getValue();
		Parameter::Instance()->avg_del = arg_delratio.getValue();
		Parameter::Instance()->avg_ins = arg_insratio.getValue();
	}

	if(!Parameter::Instance()->input_vcf.empty()){
		cout<<"\tForce calling mode enabled. Setting parameter accordingly"<<endl;
		Parameter::Instance()->min_mq=1;
	}

	//Parse IDS:
	/*std::string buffer = arg_chrs.getValue();
	 int count = 0;
	 std::string name = "";
	 for (size_t i = 0; i < buffer.size(); i++) {
	 if (buffer[i] == ',') {
	 Parameter::Instance()->chr_names[name] = true;
	 name.clear();
	 } else {
	 name += buffer[i];
	 }
	 }
	 if (!name.empty()) {
	 Parameter::Instance()->chr_names[name] = true;
	 }
	 */
	if (Parameter::Instance()->min_allelel_frequency > 0 || !Parameter::Instance()->input_vcf.empty()) {
		std::cerr << "Automatically enabling genotype mode" << std::endl;
		Parameter::Instance()->genotype = true;
	}

	if (Parameter::Instance()->tmp_file.empty()) { //TODO change to genotyper file and phasing file!
		if (Parameter::Instance()->output_bedpe.empty()) {
			Parameter::Instance()->tmp_file = Parameter::Instance()->output_vcf;
		} else {
			Parameter::Instance()->tmp_file = Parameter::Instance()->output_bedpe;
		}

		Parameter::Instance()->tmp_file += "_tmp";
	}

	Parameter::Instance()->tmp_genotyp = Parameter::Instance()->tmp_file;
	Parameter::Instance()->tmp_phasing = Parameter::Instance()->tmp_file;
	Parameter::Instance()->tmp_genotyp += "_genotype";
	Parameter::Instance()->tmp_phasing += "_phase";
	//should I check tmp file path??
}



int main(int argc, char *argv[]) {

	try {
		Parameter::Instance()->testing = true;
		//init parameter and reads user defined parameter from command line.
		read_parameters(argc, argv);

	//	Parameter::Instance()->read_name = "51bda775-ff02-44bb-918b-5193f2c0adce" ;

		//init openmp:
		omp_set_dynamic(0);
		omp_set_num_threads(Parameter::Instance()->num_threads);

		if ((!Parameter::Instance()->output_vcf.empty()) && (!Parameter::Instance()->output_bedpe.empty())) {
			std::cerr << "Please select only vcf OR bedpe output format!" << std::endl;
			exit(EXIT_FAILURE);
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
		if (Parameter::Instance()->input_vcf.empty()) {
			//regular calling
			detect_breakpoints(Parameter::Instance()->bam_files[0], printer); //we could write out all read names for each sVs
		} else {
			//force calling was selected:
			force_calling(Parameter::Instance()->bam_files[0], printer);
		}

		printer->close_file();

		//cluster the SVs together:
		if (Parameter::Instance()->phase) {
			std::cout << "Start phasing: " << std::endl;
			Cluster_SVS *cluster = new Cluster_SVS();
			cluster->update_SVs();
		}

		//determine genotypes:
		if (Parameter::Instance()->genotype) {
			std::cout << "Start genotype calling:" << std::endl;
			Genotyper * go = new Genotyper();
		//	if (!Parameter::Instance()->testing) {

				go->update_SVs2();
			//} else {
		//		go->update_SVs();
		//	}
		}

	} catch (TCLAP::ArgException &e)  // catch any exceptions
	{
		std::cerr << "Sniffles error: " << e.error() << " for arg " << e.argId() << std::endl;
		return -1;
	}
	return 0;
}
