/*
 * Detect_Breapoints.cpp
 *
 *  Created on: Jun 19, 2015
 *      Author: fsedlaze
 */

#include "Detect_Breakpoints.h"
std::string TRANS_type(char type) {
	string tmp;
	if (type & DEL) {
		tmp += "DEL";
	}
	if (type & INV) {
		if (!tmp.empty()) {
			tmp += '/';
		}
		tmp += "INV";
	}
	if (type & DUP) {
		if (!tmp.empty()) {
			tmp += '/';
		}
		tmp += "DUP";
	}
	if (type & INS) {
		if (!tmp.empty()) {
			tmp += '/';
		}
		tmp += "INS";
	}
	if (type & TRA) {
		if (!tmp.empty()) {
			tmp += '/';
		}
		tmp += "TRA";
	}

	return tmp; // should not occur!
}

long get_ref_lengths(int id, RefVector ref) {
	long length = 0;

	for (size_t i = 0; i < (size_t) id && i < ref.size(); i++) {
		length += (long) ref[i].RefLength + (long) Parameter::Instance()->max_dist;
	}
	return length;
}

bool should_be_stored(Breakpoint *& point) {
	point->calc_support(); // we need that before:
	if (point->get_SVtype() & TRA) { // we cannot make assumptions abut support yet.
		return (point->get_support() > 2); // this is needed as we take each chr independently and just look at the primary alignment
	} else if (point->get_support() > Parameter::Instance()->min_support) {
		point->predict_SV();
		return (point->get_support() > Parameter::Instance()->min_support && point->get_length() > Parameter::Instance()->min_length);
	}

	return false;
}
void polish_points(std::vector<Breakpoint *> & points) { //TODO might be usefull! but why does the tree not fully work??
	return;
	if (!points.empty()) {
		std::vector<Breakpoint *> new_points;
		points[0]->calc_support();
		new_points.push_back(points[0]);

		for (size_t i = 1; i < points.size(); i++) {
			points[i]->calc_support();
			if (should_be_stored(points[i])) {
				if (abs(points[i - 1]->get_coordinates().start.min_pos - points[i]->get_coordinates().start.min_pos) < Parameter::Instance()->min_length && abs(points[i - 1]->get_coordinates().stop.max_pos - points[i]->get_coordinates().stop.max_pos) < Parameter::Instance()->min_length) {
					//merge!
					std::cout << "HIT: " << points[i - 1]->get_coordinates().start.min_pos << " " << points[i - 1]->get_support() << " OTHER: " << points[i]->get_coordinates().start.min_pos << " " << points[i]->get_support() << std::endl;
				}
			}
		}
	}
}
void detect_breakpoints(std::string read_filename, IPrinter *& printer) {
	estimate_parameters(read_filename);
	BamParser * mapped_file = 0;
	RefVector ref;
	if (read_filename.find("bam") != string::npos) {
		mapped_file = new BamParser(read_filename);
		ref = mapped_file->get_refInfo();
	} else {
		cerr << "File Format not recognized. File must be a sorted .bam file!" << endl;
		exit(0);
	}
//Using PlaneSweep to comp coverage and iterate through reads:
//PlaneSweep * sweep = new PlaneSweep();
	std::cout << "Start parsing..." << std::endl;
//Using Interval tree to store and manage breakpoints:

	IntervallTree final;
	TNode * root_final = NULL;
	int current_RefID = 0;

	IntervallTree bst;
	TNode *root = NULL;
	//FILE * alt_allel_reads;
	FILE * ref_allel_reads;
	if (Parameter::Instance()->genotype) {

		std::string output = Parameter::Instance()->tmp_file.c_str();
		output += "ref_allele";
		ref_allel_reads = fopen(output.c_str(), "wb");

//	output = Parameter::Instance()->tmp_file.c_str();
		//output += "alt_allele";
		//alt_allel_reads = fopen(output.c_str(), "wb");
	}
	Alignment * tmp_aln = mapped_file->parseRead(Parameter::Instance()->min_mq);
	long ref_space = get_ref_lengths(tmp_aln->getRefID(), ref);
	std::string prev = "test";
	std::string curr = "wtf";
	long num_reads = 0;
	while (!tmp_aln->getQueryBases().empty()) {

		if ((tmp_aln->getAlignment()->IsPrimaryAlignment()) && (!(tmp_aln->getAlignment()->AlignmentFlag & 0x800) && tmp_aln->get_is_save())) {
			//flush_tree(local_tree, local_root, final, root_final);

			if (current_RefID != tmp_aln->getRefID()) {	// Regular scan through the SV and move those where the end point lies far behind the current pos or reads. Eg. 1MB?
				current_RefID = tmp_aln->getRefID();
				ref_space = get_ref_lengths(tmp_aln->getRefID(), ref);
				std::cout << "\tSwitch Chr " << ref[tmp_aln->getRefID()].RefName << " " << ref[tmp_aln->getRefID()].RefLength << std::endl;
				std::vector<Breakpoint *> points;
				//	clarify(points);
				bst.get_breakpoints(root, points);
				polish_points(points);
				for (int i = 0; i < points.size(); i++) {
					if (should_be_stored(points[i])) {
						if (points[i]->get_SVtype() & TRA) {
							final.insert(points[i], root_final);
						} else {
							printer->printSV(points[i]);
						}
					}
				}
				bst.makeempty(root);
			}
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
							//TODO ignore Splits that are shorter then XYbp??
							//		clock_t begin_split = clock();
							split_events = tmp_aln->getSA(ref);
							//		Parameter::Instance()->meassure_time(begin_split," Split reads ");
						}
					}
				}
				//tmp_aln->set_supports_SV(aln_event.empty() && split_events.empty());

				str_read tmp;
				tmp.SV_support = !(aln_event.empty() && split_events.empty());
				if ((Parameter::Instance()->genotype && !tmp.SV_support) && (score == -1 || score > Parameter::Instance()->score_treshold)) {
					//write read:
					tmp.chr = ref[tmp_aln->getRefID()].RefName;
					tmp.start = tmp_aln->getPosition();
					tmp.length = tmp_aln->getRefLength();
					tmp.SV_support = false;
					fwrite(&tmp, sizeof(struct str_read), 1, ref_allel_reads);
				}
				/*else { // we store the reads that support the SVs in IPrinter when writing out the SVs.
				 for (size_t i = 0; i < split_events.size(); i++) {
				 tmp.chr = ref[split_events[i].RefID].RefName;
				 tmp.start = split_events[i].pos;
				 tmp.length = split_events[i].length;
				 fwrite(&tmp, sizeof(struct str_read), 1, alt_allel_reads);
				 }
				 if (split_events.empty()) { //splits store the primary aln as well!
				 tmp.chr = ref[tmp_aln->getRefID()].RefName;
				 tmp.start = tmp_aln->getPosition();
				 tmp.length = tmp_aln->getRefLength();
				 fwrite(&tmp, sizeof(struct str_read), 1, alt_allel_reads);
				 }
				 }*/

				//	clock_t begin = clock();
				//maybe just store the extreme intervals for coverage -> If the cov doubled within Xbp or were the coverage is 0.
				if (!aln_event.empty()) {
					add_events(tmp_aln, aln_event, 0, ref_space, bst, root, num_reads);
				}
				//	Parameter::Instance()->meassure_time(begin, " add event ");
				//cout<<"event: "<<aln_event.size()<<endl;
				//	begin = clock();
				if (!split_events.empty()) {
					add_splits(tmp_aln, split_events, 1, ref, bst, root, num_reads);

				}
				//	Parameter::Instance()->meassure_time(begin, " add split ");

			}
		}
		mapped_file->parseReadFast(Parameter::Instance()->min_mq, tmp_aln);
		num_reads++;
		if (num_reads % 10000 == 0) {
			curr = tmp_aln->getName();
			cout << "\t\t# Processed reads: " << num_reads << endl;
			if (curr.size() == prev.size() && strcmp(prev.c_str(), curr.c_str()) == 0) {
				std::cerr << "Read occurred twice with primary alignment marked. Possibly no eof recognized in bam file." << std::endl;
				break;
			}
			prev = curr;
		}
	}

//	cout << "Print tree:" << std::endl;
//	bst.inorder(root);
//filter and copy results:
	std::cout << "Finalizing  .." << std::endl;
	std::vector<Breakpoint *> points;
	bst.get_breakpoints(root, points);
	polish_points(points);
	for (int i = 0; i < points.size(); i++) {
		if (should_be_stored(points[i])) {
			if (points[i]->get_SVtype() & TRA) {
				final.insert(points[i], root_final);
			} else {
				printer->printSV(points[i]);
			}
		}
	}
	bst.makeempty(root);
	if (Parameter::Instance()->genotype) {
		fclose(ref_allel_reads);
	}
//	sweep->finalyze();
	points.clear();
	final.get_breakpoints(root_final, points);
	polish_points(points);
	for (size_t i = 0; i < points.size(); i++) {
		if (points[i]->get_SVtype() & TRA) {
			points[i]->calc_support();
			points[i]->predict_SV();
		}
		if (points[i]->get_support() > Parameter::Instance()->min_support && points[i]->get_length() > Parameter::Instance()->min_length) {
			printer->printSV(points[i]);
		}
	}
}

void add_events(Alignment *& tmp, std::vector<str_event> events, short type, long ref_space, IntervallTree & bst, TNode *&root, long read_id) {

	bool flag = (strcmp(tmp->getName().c_str(), Parameter::Instance()->read_name.c_str()) == 0);
	for (size_t i = 0; i < events.size(); i++) {
		position_str svs;
		read_str read;
		read.type = 0;
		read.SV = events[i].type;

		if (flag) {
			std::cout << "ADD EVENT " << tmp->getName() << " " << tmp->getRefID() << " " << events[i].pos << " " << abs(events[i].length) << std::endl;
		}
		svs.start.min_pos = (long) events[i].pos;
		svs.start.min_pos += ref_space;
		svs.stop.max_pos = svs.start.min_pos;

		if (!(events[i].type & INS)) { //for all events but not INS!
			svs.stop.max_pos += events[i].length;
		}

		if (tmp->getStrand()) {
			read.strand.first = (tmp->getStrand());
			read.strand.second = !(tmp->getStrand());
		} else {
			read.strand.first = !(tmp->getStrand());
			read.strand.second = (tmp->getStrand());
		}
		//	start.support[0].read_start.min = events[i].read_pos;

		if (flag) {
			std::cout << tmp->getName() << " " << tmp->getRefID() << " " << svs.start.min_pos << " " << svs.stop.max_pos << " " << svs.stop.max_pos - svs.start.min_pos << std::endl;
		}

		if (svs.start.min_pos > svs.stop.max_pos) {
			//can this actually happen?
			read.coordinates.first = svs.stop.max_pos;
			read.coordinates.second = svs.start.min_pos;
		} else {
			read.coordinates.first = svs.start.min_pos;
			read.coordinates.second = svs.stop.max_pos;
		}

		svs.start.max_pos = svs.start.min_pos;
		svs.stop.min_pos = svs.stop.max_pos;

		if (svs.start.min_pos > svs.stop.max_pos) { //incase they are inverted
			svs_breakpoint_str pos = svs.start;
			svs.start = svs.stop;
			svs.stop = pos;
			pair<bool, bool> tmp = read.strand;
			read.strand.first = tmp.second;
			read.strand.second = tmp.first;
		}

		//TODO: we might not need this:
		if (svs.start.min_pos > svs.stop.max_pos) {
			read.coordinates.first = svs.stop.max_pos;
			read.coordinates.second = svs.start.min_pos;
		} else {
			read.coordinates.first = svs.start.min_pos;
			read.coordinates.second = svs.stop.max_pos;
		}
		read.id = read_id;
		svs.support[tmp->getName()] = read;
		Breakpoint * point = new Breakpoint(svs, events[i].length);
		bst.insert(point, root);
	}
}

void add_splits(Alignment *& tmp, std::vector<aln_str> events, short type, RefVector ref, IntervallTree & bst, TNode *&root, long read_id) {
	bool flag = (strcmp(tmp->getName().c_str(), Parameter::Instance()->read_name.c_str()) == 0);
	//flag = true;
	if (false) {
		cout << "SPLIT: " << std::endl;
		for (size_t i = 0; i < events.size(); i++) {
			std::cout << events[i].pos << " stop: " << events[i].pos + events[i].length << " " << events[i].RefID << " READ: " << events[i].read_pos_start << " " << events[i].read_pos_stop;
			if (events[i].strand) {
				cout << " +" << endl;
			} else {
				cout << " -" << endl;
			}
		}
	}

	for (size_t i = 1; i < events.size(); i++) {
		position_str svs;
		//position_str stop;
		read_str read;
		//read.name = tmp->getName();
		read.type = type;
		read.SV = 0;
		//stop.support.push_back(read);
		if (events[i].RefID == events[i - 1].RefID) { //IF different chr -> tra
			if (events[i - 1].strand == events[i].strand) { //IF same strand -> del/ins/dup
				if (events[i - 1].strand) {
					read.strand.first = events[i - 1].strand;
					read.strand.second = !events[i].strand;
				} else {
					read.strand.first = !events[i - 1].strand;
					read.strand.second = events[i].strand;
				}
				int len1 = 0;
				int len2 = 0;
				svs.read_start = events[i - 1].read_pos_stop; // (short) events[i - 1].read_pos_start + (short) events[i - 1].length;
				svs.read_stop = events[i].read_pos_start;
				if (events[i - 1].strand) {
					len1 = events[i - 1].pos;
					len2 = events[i].pos - events[i].length;
					svs.start.min_pos = events[i - 1].pos + events[i - 1].length + get_ref_lengths(events[i - 1].RefID, ref);
					svs.stop.max_pos = events[i].pos + get_ref_lengths(events[i].RefID, ref);
				} else {
					len1 = events[i].pos;
					len2 = events[i - 1].pos - events[i - 1].length;
					svs.start.min_pos = events[i].pos + events[i].length + get_ref_lengths(events[i].RefID, ref);
					svs.stop.max_pos = events[i - 1].pos + get_ref_lengths(events[i - 1].RefID, ref);
				}

				if (flag) {
					cout << "Debug: SV_Size: " << (svs.start.min_pos - svs.stop.max_pos) << " Ref_start: " << svs.start.min_pos - get_ref_lengths(events[i].RefID, ref) << " Ref_stop: " << svs.stop.max_pos - get_ref_lengths(events[i].RefID, ref) << " readstart: " << svs.read_start << " readstop: "
							<< svs.read_stop << "readstart+len: " << len1 << " readstop+len: " << len2 << endl;
				}

				if ((svs.stop.max_pos - svs.start.min_pos) > 0 && ((svs.stop.max_pos - svs.start.min_pos) + (Parameter::Instance()->min_cigar_event) < (svs.read_stop - svs.read_start) && (svs.read_stop - svs.read_start) > (Parameter::Instance()->min_cigar_event * 2))) {
					if (flag) {
						cout << "INS: " << endl;
					}
					//read.SV = 'n'; //TODO redefine criteria!
					read.SV |= INS;
					/*	if (events[i - 1].strand) { //TODO check f
					 svs.start.min_pos -=events[i - 1].length;
					 }else{
					 svs.start.min_pos -=events[i].length;
					 }*/
					//				cout<<"Split INS: "<<events[i - 1].length<<" "<<(svs.stop.max_pos - svs.start.min_pos)<<" "<<svs.read_stop - svs.read_start<<endl;
				} else if ((svs.start.min_pos - svs.stop.max_pos) * -1 > (svs.read_stop - svs.read_start) + (Parameter::Instance()->min_cigar_event)) {
					read.SV |= DEL;
					if (flag) {
						cout << "DEL" << endl;
					}
				} else if ((svs.start.min_pos - svs.stop.max_pos) > Parameter::Instance()->min_cigar_event && (svs.read_start - svs.read_stop)<Parameter::Instance()->min_length ) { //check with respect to the coords of reads!
					read.SV |= DUP;
					if (flag) {
						cout << "DUP: " << endl;
					}
					//TODO ADDED &&(svs.read_stop - svs.read_start) > (Parameter::Instance()->min_cigar_event)
				} else {
					if (flag) {
						cout << "N" << endl;
					}
					read.SV = 'n';
				}
			} else { // if first part of read is in a different direction as the second part-> INV
				if (flag) {
					cout << "INV" << endl;
				}
				read.strand.first = events[i - 1].strand;
				read.strand.second = !events[i].strand;
				read.SV |= INV;
				if (events[i - 1].strand) {
					svs.start.min_pos = events[i - 1].pos + events[i - 1].length + get_ref_lengths(events[i - 1].RefID, ref);
					svs.stop.max_pos = (events[i].pos + events[i].length) + get_ref_lengths(events[i].RefID, ref);
				} else {
					svs.start.min_pos = events[i - 1].pos + get_ref_lengths(events[i - 1].RefID, ref);
					svs.stop.max_pos = events[i].pos + get_ref_lengths(events[i].RefID, ref);
				}
			}

		} else { //if not on the same chr-> TRA
			read.strand.first = events[i - 1].strand;
			read.strand.second = !events[i].strand;
			if (events[i - 1].strand == events[i].strand) {
				if (events[i - 1].strand) {
					svs.start.min_pos = events[i - 1].pos + events[i - 1].length + get_ref_lengths(events[i - 1].RefID, ref);
					svs.stop.max_pos = events[i].pos + get_ref_lengths(events[i].RefID, ref);
				} else {
					svs.start.min_pos = events[i - 1].pos + get_ref_lengths(events[i - 1].RefID, ref);
					svs.stop.max_pos = events[i].pos + events[i].length + get_ref_lengths(events[i].RefID, ref);
				}
			} else {
				if (events[i - 1].strand) {
					svs.start.min_pos = events[i - 1].pos + events[i - 1].length + get_ref_lengths(events[i - 1].RefID, ref);
					svs.stop.max_pos = events[i].pos + events[i].length + get_ref_lengths(events[i].RefID, ref);
				} else {
					svs.start.min_pos = events[i - 1].pos + get_ref_lengths(events[i - 1].RefID, ref);
					svs.stop.max_pos = events[i].pos + get_ref_lengths(events[i].RefID, ref);
				}
			}
			read.SV |= TRA;
		}

		if (read.SV != 'n') {
			if (flag) {
				std::cout << "SPLIT: " << TRANS_type(read.SV) << " start: " << svs.start.min_pos << " stop: " << svs.stop.max_pos; //- get_ref_lengths(events[i].RefID, ref)
				if (events[i - 1].strand) {
					std::cout << " +";
				} else {
					std::cout << " -";
				}
				if (events[i].strand) {
					std::cout << " +";
				} else {
					std::cout << " -";
				}
				std::cout << " " << tmp->getName() << std::endl;
				std::cout<<"READ: "<<svs.read_start<<" "<<svs.read_stop<<" "<<svs.read_start - svs.read_stop<<std::endl;
			}
			//std::cout<<"split"<<std::endl;

			svs.start.max_pos = svs.start.min_pos;
			svs.stop.min_pos = svs.stop.max_pos;
			if (svs.start.min_pos > svs.stop.max_pos) {
				//maybe we have to invert the directions???
				svs_breakpoint_str pos = svs.start;
				svs.start = svs.stop;
				svs.stop = pos;

				pair<bool, bool> tmp = read.strand;

				read.strand.first = tmp.second;
				read.strand.second = tmp.first;
			}

			//TODO: we might not need this:
			if (svs.start.min_pos > svs.stop.max_pos) {
				read.coordinates.first = svs.stop.max_pos;
				read.coordinates.second = svs.start.min_pos;
			} else {
				read.coordinates.first = svs.start.min_pos;
				read.coordinates.second = svs.stop.max_pos;
			}

			//pool out?
			read.id = read_id;
			svs.support[tmp->getName()] = read;
			Breakpoint * point = new Breakpoint(svs, events[i].length);

			bst.insert(point, root);
			//	breakpoints_tmp1.push_back(pair<position_str,int> (svs, events[i].length)); TODO we could create 2 loops: 1st: collect evidence 2nd: further interpretation of complex Var
		}
	}
}

void clarify(std::vector<Breakpoint *> & points) {
//if WTF regions next to duplications-> delete!
	/*for(size_t i=0;i<points.size();i++){

	 }*/
}

void estimate_parameters(std::string read_filename) {
	cout << "Estimating parameter..." << endl;
	BamParser * mapped_file = 0;
	RefVector ref;
	if (read_filename.find("bam") != string::npos) {
		mapped_file = new BamParser(read_filename);
		ref = mapped_file->get_refInfo();
	} else {
		cerr << "File Format not recognized. File must be a sorted .bam file!" << endl;
		exit(0);
	}

	Alignment * tmp_aln = mapped_file->parseRead(Parameter::Instance()->min_mq);
	double num = 0;
	double avg_score = 0;
	double avg_mis = 0;
	double avg_indel = 0;
	double avg_diffs_perwindow = 0;
	vector<int> mis_per_window; //histogram over #differences
	vector<int> scores;
	std::string curr, prev = "";
	double avg_dist = 0;
	while (!tmp_aln->getQueryBases().empty() && num < 1000) {				//1000
		if (rand() % 100 < 20 && ((tmp_aln->getAlignment()->IsPrimaryAlignment()) && (!(tmp_aln->getAlignment()->AlignmentFlag & 0x800)))) {				//}&& tmp_aln->get_is_save()))) {
			//1. check differences in window => min_treshold for scanning!
			//2. get score ration without checking before hand! (above if!)
			double dist = 0;
			vector<int> tmp = tmp_aln->get_avg_diff(dist);
			avg_dist += dist;
			for (size_t i = 0; i < tmp.size(); i++) {
				while (tmp[i] + 1 > mis_per_window.size()) { //adjust length
					mis_per_window.push_back(0);
				}
				mis_per_window[tmp[i]]++;
			}
			//avg_diffs_perwindow+=tmp_aln->get_avg_diff();
			//get score ratio
			double score = round(tmp_aln->get_scrore_ratio());
			while (score + 1 > scores.size()) {
				scores.push_back(0);
			}
			scores[score]++;
			num++;
		}
		mapped_file->parseReadFast(Parameter::Instance()->min_mq, tmp_aln);
		curr = tmp_aln->getName();
		if (curr.size() == prev.size() && strcmp(prev.c_str(), curr.c_str()) == 0) {
			std::cerr << "Read occurred twice with primary alignment marked. Possibly no eof recognized in bam file." << std::endl;
			break;
		}
		prev = curr;
	}
	vector<int> nums;
	size_t pos = 0;
	Parameter::Instance()->max_dist_alns = floor(avg_dist / num) / 2;
	Parameter::Instance()->window_thresh = 25;
	if (!mis_per_window.empty()) {

		for (size_t i = 0; i < mis_per_window.size(); i++) {
			if (mis_per_window[i] != 0) {
				//		std::cout << i << ": " << mis_per_window[i] << std::endl;
			}
			for (size_t j = 0; j < mis_per_window[i]; j++) {
				nums.push_back(i);
			}
		}
		pos = nums.size() * 0.95; //the highest 5% cutoff
		if (pos <= nums.size()) {
			Parameter::Instance()->window_thresh = std::max(50, nums[pos]); //just in case we have too clean data! :)
		}
		nums.clear();
	}

	for (size_t i = 0; i < scores.size(); i++) {
		for (size_t j = 0; j < scores[i]; j++) {
			nums.push_back(i);
		}
	}
	pos = nums.size() * 0.05; //the lowest 5% cuttoff
	Parameter::Instance()->score_treshold = 2; //nums[pos]; //prev=2
	std::cout << "\tMax dist between aln events: " << Parameter::Instance()->max_dist_alns << std::endl;
	std::cout << "\tMax diff in window: " << Parameter::Instance()->window_thresh << std::endl;
	std::cout << "\tMin score ratio: " << Parameter::Instance()->score_treshold << std::endl;
}

