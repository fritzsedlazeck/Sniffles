/*
 * Detect_Breapoints.cpp
 *
 *  Created on: Jun 19, 2015
 *      Author: fsedlaze
 */

#include "Detect_Breakpoints.h"
#include "../tree/IntervallTree.h"
#include "../tree/TNode.h"

long get_ref_lengths(int id, RefVector ref) {
	long length = 0;

	for (size_t i = 0; i < (size_t) id && i < ref.size(); i++) {
		length += ref[i].RefLength + Parameter::Instance()->max_dist;
	}
	return length;
}

bool should_be_stored(Breakpoint *& point) {
	point->calc_support();
	if (point->get_SVtype() & TRA) {
		return (point->get_support() > 2); // this is needed as we take each chr independently and just look at the primary alignment
	} else {
		point->predict_SV();
		return (point->get_support() > Parameter::Instance()->min_support && point->get_length() > Parameter::Instance()->min_length);
	}
}
void flush_tree(IntervallTree & local_tree, TNode *& local_root, IntervallTree & final, TNode *& root_final, long pos) {
	IntervallTree tmp_tree;
	TNode *tmp_root = NULL;
	std::vector<Breakpoint *> points;
	local_tree.get_breakpoints(local_root, points);
	clarify(points);
	for (int i = 0; i < points.size(); i++) {
		if (abs(pos - points[i]->get_coordinates().start.min_pos) < 100000 || abs(pos - points[i]->get_coordinates().stop.max_pos) < 100000) { //TODO arbitrary threshold!
			tmp_tree.insert(points[i], tmp_root);
		} else if (points[i]->get_support() > Parameter::Instance()->min_support && points[i]->get_length() > Parameter::Instance()->min_length) {
			final.insert(points[i], root_final);
		}
	}
	local_tree.makeempty(local_root);
	local_tree = tmp_tree;
	local_root = tmp_root;
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
	std::cout << "start parsing..." << std::endl;
//Using Interval tree to store and manage breakpoints:

	IntervallTree final;
	TNode * root_final = NULL;
	int current_RefID = 0;

	IntervallTree bst;
	TNode *root = NULL;

	Alignment * tmp_aln = mapped_file->parseRead(Parameter::Instance()->min_mq);
	long ref_space = get_ref_lengths(tmp_aln->getRefID(), ref);
	while (!tmp_aln->getSequence().first.empty()) {

		if ((tmp_aln->getAlignment()->IsPrimaryAlignment()) && (!(tmp_aln->getAlignment()->AlignmentFlag & 0x800) && tmp_aln->get_is_save())) {
			//flush_tree(local_tree, local_root, final, root_final);
			if (current_RefID != tmp_aln->getRefID()) {	// Regular scan through the SV and move those where the end point lies far behind the current pos or reads. Eg. 1MB?
				current_RefID = tmp_aln->getRefID();
				ref_space = get_ref_lengths(tmp_aln->getRefID(), ref);
				std::cout << "Switch Chr " << ref[tmp_aln->getRefID()].RefName << " " << ref[tmp_aln->getRefID()].RefLength << std::endl;
				std::vector<Breakpoint *> points;
			//	clarify(points);
				bst.get_breakpoints(root, points);
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
			std::vector<str_event> cigar_event;
			std::vector<str_event> md_event;
			std::vector<aln_str> split_events;

			if (tmp_aln->getMappingQual() > Parameter::Instance()->min_mq) {
#pragma omp parallel // starts a new team
				{
#pragma omp sections
					{
						{
							//if (Parameter::Instance()->useMD_CIGAR) {
							cigar_event = tmp_aln->get_events_CIGAR();
							//}
						}
#pragma omp section
						{
							//if (Parameter::Instance()->useMD_CIGAR) {
							md_event = tmp_aln->get_events_MD(20);
							//}
						}
#pragma omp section
						{
							split_events = tmp_aln->getSA(ref);
						}
					}
				}
				tmp_aln->set_supports_SV((cigar_event.empty() && md_event.empty()) && split_events.empty());

				//sweep->add_read(tmp_aln);

				//maybe flush the tree after each chr.....?

				int cov = 0; //sweep->get_num_reads();

				//maybe just store the extreme intervals for coverage -> If the cov doubled within Xbp or were the coverage is 0.
				add_events(tmp_aln, cigar_event, 0, ref_space, bst, root, cov, tmp_aln->getQueryBases());
				add_events(tmp_aln, md_event, 1, ref_space, bst, root, cov, tmp_aln->getQueryBases());
				add_splits(tmp_aln, split_events, 2, ref, bst, root, cov, tmp_aln->getQueryBases());
			}
		}
		mapped_file->parseReadFast(Parameter::Instance()->min_mq, tmp_aln);
	}
	//filter and copy results:
	std::cout << "Finalizing  .." << std::endl;
	std::vector<Breakpoint *> points;
	clarify(points);
	bst.get_breakpoints(root, points);
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

//	sweep->finalyze();
	points.clear();
	final.get_breakpoints(root_final, points);
	//std::cout<<"Points: "<<points.size()<<endl;
	//clarify(points);
	for(size_t i =0;i<points.size();i++){
		if (points[i]->get_SVtype() & TRA) {
			points[i]->calc_support();
			points[i]->predict_SV();
		}
		if (points[i]->get_support() > Parameter::Instance()->min_support && points[i]->get_length() > Parameter::Instance()->min_length) {
			printer->printSV(points[i]);
		}
	}
}

void add_events(Alignment * tmp, std::vector<str_event> events, short type, long ref_space, IntervallTree & bst, TNode *&root, int cov, std::string read_seq) {
	bool flag = (strcmp(tmp->getName().c_str(), Parameter::Instance()->read_name.c_str()) == 0);
	for (size_t i = 0; i < events.size(); i++) {
		position_str svs;
		//position_str stop;
		read_str read;
		//read.name = tmp->getName();
		read.type = type;
		read.SV = 0;
		//start.support.push_back(read); //not very nice!
		//stop.support.push_back(read);
		if (flag) {
			std::cout << tmp->getName() << " " << tmp->getRefID() << " " << events[i].pos << " " << abs(events[i].length) << std::endl;
		}
		svs.start.min_pos = (long) events[i].pos;
		svs.start.min_pos += ref_space;
		svs.stop.max_pos = svs.start.min_pos;
		if (events[i].length > 0) { //length ==- length for insertions!
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

		if (type == 0 && events[i].length < 0) {
			read.SV |= INS; //insertion
		} else if (type == 0) {
			read.SV |= DEL; //deletion
		} else {
			read.SV |= DEL;
			read.SV |= INV;
		}

		if (flag) {
			std::cout << tmp->getName() << " " << tmp->getRefID() << " " << svs.start.min_pos - ref_space << " " << svs.stop.max_pos - ref_space << std::endl;
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

		if (svs.start.min_pos > svs.stop.max_pos) {
			//maybe we have to invert the directions???
			svs_breakpoint_str pos = svs.start;
			svs.start = svs.stop;
			svs.stop = pos;

			pair<bool, bool> tmp = read.strand;

			read.strand.first = tmp.second;
			read.strand.second = tmp.first;

			//read.strand.first = !tmp.first;
			//read.strand.second = !tmp.second;
		}

		//TODO: we might not need this:
		if (svs.start.min_pos > svs.stop.max_pos) {
			read.coordinates.first = svs.stop.max_pos;
			read.coordinates.second = svs.start.min_pos;
		} else {
			read.coordinates.first = svs.start.min_pos;
			read.coordinates.second = svs.stop.max_pos;
		}
		svs.support[tmp->getName()] = read;
		Breakpoint * point = new Breakpoint(svs, cov, std::abs(events[i].length));
		bst.insert(point, root);
	}
}

void add_splits(Alignment * tmp, std::vector<aln_str> events, short type, RefVector ref, IntervallTree & bst, TNode *&root, int cov, std::string read_seq) {
	/*	bool flag = false;
	 if (Parameter::Instance()->overlaps(ref[tmp->getRefID()].RefName,
	 tmp->getPosition(), tmp->getPosition() + tmp->getRefLength())) {
	 Parameter::Instance()->read_name = tmp->getName();
	 flag = true;
	 }
	 */

	bool flag = (strcmp(tmp->getName().c_str(), Parameter::Instance()->read_name.c_str()) == 0);
	for (size_t i = 1; i < events.size() && events.size() < Parameter::Instance()->max_splits; i++) {
		if (flag) {
			std::cout << "Genome pos: " << tmp->getName() << " ";
			if (events[i - 1].strand) {
				std::cout << "+";
			} else {
				std::cout << "-";
			}
			std::cout << events[i - 1].pos << " " << events[i - 1].pos + events[i - 1].length << " p2: ";

			if (events[i - 1].strand) {
				std::cout << "+";
			} else {
				std::cout << "-";
			}
			std::cout << events[i].pos << " " << events[i].pos + events[i].length << " Ref: " << events[i - 1].RefID << " " << events[i].RefID << std::endl;
		}
		position_str svs;
		//position_str stop;
		read_str read;
		//read.name = tmp->getName();
		read.type = 2;
		read.SV = 0;
		//stop.support.push_back(read);
		//they mimic paired end sequencing:
		if (events[i].RefID == events[i - 1].RefID) {

			//TODO: changed because of test:
			if (events[i - 1].strand == events[i].strand) {
				if (events[i - 1].strand) {
					read.strand.first = events[i - 1].strand;
					read.strand.second = !events[i].strand;
				} else {
					read.strand.first = !events[i - 1].strand;
					read.strand.second = events[i].strand;
				}

				svs.read_start = events[i - 1].read_pos_start + events[i - 1].length;
				svs.read_stop = events[i].read_pos_start;
				if (events[i - 1].strand) {
					svs.start.min_pos = events[i - 1].pos + events[i - 1].length + get_ref_lengths(events[i - 1].RefID, ref);
					svs.stop.max_pos = events[i].pos + get_ref_lengths(events[i].RefID, ref);
				} else {
					svs.start.min_pos = events[i].pos + events[i].length + get_ref_lengths(events[i].RefID, ref);
					svs.stop.max_pos = events[i - 1].pos + get_ref_lengths(events[i - 1].RefID, ref);
				}

				if ((svs.start.min_pos - svs.stop.max_pos) > 100) {
					read.SV |= DUP;
					//TODO ADDED &&(svs.read_stop - svs.read_start) > (Parameter::Instance()->min_cigar_event * 2)
				} else if (abs(svs.stop.max_pos - svs.start.min_pos) + (Parameter::Instance()->min_cigar_event * 2) < (svs.read_stop - svs.read_start) && (svs.read_stop - svs.read_start) > (Parameter::Instance()->min_cigar_event * 2)) {
					read.SV |= INS;
				} else if (abs(svs.stop.max_pos - svs.start.min_pos) > (svs.read_stop - svs.read_start) + (Parameter::Instance()->min_cigar_event * 2)) {
					read.SV |= DEL;
				} else {
					read.SV = 'n';
				}
			} else {					// if first part of read is in a different direction as the second part-> INV
				read.strand.first = events[i - 1].strand;
				read.strand.second = !events[i].strand;
				read.SV |= INV;
				if (events[i - 1].strand) {
					//std::cout<<events[i].pos<<"\t"<<events[i].RefID<<"\t"<<get_ref_lengths(events[i].RefID, ref)<<endl;
					svs.start.min_pos = events[i - 1].pos + events[i - 1].length + get_ref_lengths(events[i - 1].RefID, ref);
					svs.stop.max_pos = events[i].pos + events[i].length + get_ref_lengths(events[i].RefID, ref);
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

		if (flag) {
			std::cout << tmp->getName() << " start: " << svs.start.min_pos << " stop: " << svs.stop.max_pos;
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
			std::cout << std::endl;
		}
		if (read.SV != 'n') {
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

				//read.strand.first = !tmp.first;
				//read.strand.second = !tmp.second;
			}

			//TODO: we might not need this:
			if (svs.start.min_pos > svs.stop.max_pos) {
				read.coordinates.first = svs.stop.max_pos;
				read.coordinates.second = svs.start.min_pos;
			} else {
				read.coordinates.first = svs.start.min_pos;
				read.coordinates.second = svs.stop.max_pos;
			}

			svs.support[tmp->getName()] = read;
			Breakpoint * point = new Breakpoint(svs, cov, events[i].length);
			bst.insert(point, root);
		}
	}
}

void clarify(std::vector<Breakpoint *> & points) {
//if WTF regions next to duplications-> delete!
	/*for(size_t i=0;i<points.size();i++){

	 }*/
}

void estimate_parameters(std::string read_filename) {
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

	while (!tmp_aln->getSequence().first.empty() && num < 1000) {

		if ((tmp_aln->getAlignment()->IsPrimaryAlignment()) && (!(tmp_aln->getAlignment()->AlignmentFlag & 0x800) && tmp_aln->get_is_save())) {

			//get score ratio
			double score = tmp_aln->get_scrore_ratio();
			if (score != -1) {
				avg_score += score;
			} else {
				avg_score += avg_score / num;
			}
			//cout<<"Para:\t"<<score;
			//get avg mismatches
			std::string md = tmp_aln->get_md();
			if (!md.empty()) {
				avg_mis += tmp_aln->get_num_mismatches(md);
				//cout<<"\t"<<tmp_aln->get_num_mismatches(md);
			}
			//cigar threshold: (without 1!)
			avg_indel += tmp_aln->get_avg_indel_length_Cigar();
			//cout<<"\t"<<tmp_aln->get_avg_indel_length_Cigar()<<endl;
			num++;
		}
		mapped_file->parseReadFast(Parameter::Instance()->min_mq, tmp_aln);
	}
	std::cout << avg_indel / num << std::endl;

	Parameter::Instance()->min_num_mismatches = 0.3;			//(avg_mis / num) * 0.3; //previously: 0.3
	Parameter::Instance()->min_cigar_event = 40;			//(avg_indel / num) * 20;	//previously: 20
	Parameter::Instance()->score_treshold = 2;			//(avg_score / num) / 2;	//previously: 2 //2

	std::cout << "score: " << Parameter::Instance()->score_treshold << std::endl;
	std::cout << "md: " << Parameter::Instance()->min_num_mismatches << std::endl;
	std::cout << "indel: " << Parameter::Instance()->min_cigar_event << std::endl;

}
