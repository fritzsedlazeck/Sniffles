//============================================================================
// Name        : plane_sweep.cpp
// Author      : Fritz Sedlazeck
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================
/*
#include <iostream>
#include <algorithm>
#include <chrono>
#include "SamParser.h"
#include "Plane-sweep.h"

using namespace std::chrono;
int main(int argc, char *argv[]) {
	if (argc > 1) {
		IParser *parser;
		std::string filename = std::string(argv[2]);
		std::cout << filename << std::endl;
		std::string samfile = "sam";
		std::size_t found = filename.find(samfile);
		if (found != std::string::npos) {
			parser = new SamParser();
		} else {
			parser = new TableParser();
		}
		std::vector<read_str> reads = parser->parse_reads(filename);
		std::cout << "Read in reads: " << reads.size() << std::endl;

		switch (atoi(argv[1])) {
		case 1:
			if (argc == 3) {
				//exchangeable by implementing IParser:
				PlaneSweep * sweep = new PlaneSweep();
				std::cout << "start storing:" << std::endl;
				high_resolution_clock::time_point t1 = high_resolution_clock::now();
				for (size_t i = 0; i < reads.size(); i++) {
					sweep->add_read(reads[i]);
				}
				sweep->finalyze();
				high_resolution_clock::time_point t2 = high_resolution_clock::now();
				auto duration = std::chrono::duration_cast<std::chrono::milliseconds>( t2 - t1 ).count();
				std::cerr<<"time elapsed: "<<duration<<" milliseconds"<<std::endl;
			}
			break;
		case 2:
			 high_resolution_clock::time_point t1 = high_resolution_clock::now();
			int maxim = 0;
			for (size_t i = 0; i < reads.size(); i++) {
				maxim = std::max(maxim, (int) reads[i].stop);
			}
			std::cerr << maxim << std::endl;
			std::vector<int> genome;
			genome.resize(maxim + 1, 0);
			for (size_t i = 0; i < reads.size(); i++) {
				size_t j = reads[i].start;
				while (j < reads[i].stop) {
					genome[j]++;
					j++;
				}
			}
			for (size_t i = 0; i < genome.size(); i++) {
				if (i > 0 && genome[i - 1] != genome[i]) {
					std::cout << i << '\t' << genome[i] << std::endl;
				}
			}
			genome.clear();
			high_resolution_clock::time_point t2 = high_resolution_clock::now();
			auto duration = std::chrono::duration_cast<std::chrono::milliseconds>( t2 - t1 ).count();
			std::cerr<<"time elapsed: "<<duration<<" milliseconds"<<std::endl;
			break;
		}
	}else{
		std::cout<<"Choose:"<<std::endl;
		std::cout<<"1: Plain sweep"<<std::endl;
		std::cout<<"2: Brute force"<<std::endl;
	}
	return 0;
}
*/
