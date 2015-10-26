/*
 * PhaserSV.cpp
 *
 *  Created on: Sep 2, 2015
 *      Author: fsedlaze
 */

#include "PhaserSV.h"

void PhaserSV::phase(std::vector<Breakpoint *> &svs) {
	//run all against all to check read names.. Is there an easier method?
	for (size_t i = 0; i < svs.size(); i++) {
		std::map<std::string, read_str> tmp=svs[i]->get_coordinates().support;
		//std::cout<<i <<" "<<tmp.size();
		for (std::map<std::string, read_str>::iterator name = tmp.begin(); name != tmp.end(); name++) {
			//std::cout<<(*name).first<<std::endl;
			for (size_t j = 0; j < svs.size(); j++) {

				if (i!=j &&svs[j]->get_coordinates().support.find((*name).first) != svs[j]->get_coordinates().support.end()) {
					//found something!
					//std::cout<<"HIT"<<std::endl;
					svs[j]->add_grouped(svs[i]->get_id());
					svs[i]->add_grouped(svs[j]->get_id());
				}
			}

		}
	}
	std::cout<<"end phasing"<<std::endl;

	for (size_t i = 0; i < svs.size(); i++) {
		std::cout<<svs[i]->get_id()<<std::endl;
		std::cout<<"\t";
		for(size_t j=0;j<svs[i]->get_groupted().size();j++){
			std::cout<<svs[i]->get_groupted()[j]<<", ";
		}
		std::cout<<std::endl;
	}

}
