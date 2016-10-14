/*
 * Cluster_SVs.h
 *
 *  Created on: Apr 28, 2016
 *      Author: fsedlaze
 */

#ifndef CLUSTER_CLUSTER_SVS_H_
#define CLUSTER_CLUSTER_SVS_H_
#include <iostream>
#include <fstream>
#include <string.h>
#include <vector>
#include <map>
#include "../Paramer.h"
struct __attribute__((packed)) name_str{
	long read_name;
	int svs_id;
};
struct combine_str{
	int curr_id;
	int alt_id;
	int support;
};

class Cluster_SVS{
private:
	std::map<long, std::vector<int> > parse_names_ids(int & max_ID) ;
	void update_SVs( std::vector<combine_str> & ids); //just because the pass is more efficient
	void add_id(int curr_id,int new_id, std::vector<combine_str> &  ids);
	int find_id(int curr_id, std::vector<combine_str> & ids);
public:
	Cluster_SVS(){
	}
	~Cluster_SVS(){
	}
	void update_SVs();
};


#endif /* CLUSTER_CLUSTER_SVS_H_ */
