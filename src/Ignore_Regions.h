/*
 * Ignore_Regions.h
 *
 *  Created on: Feb 4, 2016
 *      Author: fsedlaze
 */

#ifndef IGNORE_REGIONS_H_
#define IGNORE_REGIONS_H_
#include "sub/Breakpoint.h"
#include "tree/Intervall_bed.h"

void ignore_regions(std::vector<Breakpoint *> & final_SV);
void initialize_bed(IntervallTree_bed &bed_tree, Leaf *&root,RefVector ref);
#endif /* IGNORE_REGIONS_H_ */
