/*
 * Force_calling.h
 *
 *  Created on: Aug 24, 2017
 *      Author: sedlazec
 */

#ifndef FORCE_CALLING_FORCE_CALLING_H_
#define FORCE_CALLING_FORCE_CALLING_H_


#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <math.h>
#include <vector>
#include <algorithm>
#include "../Paramer.h"
#include "../BamParser.h"
#include "../tree/BinTree.h"
#include "../sub/Breakpoint.h"
#include "../sub/Detect_Breakpoints.h"
#include "VCF_parser.h"


void force_calling(std::string bam_file, IPrinter *& printer);




#endif /* FORCE_CALLING_FORCE_CALLING_H_ */
