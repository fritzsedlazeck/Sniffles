/*
 * PhaserSV.h
 *
 *  Created on: Sep 2, 2015
 *      Author: fsedlaze
 */

#include "api/BamReader.h"
#include "../sub/Breakpoint.h"
#include "../tree/BinTree.h"

class PhaserSV{
private:

public:
	PhaserSV(){

	}
	~PhaserSV(){

	}
	void phase(std::vector<Breakpoint *> &svs);

};
