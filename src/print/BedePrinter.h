/*
 * BedePrinter.h
 *
 *  Created on: Aug 24, 2015
 *      Author: fsedlaze
 */

#ifndef PRINT_BEDEPRINTER_H_
#define PRINT_BEDEPRINTER_H_
#include "IPrinter.h"
class BedePrinter:public IPrinter{
private:
	void print_header();
	void print_body(std::vector<Breakpoint *>& SV, RefVector ref);
public:
	BedePrinter(){

	}
	~BedePrinter(){

	}
};



#endif /* PRINT_BEDEPRINTER_H_ */
