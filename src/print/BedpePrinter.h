/*
 * BedePrinter.h
 *
 *  Created on: Aug 24, 2015
 *      Author: fsedlaze
 */

#ifndef PRINT_BEDPEPRINTER_H_
#define PRINT_BEDPEPRINTER_H_
#include "IPrinter.h"
class BedpePrinter:public IPrinter{
private:
	void print_header();
	void print_body(Breakpoint *& SV, RefVector ref);
public:
	BedpePrinter(){

	}
	~BedpePrinter(){

	}
};



#endif /* PRINT_BEDPEPRINTER_H_ */
