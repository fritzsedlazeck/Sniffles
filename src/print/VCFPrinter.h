/*
 * VCFPrinter.h
 *
 *  Created on: Aug 24, 2015
 *      Author: fsedlaze
 */

#ifndef PRINT_VCFPRINTER_H_
#define PRINT_VCFPRINTER_H_

#include "IPrinter.h"

class VCFPrinter:public IPrinter{
private:
	void print_header();
 	void print_body(std::vector<Breakpoint *> &SV, RefVector ref);
public:
	VCFPrinter(){

	}

	~VCFPrinter(){
	}

};



#endif /* PRINT_VCFPRINTER_H_ */
