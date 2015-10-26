/*
 * NGMPrinter.h
 *
 *  Created on: Sep 23, 2015
 *      Author: fsedlaze
 */

#ifndef PRINT_NGMPRINTER_H_
#define PRINT_NGMPRINTER_H_

#include "IPrinter.h"
class NGMPrinter:public IPrinter{

private:
	void print_header();
 	void print_body(std::vector<Breakpoint *> &SV, RefVector ref);
public:
 	NGMPrinter(){

	}

	~NGMPrinter(){
	}

};




#endif /* PRINT_NGMPRINTER_H_ */
