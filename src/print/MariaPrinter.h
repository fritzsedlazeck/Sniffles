/*
 * MariaPrinter.h
 *
 *  Created on: Sep 4, 2015
 *      Author: fsedlaze
 */

#include "IPrinter.h"



class MariaPrinter: public IPrinter{
private:
	void print_header();
 	void print_body(Breakpoint * &SV, RefVector ref);

public:
	MariaPrinter(){

	}
	~MariaPrinter(){

	}
};
