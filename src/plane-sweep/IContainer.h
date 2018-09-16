/*
 * IContainer.h
 *
 *  Created on: Apr 14, 2015
 *      Author: fsedlaze
 */

#ifndef ICONTAINER_H_
#define ICONTAINER_H_
#include "Node.h"
class IContainer{
protected:

public:
	IContainer(){}
	virtual ~IContainer(){};
	virtual void push(Alignment * read)=0;
	virtual int size()=0;
	virtual Node * get_entries(int new_start)=0;
};



#endif /* ICONTAINER_H_ */
