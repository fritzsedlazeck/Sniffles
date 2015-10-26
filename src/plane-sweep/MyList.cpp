/*
 * MyList.cpp
 *
 *  Created on: Apr 14, 2015
 *      Author: fsedlaze
 */

#include "MyList.h"

void MyList::push(Alignment * read) {
	Node * new_node = new Node(read);
	//std::cout << new_node->get_read()->getPosition() << " "
	//		<< new_node->get_read()->getPosition() + new_node->get_stop()
	//		<< std::endl;

	this->num_nodes++;
	if (this->start == NULL) { //set first node:
		start = new_node;
		last = new_node;
	} else { //we have already something:
		if (this->last != NULL
				&& new_node->get_stop() > this->last->get_stop()) { //should be inserted after the last node
			last->set_next(new_node);
			last = new_node;
		} else { //insert within the list:
			if (start->get_stop() > new_node->get_stop()) {
				//insert before start;
				new_node->set_next(start);
				start = new_node;
			} else {
				Node * prev = start;
				Node * curr = start->get_next();

				while (curr != NULL && curr->get_stop() < new_node->get_stop()) {
					prev = curr;
					curr = curr->get_next();
				}
				new_node->set_next(curr);
				prev->set_next(new_node);
			}
		}
	}
}

Node* MyList::pop() {
	Node * tmp = start;
	if (tmp != NULL) {
		start = start->get_next();
		if (tmp == this->last) {
			last = NULL;
		}
		this->num_nodes--;
		tmp->set_next(NULL);
	}
	return tmp; //TODO: care about deleting the object!
}

int MyList::size() {
	return num_nodes;// + num_nodes_lowMQ;
}

Node * MyList::get_entries(int new_start) {
	if (new_start == -1
			|| (this->start != NULL && this->start->get_stop() <= new_start)) {
		return pop();
	}
	return NULL;
}
