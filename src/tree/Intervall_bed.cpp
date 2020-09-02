/*
 * Intervall_bed.cpp
 *
 *  Created on: Feb 4, 2016
 *      Author: fsedlaze
 */

#include "Intervall_bed.h"


// Inserting a node
void IntervallTree_bed::insert(long start, long stop, Leaf *&p) {
	if (p == NULL) {
		p = new Leaf(start, stop);
		if (p == NULL) {
			std::cout << "Out of Space\n" << std::endl;
		}
	} else {

		long score = p->overlap(start, stop); //comparison function

		if (score > 0) {
			insert(start, stop, p->left);
			if ((bsheight(p->left) - bsheight(p->right)) == 2) {
				score = p->left->overlap(start, stop);
				if (score > 0) {
					p = srl(p);
				} else {
					p = drl(p);
				}
			}
		} else if (score < 0) {
			insert(start, stop, p->right);
			if ((bsheight(p->right) - bsheight(p->left)) == 2) {
				score = p->right->overlap(start, stop);
				if (score < 0) {
					p = srr(p);
				} else {
					p = drr(p);
				}
			}
		} else { //overlaps!
			std::cerr << "Two regions overlap and are thus ignored:" << std::endl;
		}
	}
	int m, n, d;
	m = bsheight(p->left);
	n = bsheight(p->right);
	d = max(m, n);
	p->set_height(d + 1);
}
// Finding the Smallest
Leaf * IntervallTree_bed::findmin(Leaf * p) {
	if (p == NULL) {
		return p;
	} else {
		while (p->left != NULL) {
			p = p->left;
			//return p;
		}
		return p;
	}
}
// Finding the Largest node
Leaf * IntervallTree_bed::findmax(Leaf * p) {
	if (p == NULL) {
		return p;
	} else {
		while (p->right != NULL) {
			p = p->right;
			//return p;
		}
		return p;
	}
}
// Finding an get_value()
bool IntervallTree_bed::is_in(long position, Leaf * &p) {
	if (p == NULL) {
		return false;
	} else {
		long score = p->overlap(position);
		if (score > 0) {
			return is_in(position, p->left);
		} else if (score < 0) {
			return is_in(position, p->right);
		} else {
			return true;
		}
	}
}
// Copy a tree
void IntervallTree_bed::copy(Leaf * &p, Leaf * &p1) {
	makeempty(p1);
	p1 = nodecopy(p);
}
// Make a tree empty
void IntervallTree_bed::makeempty(Leaf * &p) {
	Leaf * d;
	if (p != NULL) {
		makeempty(p->left);
		makeempty(p->right);
		d = p;
		free(d);
		p = NULL;
	}
}
// Copy the nodes
Leaf * IntervallTree_bed::nodecopy(Leaf * &p) {
	Leaf * temp;
	if (p == NULL) {
		return p;
	} else {
		temp = new Leaf(p->get_start(), p->get_stop()); //TODO!
		temp->left = nodecopy(p->left);
		temp->right = nodecopy(p->right);
		return temp;
	}
}

// Deleting a node
void IntervallTree_bed::del(long start, long stop, Leaf * &p) {
	Leaf * d;
	if (p == NULL) {
		std::cout << "Sorry! get_value() not found\n" << std::endl;
	} else {
		long score = p->overlap(start, stop);
		if (score > 0) {
			del(start, stop, p->left);
		} else if (score < 0) {
			del(start, stop, p->right);
		} else if ((p->left == NULL) && (p->right == NULL)) {
			d = p;
			free(d);
			p = NULL;
			std::cout << "get_value() deleted successfully\n" << std::endl;
		} else if (p->left == NULL) {
			d = p;
			free(d);
			p = p->right;
			std::cout << "get_value() deleted successfully\n" << std::endl;
		} else if (p->right == NULL) {
			d = p;
			p = p->left;
			free(d);
			std::cout << "get_value() deleted successfully\n" << std::endl;
		} else {
			//p->set_value(deletemin(p->right));
		}
	}
}

int IntervallTree_bed::deletemin(Leaf * &p) {
	int c;
	std::cout << "inside deltemin\n" << std::endl;
	if (p->left == NULL) {
		//c = p->get_value();
		p = p->right;
		return c;
	} else {
		c = deletemin(p->left);
		return c;
	}
}
void IntervallTree_bed::preorder(Leaf * p) {
	if (p != NULL) {
		//std::cout << p->get_data()->to_string() << "\t";
		preorder(p->left);
		preorder(p->right);
	}
}

// Inorder Printing
void IntervallTree_bed::inorder(Leaf * p, Leaf * root) {
	if (p != NULL) {
		inorder(p->left, root);
		//std::cout << p->get_data()->to_string();
		if (p == root) {
			std::cout << "*\t";
		} else {
			std::cout << "\t";
		}
		inorder(p->right, root);
	}
}

// PostOrder Printing
void IntervallTree_bed::postorder(Leaf * p) {
	if (p != NULL) {
		postorder(p->left);
		postorder(p->right);
		std::cout << p->get_start()<<" "<<p->get_stop()<< "\t";
	}
}

int IntervallTree_bed::max(int value1, int value2) {
	return ((value1 > value2) ? value1 : value2);
}
int IntervallTree_bed::bsheight(Leaf * p) {
	int t;
	if (p == NULL) {
		return -1;
	} else {
		t = p->get_height();
		return t;
	}
}

Leaf * IntervallTree_bed::srl(Leaf * &p1) {
	Leaf * p2;
	p2 = p1->left;
	p1->left = p2->right;
	p2->right = p1;
	p1->set_height(max(bsheight(p1->left), bsheight(p1->right)) + 1);
	p2->set_height(max(bsheight(p2->left), p1->get_height()) + 1);
	return p2;
}
Leaf * IntervallTree_bed::srr(Leaf * &p1) {
	Leaf * p2;
	p2 = p1->right;
	p1->right = p2->left;
	p2->left = p1;
	p1->set_height(max(bsheight(p1->left), bsheight(p1->right)) + 1);
	p2->set_height(max(p1->get_height(), bsheight(p2->right)) + 1);
	return p2;
}
Leaf * IntervallTree_bed::drl(Leaf * &p1) {
	p1->left = srr(p1->left);
	return srl(p1);
}
Leaf * IntervallTree_bed::drr(Leaf * &p1) {
	p1->right = srl(p1->right);
	return srr(p1);
}

int IntervallTree_bed::nonodes(Leaf * p) {
	int count = 0;
	if (p != NULL) {
		nonodes(p->left);
		nonodes(p->right);
		count++;
	}
	return count;
}

