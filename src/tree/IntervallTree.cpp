/*
 * IntervallTree.cpp
 *
 *  Created on: Jun 23, 2015
 *      Author: fsedlaze
 */

#include "IntervallTree.h"

void IntervallTree::careful_screening(Breakpoint *& new_break, TNode *p) { //maybe I just need the pointer not a ref.
	if (p != NULL && !(new_break->get_coordinates().start.min_pos == -1 && new_break->get_coordinates().stop.max_pos == -1)) {
		careful_screening(new_break, p->left);
		if (p->get_data()->overlap(new_break) == 0) { //SV type
			p->get_data()->add_read(new_break);
			//cout<<"Merged: "<<endl;
			new_break->set_coordinates(-1, -1);
			return;
		}
		careful_screening(new_break, p->right);
	}
}

// Inserting a node
void IntervallTree::insert(Breakpoint * new_break, TNode *&p) {
	if (new_break->get_coordinates().start.min_pos == -1 && new_break->get_coordinates().stop.max_pos == -1) {
		return;
	}
	if (p == NULL) { // add to tree:
		p = new TNode(new_break);
		if (p == NULL) {
			std::cout << "Out of Space\n" << std::endl;
		}
	} else { // find on tree:
		long score = p->get_data()->overlap(new_break); //comparison function
		if (score == 0) { //add SV types?
			//cout<<"Merged"<<endl;
			p->get_data()->add_read(new_break);
			new_break->set_coordinates(-1, -1);
			//delete new_break;
			return;
		} else if (abs(score) < Parameter::Instance()->max_dist) { // if two or more events are too close:
			//std::cout<<"Screen"<<std::endl;
			careful_screening(new_break, p);
			if (new_break->get_coordinates().start.min_pos == -1 && new_break->get_coordinates().stop.max_pos == -1) {
				return;
			}
		}

		if (score > 0) { // go left
			insert(new_break, p->left);
			if ((bsheight(p->left) - bsheight(p->right)) == 2) {
				score = p->left->get_data()->overlap(new_break);
				if (score > 0) {
					p = srl(p);
				} else {
					p = drl(p);
				}
			}
		} else if (score < 0) { // go right
			insert(new_break, p->right);
			if ((bsheight(p->right) - bsheight(p->left)) == 2) {
				score = p->right->get_data()->overlap(new_break);
				if (score < 0) {
					p = srr(p);
				} else {
					p = drr(p);
				}
			}
		}
	}
	int m, n, d;
	m = bsheight(p->left);
	n = bsheight(p->right);
	d = max(m, n);
	p->set_height(d + 1);
}

void IntervallTree::insert_existant(Breakpoint * new_break, TNode *&p) {
	if (new_break->get_coordinates().start.min_pos == -1 && new_break->get_coordinates().stop.max_pos == -1) {
		return;
	}
	if (p == NULL) { // add to tree:
		return;
	} else { // find on tree:
		long score = p->get_data()->overlap(new_break); //comparison function
		if (score == 0) { //add SV types?
			p->get_data()->add_read(new_break);
			new_break->set_coordinates(-1, -1);
			//delete new_break;
			return;
		} else if (abs(score) < Parameter::Instance()->max_dist) { // if two or more events are too close:
			//std::cout<<"Screen"<<std::endl;
			careful_screening(new_break, p);
			if (new_break->get_coordinates().start.min_pos == -1 && new_break->get_coordinates().stop.max_pos == -1) {
				return;
			}
		}

		if (score > 0) { // go left
			insert_existant(new_break, p->left);
			if ((bsheight(p->left) - bsheight(p->right)) == 2) {
				score = p->left->get_data()->overlap(new_break);
				if (score > 0) {
					p = srl(p);
				} else {
					p = drl(p);
				}
			}
		} else if (score < 0) { // go right
			insert_existant(new_break, p->right);
			if ((bsheight(p->right) - bsheight(p->left)) == 2) {
				score = p->right->get_data()->overlap(new_break);
				if (score < 0) {
					p = srr(p);
				} else {
					p = drr(p);
				}
			}
		}
	}
	int m, n, d;
	m = bsheight(p->left);
	n = bsheight(p->right);
	d = max(m, n);
	p->set_height(d + 1);
}

bool IntervallTree::overlaps(long start, long stop, TNode *p) {
	if (p == NULL) {
		return false;
	} else {
		long score = p->get_data()->overlap_breakpoint(start,stop);
		if (score > 0) {
			return overlaps(start,stop, p->left);
		} else if (score < 0) {
			return overlaps(start,stop, p->right);
		} else {
			return true;
		}
	}
//	return false;
}


// Finding the Smallest
/*TNode * IntervallTree::findmin(TNode * p) {
	if (p == NULL) {
		std::cout << "The tree is empty\n" << std::endl;
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
TNode * IntervallTree::findmax(TNode * p) {
	if (p == NULL) {
		std::cout << "The tree is empty\n" << std::endl;
		return p;
	} else {
		while (p->right != NULL) {
			p = p->right;
			//return p;
		}
		return p;
	}
}*/


// Finding an get_value()
void IntervallTree::find(Breakpoint * point, TNode * &p) {
	if (p == NULL) {
		std::cout << "Sorry! get_value() not found\n" << std::endl;
	} else {
		long score = p->get_data()->overlap(point);
		if (score > 0) {
			find(point, p->left);
		} else if (score < 0) {
			find(point, p->right);
		} else {
			std::cout << "get_value() found!\n" << std::endl;
		}

	}
}
// Copy a tree
void IntervallTree::copy(TNode * &p, TNode * &p1) {
	clear(p1);
	p1 = nodecopy(p);
}
// Make a tree empty
void IntervallTree::clear(TNode * &p) {
	TNode * d;
	if (p != NULL) {
		clear(p->left);
		clear(p->right);
		d = p;
		free(d);
		p = NULL;
	}
}
// Copy the nodes
TNode * IntervallTree::nodecopy(TNode * &p) {
	TNode * temp;
	if (p == NULL) {
		return p;
	} else {
		temp = new TNode(p->get_data()); //TODO!
		temp->left = nodecopy(p->left);
		temp->right = nodecopy(p->right);
		return temp;
	}
}

// Deleting a node
void IntervallTree::del(Breakpoint * point, TNode * &p) {
	TNode * d;
	if (p == NULL) {
		std::cout << "Sorry! get_value() not found\n" << std::endl;
	} else {
		long score = p->get_data()->overlap(point);
		if (score > 0) {
			del(point, p->left);
		} else if (score < 0) {
			del(point, p->right);
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

int IntervallTree::deletemin(TNode * &p) {
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
void IntervallTree::preorder(TNode * p) {
	if (p != NULL) {
		//std::cout << p->get_data()->to_string() << "\t";
		preorder(p->left);
		preorder(p->right);
	}
}


void IntervallTree::get_breakpoints(TNode *p, std::vector<Breakpoint *> & points) {
	if (p != NULL) {
		get_breakpoints(p->right, points);
	//	std::cout << "( " << p->get_data()->get_coordinates().start.min_pos << "-" << p->get_data()->get_coordinates().stop.max_pos << " "<< p->get_data()->get_coordinates().support.size()<<" "<< TRANS_type2((*p->get_data()->get_coordinates().support.begin()).second.SV)<<" )"<<std::endl;
		points.push_back(p->get_data());
		get_breakpoints(p->left, points);
	}
}

// Inorder Printing
void IntervallTree::inorder(TNode * p) {
	if (p != NULL) {
		inorder(p->left);
		std::cout << p->get_data()->to_string() << endl;
		inorder(p->right);
	}
}
void IntervallTree::print(TNode *p) {
	if (p != NULL) {
		print(p->left);
		std::string msg = p->get_data()->to_string();
		if (msg.size() > 3) {
			std::cout << msg << endl;
		}
		//std::cout << "( " << p->get_data()->get_coordinates().start.min_pos << "-" << p->get_data()->get_coordinates().stop.max_pos << " "<< p->get_data()->get_coordinates().support.size()<<" )"<<std::endl;
		print(p->right);
	}
}

// PostOrder Printing
void IntervallTree::postorder(TNode * p) {
	if (p != NULL) {
		postorder(p->left);
		postorder(p->right);
		//std::cout << p->get_data()->to_string() << "\t";
	}
}

int IntervallTree::max(int value1, int value2) {
	return ((value1 > value2) ? value1 : value2);
}
int IntervallTree::bsheight(TNode * p) {
	int t;
	if (p == NULL) {
		return -1;
	} else {
		t = p->get_height();
		return t;
	}
}

TNode * IntervallTree::srl(TNode * &p1) {
	TNode * p2;
	p2 = p1->left;
	p1->left = p2->right;
	p2->right = p1;
	p1->set_height(max(bsheight(p1->left), bsheight(p1->right)) + 1);
	p2->set_height(max(bsheight(p2->left), p1->get_height()) + 1);
	return p2;
}
TNode * IntervallTree::srr(TNode * &p1) {
	TNode * p2;
	p2 = p1->right;
	p1->right = p2->left;
	p2->left = p1;
	p1->set_height(max(bsheight(p1->left), bsheight(p1->right)) + 1);
	p2->set_height(max(p1->get_height(), bsheight(p2->right)) + 1);
	return p2;
}
TNode * IntervallTree::drl(TNode * &p1) {
	p1->left = srr(p1->left);
	return srl(p1);
}
TNode * IntervallTree::drr(TNode * &p1) {
	p1->right = srl(p1->right);
	return srr(p1);
}

int IntervallTree::nonodes(TNode * p) {
	int count = 0;
	if (p != NULL) {
		nonodes(p->left);
		nonodes(p->right);
		count++;
	}
	return count;
}

void IntervallTree::collapse_intervalls(TNode *&p) {
	std::cout << "\t Collapse" << std::endl;
	TNode * new_root = NULL;
	std::vector<Breakpoint *> points;
	get_breakpoints(p, points);

	for (size_t i = 0; i < points.size(); i++) {
		if (points[i]->get_support() > Parameter::Instance()->min_support) {
			//std::cout << "\tpoints: " << points[i]->to_string(ref) << std::endl;
			this->insert(points[i], new_root);
		}
	}
	this->clear(p);
	p = new_root;
}
