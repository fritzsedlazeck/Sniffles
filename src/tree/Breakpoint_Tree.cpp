/*
 * Breakpoint_Tree.cpp
 *
 *  Created on: Mar 28, 2016
 *      Author: fsedlaze
 */

#include "Breakpoint_Tree.h"
void Breakpoint_Tree::find(int position, std::string chr, breakpoint_node *par, breakpoint_node *&loc) {
	if (par == NULL) { //not found
		loc = NULL;
		par = NULL;
		return;
	}
	if (position == par->position && strcmp(chr.c_str(), par->chr.c_str()) == 0) { //found
		loc = par;
		par = NULL;
		return;
	}
	//search goes on:
	if (position < par->position) {
		find(position, chr, par->left, loc);
	} else {
		find(position, chr, par->right, loc);
	}
}

void Breakpoint_Tree::overalps(int start, int stop, std::string chr, breakpoint_node *par, bool SV_support) {
	//start + stop: read coordinates.
	if (par == NULL) { //not found
		return;
	}
	if (par->direction) { //start
		if ((par->position-100 > start && par->position+100 < stop) && strcmp(chr.c_str(), par->chr.c_str()) == 0) { //found
			par->ref_support++;
	//		std::cout<<"start: "<<start<<" "<<stop<<std::endl;
		}
	} else { //stop coordinate
		if ((par->position > start+100 && par->position < stop-100) && strcmp(chr.c_str(), par->chr.c_str()) == 0) { //found
			par->ref_support++;
	//		std::cout<<"stop: "<< start<<" "<<stop<<std::endl;
		}
	}

	//search goes on:
	if (start < par->position) {
		overalps(start, stop, chr, par->left, SV_support);
	} else {
		overalps(start, stop, chr, par->right, SV_support);
	}
}

/*
 * Inserting Element into the Tree
 */
void Breakpoint_Tree::insert(breakpoint_node *&tree, std::string chr, int position, bool direction) {
	if (tree == NULL) {
		tree = new breakpoint_node;
		tree->position = position;
		tree->ref_support = 0;
		tree->chr = chr;
		tree->direction = direction;
		tree->left = NULL;
		tree->right = NULL;
	} else if (tree->position > position) {
		insert(tree->left, chr, position, direction);
	} else if (tree->position < position) {
		insert(tree->right, chr, position, direction);
	} else if (strcmp(chr.c_str(), tree->chr.c_str()) == 0) { // found element -> already exist!
		//std::cerr << "Element exists!" << std::endl;
		//TODO we should use this information to assess the reliability of this call!
	} else {
		insert(tree->left, chr, position,position); //think about that!
	}
}

int Breakpoint_Tree::get_ref(breakpoint_node *&tree, std::string chr, int position) {
	if (tree == NULL) {
		return -1;
	}
	if (tree->position > position) {
		return get_ref(tree->left, chr, position);
	} else if (tree->position < position) {
		return get_ref(tree->right, chr, position);
	} else if (strcmp(chr.c_str(), tree->chr.c_str()) == 0) { // found element
		return tree->ref_support;
	} else {
		return get_ref(tree->left, chr, position); //just in case.
	}
}
/*

 * Delete Element from the tree

 */

void Breakpoint_Tree::del(int position, std::string chr) {
	breakpoint_node *parent, *location;

	if (parent == NULL) {
		std::cout << "Tree empty" << std::endl;
		return;
	}

	find(position, chr, parent, location);
	if (location == NULL) {
		std::cout << "Item not present in tree" << std::endl;
		return;
	}

	if (location->left == NULL && location->right == NULL) {
		case_a(parent, location);
	}
	if (location->left != NULL && location->right == NULL) {
		case_b(parent, location);
	}
	if (location->left == NULL && location->right != NULL) {
		case_b(parent, location);
	}
	if (location->left != NULL && location->right != NULL) {
		case_c(parent, location);
	}
	delete location;
}

/*

 * Case A

 */

void Breakpoint_Tree::case_a(breakpoint_node *par, breakpoint_node *loc) {
	if (par == NULL) {
		loc = NULL;
	} else {
		if (loc == par->left) {
			par->left = NULL;
		} else {
			par->right = NULL;
		}
	}
}

/*

 * Case B

 */

void Breakpoint_Tree::case_b(breakpoint_node *par, breakpoint_node *loc) {
	breakpoint_node *child;
	if (loc->left != NULL) {
		child = loc->left;
	} else {
		child = loc->right;
	}
	if (par == NULL) {
		loc = child;
	} else {
		if (loc == par->left) {
			par->left = child;
		} else {
			par->right = child;
		}
	}
}

/*

 * Case C

 */

void Breakpoint_Tree::case_c(breakpoint_node *par, breakpoint_node *loc) {
	breakpoint_node *ptr, *ptrsave, *suc, *parsuc;
	ptrsave = loc;
	ptr = loc->right;
	while (ptr->left != NULL) {
		ptrsave = ptr;
		ptr = ptr->left;
	}
	suc = ptr;
	parsuc = ptrsave;
	if (suc->left == NULL && suc->right == NULL) {
		case_a(parsuc, suc);
	} else {
		case_b(parsuc, suc);
	}

	if (par == NULL) {
		loc = suc;
	} else {
		if (loc == par->left) {
			par->left = suc;
		} else {
			par->right = suc;
		}
	}
	suc->left = loc->left;
	suc->right = loc->right;
}

/*

 * Pre Order Traversal

 */

void Breakpoint_Tree::preorder(breakpoint_node *ptr) {
	if (ptr == NULL) {
		std::cout << "Tree is empty" << std::endl;
		return;
	}
	if (ptr != NULL) {
		std::cout << ptr->position << "  ";
		preorder(ptr->left);
		preorder(ptr->right);
	}
}

/*

 * In Order Traversal

 */

void Breakpoint_Tree::inorder(breakpoint_node *ptr) {
	if (ptr == NULL) {
		std::cout << "Tree is empty" << std::endl;
		return;
	}
	if (ptr != NULL) {
		inorder(ptr->left);
		std::cout << ptr->chr << " " << ptr->position << "  " << ptr->ref_support << std::endl;
		inorder(ptr->right);
	}
}

/*

 * Postorder Traversal

 */

void Breakpoint_Tree::postorder(breakpoint_node *ptr) {
	if (ptr == NULL) {
		return;
	}
	if (ptr != NULL) {
		postorder(ptr->left);
		postorder(ptr->right);
		std::cout << ptr->position << "  ";
	}
}

