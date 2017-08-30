/*
 * Breakpoint_Tree.h
 *
 *  Created on: Mar 28, 2016
 *      Author: fsedlaze
 */

#ifndef TREE_BREAKPOINT_TREE_H_
#define TREE_BREAKPOINT_TREE_H_
#include <vector>
#include <cstddef>
#include <iostream>
#include <string>
#include <string.h>
struct breakpoint_node {
	std::string chr;
	int position; // value to store!
	bool direction;
	int ref_support;
	breakpoint_node *left;
	breakpoint_node *right;
};

class Breakpoint_Tree {
private:

public:
	Breakpoint_Tree() {
	}
	~Breakpoint_Tree(){
	}
	void find(int position,std::string chr, breakpoint_node *par, breakpoint_node *&loc);
	void insert(breakpoint_node *&tree, std::string chr,int position,bool direction);
	void del(int position,std::string chr);
	void case_a(breakpoint_node *par, breakpoint_node *loc);
	void case_b(breakpoint_node *par, breakpoint_node *loc);
	void case_c(breakpoint_node *par, breakpoint_node *loc);
	void preorder(breakpoint_node *ptr);
	void inorder(breakpoint_node *ptr);
	void postorder(breakpoint_node *ptr);
	void display(breakpoint_node *ptr, int);
	void get_nodes(breakpoint_node *ptr, std::vector<int> & nodes);
	void overalps(int start,int stop,std::string chr, breakpoint_node *par, bool SV_support);
	int get_ref(breakpoint_node *&tree, std::string chr, int position);
};





#endif /* TREE_BREAKPOINT_TREE_H_ */
