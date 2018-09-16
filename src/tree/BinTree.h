/*
 * BinTree.h
 *
 *  Created on: Sep 3, 2015
 *      Author: fsedlaze
 */

#ifndef TREE_BINTREE_H_
#define TREE_BINTREE_H_
struct tree_node {
	int key; // value to store!
	int num;//times hit 1-> unique
	struct tree_node *left;
	struct tree_node *right;
};
#include <vector>
#include <cstddef>
#include <iostream>
class BinTree {
private:
	tree_node *root;
public:
	BinTree() {
		root = NULL;
	}
	~BinTree(){

	}
	void find(int item, tree_node **par, tree_node **loc);
	void insert(tree_node *tree, int value);
	void del(int key);
	void case_a(tree_node *par, tree_node *loc);
	void case_b(tree_node *par, tree_node *loc);
	void case_c(tree_node *par, tree_node *loc);
	void preorder(tree_node *ptr);
	void inorder(tree_node *ptr);
	void postorder(tree_node *ptr);
	void display(tree_node *ptr, int);
	void get_nodes(tree_node *ptr, std::vector<int> & nodes);
};

#endif /* TREE_BINTREE_H_ */
