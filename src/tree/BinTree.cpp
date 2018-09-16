/*
 * BinTree.cpp
 *
 *  Created on: Sep 3, 2015
 *      Author: fsedlaze
 */

#include "BinTree.h"

void BinTree::find(int item, tree_node **par, tree_node **loc) {
	tree_node *ptr, *ptrsave;
	if (root == NULL) {
		*loc = NULL;
		*par = NULL;
		return;
	}

	if (item == root->key) {
		*loc = root;
		*par = NULL;
		return;
	}

	if (item < root->key) {
		ptr = root->left;
	} else {
		ptr = root->right;
	}
	ptrsave = root;
	while (ptr != NULL) {
		if (item == ptr->key) {
			*loc = ptr;
			*par = ptrsave;
			return;
		}
		ptrsave = ptr;
		if (item < ptr->key) {
			ptr = ptr->left;
		} else {
			ptr = ptr->right;
		}
	}
	*loc = NULL;
	*par = ptrsave;
}

/*

 * Inserting Element into the Tree

 */

void BinTree::insert(tree_node *tree, int value) {
	if (root == NULL) {
		root = new tree_node;
		root->key = value;
		root->num = 1;
		root->left = NULL;
		root->right = NULL;
		std::cout << "Root tree_node is Added" << std::endl;
		return;
	}

	if (tree->key > value) {
		if (tree->left != NULL) {
			insert(tree->left, value);
		} else {

			tree->left = new tree_node;
			tree->left->key = value;
			tree->left->num = 1;
			(tree->left)->left = NULL;
			(tree->left)->right = NULL;
			std::cout << "tree_node Added To Left" << std::endl;
			return;
		}
	} else if (tree->key < value) {
		if (tree->right != NULL) {
			insert(tree->right, value);
		} else {
			tree->right = new tree_node;
			tree->right->key = value;
			tree->right->num = 1;
			(tree->right)->left = NULL;
			(tree->right)->right = NULL;
			std::cout << "tree_node Added To Right" << std::endl;
			return;
		}
	} else { // found element -> already exist!
		tree->num++;
	}
}

/*

 * Delete Element from the tree

 */

void BinTree::del(int key) {
	tree_node *parent, *location;

	if (root == NULL) {
		std::cout << "Tree empty" << std::endl;
		return;
	}

	find(key, &parent, &location);
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

void BinTree::case_a(tree_node *par, tree_node *loc) {
	if (par == NULL) {
		root = NULL;
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

void BinTree::case_b(tree_node *par, tree_node *loc) {
	tree_node *child;
	if (loc->left != NULL) {
		child = loc->left;
	} else {
		child = loc->right;
	}
	if (par == NULL) {
		root = child;
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

void BinTree::case_c(tree_node *par, tree_node *loc) {
	tree_node *ptr, *ptrsave, *suc, *parsuc;
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
		root = suc;
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

void BinTree::get_nodes(tree_node *ptr, std::vector<int> & nodes) {
	std::cout<<"get_nodes"<<std::endl;
	if (root == NULL) {
		std::cout << "Tree is empty" << std::endl;
		return;
	}
	if (ptr != NULL) {
		nodes.push_back(ptr->key);
		get_nodes(ptr->left,nodes);
		get_nodes(ptr->right,nodes);
	}
}

/*

 * Pre Order Traversal

 */

void BinTree::preorder(tree_node *ptr) {
	if (root == NULL) {
		std::cout << "Tree is empty" << std::endl;
		return;
	}
	if (ptr != NULL) {
		std::cout << ptr->key << "  ";
		preorder(ptr->left);
		preorder(ptr->right);
	}
}

/*

 * In Order Traversal

 */

void BinTree::inorder(tree_node *ptr) {
	if (root == NULL) {
		std::cout << "Tree is empty" << std::endl;
		return;
	}
	if (ptr != NULL) {
		inorder(ptr->left);
		std::cout << ptr->key << "  ";
		inorder(ptr->right);
	}
}

/*

 * Postorder Traversal

 */

void BinTree::postorder(tree_node *ptr) {
	if (root == NULL) {
		std::cout << "Tree is empty" << std::endl;
		return;
	}
	if (ptr != NULL) {
		postorder(ptr->left);
		postorder(ptr->right);
		std::cout << ptr->key << "  ";
	}
}

/*

 * Display Tree Structure

 */
void BinTree::display(tree_node *ptr, int level) {
	int i;
	if (ptr != NULL) {
		display(ptr->right, level + 1);
		std::cout << std::endl;
		if (ptr == root)
			std::cout << "Root->:  ";
		else {
			for (i = 0; i < level; i++) {
				std::cout << "       ";
			}
		}
		std::cout << ptr->key;
		display(ptr->left, level + 1);
	}
}
