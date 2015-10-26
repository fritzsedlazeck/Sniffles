/*
 * Scapegoat.h
 *
 *  Created on: Jun 23, 2015
 *      Author: fsedlaze
 */

#ifndef TREE_SCAPEGOAT_H_
#define TREE_SCAPEGOAT_H_
#include "TNode.h"
/*
struct TNode {
	TNode * parent;
	TNode * left;
	TNode * right;
	Breakpoint * data;
	int value;
};*/

/*
class Scapegoat {
private:
	TNode * root;
	int n, q;
public:
	Scapegoat() {
		root = NULL;
		n = 0;
	}

	bool isEmpty() {
		return root == NULL;
	}

	void makeEmpty() {
		root = NULL;
		n = 0;
	}

	int size(TNode *r) {
		if (r == NULL)
			return 0;
		else {
			int l = 1;
			l += size(r->left);
			l += size(r->right);
			return l;
		}
	}

	bool search(int value, Breakpoint * point) {
		return search(root,value, point);
	}
	bool search(TNode *r, int value, Breakpoint * point) {
		bool found = false;
		while ((r != NULL) && !found) {

			//int score=r->get_data()->overlap(point);
			//if (score>0) {

			if(r->get_value() > value){
				r = r->left;
		//	} else if (score<0) {
			}else if (r->get_value() < value){
				r = r->right;
			} else {
				found = true;
				break;
			}
			found = search(r, value,point);
		}
		return found;
	}

	int size() {
		return n;
	}
	void inorder() {
		inorder(root);
	}

	void inorder(TNode *r) {
		if (r != NULL) {
			inorder(r->left);
			std::cout << r->get_value() ;//->get_data()->to_string() << "   ";
			if(r==this->root){
				std::cout<<"* ";
			}else{
				std::cout<<" ";
			}
			inorder(r->right);
		}
	}


	void preorder() {
		preorder(root);
	}

	void preorder(TNode *r) {

		if (r != NULL) {
			std::cout << r->get_value() <<" ";//get_data()->to_string() << "   ";
			preorder(r->left);
			preorder(r->right);
		}

	}

	void postorder() {
		postorder(root);
	}

	void postorder(TNode *r) {
		if (r != NULL) {
			postorder(r->left);
			postorder(r->right);
			std::cout << r->get_value() <<" ";//->get_data()->to_string() << "   ";
		}
	}

	static int const log32(int q) {
		double const log23 = 2.4663034623764317;
		return (int) ceil(log23 * log(q));
	}

	bool add(int value,Breakpoint * point) {
		TNode *u = new TNode(value,point);
		int d = addWithDepth(u);
		if (d > log32(q)) {
			TNode *w = u->parent;
			while (3 * size(w) <= 2 * size(w->parent)) {
				w = w->parent;
			}
			rebuild(w->parent);
		}
		return d >= 0;
	}


	void rebuild(TNode *u) {
		int ns = size(u);
		TNode *p = u->parent;
		TNode **a = new TNode*[ns];
		packIntoArray(u, a, 0);
		if (p == NULL) {
			root = buildBalanced(a, 0, ns);
			root->parent = NULL;
		} else if (p->right == u) {
			p->right = buildBalanced(a, 0, ns);
			p->right->parent = p;
		} else {
			p->left = buildBalanced(a, 0, ns);
			p->left->parent = p;
		}
	}


	int packIntoArray(TNode *u, TNode *a[], int i) {

		if (u == NULL) {
			return i;
		}
		i = packIntoArray(u->left, a, i);
		a[i++] = u;
		return packIntoArray(u->right, a, i);
	}


	TNode *buildBalanced(TNode **a, int i, int ns)

	{
		if (ns == 0) {
			return NULL;
		}
		int m = ns / 2;
		a[i + m]->left = buildBalanced(a, i, m);
		if (a[i + m]->left != NULL) {
			a[i + m]->left->parent = a[i + m];
		}
		a[i + m]->right = buildBalanced(a, i + m + 1, ns - m - 1);\
		if (a[i + m]->right != NULL) {
			a[i + m]->right->parent = a[i + m];
		}
		return a[i + m];
	}


	int addWithDepth(TNode *u) {
		TNode *w = root;
		if (w == NULL) {
			root = u;
			n++;
			q++;
			return 0;
		}
		bool done = false;
		int d = 0;
		do {
			//int score=u->get_data()->overlap(w->get_data());
			//if (score>0) {

			if(w->get_value() > u->get_value()){
				if (w->left == NULL) {
					w->left = u;
					u->parent = w;
					done = true;
				} else {
					w = w->left;
				}
			//} else if (score<0) {
			}else if(w->get_value() < u->get_value()){
				if (w->right == NULL) {
					w->right = u;
					u->parent = w;
					done = true;
				} else {

					w = w->right;
				}
			} else {
				//are equal!
				return -1;
			}
			d++;
		} while (!done);
		n++;
		q++;
		return d;
	}

};
*/
#endif /* TREE_SCAPEGOAT_H_ */
