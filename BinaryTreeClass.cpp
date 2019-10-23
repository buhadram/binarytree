/* Binary TreeNode */

#include <iostream>
#include <deque>
#include <climits>
#include <vector>
#include <algorithm>
#include <utility>
#include <stack>
#include <array>
#include <random>
#include <memory>

using namespace std;

template <typename T>
class Tree
{
public:
	struct TreeNode
	{
		TreeNode(T val, TreeNode* parent=nullptr) :
			m_data(val),
			m_left(nullptr),
			m_right(nullptr),
			m_parent(parent)
		{}

		T 			 m_data;
		TreeNode 	*m_left;
		TreeNode 	*m_right;  
		TreeNode 	*m_parent;
	};

	// default ctor
    Tree() : m_root(nullptr)
	{}

	// copy ctor
	Tree(Tree const & rhs)
	{
		copy(rhs.m_root);
	}

	// dtor
    virtual ~Tree()
	{
		deleteNodes(m_root);
		m_root = nullptr;
	}

    
    TreeNode *leftMost()
    {
        TreeNode* node = m_root;
		while (node->m_right)
		{
			if (node->m_left)
			{
				node = node->m_left;
				break;
			}
		}
        while (node->m_left != nullptr)
		{
            node = node->m_left;
		}
        return node;
    }

    TreeNode *insert(T data)
    {
        TreeNode *newNode = new TreeNode(data);
        if (m_root == nullptr)
        {
			m_root = newNode;
        }
		else
		{
			bool left = false;
			TreeNode* parent = nullptr;
			TreeNode* searchNode = m_root;
			while (searchNode != nullptr)
			{
				parent = searchNode;
				if (newNode->m_data <= searchNode->m_data)
				{
					searchNode = searchNode->m_left;
					left = true;
				}
				else
				{
					searchNode = searchNode->m_right;
					left = false;
				}
			}
			newNode->m_parent = parent;
			if (left)
			{
				parent->m_left = newNode;
			}
			else
			{
				parent->m_right = newNode;
			}
		}
        return newNode;

    }

	void deleteNodes(TreeNode *node)
	{
		if (node == nullptr)
		{
			return;
		}
		TreeNode* parent = node->m_parent;

		// delete children
		deleteNodes(node->m_right);
		deleteNodes(node->m_left);

		if (parent != nullptr)
		{
			if (parent->m_left == node)
			{
				parent->m_left == nullptr;
			}
			if (parent->m_right == node)
			{
				parent->m_right == nullptr;
			}
		}

		// delete the node
		/*
		std::cout << "Deleting node" << std::hex << node 
					<< " ('" << node->m_data << "')" << std::endl;
		*/
		delete node;
	}

	
	TreeNode *find(T key) { return search(m_root, key); }

	bool isBST() const
	{ 
		//std::cout << "m_root = " << std::hex << m_root << std::endl;
		return (isBST(m_root) == 1);
	}

	int size() const { return size(m_root); }
	int size() { return size(m_root); }

	int maxDepth() 		const { return maxDepth(m_root); }
	int minDepth() 		const { return minDepth(m_root); }
	int maxDepth() 	          { return maxDepth(m_root); }
	int minDepth() 	          { return minDepth(m_root); }
		
	TreeNode* maxTree() const { return maxTree(m_root); }
	TreeNode* minTree() const { return minTree(m_root); }
	TreeNode* maxTree()       { return maxTree(m_root); }
	TreeNode* minTree()       { return minTree(m_root); }

	TreeNode* getRoot() { return m_root; }

    // swap two subtrees/nodes
    void swap(TreeNode** n1, TreeNode** n2)
    {
        if (n1 != nullptr && n2 != nullptr)
        {
            std::cout << "swapping n1=" << std::hex << n1 << " with n2=" << std::hex << n2 << std::endl;
            TreeNode* tmp = *n1;
            *n1 = *n2;
            *n2 = tmp;
        }
    }

    // swap left & right subtree
    void swap(TreeNode* node)
    {
        if (node)
        {
            swap(&node->m_left, &node->m_right);
        }
    }

    // swap tree's root
    void swap(Tree& rhs)
    {
        std::cout << "swapping trees\n";
        swap(&m_root, &rhs.m_root);
    }


	template<class InputIt, class UnaryFunction>
	UnaryFunction for_each(InputIt first, UnaryFunction func)
	{
		// using non-recursive nor stack (Morris Traversal)
		TreeNode* current;

		current = first;
		while (current != nullptr)
		{
			if (current->m_left == nullptr)
			{
				func(current->m_data);
				current = current->m_right;
			}
			else
			{
				// find the inorder predecessor of current
				TreeNode* pre = current->m_left;
				while (pre->m_right != nullptr && pre->m_right != current)
				{
					pre = pre->m_right;
				}
				// make current as right child of its inorder predecesor
				if (pre->m_right == nullptr)
				{
					pre->m_right = current;
					current = current->m_left;
				}
				else
				{
					// revert the change made in if part to restore the original tree
					pre->m_right = nullptr;
					func(current->m_data);
					current = current->m_right;
				}
			} // while
		}

		return func;
	}

	bool isBalanced()
	{
		if(maxDepth(m_root)-minDepth(m_root) <= 1) 
			return true;
		else
			return false;
	}

	// operators
	Tree& operator=(Tree &other)
	{
		std::cout << __FUNCTION__ << "(): other=" << std::hex << other.getRoot() << std::endl;

		// check for self-assignment
        if (&other == this)
            return *this;

		copy(other.getRoot());
		
        return *this;
	}

	bool operator==(Tree& rhs)
	{
		// TODO
		//FIXME
		if (m_root == nullptr || rhs.m_root == nullptr)
		{
			return false;
		}
		// same root means both are identical
		if (rhs.m_root == m_root)
		{
			return true;
		}

		if (size() != rhs.size())
		{
			return false;
		}
	}

	// misc
	friend std::ostream& operator<<(std::ostream& os, TreeNode* node)
	{
		if (node != nullptr)
		{
			os << node->m_data;
		}
		else
		{
			os << "null";
		}

		return os;
	}

	std::ostream& operator<<(std::ostream& os)
	{
		//TODO
		return os;
	}

	friend std::ostream& operator<<(std::ostream& os, Tree& t)
	{
		t.for_each( t.getRoot(), 
				    [](T const &data){ 
						std::cout << ", " << data;
					} );
		return os;
	}

protected:

	int size(TreeNode *node)
	{
		if(node == nullptr) 
		{
			return 0;
		}
		else
		{
			return size(node->m_left) + 1 + size(node->m_right);
		}
	}


    void shallowCopy(Tree& rhs)
    {
        m_root = rhs.getRoot();
    }


	void copy(TreeNode *root)
	{
		std::cout << __FUNCTION__ << "(): root=" << std::hex << root << std::endl;

		if (root == nullptr || root == m_root)
		{
			// don't copy from empty tree or identical tree
			return;
		}

		#if USE_RECURSIVE

		// just using recursive preorder copy/insert
		insert(root->m_data);
		copy(root->m_left);
		copy(root->m_right);

		#else

		#if 1 
		// using non-recursive nor stack (Morris Traversal)
		for_each(root, [&](T data) { insert(data); });		
		#else

		// use in-order traversal without recursion        
		TreeNode* node = root;
		std::stack<TreeNode*> s;

		while (node != nullptr || !s.empty())
		{
			// reach the left most node of the current node
			while (node != nullptr)
			{
				//save the node into stack
				s.push(node);
				node = node->m_left;
			}

			// current node must be null at this point
			node = s.top();
			s.pop();
			insert(node->m_data);

			//we've visited the node and its subtree.
			// now it's right subtree's turn
			node = node->m_right;
		}
		#endif // Morris-Nostack

		#endif // recursive
	}

	TreeNode *search(TreeNode* root, T key)
    {
		#ifdef USE_RECURSIVE
        if ((root == nullptr) ||
            (root->m_data == key))
		{
            return root;
		}
        else
        {
            if (root->m_data < key)
                return lookUp(root->m_right, key);
            else
                return lookUp(root->m_left, key);
			
        }
		#else
		TreeNode *foundNode = nullptr;
        if ((root == nullptr) ||
            (root->m_data == key))
		{
            return root;
		}
		else
		{
			// Traverse untill root reaches to dead end 
			while (root != nullptr) 
			{ 
				// pass right subtree as new tree 
				if (key > root->m_data) 
					root = root->m_right; 

				// pass left subtree as new tree 
				else if (key < root->m_data) 
				{
					root = root->m_left;
				}
				else 
				{
					foundNode = root;//the key is found
				}
			} 
		}
		return foundNode;
		#endif
    }


	int isBST(TreeNode *node) const
	{
		/*
		A binary search tree (BST) is a node based binary tree data
		structure which has the following properties.
	       	a) The left subtree of a node contains only nodes 
		       with keys less than the node’s key.
		   	b) The right subtree of a node contains only nodes 
		       with keys greater than the node’s key.
	    	c) Both the left and right subtrees must also be binary
			   search trees.
		*/
		if(node == nullptr) 
		{
			//std::cout << std::endl << "   Node is NPTR" << std::endl;
			return 0;
		}
		
		// false if left is > than node
		if (node->m_left != nullptr &&
		    node->m_left->m_data > node->m_data)
		{
			std::cout << "Left Node > This Node" << std::endl;
			return -1;
		}

		// false if, recursively, the left or right is not a BST
		if ( isBST(node->m_left)==-1 || isBST(node->m_right)==-1 )
		{
			return -1 ;
		}

		return 1;
	}


	int maxDepth(TreeNode *node) const
	{
		if (node == nullptr || 
		   (node->m_left == nullptr && node->m_right == nullptr))
			return 0;

		int leftDepth = maxDepth(node->m_left);
		int rightDepth = maxDepth(node->m_right);

		return leftDepth > rightDepth ? 
					leftDepth + 1 : rightDepth + 1;
	}

	int minDepth(TreeNode *node) const
	{
		if (node == nullptr || 
		   (node->m_left == nullptr && node->m_right == nullptr))
			return 0;

		int leftDepth = minDepth(node->m_left);
		int rightDepth = minDepth(node->m_right);

		return leftDepth < rightDepth ? 
					leftDepth + 1 : rightDepth + 1;
	}
	

	/* Tree Minimum */
	TreeNode* minTree(TreeNode *node)
	{
		if (node == nullptr) return nullptr;
		while (node->m_left) 
			node = node->m_left;
		return node;
	}

	/* Tree Maximum */
	TreeNode* maxTree(TreeNode *node)
	{
		if (node == nullptr) return nullptr;

		while(node->m_right) 
			node = node->m_right;
		return node;
	}

private:
	//std::shared_ptr<TreeNode>		m_root = nullptr;
    TreeNode                       *m_root = nullptr;
};



#if 0




/* In Order Successor - a node which has the next higher key */ 
TreeNode *succesorInOrder(TreeNode *node)
{
	/* if the node has right child, seccessor is Tree-Minimum */
	if(node->m_right != nullptr) return minTree(node->m_right);

	TreeNode *y = node->parent;
	while(y != nullptr && node == y->m_right) 
    {
	    node = y;
	    y = y->parent;
	}
	return y;
}

/* In Order Predecessor - a node which has the next lower key */
TreeNode *predecessorInOrder(TreeNode *node)
{
	/* if the node has left child, predecessor is Tree-Maximum */
	if(node->m_left != nullptr) return maxTree(node->m_left);

	TreeNode *y = node->parent;
	/* if it does not have a left child, 
	predecessor is its first left ancestor */
	while(y != nullptr && node == y->m_left) {
	    node = y;
	    y = y->parent;
	}
	return y;
}

void reverseOrderPrint(TreeNode *node)
{
	if(node == nullptr) return;
	if(node->m_left == nullptr && node->m_right == nullptr) {
		cout << node->data << " ";
		return;
	}
	
	reverseOrderPrint(node->m_right);
	cout << node->data << " ";
	reverseOrderPrint(node->m_left);
}
 
TreeNode *lowestCommonAncestor(TreeNode *node, TreeNode *p, TreeNode *q) 
{
	TreeNode *left, *right;
	if(node == nullptr) return nullptr;
	if(node->m_left == p || node->m_left == q
		|| node->m_right == p || node->m_right == q) return node;
	
	left = lowestCommonAncestor(node->m_left,p,q);
	right = lowestCommonAncestor(node->m_right, p,q);
	if(left && right) 
	    return node;
	else 
	    return (left) ? left : right;	
}

/* print tree in order */
/* 1. Traverse the left subtree. 
   2. Visit the root. 
   3. Traverse the right subtree. 
*/

void printTreeInOrder(TreeNode *node)
{
	if(node == nullptr) return;

	printTreeInOrder(node->m_left);
	cout << node->data << " ";
	printTreeInOrder(node->m_right);
}

/* print tree in postorder*/
/* 1. Traverse the left subtree. 
   2. Traverse the right subtree. 
   3. Visit the root. 
*/
void printTreePostOrder(TreeNode *node)
{
	if(node == nullptr) return;

	printTreePostOrder(node->m_left);
	printTreePostOrder(node->m_right);
	cout << node->data << " ";
}

/* print in preorder */
/* 1. Visit the root. 
   2. Traverse the left subtree. 
   3. Traverse the right subtree. 
*/
void printTreePreOrder(TreeNode *node)
{
	if(node == nullptr) return;

	cout << node->data << " ";
	printTreePreOrder(node->m_left);
	printTreePreOrder(node->m_right);
}

/* In reverse of printTreeInOrder() */
void printTreeReverseOrder(TreeNode *node)
{
	if(node == nullptr) return;
	if(node->m_left == nullptr && node->m_right == nullptr) {
	    cout << node->data << " ";
	    return;
	}
	
	printTreeReverseOrder(node->m_right);
	cout << node->data << " ";
	printTreeReverseOrder(node->m_left);
}
/* recursion routine to find path */
void pathFinder(TreeNode *node, int path[], int level)
{
	if(node == nullptr) return;
        // save leaf node
	if(node->m_left == nullptr && node->m_right == nullptr) {
	    path[level] = node->data;
	    for(int i = 0; i <= level; i++) {
		cout << (char)path[i];
	    }
	    cout << endl;
	    return;
	}
        // save parent node
	path[level] = node->data;
	pathFinder(node->m_left, path, level+1);
	pathFinder(node->m_right, path, level+1);
}

bool matchTree(TreeNode *r1, TreeNode *r2)
{
	/* Nothing left in the subTreeNode */
	if(r1 == nullptr && r2 == nullptr)
	    return true;
	/* Big tree empty and subtree not found */
	if(r1 == nullptr || r2 == nullptr)
	    return false;
	/* Not matching */
	if(r1->data != r2->data)
	    return false;
	return (matchTree(r1->m_left, r2->m_left) && 
			matchTree(r1->m_right, r2->m_right));
}

bool subTree(TreeNode *r1, TreeNode *r2)
{
	/*Big tree empty and subtree not found */
	if(r1 == nullptr)
	    return false;
	if(r1->data == r2->data)
	    if(matchTree(r1, r2)) return true;
	return 
            (subTree(r1->m_left, r2) || subTree(r1->m_right, r2));
}

bool isSubTree(TreeNode *r1, TreeNode *r2)
{
	/* Empty tree is subTreeNode */
	if(r2 == nullptr) 
	    return true;
	else
	    return subTree(r1, r2);
}

/* change a tree so that the roles of the left 
and right hand pointers are swapped at every node */
void mirror(TreeNode *r)
{
	if(r == nullptr) return;
	
	mirror(r->m_left);
	mirror(r->m_right);

	/* swap pointers */
	std::swap(r->m_left, r->m_right);
}

/* create a new tree from a sorted array */
TreeNode *addToBST(char arr[], int start, int end)
{
	if(end < start) return nullptr;
	int mid = (start + end)/2;

	TreeNode *r = new Tree;
	r->data = arr[mid];
	r->m_left = addToBST(arr, start, mid-1);
	r->m_right = addToBST(arr, mid+1, end);
	return r;
}

TreeNode *createMinimalBST(char arr[], int size)
{
	return addToBST(arr,0,size-1);
}

/* Breadth first traversal using queue */
void BreadthFirstTraversal(TreeNode *root)
{
	if (root == nullptr) return;
	deque <TreeNode *> queue;
	queue.push_back(root);

	while (!queue.empty()) {
	    TreeNode *p = queue.front();
	    cout << p->data << " ";
	    queue.pop_front();

	    if (p->m_left != nullptr)
		queue.push_back(p->m_left);
	    if (p->m_right != nullptr)
		queue.push_back(p->m_right);
	}
	cout << endl;
} 

/* get the level of a node element: root level = 0 */
int getLevel(TreeNode *node, int elm, int level)
{
	if(node == nullptr) return 0;
	if(elm == node->data) 
	    return level;
	else if(elm < node->data) 
	    return getLevel(node->m_left, elm, level+1);
	else 
	    return getLevel(node->m_right, elm, level+1);
}

/* This code prints out all nodes at the same depth (level) */
void BreadthFirst_LevelElement_Print
               (TreeNode *root, vector<vector<int> > &v)
{
	if(root == nullptr) return;
	deque<TreeNode *> q;
	q.push_back(root);
	while(!q.empty()) {
	    TreeNode *p =  q.front();
	    int lev = getLevel(root, p->data, 0);
	    v[lev].push_back(p->data);
	    q.pop_front();
	    if(p->m_left) q.push_back(p->m_left);
	    if(p->m_right)q.push_back(p->m_right);
	}
	return;
}

/* levelPrint()
prints nodes at the same level
This is simpler than the BreadthFirstTraversal(root) above
It takes 2D vector with the same size of level (= MaxDepth+1)
and fills elements as we traverse (preOrder) */

void levelPrint(TreeNode *node, vector<vector<char> >&elm, int level)
{
	if(node == nullptr) return;
	// leaf nodes
	if(node->m_left == nullptr && node->m_right == nullptr) {
	    elm[level].push_back(node->data);
	    return;
	}
	// other nodes
	elm[level++].push_back(node->data);
	levelPrint(node->m_left, elm, level);
	levelPrint(node->m_right, elm, level);
}

/* find n-th max node from a TreeNode */
void NthMax(TreeNode* t)
{        
	static int n_th_max = 5;
	static int num = 0;
	if(t == nullptr) return;        
	NthMax(t->m_right);        
	num++;        
	if(num == n_th_max)                
	    cout << n_th_max << "-th maximum data is " 
                 << t->data << endl;        
	NthMax(t->m_left);
}

/* Converting a BST into an Array */ 
void TreeToArray(TreeNode *node, int a[]){ 
	static int pos = 0; 
  
	if(node){ 
	    TreeToArray(node->m_left,a); 
	    a[pos++] = node->data; 
	    TreeToArray(node->m_right,a); 
       } 
} 

    /* Separate even/odd level elements */
    /* This function is using BFS */
    void level_even_odd(TreeNode *node)
    {
        vector<char> evenVec, oddVec;
        if (node == nullptr) return;
        deque<TreeNode*> que;
        que.push_back(node);

        while(!que.empty()) 
        {
            TreeNode *p = que.front();
            int level = getLevel(node, p->data, 0) ;
            // even level
            if (level % 2 == 0) 
            evenVec.push_back(p->data);
            else
            oddVec.push_back(p->data);
            que.pop_front();
            if(p->m_left)  que.push_back(p->m_left);
            if(p->m_right) que.push_back(p->m_right);
        }
        
        cout << "even level elements : ";
        for(int i = 0; i < evenVec.size(); i++) 
                cout << evenVec[i] << " ";
        cout << endl << "odd level elements : ";
            for(int i = 0; i < oddVec.size(); i++) 
                cout << oddVec[i] << " ";
        cout << endl;
    }
#endif

int main(int argc, char **argv)
{
	char ch, ch1, ch2;
	Tree<char> *found;
	Tree<char> *succ;
	Tree<char> *pred;
	Tree<char> *ancestor;
	char charArr[9] 
	    = {'A','B','C','D','E','F','G','H','I'};

	// prepare random seed
	std::random_device rd;
	std::seed_seq seed {rd(), rd(), rd(), rd(), rd(), rd(), rd(), rd()};
	std::mt19937 rng(seed);

	Tree<char> t;

	for(auto c : charArr)
		t.insert(c);

	t.insert('I');
	t.insert('B');
	t.insert('G');
	t.insert('H');
	t.insert('A');  
	t.insert('F');
	t.insert('E');
	t.insert('D');
	t.insert('C');

	/* is the tree BST? */
	std::cout << "Tree is " << (t.isBST() ? "" : "not") << " BST" << std::endl;
	std::cout << "size of Tree t = " << t.size() << std::endl;

	// copy 
	Tree<char> dupt = t;

	// operator=
	Tree<char> dupt2;
	dupt2 = t;

	std::cout << "dupt2 = " << dupt2 << std::endl;

	const int n = 100;
	std::uniform_int_distribution<int> dist(0,1000);

	Tree<int> intTree;
	for(auto i=0; i<n; ++i)
	{
		intTree.insert(dist(rng));
	}

	std::cout << "intTree = " << std::dec << intTree << std::endl;

	/* size of TreeNode */
	std::cout << "size = " << intTree.size() << std::endl;

	/* max depth */
	std::cout << "max depth = " << intTree.maxDepth() << endl;

	/* min depth */
	std::cout << "min depth = " << intTree.minDepth() << endl;

	/* balanced tree? */
	std::cout << "This tree is " << (intTree.isBalanced() ? "" : "not ") 
              << "balanced!" << std::endl;
	
	/* min value of the tree*/
	Tree<int>::TreeNode *tree = intTree.minTree();
	if (tree != nullptr)
		std::cout << "Min value = " << tree->m_data << std::endl;

	/* max value of the tree*/
	tree = intTree.maxTree();
	if (tree != nullptr)
		std::cout << "Max value = " << tree->m_data << std::endl;

    Tree<int> anotherTree;

    anotherTree.swap(intTree);
	std::cout << "intTree = " << std::dec << intTree << std::endl;
    std::cout << "anotherTree = " << std::dec << anotherTree << std::endl;

	#if 0
	/* get the level of a data: root level = 0 */
	std::cout << "Node B is at level: " << getLevel(root, 'B', 0) << endl;
	std::cout << "Node H is at level: " << getLevel(root, 'H', 0) << endl;
	std::cout << "Node F is at level: " << getLevel(root, 'F', 0) << endl;

    /* separate even/odd level elements */
    level_even_odd(root);

	ch = 'B';
	found = lookUp(root,ch);
	if(found) {
	    std::cout << "Min value of subtree " << ch << " as a root is "
			 << minTree(found)->data << endl;
	    std::cout << "Max value of subtree " << ch << " as a root is "
			 << maxTree(found)->data << endl;
	}

	ch = 'B';
	found = lookUp(root,ch);
	if(found) {
	    succ = succesorInOrder(found);
	    if(succ)
		cout << "In Order Successor of " << ch << " is "
		     << succesorInOrder(found)->data << endl;
	    else 
		cout << "In Order Successor of " << ch << " is None\n";
	}

	ch = 'E';
	found = lookUp(root,ch);
	if(found) {
	    succ = succesorInOrder(found);
	    if(succ)
		cout << "In Order Successor of " << ch << " is "
			 << succesorInOrder(found)->data << endl;
	    else 
		cout << "In Order Successor of " << ch << " is None\n";
	}

	ch = 'I';
	found = lookUp(root,ch);
	if(found) {
	    succ = succesorInOrder(found);
	    if(succ)
		cout << "In Order Successor of " << ch << " is "
			 << succesorInOrder(found)->data << endl;
	    else 
		cout << "In Order Successor of " << ch << " is None\n";
	}

	ch = 'B';
	found = lookUp(root,ch);
	if(found) {
	    pred = predecessorInOrder(found);
	    if(pred)
		cout << "In Order Predecessor of " << ch << " is "
			 << predecessorInOrder(found)->data << endl;
	    else 
		cout << "In Order Predecessor of " << ch << " is None\n";
	}

	ch = 'E';
	found = lookUp(root,ch);
	if(found) {
	    pred = predecessorInOrder(found);
	    if(pred)
		cout << "In Order Predecessor of " << ch << " is "
			 << predecessorInOrder(found)->data << endl;
	    else 
		cout << "In Order Predecessor of " << ch << " is None\n";
	}

	ch = 'I';
	found = lookUp(root,ch);
	if(found) {
	    pred = predecessorInOrder(found);
	    if(pred)
		cout << "In Order Predecessor of " << ch << " is "
			 << predecessorInOrder(found)->data << endl;
	    else 
		cout << "In Order Predecessor of " << ch << " is None\n";
	}

	/* Lowest Common Ancestor */
	ch1 = 'A';
	ch2 = 'C';
	ancestor = 
	    lowestCommonAncestor(root, 
			lookUp(root,ch1), lookUp(root,ch2));
	if(ancestor) 
	    cout << "The lowest common ancestor of " << ch1 << " and "
		<< ch2 << " is " << ancestor->data << endl;

	ch1 = 'E';
	ch2 = 'H';
	ancestor = 
	    lowestCommonAncestor(root, 
			lookUp(root,ch1), lookUp(root,ch2));
	if(ancestor) 
	    cout << "The lowest common ancestor of " << ch1 << " and "
		<< ch2 << " is " << ancestor->data << endl;

	ch1 = 'D';
	ch2 = 'E';
	ancestor = 
	    lowestCommonAncestor(root, 
			lookUp(root,ch1), lookUp(root,ch2));
	if(ancestor) 
	    cout << "The lowest common ancestor of " << ch1 << " and "
		<< ch2 << " is " << ancestor->data << endl;

	ch1 = 'G';
	ch2 = 'I';
	ancestor = 
	    lowestCommonAncestor(root, 
			lookUp(root,ch1), lookUp(root,ch2));
	if(ancestor) 
	    cout << "The lowest common ancestor of " << ch1 << " and "
		<< ch2 << " is " << ancestor->data << endl;

	ch1 = 'H';
	ch2 = 'I';
	ancestor = 
	    lowestCommonAncestor(root, 
			lookUp(root,ch1), lookUp(root,ch2));
	if(ancestor) 
	    cout << "The lowest common ancestor of " << ch1 << " and "
		<< ch2 << " is " << ancestor->data << endl;

	/* print tree in order */
	cout << "increasing sort order\n";
	printTreeInOrder(root);
	cout << endl;

	/* print tree in postorder*/
	cout << "post order \n";
	printTreePostOrder(root);
	cout << endl;

	/* print tree in preorder*/
	cout << "pre order \n";
	printTreePreOrder(root);
	cout << endl;

	/* print tree in reverse order*/
	cout << "reverse order \n";
	printTreeReverseOrder(root);
	cout << endl;

	/* lookUp */
	ch = 'D';
	found = lookUp(root,ch);
	if(found) 
	    cout << found->data << " is in the tree\n";
	else
	    cout << ch << " is not in the tree\n";

	/* lookUp */
	ch = 'M';
	found = lookUp(root,ch);
	if(found) 
	    cout << found->data << " is in the tree\n";
	else
	    cout << ch << " is not in the tree\n";

	/* printing all paths :
	Given a binary tree, print out all of its root-to-leaf 
	paths, one per line. Uses a recursive helper to do the work. */
	cout << "printing paths ..." << endl;
	int path[10];
	pathFinder(root, path, 0);

	/* find n-th maximum node */
	NthMax(root);

	/* Traversing level-order. 
	We visit every node on a level before going to a lower level. 
	This is also called Breadth-first traversal.*/
	cout << "printing with Breadth-first traversal" << endl;
	BreadthFirstTraversal(root);

	/* Prints all element at the same depth (level) */
	int row = maxDepth(root);
	vector<vector<int> > levVec(row+1);
	BreadthFirst_LevelElement_Print(root, levVec);
	for(int m = 0; m < levVec.size(); m++) {
	    cout << "Level at " << m << ": ";
	    for(int n = 0; n < levVec[m].size(); n++) 
		cout << (char)levVec[m][n] << " ";
	    cout << endl;
	}

	/* levelPrint()
	prints nodes at the same level
	This is simpler than the BreadthFirstTraversal(root) above
	It takes 2D vector (elm) with the same size of level (= MaxDepth+1)
	and fills elements as we traverse (preOrder) */
	vector<vector<char> > elm;
	int mxDepth = maxDepth(root);
	elm.resize(mxDepth+1);
	int level = 0;
	levelPrint(root, elm, level);
	cout << "levelPrint() " << endl;
	for(int i = 0; i <= mxDepth; i++) {
	    cout << "level " << i << ": " ;
	    for(int j = 0; j < elm[i].size(); j++) 
		cout << elm[i][j] << " ";
	    cout << endl;
	}

	/* convert the tree into an array */
	int treeSz = treeSize(root);
	int *array = new int[treeSz];
	TreeToArray(root,array);
	cout << "New array: ";
	for (int i = 0; i < treeSz; i++)
	    cout << (char)array[i] << ' ';
	cout << endl;
	delete [] array;

	/* subTreeNode */
	TreeNode *root2 = newTreeNode('D');
	insertTreeNode(root2,'C');  
	insertTreeNode(root2,'E');
	cout << "1-2 subtree: " << isSubTree(root, root2) << endl;

	TreeNode *root3 = newTreeNode('B');
	insertTreeNode(root3,'A');  
	insertTreeNode(root3,'D');
	insertTreeNode(root3,'C');  
	insertTreeNode(root3,'E');
	cout << "1-3 subtree: " << isSubTree(root, root3) << endl;

	TreeNode *root4 = newTreeNode('B'); 
	insertTreeNode(root4,'D');
	insertTreeNode(root4,'C');  
	insertTreeNode(root4,'E');
	cout << "1-4 subtree: " << isSubTree(root, root4) << endl;

	cout << "2-3 subtree: " << isSubTree(root2, root3) << endl;
	cout << "3-2 subtree: " << isSubTree(root3, root2) << endl;

	/* swap left and right */
	mirror(root);

	/* deleting a TreeNode */
	clear(root);

	/* make a new tree with minimal depth */
	TreeNode *newRoot = createMinimalBST(charArr,9);
#endif

	(void)argc;
	(void)argv;

	return 0;
}
