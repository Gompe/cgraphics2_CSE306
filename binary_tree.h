#include<iostream>

namespace binary_tree {

    template <typename T>
    class Node
    {
    public:
        T data;
        Node<T> *left = nullptr;
        Node<T> *right = nullptr;

        Node() = default;
        Node(T data) : data(data) {}

        bool is_leaf() const {
            return (left == nullptr) && (right == nullptr);
        }
    };


    template <typename T>
    class BinaryTree 
    {
    public:
        Node<T> *root = nullptr; 

        BinaryTree() = default;  
        BinaryTree(Node<T> *root) : root(root) {}

        BinaryTree(const BinaryTree& other) = delete;
        BinaryTree& operator=(const BinaryTree& other) = delete;

        void reset() {
            destroyRecursive(root);
        }

        void destroyRecursive(Node<T> *node) {
            if (node) {
                destroyRecursive(node->left);
                destroyRecursive(node->right);
                delete node;
            }
        }

        ~BinaryTree(){
            reset();
        }
    };
}