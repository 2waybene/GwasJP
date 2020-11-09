'''
Credit: https://www.geeksforgeeks.org/iterative-preorder-traversal/
'''
# Python program to perform iterative preorder traversal

# A binary tree node
class Node:

    # Constructor to create a new node
    def __init__(self, data):
        self.data = data
        self.left = None
        self.right = None


# An iterative process to print preorder traveral of BT
def iterativePreorder(root):

    # Base CAse
    if root is None:
        return

    # create an empty stack and push root to it
    nodeStack = []
    nodeStack.append(root)

    # Pop all items one by one. Do following for every popped item
    # a) print it
    # b) push its right child
    # c) push its left child
    # Note that right child is pushed first so that left
    # is processed first */
    while(len(nodeStack) > 0):

        # Pop the top item from stack and print it
        node = nodeStack.pop()
        print (node.data)

        # Push right and left children of the popped node
        # to stack
        if node.right is not None:
            nodeStack.append(node.right)
        if node.left is not None:
            nodeStack.append(node.left)

# Iterative function for inorder tree traversal
def inOrder(root):

    # Set current to root of binary tree
    current = root
    stack = [] # initialize stack
    done = 0

    while True:

        # Reach the left most Node of the current Node
        if current is not None:

            # Place pointer to a tree node on the stack
            # before traversing the node's left subtree
            stack.append(current)

            current = current.left


        # BackTrack from the empty subtree and visit the Node
        # at the top of the stack; however, if the stack is
        # empty you are done
        elif(stack):
            current = stack.pop()
            print(current.data, end=" ") # Python 3 printing

            # We have visited the node and its left
            # subtree. Now, it's right subtree's turn
            current = current.right

        else:
            break

    print()

# Driver program to test above function




def main():
    # Driver program to test above function
    root = Node(10)
    root.left = Node(8)
    root.right = Node(2)
    root.left.left = Node(3)
    root.left.right = Node(5)
    root.right.left = Node(2)
    iterativePreorder(root)

    inOrder(root)

    """ Constructed binary tree is 
            1 
          /   \ 
         2     3 
       /  \ 
      4    5   
    """

    root = Node(1)
    root.left = Node(2)
    root.right = Node(3)
    root.left.left = Node(4)
    root.left.right = Node(5)

    inOrder(root)
    iterativePreorder(root)

# This code is contributed by Nikhil Kumar Singh(nickzuck_007)

if __name__ == "__main__":
    main()
