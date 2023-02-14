#!/usr/bin/env python
# -*- coding: UTF-8 -*-

#There are N students in a class. Some of them are friends, while some are not. Their friendship is transitive in nature, i.e., if A is friend of B and B is friend of C, then A is also friend of C. A friend circle is a group of students who are directly or indirectly friends.

#You are given a N×N − matrix M which consists of characters Y or N. If M[i][j]=Y, then ith and jth students are friends with each other, otherwise not. You have to print the total number of friend circles in the class.

#Each element of matrix friends will be Y or N.
#Number of rows and columns will be equal in the matrix.

#Example 1
#Input:
    #YYNN
    #YYYN
    #NYYN
    #NNNY
#Output: 2

import sys

### Union Find ######
class UnionFind:
    ## initialization of all the nodes in the graph
    # all the node has its unique label
    # n circles in the graph
    def __init__(self, n):
        self.label = range(n)
        self.n_circle = n

    ## find the label for a node
    def find(self, node):
        if self.label[node] == node:
            return node
        self.label[node] = self.find(self.label[node])
        return self.label[node]

    ## group two nodes with edges in between into a single circle
    ## meanwhile, keep track of the number of circles in the graph
    def union(self, node1, node2):
        if self.find(node1) != self.find(node2):
            self.label[self.find(node1)] = self.find(node2)
            self.n_circle -= 1
        return

    ## get the number of circles in the graph
    def get_n_circles(self):
        return self.n_circle

    ## check if two nodes are in the same circle
    def is_connected(self, node1, node2):
        return self.label[node1] == self.label[node2]

def friend_circles(friends):
    if not friends or len(friends) == 0:
        return 0

    n_friend = len(friends)
    friend_union = UnionFind(n_friend)

    for i in range(n_friend):
        if len(friends[i]) != n_friend:
            raise ValueError('The size of friends matrix is not valid!')
        for j in range(i, n_friend):
            if friends[i][j] == 'Y' and not friend_union.is_connected(i, j):
                friend_union.union(i, j)

    return friend_union.get_n_circles()

if __name__ == "__main__":
    friends = ["YNNNN", "NYNNN", "NNYNN", "NNNYN", "NNNNY"]
    print friend_circles(friends)
