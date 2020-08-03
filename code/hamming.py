#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
preamble
"""
import numpy as np

#%%
"""
Auxiliary functions and class 
"""

# hamming distance between two strings of equal length
def hamming_distance(string1, string2): 
    return sum([string1[i] != string2[i] for i in range(len(string1))])

# compute hamming ball of radius 2 around a barcode
def hamming_ball(word):
    n = len(word)
    ball = set()
    for i in range(n-1):
        for j in range(i+1, n):
            for c1 in 'ACGT':
                for c2 in 'ACGT':
                    tmp = word[:i] + c1 + word[(i+1):j] + c2 + word[(j+1):]
                    ball.add(tmp)
    return ball

# Nucleotide is a class that represents the four nucleotides 'A' 'C' 'G' 'T'
# as the field with four elements. Specifically, 'A' is the additive identity,
# 'C' the multiplicative identity.
class Nucleotide:
    additionMatrix = [[0, 1, 2, 3], 
                      [1, 0, 3, 2],
                      [2, 3, 0, 1],
                      [3, 2, 1, 0]]
    multiplicationMatrix = [[0, 0, 0, 0],
                            [0, 1, 2, 3],
                            [0, 2, 3, 1],
                            [0, 3, 1, 2]]
    divisionMatrix = [[5, 0, 0, 0],  # 4 ('-') codes for Inf, 5 ('*') codes for NaN
                     [4, 1, 2, 3],
                     [4, 3, 1, 2],
                     [4, 2, 3, 1]]
    def __init__(self, ind):
        self.ind = ind
    def __str__(self):
        return ['A', 'C', 'G', 'T', '-', '*'][self.ind]
    def __add__(self, other):
        return Nucleotide(Nucleotide.additionMatrix[self.ind][other.ind])
    def __mul__(self, other):
        return Nucleotide(Nucleotide.multiplicationMatrix[self.ind][other.ind])
    def __truediv__(self, other):
        return Nucleotide(Nucleotide.divisionMatrix[self.ind][other.ind])

# display number representation for Nucleotides, mx is either a matrix or a 
# vector
def show(mx):
    if len(mx.shape) == 2: # if mx is a matrix
        return [[letter.ind for letter in row] for row in mx]
    else: # if mx is a vector
        return [letter.ind for letter in mx]


#%%
"""
Method 1: Generated barcodes with pairwise hamming distance at least 3

Pick random unused barcode that does not lie in any hamming balls
of radius 2 from previously selected barcodes. This can achieve the Gilbert-
Varshamov lower bound of 4^n / (1 + 3n + 9n(n-1)/2). 
"""
def greedy_generation(n):
    all_words = ['']
    for i in range(n):
        all_words = [word + c for word in all_words for c in 'ACGT']
    unused = set()
    codes = set()
    
    # add all words to the set 'unused', for faster deletion operation
    for word in all_words:
        unused.add(word)
        
    while len(unused) > 0:
        next_word = unused.pop()
        codes.add(next_word)
        for word in hamming_ball(next_word):
            if word in unused:
                unused.remove(word)
    
    return codes 

codes = greedy_generation(10)
print(len(codes))  # 10958

codes = greedy_generation(12)
print(len(codes))  # 132832


#%%
"""
Method 2: generate a linear code of code-length n. 

First generate a parity matrix and the corresponding generator matrix. Then
generate barcodes as the linear span of rows of the generator matrix. This 
method achieves the strengthed Gilbert-Varshamov bound of q^k where k is the
largest integer satisfying q^(n-k) > 3n - 2. 

The code currently can handle 4 <= n <= 21
"""
def linear_generation(n):
    r = 3
    k = n - r
    # A list 21 unparallel vectors in F_4^3
    A = np.array([[0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
                  [0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3],
                  [1, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3]], dtype='int')
    # H is an r x n matrix using columns from A and having final three columns 
    # being identity. H is used to generate the parity matrix.
    H = np.concatenate((A[:, [i for i in np.arange(21) if i not in [5 ,1, 0]][-k:]],
                             np.identity(r, dtype='int')), axis=1)
    # parity = np.array([[Nucleotide(H[i,j]) for j in range(H.shape[1])] for i in range(H.shape[0])])
    
    # G is a k x n matrix obtained from H
    # G is used to construct the generator matrix
    G = np.concatenate((np.identity(k, dtype='int'), np.transpose(H[:, :k])), axis=1)
    generator = np.array([[Nucleotide(G[i,j]) for j in range(G.shape[1])] for i in range(G.shape[0])])
    
    # check generator and parity matrix have orthogonal rows
    # show(generator.dot(parity.T))
    
    # W is a 4^k x k matrix of all possible words using 0, 1, 2, 3
    # It is used to construct the weights matrix.
    W = [[]]
    for i in range(k):
        W = [row + [c] for row in W for c in [0,1,2,3]]
    weights = np.array([[Nucleotide(a) for a in vec] for vec in W])
    
    # Generate desired barcode set
    codes = weights.dot(generator)
    codes_str = [''.join([letter.__str__() for letter in code]) for code in codes]
    
    return(codes_str)

codes = linear_generation(10) 
len(codes)  # 16384
codes = linear_generation(12)
len(codes)  # 262144





