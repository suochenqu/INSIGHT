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
Generate a linear code of code-length n.
(for left and right barcodes each of length 5, take n = 10)

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

codes = linear_generation(10) # `codes` contains 5-nt left and right barcode pairs
len(codes)  # 16384
codes = linear_generation(12) # `codes` contains 6-nt left and right barcode pairs
len(codes)  # 262144





