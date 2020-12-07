#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script can be used to generate a large number of Hamming-distance d
separated 4-ary codes (e.g. DNA barcodes) of length n. Code constructed via
linear block code. 

usage example: hamming.py -n 12 -d 3 -o 'codes.txt'
"""
#%%
"""
preamble
"""
import numpy as np
import scipy.special
import math
import itertools
import sys, getopt

#%%
"""
Define a class `Nucleotide` that encodes 'A C G T' by a finite field of 4 
elements GF(4). 
"""
# define arithmetics over the field GF(4)
additionMatrix          = [[0, 1, 2, 3], 
                           [1, 0, 3, 2],
                           [2, 3, 0, 1],
                           [3, 2, 1, 0]]
subtractionMatrix       = [[0, 1, 2, 3], 
                           [1, 0, 3, 2],
                           [2, 3, 0, 1],
                           [3, 2, 1, 0]]
multiplicationMatrix    = [[0, 0, 0, 0],
                           [0, 1, 2, 3],
                           [0, 2, 3, 1],
                           [0, 3, 1, 2]]
# divisionMatrix[i,j] is equal to i / j in GF(4)
# 4 ('-') codes for Inf, 5 ('*') codes for NaN
divisionMatrix          = [[5, 0, 0, 0],  
                           [4, 1, 3, 2],
                           [4, 2, 1, 3],
                           [4, 3, 2, 1]]

# Nucleotide is a class that represents the four nucleotides 'A' 'C' 'G' 'T'
# as the field with four elements. Specifically, 'A' is the additive identity,
# 'C' the multiplicative identity.
class Nucleotide:
    def __init__(self, ind):
        self.ind = ind
    def __add__(self, other):
        return Nucleotide(additionMatrix[self.ind][other.ind])
    def __sub__(self, other):
        return Nucleotide(subtractionMatrix[self.ind][other.ind])
    def __mul__(self, other):
        return Nucleotide(multiplicationMatrix[self.ind][other.ind])
    def __truediv__(self, other):
        return Nucleotide(divisionMatrix[self.ind][other.ind])
    def __eq__(self, other): 
        return self.ind == other.ind
    def __str__(self):
        return ['A', 'C', 'G', 'T', '-', '*'][self.ind]
    def __repr__(self):
        return str(self.ind)
        
# convert from a matrix/vector of GF(4) elements to the corresponding
# matrix/vector of {0,1,2,3}.
def nuc_to_ind(mx):
    if type(mx) is list:
        mx = np.array(mx)
    if len(mx.shape) == 2: # if mx is a matrix
        return np.array([[letter.ind for letter in row] for row in mx])
    else: # if mx is a vector
        return np.array([letter.ind for letter in mx])

# convert from a matrix/vector of indices to the corresponding matrix/vector
# of GF(4) elements
def ind_to_nuc(A_ind):
    if type(A_ind) is list: 
        A_ind = np.array(A_ind)
    if len(A_ind.shape) == 2:
        return np.array([[Nucleotide(A_ind[i,j]) for j in range(A_ind.shape[1])] for i in range(A_ind.shape[0])])
    else:
        return np.array([Nucleotide(ind) for ind in A_ind])


#%%
"""
Row reduced echelon form of a matrix over GF(4)
"""
def RREF(M):
    zero = M[0,0] - M[0,0]
    A = M.copy()
    n, p = A.shape
    # transform into a short and wide matrix
    if n > p: 
        A = A.T
        n, p = A.shape
        
    lead = 0
    for r in range(n):
        tf = A[r:n, lead:p].flatten('F') != zero
        if any(tf):
            loc = np.min(np.where(tf))
            loc_row = loc % (n-r) + r # row index of first nonzero entry in bottomright block
            loc_col = loc // (n-r) + lead # col index of first nonzero entry in bottomright block
            
            if loc_row != r: # swap rows if necessary to make sure i,loc_col-th entry is nonzero
                 A[[r, loc_row], :] = A[[loc_row, r], :] 
    
            # normalise row r and do row operation in subsequent rows
            A[r, :] = A[r, :] / A[r, loc_col] 
            if r == n-1:
                rank = n
                break
            else:
                A[range(r+1, n), :] -= A[range(r+1, n), loc_col:(loc_col+1)].dot(A[r:(r+1), :])
                lead = loc_col + 1
        else:
            rank = r
            break

    return (A, rank)

#%% 
"""
This function creates the generator matrix for a 4-ary linear code of word
length n and hamming distance >= d. The k x n generator matrix is obtained from
the r x n check matrix over GF(4). A matrix is a check matrix if no d-1 columns
are linearly dependent. As we are unaware of any constructive method the check
matrix, we generate an r x k random matrix over GF(4) and concatenate with I_r,
and then verify that all d-1 columns are independent. We first compute the 
smallest possible redundancy r based on Gilbert Varshamov bound. Then try
max_trials number of random check matrices. If none satisfy the check matrix
requirement, we increment r (which decreases the number of codes by 4 fold). 

Input: n - code length, d - min Hamming separation between codes
Output: G_nuc - generator matrix, H_nuc - check matrix, 
    k - linear code dimension, r - linear code redundancy
"""
def generator(n, d, max_trials=10000):
    # The distance d of a linear code C also equals the minimum number of 
    # linearly dependent columns of the check matrix H. 
    r = math.ceil(math.log(sum([scipy.special.binom(n-1, j) * 3**j for j in range(d-1)]), 4))
    print('Creating generator matrix for k =', n - r)
    while True:
        k = n - r
        for trial in range(max_trials):
            # A_ind is a r x k matrix
            # H_nuc is an r x n check matrix [A_nuc, I_r]
            # we want to make sure no d-1 columns of H_nuc are dependent
            A_ind = np.random.randint(0, 4, size=(r, k)) # randomly generate A
            H_ind = np.concatenate((A_ind, np.identity(r, dtype='int')), axis=1)
            H_nuc = ind_to_nuc(H_ind)
            
            # look through all subset of cardinality d-1 of {0, ..., n-1}, and 
            # check submatrix has full rank
            all_full_rank = True
            for S in itertools.combinations(range(n), d-1):
                if RREF(H_nuc[:, S])[1] != d - 1:
                    all_full_rank = False
                    break
            if all_full_rank:
                # G_nuc is a k x n generator matrix [I_k, A_nuc^T]
                G_ind = np.concatenate((np.identity(k, dtype='int'), A_ind.T), axis=1)
                G_nuc = ind_to_nuc(G_ind)
                return (G_nuc, H_nuc, k, r)
        r = r + 1  # increase redundancy to make linear independence more likely
        print('Max_trials reached... creating generator matrix for k =', n - r)
        
#%%
"""
Generate a 4-ary linear code of word length n and pairwise hamming distance >=d

We first generate a k x n generator matrix G_nuc, with k as large as possible.
We then generate barcodes as all possible linear combinations of rows of the 
generator matrix (4^k in total). This method can achieve the strengthed 
Gilbert-Varshamov lower bound. 
"""
def linear_code(n, d, outFile = 'codes.txt'):
    # The distance d of a linear code C also equals the minimum number of 
    # linearly dependent columns of the check matrix H. 
    G_nuc, _,  k, _ = generator(n, d)
    print('k =', k)
    
    # generate all codes by taking all possible linear combinations (over GF(4))
    # of rows in G_nuc. 
    all_nuc = ind_to_nuc([0,1,2,3])
    G_nuc_transpose = G_nuc.T
    textFile = open(outFile, 'w')
    counter = 0
    for w in itertools.product(all_nuc, repeat=k):
        counter += 1
        if counter % 10000 == 0: print(counter, '/', 4**k) # print progress
        code = G_nuc_transpose.dot(np.array(w))
        s = ''.join([letter.__str__() for letter in code])
        textFile.write('%s\n' % s)
        
    textFile.close()

#%%
"""
main function for running from command line
usage example: hamming.py -n 12 -d 3 -o 'codes.txt'
"""
def main(argv):
    # default outFile name
    outFile = 'codes.txt'
    
    # get commandline input arguments
    try:
        opts, args = getopt.getopt(argv,'hn:d:o:',['outfile='])
    except getopt.GetoptError:
        print('hamming.py -n <code_length> -d <hamming_dist> -o <out_file>')
        sys.exit(2)
    
    # parse input arguments
    for opt, arg in opts:
        if opt == '-h':
            print('hamming.py -n <code_length> -d <hamming_dist> -o <out_file>')
            sys.exit()
        elif opt in ('-n'):
            n = int(arg)
        elif opt in ('-d'):
            d = int(arg)
        elif opt in ('-o', '--outfile'):
            outFile = str(arg)
    
    # run linear_code function
    linear_code(n, d, outFile)
    
# run main function if directly executed from commandline
if __name__ == '__main__':
   main(sys.argv[1:])


