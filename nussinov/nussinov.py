#! usr/bin/python

#=========================================================
# Michael Ting
# nussinov.py
# Nussinov Algorithm for RNA Secondary Structure Prediction
# 10/26/12
#
# This program reads command-line arguments of the form
# $ nussinov.py AAGGUCUCAGA
# 
# Output of the program includes the input sequence length,
# maximum number of Watson-Crick base pairs, and the total
# calculation time.
#=========================================================

import sys
import re
import time

# scoring dictionary for WC base pairs
scoredict = {   'AU':1,
                'UA':1,
                'CG':1,
                'GC':1,
                'AA':0,
                'AC':0,
		'AG':0,
                'CA':0,
		'CC':0,
		'CU':0,
		'GA':0,
		'GG':0,
		'GU':0,
		'UC':0,
		'UG':0,
		'UU':0 }

def nussinov_fill(seq):

    maxpairs = 0
    seqlen = len(seq)
	
    # initialize the 2D square array (table), with axes corresponding to the sequence
    matrix = []
    for index in range(seqlen):
	matrix.append([0])                            
    for index in range(seqlen):
	for jndex in range(seqlen):
	    matrix[index].append(0)
	
    # bases cannot pair with themselves, set score to 0
    for i in range(1,seqlen):            
	matrix[i].insert(i,0)
	matrix[i].insert(i-1,0)      
	
    # calculate the maximum number of base pairs in the input sequence
    for h in range(1,seqlen):           # travel down the diagonal of the table
	i = 0   
	for j in range(h,seqlen):       # travel across each element of diagonal in table
	    matrix[i].insert(j,matrix[i+1][j])    # base i is unpaired
	    maxpairs = matrix[i][j-1]        # base j is unpaired
	    if maxpairs > matrix[i][j]:
		matrix[i][j] = maxpairs
			
	    pair = seq[i] + seq[j]   			    # select bases for WC-pair check
	    maxpairs = matrix[i+1][j-1] + scoredict.get(pair)    # bases i and j are WC pairs
	    if maxpairs > matrix[i][j]:                	    # return score from dictionary
		matrix[i].insert(j,maxpairs)
			
	    # bifurcation point in the RNA structure
	    for k in range(i+1,j):             
		maxpairs = matrix[i][k] + matrix[k+1][j]
		if maxpairs > matrix[i][j]: 
		    matrix[i].insert(j,maxpairs)
			
	    # increment for the j loop
	    i+= 1

    # maximum number of WC base pairs
    return matrix[0][seqlen-1]

def main():

    # read command-line arguments of the form $ nussinov.py AAGGUCUCAGA
    if (len(sys.argv) == 1):
        raise TypeError("Please enter a sequence.")
    elif (len(sys.argv) == 2):
        sequence = sys.argv[1]
    else:
        raise TypeError("Too many arguments!")

    # check that sequence only has A,C,G,U
    sequence = sequence.strip().upper()
    if re.search(r"[^ACGU]", sequence):
        raise TypeError("Incorrect sequence format!")

    # Calculate maximum number of base pairs with Nussinov
    starttime = time.time()    
    maxbp = nussinov_fill(sequence)
    endtime = time.time()-starttime
    print "Input sequence length: %d" % len(sequence)
    print "Maximum number of WC base pairs: %d" % maxbp
    print "Total calculation time: %f seconds" % endtime
    
if __name__ == "__main__":
    main()
