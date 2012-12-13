#! /usr/bin/python

#=================================================
# Michael Ting
# revcomp.py
# 9/21/12
#
# Outputs the reverse complement of all nucleotide
# sequences in an input FASTA file.
#=================================================

import sys

# Processes a FASTA file by outputting the reverse complement with
#	max 80bp per line format
# Input:
#	infile - FASTA file of sequence(s)
#	nuc - "dna" for DNA, "rna" for RNA
# Output:
#	FASTA file of reverse complemented sequences, with length and 80bp per line	
def process_file(infile, nuc):
	for name, seq in process_seq(infile):
		final_seq = format_nuc(seq, nuc)
		print_FASTA(name, final_seq)	

# Generator that finds all sequences in a FASTA file
# Input:
#	infile - FASTA file of sequence(s)
# Output:
#	name - the name of the sequence following the ">"
#	seq - the nucleotide sequence
def process_seq(infile):
    name, seq = None, []
    for line in infile:
        if line.startswith(">"):
            if name: 
				yield (name, ''.join(seq))
            name, seq = line.strip(), []
        else:
            seq.append(line.strip())
    if name: 
		yield (name, ''.join(seq))		
	
# Prints a formatted sequence in FASTA format
# Input:
#	name - comment with sequence name
#	seq - A List of the nucleotide sequence
# Output:
#	FASTA format of the sequence including length
def print_FASTA(name, seq):
	seq_length = len(linearize(seq))
	full_name = name + ", " + "%d " % seq_length + "bp"
	print full_name
	for row in seq:
		print row
	
# Takes in a list of nucleotides, performs the reverse complement,
#	outputs in acceptable format (80 bases per line)
# Input:
#	seqs - A list of Strings corresponding to nucleotide sequences
#	nuc - "dna" for DNA, "rna" for RNA, else rev_comp() raises TypeError
# Output:
#	A list of the reverse complemented sequence formatted 80 bases per line
def format_nuc(seqs, nuc):
	reversed = rev_no_mod(seqs)			# reverse line order
	rev_comp_seq = []
	for phrase in reversed:				# reverse complement the sequence
		rev_comp_seq.append(rev_comp(phrase, nuc))
	single_seq = linearize(rev_comp_seq)
	formatted_seq = split80(single_seq)
	return formatted_seq
	
# Test case for format_nuc()
def test_format_nuc():
	seq1 = ["ATCG", "GGTA"]
	final_seq1_dna = ["TACCCGAT"]
	final_seq1_rna = ["UACCCGAU"]
	
	dna_format_1 = format_nuc(seq1, "dna")
	rna_format_1 = format_nuc(seq1, "rna")
	assert dna_format_1 == final_seq1_dna, "dna test 1 failed"
	assert rna_format_1 == final_seq1_rna, "rna test 1 failed"
	
	seq2 = ["AAAAAGGGGGAAAAAGGGGGAAAAAGGGGGAAAAAGGGGGAAAAAGGGGGAAAAAGGGGGAAAAAGGGGGAAAAAGGGGGAAAAAGGGGG"]
	final_seq2_dna = ["CCCCCTTTTTCCCCCTTTTTCCCCCTTTTTCCCCCTTTTTCCCCCTTTTTCCCCCTTTTTCCCCCTTTTTCCCCCTTTTT", "CCCCCTTTTT"]
	final_seq2_rna = ["CCCCCUUUUUCCCCCUUUUUCCCCCUUUUUCCCCCUUUUUCCCCCUUUUUCCCCCUUUUUCCCCCUUUUUCCCCCUUUUU", "CCCCCUUUUU"]
	
	dna_format_2 = format_nuc(seq2, "dna")
	rna_format_2 = format_nuc(seq2, "rna")
	assert dna_format_2 == final_seq2_dna, "dna test 2 failed"
	assert rna_format_2 == final_seq2_rna, "rna test 2 failed"
	print "test_format_nuc() passed"
	
# Concatenates the nucleotide sequence into one linear sequence.
# Input: 
#	seqs - A List of String sequence lines
# Output: A single linear String of the nucleotide sequence	
def linearize(seqs):
	linear_strand = ""
	for strand in seqs: #range(len(seqs))
		linear_strand = linear_strand + strand
	return linear_strand

# Test case for linearize()
def test_linearize():
	a = ["abcd", "efgh"]
	assert linearize(a) == "abcdefgh", "linearize test a failed"
	a.append("ijkl")
	assert linearize(a) == "abcdefghijkl", "linearize test b failed"
	print "test_linearize() passed"

# Splits a linear sequence into a list of sequences no more than 80
# 	characters long. Intended for use in combination with linearize().
# Input: 
#	seq - A single linear String of the nucleotide sequence
# Output: A List of Strings, each String no more than 80 characters long
def split80(seq):
	remaining = seq
	list_of_seq = []
	while (len(remaining) >= 80):
		list_of_seq.append(remaining[:80])
		remaining = remaining[80:]
	list_of_seq.append(remaining)
	return list_of_seq
	
# Test case for split80()
def test_split80():
	x = "ABCDEFGHIJABCDEFGHIJABCDEFGHIJABCDEFGHIJABCDEFGHIJABCDEFGHIJABCDEFGHIJABCDEFGHIJABCDEFGHIJ"
	result = 0
	split1 = "ABCDEFGHIJABCDEFGHIJABCDEFGHIJABCDEFGHIJABCDEFGHIJABCDEFGHIJABCDEFGHIJABCDEFGHIJ"
	split2 = "ABCDEFGHIJ"
	
	result = split80(x)
	assert result != 0, "split80 not storing in result"
	assert result[0] == split1, "split80 test a failed"
	assert result[1] == split2, "split80 test b failed"
	
	y = "ABCDEFGHIJABCDEFGHIJABCDEFGHIJABCDEFGHIJABCDEFGHIJABCDEFGHIJABCDEFGHIJABCDEFGHIJABCDEFGHIJABCDEFGHIJABCDEFGHIJABCDEFGHIJABCDEFGHIJABCDEFGHIJABCDEFGHIJABCDEFGHIJABCDEFGHIJ"
	result2 = 0
	split3 = "ABCDEFGHIJABCDEFGHIJABCDEFGHIJABCDEFGHIJABCDEFGHIJABCDEFGHIJABCDEFGHIJABCDEFGHIJ"
	split4 = "ABCDEFGHIJABCDEFGHIJABCDEFGHIJABCDEFGHIJABCDEFGHIJABCDEFGHIJABCDEFGHIJABCDEFGHIJ"
	split5 = "ABCDEFGHIJ"
	
	result2 = split80(y)
	assert result2 != 0, "split80 not storing in result2"
	assert result2[0] == split3, "split80 test c failed"
	assert result2[1] == split4, "split80 test d failed"
	assert result2[2] == split5, "split80 test e failed"
	print "test_split80() passed"
	
# Reverse the order of lines in a List of sequences
# Note: This differs from reverse() in that it DOES NOT modify
#		the original List
# Input: 
#	seqs - A List of sequences
# Output: A List of sequences in the reverse order	
def rev_no_mod(seqs):
	rev_lines = seqs[::-1]
	return rev_lines
	
# Test case for rev_seq_lines()
def test_rev_no_mod():
	s1 = "abcd"
	s2 = "efgh"
	s3 = "ijkl"
	s4 = "mnop"
	
	seq_list = [s1, s2, s3, s4]
	rev_seq_list = rev_no_mod(seq_list)
	
	for i in range(len(rev_seq_list)):
		assert rev_seq_list[i] == seq_list[len(rev_seq_list)-i-1], "test a %d failed" % i
	
	rev_seq_list[2] = "acac"
	assert rev_seq_list[2] != seq_list[len(rev_seq_list)-2-1], "test b failed"
	
	print "test_rev_no_mod() passed"

# Calculates the reverse complement of a single nucleotide sequence
# Input: 
#	seq - A String corresponding to a single nucleotide sequence
#	nuc - "dna" for DNA, "rna" for RNA, else raise TypeError
# Output: A String corresponding to the reverse complement of the input DNA sequence
def rev_comp(seq, nuc):
	
	# DNA complementary base pairs
	dna_comp_pairs = {
		'A':'T',
		'a':'t',
		'T':'A',
		't':'a',
		'G':'C',
		'g':'c',
		'C':'G',
		'c':'g'
		}
	# RNA complementary base pairs
	rna_comp_pairs = {
		'A':'U',
		'a':'u',
		'T':'A',
		't':'a',
		'G':'C',
		'g':'c',
		'C':'G',
		'c':'g'
		}
		
	# Choose type of complement
	base_dict = {}
	if nuc.lower() == "rna":
		base_dict = rna_comp_pairs
	elif nuc.lower() == "dna":
		base_dict = dna_comp_pairs
	else:
		raise TypeError("Invalid nucleotide option: %s" % nuc)
	
	# Reverse the sequence
	rev_seq = rev_no_mod(seq)
	
	# Substitute complementary bases
	rcomp = ""
	for base in rev_seq:
		rcomp = rcomp + base_dict.get(base)
	
	return rcomp

# Test case for rev_comp()
def test_rev_comp():
	seq1 = "ATCGGCTA"
	rev_seq1 = "TAGCCGAT"
	rna_rev_seq1 = "UAGCCGAU"
	assert rev_comp(seq1, "dna") == rev_seq1
	assert rev_comp(seq1, "rna") == rna_rev_seq1
	
	seq2 = "GGTCTA"
	rev_seq2 = "TAGACC"
	rna_rev_seq2 = "UAGACC"
	assert rev_comp(seq2, "dna") == rev_seq2
	assert rev_comp(seq2, "rna") == rna_rev_seq2
	
	print "test_rev_comp() passed"

# Runs all test functions
def test():
	test_format_nuc()
	test_linearize()
	test_split80()
	test_rev_comp()
	test_rev_no_mod()

	
# ====================================================
# Running the program by reading from the command line
#=====================================================

rna_check = ""
if (len(sys.argv) == 1):
	raise TypeError("Please enter a filename.")
elif (len(sys.argv) == 2):
	rna_check = "dna"
elif (len(sys.argv) == 3):
	rna_check = sys.argv[2]
else:
	raise TypeError("Too many arguments!")

# Check that the file can be opened
try:
	open(sys.argv[1])
except IOError:
	print "Could not open file!"

infile = open(sys.argv[1])
	
# Checker for processing lines
if rna_check == "rna":
	process_file(infile, "rna") # process as RNA
elif rna_check == "dna":
	process_file(infile, "dna") # process as DNA
else:
	raise TypeError("Incorrect option. For rna use \"rna\" and for dna omit argument.")

## Enable this for running tests
# test()	
	
infile.close()