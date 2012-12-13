#! /usr/bin/python

#=================================================================
# Michael Ting
# simulator.py
# 9/28/12
#
# Generates a random nucleotide sequence
# - no additional arguments: random 1000bp dna sequence
# - load length and composition statistics from a given FASTA file
# - load statistics from a given params file
# - save statistics from the generated sequence into a new file
#=================================================================

import sys
import optparse
import random

# Constants
SEQ_COMP = 0	# index of comp is 0
SEQ_LEN = 1		# index of length is 1

# default stats
DEFAULT_LEN = 1000
DEFAULT_COMP = {
	'a':0.25,
	't':0.25,
	'g':0.25,
	'c':0.25
	}
	
# Processes a FASTA file
# Input:
#	infile - FASTA file of sequence(s)
# Output:
#	List with average composition (0) and average length (1) of
#	sequences from infile
def process_file(infile):
	# list of (comp, len)
	seq_list = []
	num_seqs = 0
	comb_len = 0
	avg_comp = {
		'a':0.0,
		't':0.0,
		'g':0.0,
		'c':0.0
		}
	avg_len = 0.0
	
	# store composition and length of each seq in file
	for name, seq in process_seq(infile):
		seq_list.append([get_comp(seq), get_length(linearize(seq))])

	# total number of sequences used to calculate averages
	num_seqs = len(seq_list)	
	
	# store total combined lengths and compositions
	for item in seq_list:
		comb_len += item[SEQ_LEN]
		for base in item[SEQ_COMP]:
			avg_comp[base] = avg_comp.get(base) + item[SEQ_COMP].get(base)
	
	# compute averages
	avg_len = comb_len*1.0/num_seqs
	for base in avg_comp:
		avg_comp[base] = avg_comp.get(base)*1.0/num_seqs
	
	return [avg_comp, avg_len]	

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

# Loads the length and composition from a stats file
# Input:
#	infile - formatted file with length and composition stats
# Output:
#	List of composition (0) and length (1)
def load_stats(infile):
	seq_len = 0
	seq_comp = {}
	for line in infile:
		line = line.strip()
		if line.startswith("len:"):
			seq_len = float(line[4:].strip())
		elif line.startswith("%c:" % line[:1]):
			seq_comp[line[:1]] = float(line[2:].strip())
		else:
			raise TypeError("Incorrect file format!")
	return [seq_comp, seq_len]
		
# Calculates nucleotide composition normalized to the total number of bases
# Input:
#	seq - a linear DNA sequence
# Output:
#	Dictionary of bases with values corresponding to base composition
def get_comp(seq):
	seq = seq.strip().lower()
	base_count = {	
		'a':0,
		't':0,
		'g':0,
		'c':0
		}
	for base in seq:
		base_count[base] = base_count.get(base) + 1
	total_bases = 0
	for base in base_count:
		total_bases += base_count.get(base)
	for base in base_count:
		base_count[base] = base_count.get(base)*1.0/total_bases
	return base_count

# Test case for get_comp()
def test_get_comp():
	seq1 = "agactctcca"
	comp1 = {
		'a':0.3,
		't':0.2,
		'g':0.1,
		'c':0.4
		}
	test1 = get_comp(seq1)
	assert test1 == comp1, "get_comp test 1 failed"
	print "test_get_comp() passed"
	
# Calculates the length of an input DNA sequence
# Input:
#	seq - linearized sequence of DNA
# Output:
#	int - length of DNA sequence
def get_length(seq):
	seq = seq.strip()
	return len(seq)

# Test case for get_length()
def test_get_length():
	seq1 = "AGAGAGAG"
	len1 = 8
	test1 = get_length(seq1)
	assert test1 == len1, "get_length test 1 failed"
	
	seq2 = "   atcgagagt"
	len2 = 9
	test2 = get_length(seq2)
	assert test2 == len2, "get_length test 2 failed"
	
	print "test_get_length() passed"		

# Produces the next random base dependent on provided bounds
# Input:
#	atb - A-T bound
#	tgb - T-G bound
#	gcb - G-C bound
# Output:
#	A nucleotide base dependent on random number generated
def next_base(atb, tgb, gcb):
	val = random.random()	# a random number between 0.0 and 1.0
	if val < atb:
		return 'a'
	elif (val > atb) and (val < tgb):
		return 't'
	elif (val > tgb) and (val < gcb):
		return 'g'
	else:
		return 'c'

# Replaces stop codons in the sequence
# Input:
#	seq - randomly generated nucleotide sequence
# Output:
#	Nucleotide sequence with stop codons replaced with a non-stop
#	codon with equivalent bases
def replace_stop(seq):
	for pos in range(len(seq)):
		if (pos == 0) or ((pos % 3) == 0):
			codon = seq[pos:pos+3]
			if codon == "taa":
				seq = seq[:pos] + "ata" + seq[pos+3:]
			elif codon == "tag":
				seq = seq[:pos] + "gta" + seq[pos+3:]
			elif codon == "tga":
				seq = seq[:pos] + "gat" + seq[pos+3:]
	return seq
	
# Test case for replace_stop()
def test_replace_stop():
	seq1 = "atgcgc"
	fix1 = "atgcgc"
	test1 = replace_stop(seq1)
	assert test1 == fix1, "replace_stop test 1 failed"
	
	seq2 = "atgtaacgc"
	fix2 = "atgatacgc"
	test2 = replace_stop(seq2)
	assert test2 == fix2, "replace_stop test 2 failed"
	
	seq3 = "atgtaatgactagcgc"
	fix3 = "atgatagatctagcgc"
	test3 = replace_stop(seq3)
	assert test3 == fix3, "replace_stop test 3 failed"
	
	print "test_replace_stop() passed"

# Concatenates the nucleotide sequence into one linear sequence.
# Input: 
#	seqs - A List of String sequence lines
# Output: A single linear String of the nucleotide sequence	
def linearize(seqs):
	linear_strand = ""
	for strand in seqs: 
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

# Runs all test functions
def test():
	test_get_comp()
	test_get_length()
	test_replace_stop()
	test_linearize()
	test_split80()
	
def main():
	
	myparser = optparse.OptionParser()
	myparser.add_option("-c", "--calc",
						type="string",
						dest="seqfile", 
						help="use len and comp of specified file in seq gen")
	myparser.add_option("-l", "--load", 
						type="string", 
						dest="loadfile", 
						help="loads params from specified file")
	myparser.add_option("-s", "--save", 
						type="string", 
						dest="savefile", 
						help="saves params into new file")

	(options, args) = myparser.parse_args()

	# set stats to default
	length = DEFAULT_LEN
	composition = DEFAULT_COMP

	# ensures that user does not try to use two separate composition files
	if (options.seqfile != None) and (options.loadfile != None):
		raise TypeError("Cannot use --calc and --load simultaneously")
	
	# calculate and set new composition and length from file
	if options.seqfile != None:
		sefile = open("%s" % options.seqfile)
		averaged = process_file(sefile)
		length = averaged[SEQ_LEN]
		composition = averaged[SEQ_COMP]
		sefile.close()
		
	# load length and composition statistics from file	
	if options.loadfile != None:
		lofile = open("%s" % options.loadfile)
		stat_list = load_stats(lofile)
		length = stat_list[SEQ_LEN]
		composition = stat_list[SEQ_COMP]
		lofile.close()
		
	# save length and composition statistics in new file	
	if options.savefile != None:
		outfile = open("%s" % options.savefile, "w")
		outfile.write("len: %d\n" % length)
		for base in composition:
			outfile.write("%c: %f\n" % (base, composition[base]))
		outfile.close()
	
	# set bounds for random bases: a,t,g,c
	at_bound = composition['a']
	tg_bound = at_bound + composition['t']
	gc_bound = tg_bound + composition['g']
	
	final_seq = "atg"
	final_len = 3
	
	# randomly generate sequence based on composition and length stats
	while final_len < length:
		final_seq += next_base(at_bound, tg_bound, gc_bound)
		final_len += 1
	
	# check for stop codons and replace
	final_seq = replace_stop(final_seq)
	
	# print FASTA format
	final_seq = split80(final_seq)
	print_FASTA("> random sequence", final_seq)

if __name__ == "__main__":
	main()