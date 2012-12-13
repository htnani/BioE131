#! usr/bin/python

#=======================================================================#
# Michael Ting                                                          #
# upgmay                                                                #
# UPGMA algorithm implementation                                        #
#                                                                       #
# This program outputs a phylogenetic tree in Newick format using the   #
# UPGMA algorithm, given a distant matrix .distmat file.                #
# The phylogenetic tree is also written to an output file tree.txt.     #
#=======================================================================#

import sys
import re
        
# Reads a distmat file and extracts the names and distances into
# separate lists. The distances are placed into a 2D matrix where [0][1]
# corresponds to the distance between sequences in the names list,
# names[0] and names[1].

def readDMatFile(dmatfile):
    names = []
    distmat = []
    for line in dmatfile:
        items = line.split()
        names.append(items.pop(0))  # extract names
        distmat.append(items)       # generate distance-only matrix
    nodenum = len(names)
    for i in range(len(distmat)):
        for j in range(len(distmat[i])):
            distmat[i][j] = float(distmat[i][j])
    return [names, distmat, nodenum]

# Finds the minimum distance in a symmetric distance matrix
# Also returns the indices of the minimum distance
# Input:
#   distmat - the distance matrix, with names already extracted
# Output:
#   mindist - the minimum distance value in the matrix
#   xpos - the row index of mindist
#   ypos - the column index of mindist

def findMinDist(distmat):
    mindist = float("inf")
    xpos = 0
    ypos = 0
    for i in range(len(distmat)):
        for j in range(len(distmat[i])):
            if (i < j) and (distmat[i][j] < mindist):
                mindist = distmat[i][j]
                xpos = i
                ypos = j
    return [xpos, ypos, mindist]

def testFMD():
    testmat = [ [0, 0.9, 0.4, 0.3],
                [0.9, 0, 0.2, 0.5],
                [0.4, 0.2, 0, 0.7],
                [0.3, 0.5, 0.7, 0] ]
    testresult = findMinDist(testmat)
    assert testresult[0] == 1, "xpos failed!"
    assert testresult[1] == 2, "ypos failed!"
    assert testresult[2] == 0.2, "mindist failed!"

    print "all tests passed!"


# Returns an ultrametric tree in Newick format using the UPGMA algorithm.
# Input:
#   names - a list of names corresponding to all genes in the matrix
#   distances - a 2D distance matrix corresponding to distances between genes
# Output:
#   newick - The ultrametric tree in Newick format.

def UPGMA(names, distances):
    nodes = names
    distmat = distances
    heightlist = []
    # initialize heights to 0 for all nodes
    for i in range(len(nodes)):
        heightlist.append(0)
    while len(nodes) >= 2:
        # find minimum distance and derive the nodes that the distance corresponds to
        minlist = findMinDist(distmat)
        xpos    = minlist[0]
        ypos    = minlist[1]
        mindist = minlist[2]
            
        # create the parent node K, calculate heights of parent and branches
        iheight = heightlist[xpos]
        jheight = heightlist[ypos]
        parheight = 0.5*( iheight + jheight + mindist )
        kibranch = parheight - iheight
        kjbranch = parheight - jheight
        iname = nodes[xpos]
        jname = nodes[ypos]
        newname = "(" + iname + ":" + str(kibranch) + "," + jname + ":" + str(kjbranch) + ")"
        
        # Define distances from K to all other nodes:
        # D(K,L) = ( D(I,L) + D(J,L) ) * 0.5 for all nodes L in N
        kdlist = []
        for loc in range(len(distmat)):
            distil = distmat[xpos][loc]
            distjl = distmat[ypos][loc]
            newdist = 0.5*(distil + distjl)
            kdlist.append(newdist)
            distmat[loc].append(newdist)    # put k,l distance into last column
        kdlist.append(0.0)  # Needed for diagonal
        distmat.append(kdlist)
        
        # insert K into N
        nodes.append(newname)
        heightlist.append(parheight)
 
        # remove the closest two nodes from N
        # calling pop() shifts indices over by 1, so need to account for
        # multiple calls of pop() by decrementing
        nodes.pop(xpos)
        nodes.pop(ypos-1)
        for distlist in distmat:
            distlist.pop(xpos)
            distlist.pop(ypos-1)
        distmat.pop(xpos)
        distmat.pop(ypos-1)
        heightlist.pop(xpos)
        heightlist.pop(ypos-1)

    newick = nodes[0]   # Extract string from list
    return newick

# Runs all test functions
def testall():
    testFMD()

# Reads from the command line with the format:
# python hw9.py DISTMATFILE

def main():

    if len(sys.argv) == 1:
        raise IOError("Please enter a distance matrix as an argument")
    elif len(sys.argv) == 2:
        dmatfile = sys.argv[1]
    else:
        raise IOError("Too many arguments!")

    try:
        dmatfile = open(dmatfile)
    except IOError:
        print "Invalid file name!"
        sys.exit()

    fileinfo = readDMatFile(dmatfile)

    dmatfile.close()

    NAMEINDEX = 0
    DISTINDEX = 1
    NODEINDEX = 2

    names   = fileinfo[NAMEINDEX]
    distmat = fileinfo[DISTINDEX]
    nodes   = fileinfo[NODEINDEX]

    newick = UPGMA(names, distmat) + ";"
    newick = re.sub(r'/\d+-\d+:',':',newick)

    try:
        newfile = open("tree.txt", "w")
        try:
            newfile.write(newick)
        finally:
            newfile.close()
    except IOError:
        print "Error with writing to output file!"
        sys.exit()

    print "Newick format output written to tree.txt"
    print ""
    print newick

if __name__ == "__main__":
    main()
