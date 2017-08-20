# Kozak.py
# Eric Johannsen    Apr 2011
# 
#
# Given a PSSM for a pattern (log-odds assumed, but not absolutely required to work), will output the score for a sliding window across a sequence (e.g., genome)
# Calculates the score for each strand and writes the score for each to two separte files for each position in the genome (actually genome_length minus pattern_length + 1)  
# op = "output" file contains the output scores with a wiggle header for ease of uploading to UCSC browser


import sys
import string
import math


def normalize(matrix):
    nrow = len(matrix)
    ncol = len(matrix[0])
    pssm = [[0.0 for col in xrange(ncol)] for row in xrange(nrow)] 
    for c in xrange(ncol):
        colsum = 0.0
        for r in xrange(nrow):
            colsum += float(matrix[r][c])
        for r in xrange(nrow):
            if (matrix[r][c] == 0):
                pssm[r][c] = -100.0
            else:
                pssm[r][c] = math.log((float(matrix[r][c]) / colsum / 0.25),2)
    return pssm

##########
# main():
##########

if (len(sys.argv) < 3):
    print "Usage:   python", sys.argv[0], "<matrix-filename> <genome-filename>"
    exit()
else:
    print "\n\nmatrix filename = ", sys.argv[1], "\ngenome filename =", sys.argv[2]
    mp = open(sys.argv[1])
    gp = open(sys.argv[2])

###############
# read inputs #
###############

genomename = gp.readline().replace("\n", "").replace(">", "")
genseq = gp.readline().replace("\n", "").upper()
genlen = len(genseq)
codedseq = genseq.translate(string.maketrans("ACGT", "0123"))
matrix = [[] for row in xrange(4)]

tfname = mp.readline().replace("\n", "").replace(">", "")
text = mp.readline().replace("\n","")
while text:
    line = text.split(" ")
    for row in xrange(4):
        matrix[row].append(int(line[row+1]))
    text = mp.readline().replace("\n","")

patternlen = len(matrix[0])
pssm = normalize(matrix)
rcpssm = [[0.0 for col in xrange(patternlen)] for row in xrange(4)]
for c in xrange(patternlen):
    for r in xrange(4):
        rcpssm[r][c] = pssm[3-r][patternlen-1-c]
print pssm
print rcpssm

################
# open outputs #
################

outname = tfname + "-F.out"
op = open(outname, 'w')
outname = tfname + "-R.out"
op2 = open(outname, 'w')

#print header row
outname = tfname + "-F"
head = 'track type=wiggle_0 name="' + outname + '" description="' + tfname + '-predicted" visibility=2 colorByStrand="255,0,0 0,0,255"'
print >> op, head
outname = tfname + "-R"
head = 'track type=wiggle_0 name="' + outname + '" description="' + tfname + '-predicted" visibility=2 colorByStrand="255,0,0 0,0,255"'
print >> op2, head
head = 'fixedStep chrom=' + genomename + ' start=1 step=1'
print >> op, head
print >> op2, head


for pos in xrange(genlen - patternlen + 1):
    score = 0.0
    revscore = 0.0
    for c in xrange(patternlen):
        score += pssm[int(codedseq[pos+c])][c]
        revscore += rcpssm[int(codedseq[pos+c])][c]
    print >>op, '%.1f' % (2**score)
    print >>op2, '%.1f' % (2**revscore)
exit()

