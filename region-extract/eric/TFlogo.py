# MapSeq.py
# Eric Johannsen    Apr 2011
# 
#
# Given a PSSM for a pattern, will output a flat file that can be pasted into WebLogo to make a log
# Only works if each position in the matrix adds up to the same number (i.e., A + C + G + T = constant)



import sys
import string
import math


def order(profile):
    ord = 0.0
    h = 0.0
    for c in xrange(len(profile[0])):
        colsum = 0
        for r in xrange(4):
            colsum += profile[r][c]
        for r in xrange(4):
            if (profile[r][c] != 0):
                odds = float(profile[r][c]) / float(colsum)
                value = odds * math.log(odds,2)
                h += value
            else:
                value = 0.0
            ord += (0.5 + value)
    h *= -1
    return ord, h


def entropy(profile):
    h = 0.0
    for r in xrange(len(profile)):
        for c in xrange(len(profile[0])):
            if (profile[r][c] != 0):
                h += profile[r][c] * math.log(profile[r][c],2)
    h = h * -1.0
    return h

##########
# main():
##########

if (len(sys.argv) < 2):
    print "Usage:   python", sys.argv[0], "<matrix-filename>"
    exit()
else:
    print "\n\nmatrix filename = ", sys.argv[1]
    mp = open(sys.argv[1])

###############
# read inputs #
###############

matrix = [[] for row in xrange(4)]


tfname = mp.readline().replace("\n", "").replace(">", "")
text = mp.readline().replace("\n","")
while text:
    line = text.split(" ")
    for row in xrange(4):
        matrix[row].append(int(line[row+1]))
    text = mp.readline().replace("\n","")

patternlen = len(matrix[0])
mh = order(matrix)
print "Matrix order = ", '%.1f' % (mh[0])
print "Matrix entropy = ", '%.1f' % (mh[1])


################
# open outputs #
################

outname = tfname + ".logo"
op = open(outname, 'w')

seqnum = 0
for r in xrange(4):
    seqnum += matrix[r][0]

for seqs in xrange(seqnum):
    text = ""
    for c in xrange(patternlen):
        for r in xrange(4):
            if matrix[r][c] > 0:
                text += str(r)
                matrix[r][c] -= 1
                break

    decode = text.translate(string.maketrans("0123", "ACGT"))
    print >>op, decode


exit()
