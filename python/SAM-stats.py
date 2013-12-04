#!/usr/bin/python
# SAM-stats.py
# Eric Johannsen    Dec 2013
#
#
# Take a bowtie2 output file in SAM format, parse each line
# counting number reads that match uniquely by chromosome number.

import sys
import re
import string

if (len(sys.argv) < 2):
    print "Usage:   python", sys.argv[0], "<aligned_reads_filename> [<out_filename>]"
    exit()
else:
    print "\n\nBowtie2 output filename = ", sys.argv[1]
    fp = open(sys.argv[1])

    if (len(sys.argv) > 2):
       op = open(sys.argv[2], 'w')
    else:
       op = open('output', 'w')

symbolCounts = {}
key = []

x = -1

text = fp.readline().replace("\n","")
while text:
    x += 1
    line = text.split("\t")

#    print "This is line ", x+1

    if (line[2] in symbolCounts):
        symbolCounts[line[2]] = symbolCounts[line[2]] + 1
    else:
        symbolCounts[line[2]] = 1

    if (x % 1000000 == 1):
        print '{0:3d}'.format((x-1)/1000000),"million reads processed"

    text = fp.readline().replace("\n","")

print "\n\nsequence filename = ", sys.argv[1]
total = (x+1)

print total, "reads counted\n"


symbols = ['chrX', 'chrY', 'chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrM', 'EBV']


y = 0
for ch in symbols:
    if (ch in symbolCounts):
        print '{0:8s}'.format(ch),'{0:8d}'.format(symbolCounts[ch])
        y = y + symbolCounts[ch]
    else:
        print '{0:8s}'.format(ch),'{0:8d}'.format(0)

print "adds up to ", y, "\n\n\n"
print symbolCounts


print >>op, "sequence filename = ", sys.argv[1], "\n"
print >>op, "counted", total, "\n"

for ch in symbols:
    if (ch in symbolCounts):
        print >>op,'{0:8s}'.format(ch),'{0:8d}'.format(symbolCounts[ch])
    else:
        print >>op,'{0:8s}'.format(ch),'{0:8d}'.format(0)



exit()





