# Sam2Wig-BFR.py
# Eric Johannsen    Jan 2013
#
#
# Take a bwa .sam output file (sorted with unmapped reads removed), parse each line, and write it as new file in wiggle format.
#
# The overall strategy here is to step through the reads mapped by bwa (which are sorted in order along the chromosome) and store the "depth"
# of each read in a list (Fdepth for forward reads and Rdepth for reverse reads).  These lists serve as a window that we slide along the
# chromosome and at each position we add the number of forward reads at position to every member of Fdepth and the number of reverse reads to Rdepth
# finally the end of the list is popped off and written to the output file.  This effectively "stacks up" the reads on each other and measures
# the depth at every position on the chromosome.
# I included the sfor and srev flags in avoid writting to the file until at least one read had been encountered (this avoids a lot of leading
# zeros in sparse data.


import sys
import re
import string

##########
# main():
##########

if (len(sys.argv) < 3):
    print "Usage: ", sys.argv[0], "<reads-filename> <genome-fastafile>"
    exit()
else:
    print "\nMaq mapview outfile name = ", sys.argv[1]
    fp = open(sys.argv[1])
    gp = open(sys.argv[2])
    outname = sys.argv[1]
    outname = re.sub(".sam", "", outname)
    outf = outname + "-F.wig"
    outr = outname + "-R.wig"
    outb = outname + "-B.wig"
    opf = open(outf, 'w')
    opr = open(outr, 'w')
    opb = open(outb, 'w')

genomename = gp.readline().replace("\n", "").replace(">", "")
genseq = gp.readline().replace("\n", "").upper()
genomelen = len(genseq)


z = 0
text = fp.readline().replace("\n","")
line = text.split("\t")
if (line[2] != genomename):
    print "mismatch between reference genome =", genomename, "\nand mapfile genome =", line[1]
    exit()

readlen = len(line[9])
print "Using read length of", readlen

Fdepth = [0]*readlen
Rdepth = [0]*readlen
Bdepth = [0]*readlen
index = 1
nfor = 0    # count number of forward reads at current position
nrev = 0
tfor = 0    # count total forward reads
trev = 0
tboth = 0
sfor = 0    # flag to keep track of whether any forward reads have been encountered yet
srev = 0
sboth = 0

# Loop through the .sam output file counting the number of forward (nfor) and reverse (nrev) at each position on the chromosome (index)
# once we encounter a read that is beyond the current index (start > index) we add our accumulated reads to the Fdepth and Rdepth lists
# and write the read depth up to the new read position (i.e, index = start) and begin reading sequences from the file anew.
# note that the second column of a sam file contains a bitwise flag that includes the strand (16 position).
# this can be extracted by line[2]&16  which will evaluate to 0 if set (sense) or 16 if set (reverse strand)

while text:
    z += 1
    if (z % 100000 == 0):
        print '{0:3d}'.format(z/1000000),'{0:s}'.format("."), '{0:1d}'.format((z/100000)%10), "million reads processed"
    line = text.split("\t")
    start = int(line[3])
    alnflag = int(line[1])
    if (start < index):               
        print "Error:  input samfile is not sorted!"
        exit()
    elif (start > genomelen):
        print "Warning: read mapped beyond end of genome -> discarded"
        continue
    elif (start > index):
        for x in xrange(readlen):
            Fdepth[x] += nfor
            Rdepth[x] += nrev
            Bdepth[x] += (nfor+nrev)
        tfor += nfor
        trev += nrev
        tboth += (nfor+nrev)
        if (tfor > 0):
            if (sfor == 0):
                head = 'track type=wiggle_0 name="' + outname + '-F" description="' + outname + '-F" visibility=2 color="0,0,255"'
                print >> opf, head
                head = 'fixedStep chrom=' + genomename + ' start=' + str(index + 1) + ' step=1'
                print >> opf, head
                sfor = 1
            for y in xrange (index, start):
                print >> opf, Fdepth.pop(0)
                Fdepth.append(0)
        if (trev > 0):
            if (srev == 0):
                head = 'track type=wiggle_0 name="' + outname + '-R" description="' + outname + '-R" visibility=2 color="255,0,0"'
                print >> opr, head
                head = 'fixedStep chrom=' + genomename + ' start=' + str(index + 1) + ' step=1'
                print >> opr, head
                srev = 1
            for y in xrange (index, start):
                print >> opr, Rdepth.pop(0)
                Rdepth.append(0)
        if (tboth > 0):
            if (sboth == 0):
                head = 'track type=wiggle_0 name="' + outname + '-B" description="' + outname + '-B" visibility=2 color="0,0,0"'
                print >> opb, head
                head = 'fixedStep chrom=' + genomename + ' start=' + str(index + 1) + ' step=1'
                print >> opb, head
                sboth = 1
            for y in xrange (index, start):
                print >> opb, Bdepth.pop(0)
                Bdepth.append(0)
        index = start
        if (alnflag&16):    # (16 is set if reverse strand)
            nfor = 0
            nrev = 1
        else:
            nfor = 1
            nrev = 0
    elif (alnflag&16 == 0):  #  if we get to this statment then start == index and we need to increment the read depth (either forward or reverse)
        nfor += 1           #    for this index and look for more reads at this position before writing
    elif (alnflag&16):
        nrev += 1
            
    text = fp.readline().replace("\n","")

tfor += nfor
trev += nrev
tboth += (nfor+nrev)
for x in xrange(readlen):
    Fdepth[x] += nfor
    Rdepth[x] += nrev
    Bdepth[x] += (nfor+nrev)
for y in xrange(readlen):
    print >> opf, Fdepth.pop(0)
    print >> opr, Rdepth.pop(0)
    print >> opb, Bdepth.pop(0)

print "Total reads = ", z
print "Forward reads = ", tfor
print "Reverse reads = ", trev
print "F+R reads = ", tfor + trev, tboth

exit()



