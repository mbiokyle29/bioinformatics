# Eland-fastq.py
# Eric Johannsen    Feb 2011
#
#
# Take and eland output file (~~export.txt), parse each line, and write all the nonmatched reads (NM), that passed purity (line[21] == y) to an output file, counting as we go.

import sys
import string


if (len(sys.argv) < 1):
    print "Usage:   python", sys.argv[0], "<Eland-output> "
else:
    print "\n\nEland output filename = ", sys.argv[1]
    fp = open(sys.argv[1])
    op = open('output', 'w')

x = -1
y = 0

text = fp.readline().replace("\n","")
print text
while text:
    x += 1
    line = text.split("\t")
    print >>op, "@", line[0]
    print >>op, line[8]
    print >>op, "+"
    print >>op, line [9]
    
    if (x % 1000000 == 1):
        print '{0:3d}'.format((x-1)/1000000),"million reads processed"

    text = fp.readline().replace("\n","")

total = (x+1)
print total, "reads counted\n"
exit()



