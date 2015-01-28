#!/usr/bin/python
import MySQLdb
import MySQLdb.cursors
from graph_tool.all import *

# Peak Class
class Peak(object):
	def __init__(self, id, chromosome, start, stop, 
		summit, sequence, ebna3a, ebna2, ebna3b, ebna3c,
		rbpj92, rbpj234):
		self.id         = id
		self.chromosome = chromosome
		self.start      = start
		self.stop       = stop
		self.summit     = summit
		self.sequence   = sequence
		self.ebna3a     = ebna3a
		self.ebna2 	    = ebna2
		self.ebna3b     = ebna3b
		self.ebna3c     = ebna3c
		self.rbpj92     = rbpj92
		self.rbpj234    = rbpj234


db = MySQLdb.connect(
	host="localhost",
	user="peaker",
	passwd="11peaks12",
	db="Peaks",
	cursorclass=MySQLdb.cursors.DictCursor
)


# Set up db handle
dbh = db.cursor()

# First get the list of chromosomes
dbh.execute("SELECT DISTINCT chromosome FROM peaks")

# chromosome list
chromosomes = []

for row in dbh:
	chromosomes.append(row['chromosome'])

for chromosome in chromosomes:
	peaks = []

	dbh.execute("SELECT * FROM peaks WHERE chromosome=%s", chromosome)

	for row in dbh:
		peak = Peak (
			id         = row['id'],
			chromosome = row['chromosome'],
			start      = row['start'],
			stop       = row['stop'],
			summit     = row['summit'],
			sequence   = row['sequence'],
			ebna3a     = row['ebna3a'],
			ebna2 	   = row['ebna2'],
			ebna3b     = row['ebna3b'],
			ebna3c     = row['ebna3c'],
			rbpj92     = row['rbpj92'],
			rbpj234    = row['rbpj234']
		)
		peaks.append(peak)

	graph = Graph(directed=False)
	for peak in peaks:
		graph.add_vertex(peak)