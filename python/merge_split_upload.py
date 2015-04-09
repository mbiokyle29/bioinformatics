#!/usr/bin/env python
"""
Kyle McChesney
Merge paired end read files with PEAR, then split them into files under 5GB and upload to one codex
"""
import argparse
import os
from math import ceil
from path import path
from subprocess import call

if __name__ == "__main__":
	parser = argparse.ArgumentParser(
		description = (" Upload to large fastq files to One Codex, merge them with PEAR and split with fastqsplitter")
	)
	parser.add_argument("--fastq_for", help="Forward - /1 fastq file ")
	parser.add_argument("--fastq_rev", help="Reverse - /2 fastq file ")
	parser.add_argument("--output",    help="Output file name and/or RELATIVE path")
	args = parser.parse_args()

	codex_max_gb = 5
	byte_to_gb = 1000000000
	forward = path(args.fastq_for).abspath()
	reverse = path(args.fastq_rev).abspath()
	output  = args.output

	# call PEAR on them
	# using default RAM Threads for now
	call(["pear", "-y", "2G", "-j", "5" ,"-f", forward, "-r", reverse, "-o", output])
	output = path(output).abspath()

	# call fastq-splitter.pl  -- needs to be in path
	# first get the file size
	gb_size = os.path.getsize(forward) / byte_to_gb
	parts_needed = int(ceil(float(gb_size) / float(codex_max_gb)))

	call(["fastq-splitter.pl", "--n-parts", parts_needed, output])