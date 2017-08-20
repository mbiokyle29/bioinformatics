#!/usr/bin/env python
"""
cms.py
Kyle McChesney

cDNA Motif Searcher!

"""

import logging, argparse, os, re
from GeneRecord import GeneRecord
from itertools import islice

# loggin
logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)

def main():
    # args
    parser = argparse.ArgumentParser(
        description = (" Scan a list of cDNA entries for a given motif "),
    )

    parser.add_argument("--motif", help="the k-mer to search for", required=True)
    parser.add_argument("--mfasta", help="fullpath to multi fasta", required=True)
    parser.add_argument("--count", help="Number of hits needed", type=int, required=True)
    
    args = parser.parse_args()

    log.info("Starting motif scan on: %s with %s", args.mfasta, args.motif)
    log.info("Cutoff is: %i", args.count)
    
    recs = scan(args.mfasta, args.motif, args.count)
    log.info("Found %i matches for %s", len(recs), args.motif)
    for rec in recs:
        print rec

def scan(file, motif, cutoff):

    # check if exist
    if not os.path.isfile(file):
        log.warn(" %s is not a file!", file)
        raise SystemExit

    # compile the re's
    patt = re.compile(motif)

    gene = ".*gene:(ENSG\d+).*"
    gene_re = re.compile(gene)

    # for results
    matches = []
    genes = {}

    # read er
    line_count = 0
    header = None
    seq = ""
    with open(file, "r") as fh:
        for line in fh:

            if line.startswith(">"):

                if header is not None:

                    # check if it is a hit
                    motifs = len(patt.findall(seq))
                    if motifs >= cutoff:
                        gene_id = parse_gene_from_header(header, gene_re)

                        if gene_id not in genes:
                            header = ">{}".format(gene_id)
                            matches.append(GeneRecord(header,seq,motifs))
                            seq = ""
                            genes[gene_id] = 1

                        else:
                            genes[gene_id] += 1

                header = line

            else:
                seq = seq+line.rstrip()
            line_count += 1

    log.info("Processed %i lines", line_count)
    dump_genes(genes)
    return matches

# we just want the gene id or gene:ENSG00000241255
def parse_gene_from_header(header, regex):
    match = re.match(regex, header)
    return match.group(1)


def dump_genes(genes):
    log.info("Dumping unqiue gene hash")

    with open("genes.tsv", "w+") as fh:
        for key in genes:
            fh.write("{}\t{}\n".format(key, str(genes[key])))

if __name__ == "__main__":
    main()