#!/usr/bin/env python
"""
venn.py
Kyle McChesney

Generate a venn diagram based on gene sets
--A         first list as file
--B         second list as file
--C         third list as file
--interest  genes to give points
"""

import logging, argparse, os
from matplotlib import pyplot as plt
from matplotlib_venn import venn2
from sets import Set
import random as rand
import time

# loggin
logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)

def main():
    # args
    parser = argparse.ArgumentParser(
        description = (" Generate a venn diagram from Gene Lists "),
    )

    parser.add_argument("--a", help="first list file", required=True)
    parser.add_argument("--b", help="second list file", required=True)
    parser.add_argument("--labelA", help="first list label", default="A")
    parser.add_argument("--labelB", help="second list label", default="B")
    parser.add_argument("--c", help="third optional list file")
    parser.add_argument("--interest", help=" interested genes to be plotted as points")
    parser.add_argument("--title", help=" Text for plot title", default="A vs B Venn Diagram Plot")
    args = parser.parse_args()

    if(args.c):
        log.warn("Three circle not supported yet")
        raise SystemExit

    firstGenes = load_gene_list(args.a)
    secondGenes = load_gene_list(args.b)

    plt.figure(figsize=(10,10))
    v = venn2([firstGenes, secondGenes],(args.labelA,args.labelB))
    plt.title(args.title)


    if(args.interest):
        interest_genes = load_gene_list(args.interest)
        plot_point_on_circle(v, 0, interest_genes, plt)   

    ts = str(time.time()).replace(".","")
    plt.savefig(args.labelA+"-"+args.labelB+"-"+ts+".png")


def load_gene_list(file):

    # does it exist
    if not os.path.isfile(file):
        log.warn("%s does not exist", file)
        raise SystemExit

    # read er
    list = [ line.rstrip('\n') for line in open(file)]
    return set(list)

def plot_point_on_circle(venn, index, genes, plt):

    # get center and radius
    center = venn.get_circle_center(index)
    radius = venn.get_circle_radius(index)

    x_points = {}
    y_points = {}

    for gene in genes:

        # random float from 0 --> 1
        found = False
        while not found:
            x = rand.uniform(center[0] - radius, center[0] + radius)
            y = rand.uniform(center[1] - radius, center[1] + radius)

            # hasn't been used
            # this is nasty
            if not x in x_points and not y in y_points:
                plt.text(x,y,gene)
                found = True
                x_points[x] = 1
                y_points[y] = 1

    return plt
if __name__ == "__main__":
    main()