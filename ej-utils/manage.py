#!/usr/bin/python
from ej import app,db,session
from ej.models import GeneDesc
import argparse
import pandas as pd

def create_db():
    db.create_all()

def drop_db():
    db.drop_all()

def load_gene_desc():
    count = 0
    genes_descriptions = pd.read_table("./ej/data/rosseta.tsv",
        names = ["id","name","desc"]
    )
    genes = []
    for index, row in genes_descriptions.iterrows():
        gene = GeneDesc(
            id = row['id'],
            uni_name = row['name'],
            description = row['desc']
        )
        genes.append(gene)
        count += 1
    session.add_all(genes)
    session.commit()
    return count

def main():
    parser = argparse.ArgumentParser(description='Manage this Flask application.')
    parser.add_argument('command', help='the name of the command you want to run')
    args = parser.parse_args()

    if args.command == 'create_db':
        create_db()
        print "DB created!"

    elif args.command == 'delete_db':
        drop_db()
        print "DB deleted!"

    elif args.command == 'genes':
        count = load_gene_desc()
        print("loaded %d gene descriptions into db" % (count))

if __name__ == '__main__':
    main()