# ej-utils
Platform for deploying in-house lab services

# setting up the application
Use virtualenv ( at the project root)
```bash
virtualenv venv
. venv/bin/activate
venv/bin/pip install -r requirements.txt
```

# setting up the gene description records in the db
Download the gene info file from ncbi: http://ftp://ftp.ncbi.nlm.nih.gov/gene/data/gene_info.gz
```bash
mv gene_info.gz $PROJECTROOT/ej/data/
cd ej/data/
gunzip gene_info.gz
awk 'BEGIN {FS = "\t"}; $1==9606 {print $2"\t"$3"\t"$9}' ej/data/gene_info > rosetta.tsv
```

Then run manage.py
```bash
python manage.py create_db
python manage.py genes
```

ready to go!
```bash
python runserver.py
```
