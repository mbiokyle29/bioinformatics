Release=80
from cogent.db.ensembl import Species, Genome
human = Genome(Species='human', Release=Release, account=None)
gene = human.getGeneByStableId(StableId='ENSG00000205274')
print gene.Symbol