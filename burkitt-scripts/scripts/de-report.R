library("cummeRbund")
cuff <- readCufflinks()
cuff
pdf("qc.pdf")
d<-dispersionPlot(genes(cuff))
d
pBoxRep<-csBoxplot(genes(cuff),replicates=T)
pBoxRep
pDendro<-csDendro(genes(cuff),replicates=T)
pDendro
pBox<-csBoxplot(genes(cuff))
pBox
sigGeneIds<-getSig(cuff,alpha=0.05,level="genes")
length(sigGeneIds)
p<-PCAplot(genes(cuff),x="PC1",y="PC2",replicates=TRUE)
p
gene_diff_data <- diffData(genes(cuff))
sig_gene_data  <- subset(gene_diff_data, (significant ==  'yes'))
up_gene_data  <- subset(sig_gene_data, (log2_fold_change > 1))
geneids <- c(up_gene_data$gene_id)
myGenes <- getGenes(cuff, geneids)
csHeatmap(myGenes, cluster="both")
myDistHeat<-csDistHeat(genes(cuff))
myRepDistHeat<-csDistHeat(genes(cuff),replicates=T)
myDistHeat
myRepDistHeat
genes.PCA<-PCAplot(genes(cuff),"PC1","PC2")
genes.PCA
dev.off();
quit();