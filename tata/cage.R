library(CAGEr)
library(BSgenome)
library(BSgenome.AKATAgfp.UW.v1)

# set each bam file
bams <- c('./var/CAGEseq-CNhi10613-AKATA-GFP.05.bam', './var/CAGEseq-CNhi10613-AKATA-GFP.04.bam')

# create obj
CAGEset <- new("CAGEset", genomeName="BSgenome.AKATAgfp.UW.v1", inputFiles = bams,
                inputFilesType="bam", sampleLabels = c("CAGEseq-04", "CAGEseq-05"))

# read the bams and detect TSS's
getCTSS(CAGEset)
ctssTable <- CTSStagCount(CAGEset)

# check if there is good correlation btwn the two
corr <- plotCorrelation(CAGEset, samples="all", method="pearson")
plotReverseCumulatives(myCAGEset, fitInRange = c(5, 10000), onePlot = TRUE)

# merge and normalize assuming corr was good
mergeSamples(CAGEset, mergeIndex= c(1,1), mergedSampleLabels = "CAGEseq-AKATA_GFP")
normalizeTagCount(CAGEset, method="powerLaw", fitInRange = c(5, 10000), alpha=1.05, T=5*10^4)

# we can look promoter width!
clusterCTSS(object = CAGEset, threshold = 1, thresholdIsTpm = TRUE,
            nrPassThreshold = 1, method = "distclu", maxDist = 20)

tagClusters <- tagClusters(CAGEset, sample="CAGEseq-AKATA_GFP")
cumulativeCTSSdistribution(CAGEset, clusters = "tagClusters")
quantilePositions(CAGEset, clusters = "tagClusters", qLow = 0.1, qUp = 0.9)
exportToBed(object = CAGEset, what = "tagClusters", qLow = 0.1, qUp = 0.9, oneFile = TRUE)