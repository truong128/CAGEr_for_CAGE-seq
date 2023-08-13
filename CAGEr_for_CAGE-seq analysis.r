#CAGEr for CAGE-seq analysis
---
title: "CAGEr for CAGE-seq analysis"
author: "Truong"
date: "2023-06-29"
output: html_document
---
```{r setup, include=TRUE}
knitr::opts_chunk$set(echo = TRUE)

library(CAGEr)
Truong <- CAGEexp( genomeName     = "BSgenome.Hsapiens.UCSC.hg38", inputFiles     = c("/Users/truong/Documents/CAGE-seq/BRRF1.1.fastq.gz.hg38plusAkata_inverted-Aligned.sortedByCoord.out.bam", "/Users/truong/Documents/CAGE-seq/BRRF1.2.fastq.gz.hg38plusAkata_inverted-Aligned.sortedByCoord.out.bam", "/Users/truong/Documents/CAGE-seq/BRRF1.3.fastq.gz.hg38plusAkata_inverted-Aligned.sortedByCoord.out.bam", "/Users/truong/Documents/CAGE-seq/C1.2.fastq.gz.hg38plusAkata_inverted-Aligned.sortedByCoord.out.bam", "/Users/truong/Documents/CAGE-seq/C2.2.fastq.gz.hg38plusAkata_inverted-Aligned.sortedByCoord.out.bam", "/Users/truong/Documents/CAGE-seq/C3.2.fastq.gz.hg38plusAkata_inverted-Aligned.sortedByCoord.out.bam") , inputFilesType = "bam", sampleLabels   = c("BR1", "BR2", "BR3", "C1", "C2", "C3"))

Truong

colData(Truong)

TSS <- getCTSS(Truong)

#output_file <- "/Users/truong/Documents/Trang_ce.rds"
#saveRDS(cee, file = output_file)
# Load the saved cee object from the output file
#loaded_cee <- readRDS(file = output_file)

TSS
CTSStagCountSE(TSS)
CTSScoordinatesGR(TSS)
CTSStagCountDF(TSS)
CTSStagCountGR(TSS, 1)  # GRanges for one sample with expression count.
sampleLabels(TSS)
#install.packages("rtracklayer")
library(rtracklayer)
#path of gtf = /Users/truong/Documents/RPMS1_YCCEL1/hg38_plus_Akata_inverted.gtf
gtf <- import("/Users/truong/Documents/RPMS1_YCCEL1/hg38_plus_Akata_inverted.gtf")
TSSano <- annotateCTSS(TSS, gtf)
colData(TSSano)[,c("librarySizes", "promoter", "exon", "intron", "unknown")]
#CTSScoordinatesGR(TRang)


#This function is to Merge and rearrange all the samples: In this case, I merge 3 samples (first three 1,2,3) together and assign it into Treatment, and merge 3 samples (last samples three 4,5,6) together and assign it into Control
#cee1 = mergeSamples(cee, mergeIndex = c(1,1,1,2,2,2), mergedSampleLabels = c("Treatment", "Control"))
#Trang <- annotateCTSS(Trang, gtf)

corr.m <- plotCorrelation2(Trang, samples = "all"
                            , tagCountThreshold = 1, applyThresholdBoth = FALSE
                            , method = "pearson")
librarySizes(TSSano)
plotReverseCumulatives(TSSano, fitInRange = c(5, 1000), onePlot = TRUE)
Truongnor <- normalizeTagCount(TSSano, method = "powerLaw", fitInRange = c(5, 1000), alpha = 1.2, T = 10^6)
Truongnor[["tagCountMatrix"]]
Truongcluster <- clusterCTSS( Truongnor
                   , threshold = 1
                   , thresholdIsTpm = TRUE
                   , nrPassThreshold = 1
                   , method = "distclu"
                   , maxDist = 20
                   , removeSingletons = TRUE
                   , keepSingletonsAbove = 5)

BR1<- tagClustersGR(Trang, sample = "BR1")
BR2<- tagClustersGR(Trang, sample = "BR2")
BR3<- tagClustersGR(Trang, sample = "BR3")
C1<- tagClustersGR(Trang, sample = "C1")
C2<- tagClustersGR(Trang, sample = "C2")
C3<- tagClustersGR(Trang, sample = "C3")
write.csv(BR1, file = "/Users/truong/Documents/CAGE-seq/CTSSs_BR1.csv", row.names = TRUE)
write.csv(BR2, file = "/Users/truong/Documents/CAGE-seq/CTSSs_BR2.csv", row.names = TRUE)
write.csv(BR3, file = "/Users/truong/Documents/CAGE-seq/CTSSs_BR3.csv", row.names = TRUE)
write.csv(C1, file = "/Users/truong/Documents/CAGE-seq/CTSSs_C1.csv", row.names = TRUE)
write.csv(C2, file = "/Users/truong/Documents/CAGE-seq/CTSSs_C2.csv", row.names = TRUE)
write.csv(C3, file = "/Users/truong/Documents/CAGE-seq/CTSSs_C3.csv", row.names = TRUE)
Truongcum <- cumulativeCTSSdistribution(Truongcluster, clusters = "tagClusters", useMulticore = T)
Truongquan <- quantilePositions(Truongcum, clusters = "tagClusters", qLow = 0.1, qUp = 0.9)

tagClustersGR( Truongquan, "BR1", returnInterquantileWidth = TRUE,  qLow = 0.1, qUp = 0.9)

plotInterquantileWidth(Truongquan, clusters = "tagClusters", tpmThreshold = 5, qLow = 0.1, qUp = 0.9)

Truongagg <- aggregateTagClusters(Truongquan, tpmThreshold = 5, qLow = 0.1, qUp = 0.9, maxDist = 100)

Truongagg <- aggregateTagClusters(Truongquan, tpmThreshold = 5, excludeSignalBelowThreshold = TRUE, qLow = NULL, qUp = NULL, maxDist = 100)

Trang$outOfClusters / Trang$librarySizes *100

Truongcon <- consensusClustersGR(Truongagg)

Truongconsen <- annotateConsensusClusters(Truongcon, gtf)

Truongcumu <- cumulativeCTSSdistribution(Truongagg, clusters = "consensusClusters", useMulticore = TRUE)

Truongbest <- quantilePositions(Truongcumu, clusters = "consensusClusters", qLow = 0.1, qUp = 0.9, useMulticore = TRUE)

consensusClustersGR( cee, sample = "R1", returnInterquantileWidth = TRUE,  qLow = 0.1, qUp = 0.9)

trk <- exportToTrack(CTSSnormalizedTpmGR(cee, "R1"))
cee |> CTSSnormalizedTpmGR("all") |> exportToTrack(cee, oneTrack = FALSE)

iqtrack <- exportToTrack(cee, what = "tagClusters", qLow = 0.1, qUp = 0.9, oneTrack = FALSE)

rtracklayer::export.bed(iqtrack,"/Users/truong/Documents/HiC-Analysis/HiC_files_4lanes/Galaxy_Output/outputFileName_CAGEr.bed")
#BiocFileList of length 2

cee <- getExpressionProfiles(cee, what = "consensusClusters", tpmThreshold = 10, 
                             nrPassThreshold = 1, method = "som", xDim = 4, yDim = 2)

#consensusClustersGR(Trang)$exprClass |> table(useNA = "always")

#plotExpressionProfiles(cee, what = "consensusClusters")

#consensusClustersGR(cee) |> subset(consensusClustersGR(cee)$exprClass ==  "0_1")


Truongbest$group <- factor(c("Test", "Test", "Test", "Control", "Control", "Control"))
tada <- consensusClustersDESeq2(Truongbest, ~group)
tada
library(DESeq2)
tada = DESeq(tada)

resultsNames(tada)
[1] "Intercept"    "group_b_vs_a"
ohyeah = results(tada)
head(results(tada, tidy =TRUE))