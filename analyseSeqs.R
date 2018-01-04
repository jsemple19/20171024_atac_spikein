# 2018-01-04 script to align spikein atac data for C. elegans

library(QuasR)

sampleFile<-"./sampleFile.txt"
genomeFile<-"/Users/semple/Documents/MeisterLab/GenomeVer/sequence/c_elegans.PRJNA13758.WS250.genomic.fa"

dir.create("./aln",showWarnings=FALSE)
proj<-qAlign(sampleFile,genomeFile,paired="fr",alignmentsDir="./aln")
