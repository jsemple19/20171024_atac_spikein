# 2018-01-04 script to align spikein atac data for C. elegans

library(QuasR)

sampleFile<-"./sampleFile.txt"
genomeFile<-"/Users/semple/Documents/MeisterLab/GenomeVer/sequence/c_elegans.PRJNA13758.WS250.genomic.fa"

#removing adaptors
#create new directory for trimmed seqs
dir.create("./trimmed",showWarnings=F)
samples<-read.delim(sampleFile,header=T,as.is=T)
samples1<-sub("rawData/", "trimmed/", samples)
names(samples1)<-names(samples)
samples1<-as.data.frame(t(samples1),stringsAsFactors=F)

#read in adaptor file
# atacAdapt<-read.delim("./atacPrimers.txt",header=T,as.is=T)
# sampleAdapt<-c("Ad1_noMX","Ad2.21_TGGGTTTC")
# atacAdapt$Sequence[match(sampleAdapt,atacAdapt$Name)]

#scroll through samples in name tables (fileName1 and fileName2 in each row) to do trimming
nexteraAdapt<-c("CTGTCTCTTATA")
for (i in 1:dim(samples)[1]) {
  preprocessReads(samples[i,1],samples1[i,1],Rpattern=nexteraAdapt)
  preprocessReads(samples[i,2],samples1[i,2],Rpattern=nexteraAdapt)
}

#create a new sampleFile so alignment will be done with trimmed seqs
write.table(samples1,file="./sampleFile_trimmed.txt",quote=F,sep="\t",row.names=F)
sampleFile<-"./sampleFile_trimmed.txt"

# do alignments
dir.create("./aln",showWarnings=FALSE)
proj<-qAlign(sampleFile,genomeFile,paired="fr",alignmentsDir="./aln")

# get alignment QC
dir.create("./plots",showWarnings=FALSE)
qQCReport(proj, "./plots/qc_report.pdf")
