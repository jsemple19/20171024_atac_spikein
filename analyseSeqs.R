# 2018-01-04 script to align spikein atac data for C. elegans

library(QuasR)

sampleFile<-"txt/sampleFile.txt"
genomeFile<-"/Users/semple/Documents/MeisterLab/GenomeVer/sequence/c_elegans.PRJNA13758.WS250.genomic.fa"

#removing adaptors
#create new directory for trimmed seqs
# dir.create("./trimmed",showWarnings=F)
# samples<-read.delim(sampleFile,header=T,as.is=T)
# samples1<-sub("rawData/", "trimmed/", samples)
# names(samples1)<-names(samples)
# samples1<-as.data.frame(t(samples1),stringsAsFactors=F)

#read in adaptor file
# atacAdapt<-read.delim("./atacPrimers.txt",header=T,as.is=T)
# sampleAdapt<-c("Ad1_noMX","Ad2.21_TGGGTTTC")
# atacAdapt$Sequence[match(sampleAdapt,atacAdapt$Name)]

#scroll through samples in name tables (fileName1 and fileName2 in each row) to do trimming
# nexteraAdapt<-c("CTGTCTCTTATA")
# for (i in 1:dim(samples)[1]) {
#   res<-preprocessReads(samples[i,1],samples1[i,1],Rpattern=nexteraAdapt,Lpattern=nexteraAdapt)
#   res
#   res<-preprocessReads(samples[i,2],samples1[i,2],Rpattern=nexteraAdapt, Lpattern=nexteraAdapt)
#   res
# }
#
# #create a new sampleFile so alignment will be done with trimmed seqs
# write.table(samples1,file="./txt/sampleFile_trimmed.txt",quote=F,sep="\t",row.names=F)
# sampleFile<-"./txt/sampleFile_trimmed.txt"

# #Align to E coli genome as auxiliary file
# library(BSgenome.Ecoli.NCBI.20080805)
#
# #write genome to fasta file
# dir.create("./tmp",showWarnings=F)
# export(Ecoli,"./tmp/Ecoli.fasta")
#
# #create auxiliary file
# AuxInput=as.data.frame(cbind(
#   FileName="./tmp/Ecoli.fasta",
#   AuxName="Ecoli"))
# write.table(AuxInput,'./txt/auxiliaryFile.txt',quote=F,row.names=F,sep='\t')

# do alignments
dir.create("./aln",showWarnings=FALSE)
QuasRdef='-k 2 --best --strata'
#cluObj=makeCluster(2)

proj<-qAlign(sampleFile,
             genomeFile,
             paired="fr",
             #auxiliaryFile="./txt/auxiliaryFile.txt",
             alignmentsDir="aln/",
             alignmentParameter=paste('-e 150 -X 2000 ',QuasRdef,sep=''))
             # bowtie_usage() for parameters:
             #"  -k <int>           report up to <int> good alignments per read (default: 1)"
             #"  --best             hits guaranteed best stratum; ties broken by quality"
             #"  --strata           hits in sub-optimal strata aren't reported (requires --best)"
             #"  -e/--maqerr <int>  max sum of mismatch quals across alignment for -n (def: 70)"
             #"  -X/--maxins <int>  maximum insert size for paired-end alignment (default: 250)")

# get alignment QC
dir.create("./plots",showWarnings=FALSE)
qQCReport(proj, "./plots/qc_report.pdf")


