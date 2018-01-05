#! /bin/bash
#cutadapt -a  CTGTCTCTTATA -g ^CTGTCTCTTATA -q 10 -m 20 

fqFilesR1=(`ls ./rawData/*_R1_*fastq.gz`)
fqFilesR2=(`ls ./rawData/*_R2_*fastq.gz`)

mkdir -p ./cutadapt

#length of array
length=`echo ${#fqFilesR1[@]}`
let length=length-1
for i in `seq 0 $length`; do
	r1=`echo ${fqFilesR1[$i]}`
	r2=`echo ${fqFilesR2[$i]}`
	#put output files in cutadapt folder with same filename
	out1=`echo $r1 | sed "s@rawData@cutadapt@"`
	out2=`echo $r2 | sed "s@rawData@cutadapt@"`
	#cut 3' adaptors and anchored 5' adaptors with nextera sequence in paired end mode
	cutadapt -a CTGTCTCTTATA -A CTGTCTCTTATA -q 10 -m 25 -o $out1 -p $out2 $r1 $r2 >cutadapt/report.txt
done

#align reads with bwa
fqFilesR1=(`ls ./cutadapt/*_R1_*fastq.gz`)
fqFilesR2=(`ls ./cutadapt/*_R2_*fastq.gz`)

genomefile="/Users/semple/Documents/MeisterLab/GenomeVer/sequence/c_elegans.PRJNA13758.WS250.genomic.fa"

#length of array
length=`echo ${#fqFilesR1[@]}`
let length=length-1
for i in `seq 0 $length`; do
	#r1=`echo ${fqFilesR1[$i]%.fastq.gz} | cut -f1 -d'/'`
	r1=`echo ${fqFilesR1[$i]}`
	r2=`echo ${fqFilesR2[$i]}`
	#change output directory to ./aln
	out=`echo $r1 | sed "s@cutadapt@aln@"`
	#change file extension to .sam
	out=`echo $out | sed "s@.fastq.gz@.sam@"`
	#perform alignment with 2 cores and default parameter
	bwa mem -t 2 $genomefile ${fqFilesR1[$i]} ${fqFilesR2[$i]} > $out 		
done

# create arrays of file names
samfiles=(`ls ./aln/*.sam`)

# ordinal number of samFile input from command line arg
length=`echo ${#samfiles[@]}`
let length=length-1
for i in `seq 0 $length`; do
	insam=`echo ${samfiles[$i]}`
	outbam=`echo $insam | sed "s@.sam@.bam@"`
	##convert samfiles to bamfiles and sort:
	samtools view -u ${samfiles[$i]} | samtools sort -o $outbam 
	#get alignment stats
	samtools flagstats $outbam > ./txt/report_flagstats.txt
	samtools stats $outbam > ./txt/report_stats.txt
done


# 
# My 2p, not really answering the question about the 30bp trimming. I had good results (=sensible) processing ATAC-Seq as pretty much a standard ChIP-Seq. This is my pipeline for reads of 75 bp:
# 
# Trim adapters with cutadapt (or similar)
# 
# Align with bwa mem with default settings
# 
# Discard alignments with flag 2304 (not primary or supplementary alignment) and mapq < 10
# 
# If a "blacklist" of ChIP-seq genomic regions is available, remove reads in these regions
# 
# Mark duplicates with picard
# 
# Remove duplicates and reads on mitochondrial genome (ATAC-Seq maps a lot of reads on chrM, up to 40-60% and I prefer to remove them before peak calling)
# 
# Call peaks with macs2. For human samples I find in tens to hundreds of thousands of peaks looking pretty sharp, often in proximity of TSS.
# 
# More detail here https://github.com/sblab-bioinformatics/dna-secondary-struct-chrom-lands/blob/master/Methods.md where I put this as supplementary material to a paper (1) we recently published using ATAC-Seq (and other stuff).

