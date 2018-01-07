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
	samtools view -q 10 -F2304 -u ${samfiles[$i]} | samtools sort -o $outbam 
	#get alignment stats
	samtools flagstat $outbam > ./txt/report_flagstats.txt
	samtools stats $outbam > ./txt/report_stats.txt
	#use picard to remove duplicates:
	cleanbam=`echo $outbam | sed "s@.bam@.clean.bam@"`
	java -Xmx5g -jar $PICARD MarkDuplicates I=$outbam O=$cleanbam M=./txt/picard_out.txt
	#to remove mitochondrial reads
	#samtools view -h -F1024 $cleanbam | grep -v -e '\tMtDNA\t' | samtools view -b -> $tmpbam
done



# create arrays of file names
cleanfiles=(`ls ./aln/*.clean.bam`)
mkdir -p ./macs2peaks

source activate py27
# ordinal number of samFile input from command line arg
length=`echo ${#cleanfiles[@]}`
let length=length-1
for i in `seq 0 $length`; do
	cleanbam=`echo ${cleanfiles[$i]}`
	bname=`echo $cleanbam | sed "s@.clean.bam@@"`
	bname=`echo $bname | sed "s@./aln@./macs2peaks@"`
	macs2 callpeak --keep-dup all -t $cleanbam -n ${bname}
done	

source deactivate


