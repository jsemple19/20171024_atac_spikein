#!/usr/bin/make -f


all: aln/atac1f2_S1_L008.noDup.bam txt/report_atac1f2_S1_L008_flatstats.txt txt/report_atac1f2_S1_L008_stats.txt
#all: ./macs2peaks/%_model.r ./macs2peaks/%_peaks.narrowPeak ./macs2peaks/%_peaks.xls ./macs2peaks/%_summits.bed

#clean: 
#	rm  -f ./aln/%.sam ./aln/%.bam ./cutadapt/%_*.fastq.gz 

clean: 
	
	#rm -f ./aln/atac1f2_S1_L008.sam
	#rm -f ./cutadapt/*.fastq.gz cutadapt/report.txt
	
.PHONY: all clean

########### VARIABLES #########
genomefile:=/Users/semple/Documents/MeisterLab/GenomeVer/sequence/c_elegans.PRJNA13758.WS250.genomic.fa


########### RULES  ############

# use cutadapt to trim nextera sequence in paired end mode
cutadapt/%_R1_001.fastq.gz cutadapt/%_R2_001.fastq.gz: rawData/%_R1_001.fastq.gz rawData/%_R2_001.fastq.gz
	mkdir -p cutadapt
	cutadapt -a CTGTCTCTTATA -A CTGTCTCTTATA -q 10 -m 25 -o cutadapt/$*_R1_001.fastq.gz -p cutadapt/$*_R2_001.fastq.gz \
	rawData/$*_R1_001.fastq.gz rawData/$*_R2_001.fastq.gz >cutadapt/report.txt

# use bwa-mem to align sequences
aln/%.sam: $(genomefile) cutadapt/%_R1_001.fastq.gz cutadapt/%_R2_001.fastq.gz 
	mkdir -p aln
	bwa mem -t 2 $^ > $@ 

# use samtools to convert to bam, sort and clean
aln/%.sorted.bam: aln/%.sam
	samtools view -q 10 -F2304 -u $^ | samtools sort -o $@

# 	get alignment stats
txt/report_%_flagstats.txt: aln/%.sorted.bam
	samtools flagstat $^ > $@

txt/report_%_stats.txt: aln/%.sorted.bam
	samtools stats $^ > $@

# find duplicates with picard	
aln/%.noDup.bam txt/picard_%.txt: aln/%.sorted.bam $(PICARD)
	java -Xmx5g -jar $(PICARD) MarkDuplicates I=aln/$*.sorted.bam O=aln/$*.noDup.bam M=txt/picard_$*.txt


# 	samtools stats $outbam > ./txt/report_stats.txt
# 	use picard to remove duplicates:
# 	cleanbam=`echo $outbam | sed "s@.bam@.clean.bam@"`
# 	java -Xmx5g -jar $PICARD MarkDuplicates I=$outbam O=$cleanbam M=./txt/picard_out.txt
# 	to remove mitochondrial reads
# 	samtools view -h -F1024 $cleanbam | grep -v -e '\tMtDNA\t' | samtools view -b -> $tmpbam

# ordinal number of samFile input from command line arg
# length=`echo ${#samfiles[@]}`
# let length=length-1
# for i in `seq 0 $length`; do
# 	insam=`echo ${samfiles[$i]}`
# 	outbam=`echo $insam | sed "s@.sam@.bam@"`
# 	#convert samfiles to bamfiles and sort:
# 	samtools view -q 10 -F2304 -u ${samfiles[$i]} | samtools sort -o $outbam 
# 	get alignment stats
# 	samtools flagstat $outbam > ./txt/report_flagstats.txt
# 	samtools stats $outbam > ./txt/report_stats.txt
# 	use picard to remove duplicates:
# 	cleanbam=`echo $outbam | sed "s@.bam@.clean.bam@"`
# 	java -Xmx5g -jar $PICARD MarkDuplicates I=$outbam O=$cleanbam M=./txt/picard_out.txt
# 	to remove mitochondrial reads
# 	samtools view -h -F1024 $cleanbam | grep -v -e '\tMtDNA\t' | samtools view -b -> $tmpbam
# done



