#!/usr/bin/make -f

<<<<<<< HEAD
OBJDIR := macs2peaks
bname :=  atac1f2_S1_L008
objects := $(addprefix $(OBJDIR)/$(bname), _model.r _peaks.narrowPeak \
	_peaks.xls _summits.bed)

#objects= macs2peaks/atac1f2_S1_L008_model.r macs2peaks/atac1f2_S1_L008_peaks.narrowPeak \
#	macs2peaks/atac1f2_S1_L008_peaks.xls macs2peaks/atac1f2_S1_L008_summits.bed

statsObjects := txt/report_$(bname)_flagstats.txt txt/report_$(bname)_stats.txt

all: $(objects) $(statsObjects)

#all: macs2peaks/atac1f2_S1_L008_model.r macs2peaks/atac1f2_S1_L008_peaks.narrowPeak macs2peaks/atac1f2_S1_L008_peaks.xls macs2peaks/atac1f2_S1_L008_summits.bed

#aln/atac1f2_S1_L008.noMito.bam txt/report_atac1f2_S1_L008_flagstats.txt txt/report_atac1f2_S1_L008_stats.txt
=======

all: aln/atac1f2_S1_L008.noDup.bam txt/report_atac1f2_S1_L008_flatstats.txt txt/report_atac1f2_S1_L008_stats.txt
>>>>>>> c304d89d5cb8de7f706a3a18bd5d0da3cb6748f1
#all: ./macs2peaks/%_model.r ./macs2peaks/%_peaks.narrowPeak ./macs2peaks/%_peaks.xls ./macs2peaks/%_summits.bed

#clean: 
#	rm  -f ./aln/%.sam ./aln/%.bam ./cutadapt/%_*.fastq.gz 

clean: 
<<<<<<< HEAD
	rm -f $(objects)
	rm -f $(statsObjects)
	
.PHONY: all clean
.INTERMEDIATE: aln/*.sam aln/*.noDup.bam cutadapt/*.fastq.gz
=======
	
	#rm -f ./aln/atac1f2_S1_L008.sam
	#rm -f ./cutadapt/*.fastq.gz cutadapt/report.txt
	
.PHONY: all clean
>>>>>>> c304d89d5cb8de7f706a3a18bd5d0da3cb6748f1

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

<<<<<<< HEAD
# find duplicates with picard  (path to picard should be set in $PICARD variable in .bash_profile or in session)
aln/%.noDup.bam txt/picard_%.txt: aln/%.sorted.bam $(PICARD)
	java -Xmx5g -jar $(PICARD) MarkDuplicates I=aln/$*.sorted.bam O=aln/$*.noDup.bam M=txt/picard_$*.txt

# 	remove mitochondrial reads
aln/%.noMito.bam: aln/%.noDup.bam
	samtools view -h -F1024 $^ | grep -v -e '\tMtDNA\t' | samtools view -b -> $@

# call peaks with macs2peaks (need to activate conda python 2.7 emvironment. you need to invoke bash shell and )
macs2peaks/%_model.r macs2peaks/%_peaks.narrowPeak macs2peaks/%_peaks.xls macs2peaks/%_summits.bed: aln/%.noMito.bam
	mkdir -p ./macs2peaks
	( bash -c "source ${HOME}/anaconda/bin/activate py27; \
		macs2 callpeak --keep-dup all -t aln/$*.noMito.bam -n macs2peaks/$* ")
	
=======
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
>>>>>>> c304d89d5cb8de7f706a3a18bd5d0da3cb6748f1



