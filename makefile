#!/usr/bin/make -f
## required software: fastqc, cutadapt, bwa, samtools, picard, qualimap, bedtools, macs2peaks

########### VARIABLES #########

genomefile:=/Users/semple/Documents/MeisterLab/GenomeVer/sequence/c_elegans.PRJNA13758.WS250.genomic.fa
OBJDIR := macs2peaks
bname :=  atac1f2_S1_L008
longbname := $(addprefix $(bname), _R1_001 _R2_001)

#list of the final macs2peaks output files
objects := $(addprefix $(OBJDIR)/$(bname), _model.r _peaks.narrowPeak \
	_peaks.xls _summits.bed)

#list of the various reports and stats produced during analysis
statsObjects := aln/report_$(bname)_flagstats.txt aln/report_$(bname)_stats.txt \
	aln/report_$(bname)_picard.txt cutadapt/report_$(bname)_cutadapt.txt \
	$(addsuffix _fastqc.html, $(addprefix fastQC/, $(longbname))) \
	$(addsuffix _picard_insert_size_metrics.txt, $(addprefix aln/, $(bname))) \
	$(addsuffix _picard_insert_size_histogram.pdf, $(addprefix aln/, $(bname))) \
	$(addsuffix _report_qualimap.pdf, $(addprefix aln/, $(bname)))
	
#list of files to delete at the end of the analysis
intermediateFiles := aln/$(bname).sam aln/$(bname).noDup.bam aln/$(bname).sorted.bam \
	$(addsuffix .fastq.gz, $(addprefix cutadapt/, $(longbname)))

#list of intermediate files to keep
secondaryFiles := $(addsuffix .noMito.bam, aln/$(bname))


########### RULES  ############

all: $(objects) $(statsObjects)

#use cleanall when you want to force rerunning all the analysis
cleanall:
	rm -f $(objects)
	rm -f $(statsObjects)
	rm -f $(secondaryFiles)
	
#use clean if the intermediate files are not properly removed (should not be required)
clean:
	rm -f $(intermediateFiles)

.PHONY: all clean cleanall
.INTERMEDIATE: $(intermediateFiles)
.SECONDARY: $(secondaryFiles)

#run fastqc on downloaded sequences
fastQC/%_fastqc.html: rawData/%.fastq.gz
	mkdir -p fastQC
	fastqc $^ -o fastQC 

# use cutadapt to trim nextera sequence in paired end mode
cutadapt/%_R1_001.fastq.gz cutadapt/%_R2_001.fastq.gz cutadapt/report_%_cutadapt.txt: \
 rawData/%_R1_001.fastq.gz rawData/%_R2_001.fastq.gz
	mkdir -p cutadapt
	cutadapt -a CTGTCTCTTATA -A CTGTCTCTTATA -q 10 -m 25 -o cutadapt/$*_R1_001.fastq.gz -p cutadapt/$*_R2_001.fastq.gz \
	rawData/$*_R1_001.fastq.gz rawData/$*_R2_001.fastq.gz >cutadapt/report_$*_cutadapt.txt

# use bwa-mem to align sequences
aln/%.sam: $(genomefile) cutadapt/%_R1_001.fastq.gz cutadapt/%_R2_001.fastq.gz
	mkdir -p aln
	bwa mem -t 2 $^ > $@

# use samtools to convert to bam, sort and clean
aln/%.sorted.bam: aln/%.sam
	samtools view -q 30 -F 1804 -u $^ | samtools sort -o $@


# 	get alignment stats
aln/report_%_flagstats.txt: aln/%.sorted.bam
	samtools flagstat $^ > $@

aln/report_%_stats.txt: aln/%.sorted.bam
	samtools stats $^ > $@

# mark duplicates with picard (path to picard should be set in $PICARD variable in .bash_profile or in session)
aln/%.noDup.bam aln/report_%_picard.txt: aln/%.sorted.bam ${PICARD}
	java -Xmx5g -jar ${PICARD} MarkDuplicates I=aln/$*.sorted.bam O=aln/$*.noDup.bam M=aln/report_$*_picard.txt

	
# Get insert size statistics and plots with picard and qualimap (set path to $QUALIMAP in .bash_profile)
aln/%_picard_insert_size_metrics.txt aln/%_picard_insert_size_histogram.pdf \
aln/%_report_qualimap.pdf: aln/%.noDup.bam ${PICARD} ${QUALIMAP}
	java -Xms1g -Xmx5g -jar ${PICARD} CollectInsertSizeMetrics I=aln/$*.noDup.bam \
       O=aln/$*_picard_insert_size_metrics.txt \
       H=aln/$*_picard_insert_size_histogram.pdf
	${QUALIMAP} bamqc -bam aln/$*.noDup.bam -c -outdir aln -outfile $*_report_qualimap.pdf -outformat PDF

      
# 	remove mitochondrial reads
aln/%.noMito.bam: aln/%.noDup.bam
	samtools view -h -F 1024 $^ | grep -v -e '\tMtDNA\t' | samtools view -b -> $@
	
# convert 


# call peaks with macs2peaks (need to activate conda python 2.7 emvironment. you need to invoke bash shell and )
macs2peaks/%_model.r macs2peaks/%_peaks.narrowPeak macs2peaks/%_peaks.xls macs2peaks/%_summits.bed: aln/%.noMito.bam
	mkdir -p ./macs2peaks
	( bash -c "source ${HOME}/anaconda/bin/activate py27; \
		macs2 callpeak -t aln/$*.noMito.bam -f BAM -n macs2peaks/$* -g 9e7 -p 5e-2 \
		--nomodel --extsize 150 --shift -75 -B --keep-dup all --call-summits" )
	#( bash -c "source ${HOME}/anaconda/bin/activate py27; \
	#	macs2 callpeak --format BAMPE --keep-dup all -t aln/$*.noMito.bam -n macs2peaks/$* ")


