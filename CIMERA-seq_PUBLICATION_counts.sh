#! /bin/sh

#Structure of files
#Directory (folder) for every sample
#Each sample directory contains 2 FASTQ files, one for each end of paired-end sequencing

#    path
#    │   
#    │       
#    │
#    └───to
#        │   
#        │       
#        │
#        └───files
#             │   
#             │       
#             │
#             └───sample1
#             │       sample1_R1.fastq.gz
#             │       sample1_R2.fastq.gz
#             │
#             └───sample2
#             │       sample2_R1.fastq.gz
#             │       sample2_R2.fastq.gz
#             │
#             └───sample3
#                     sample3_R1.fastq.gz
#                     sample3_R2.fastq.gz


#Create an environment variable for the path to the directory where the sample directories are stored
#Make sure the path ends with a forward slash

#Format
#directory="/path/to/files/"

#Example
directory="/media/sf_UbuntuSharing/files/"

#Create an environment variable for the path to the directory where the reference files are stored

#Format
#reference_directory="/path/to/references/"

#Example
reference_directory="/media/sf_UbuntuSharing/mouse/"

#FASTQ file names (suffix .fastq not .fq)
#Assuming paired-end reads are in format "file1_R1.fastq" and "file1_R2.fastq"

#Format
#samples="
#sample1
#sample2
#sample3
#"

#Example
samples="
CS1
CS2
"

##############################################################################
#Initialize packages used

#FastQC
PATH=$PATH:FastQC/
#FLASh
PATH=$PATH:FLASH-1.2.11-Linux-x86_64/
#Cutadapt
PATH=$PATH:/home/william/.local/bin/
#BLAST
PATH=$PATH:ncbi-blast-2.12.0+/bin/


#for sample in $samples
#do

###Quality Control (FastQC)###

#Make output directory for FastQC reports

#	mkdir ${directory}FastQC_Reports

#Run FastQC

#	fastqc \
#	-o ${directory}FastQC_Reports \
#	${directory}${sample}/${sample}_R1.fastq.gz

#	fastqc \
#	-o ${directory}FastQC_Reports \
#	${directory}${sample}/${sample}_R2.fastq.gz

#done


###Quality Control (MultiQC)###

#Make output directory for MultiQC report

#	mkdir ${directory}MultiQC_Report

#Run MultiQC
#For further usage details: https://multiqc.info/docs/

#	multiqc ${directory}FastQC_Reports -o ${directory}MultiQC_Report


for sample in $samples
do

##########
	[ -e ${directory}${sample}/${sample}.summary.txt ] && rm ${directory}${sample}/${sample}.summary.txt

	echo "$sample" >> ${directory}${sample}/${sample}.summary.txt

	echo "Start:" >> ${directory}${sample}/${sample}.summary.txt
	date >> ${directory}${sample}/${sample}.summary.txt

	gunzip -k ${directory}${sample}/${sample}_R1.fastq.gz
	echo "${sample} R1 Reads: $(( $(wc -l ${directory}${sample}/${sample}_R1.fastq | awk '{print $1}') / 4 ))" >> ${directory}${sample}/${sample}.summary.txt
	rm ${directory}${sample}/${sample}_R1.fastq

	gunzip -k ${directory}${sample}/${sample}_R2.fastq.gz
	echo "${sample} R2 Reads: $(( $(wc -l ${directory}${sample}/${sample}_R2.fastq | awk '{print $1}') / 4 ))" >> ${directory}${sample}/${sample}.summary.txt
	rm ${directory}${sample}/${sample}_R2.fastq
##########	

###Combine Paired-End Reads (FLASh)###

#Make directory for FLASh outpout

	mkdir ${directory}${sample}/${sample}_FLASh

#Run FLASh
#For further usage details: http://gensoft.pasteur.fr/docs/FLASH/1.2.11/flash

	flash \
	--allow-outies \
	--output-directory=${directory}${sample}/${sample}_FLASh/ \
	--output-prefix=${sample} \
	--max-overlap=150 \
	--min-overlap=6 \
	--compress \
	${directory}${sample}/${sample}_R1.fastq.gz \
	${directory}${sample}/${sample}_R2.fastq.gz \
	2>&1 | tee ${directory}${sample}/${sample}_FLASh/FLASh_${sample}.log


##########
	gunzip -k ${directory}${sample}/${sample}_FLASh/${sample}.extendedFrags.fastq.gz
	echo "${sample} FLASh Reads: $(( $(wc -l ${directory}${sample}/${sample}_FLASh/${sample}.extendedFrags.fastq | awk '{print $1}') / 4 ))" >> ${directory}${sample}/${sample}.summary.txt
	rm ${directory}${sample}/${sample}_FLASh/${sample}.extendedFrags.fastq
##########


###Remove Adapters, Qualty Trimming, Minimum Length Threshold (Cutadapt)###

#Run Cutadapt
#For further usage details: https://cutadapt.readthedocs.io/en/stable/guide.html

	cutadapt \
	-g CTACAGTCCGACGATC...TGGAATTCTCGGGTGCCAAGG \
	-q 30 \
	-m 30 \
	-o ${directory}${sample}/${sample}.cutadapt.fastq.gz \
	${directory}${sample}/${sample}_FLASh/${sample}.extendedFrags.fastq.gz

	rm -r ${directory}${sample}/${sample}_FLASh


###Deduplicate Reads###

#Unzip FASTQ file from cutadapt

	gunzip ${directory}${sample}/${sample}.cutadapt.fastq.gz

##########
	echo "${sample} Cutadapt Reads: $(( $(wc -l ${directory}${sample}/${sample}.cutadapt.fastq | awk '{print $1}') / 4 ))" >> ${directory}${sample}/${sample}.summary.txt
##########


#Remove duplicate reads
#Read ID now contains rank (by duplicate number) and number of duplicates

	awk 'NR%4==2' ${directory}${sample}/${sample}.cutadapt.fastq | \
	sort -T ${directory}${sample} | \
	uniq -c | \
	sort -k1,1nr -T ${directory}${sample} | \
	awk '{print $0,NR}' | \
	awk '{print ">"$3"-"$1"\n"$2}' \
	> ${directory}${sample}/${sample}.cutadapt.deduped.fasta

	rm ${directory}${sample}/${sample}.cutadapt.fastq

##########
	echo "${sample} Deduplicated Reads: $(( $(wc -l ${directory}${sample}/${sample}.cutadapt.deduped.fasta | awk '{print $1}') / 2 ))" >> ${directory}${sample}/${sample}.summary.txt
##########


###Barcode Reads###

	awk 'BEGIN{RS=">";OFS="\t"}NR>1{print ">"$1,$2}' ${directory}${sample}/${sample}.cutadapt.deduped.fasta | \
	awk '{print $1"\t"$2"\t"substr($2,1,4)""substr($2,length($2)-1,2)}' | \
	awk '{$2 = substr($2, 5); print }' | \
	awk '{$2 = substr($2, 1, length($2)-2); print }' | \
	awk '{print $1"_"$3"\n"$2}' \
	> ${directory}${sample}/${sample}.cutadapt.deduped.barcoded.fasta

	rm ${directory}${sample}/${sample}.cutadapt.deduped.fasta


###Identify sncRNA (BLAST)###
			#Describe how to make reference list
			#Describe how to convert to db

#Run BLAST
#For further usage details: https://www.ncbi.nlm.nih.gov/books/NBK279690/pdf/Bookshelf_NBK279690.pdf

	blastn \
	-db ${reference_directory}sncRNA.fasta \
	-query ${directory}${sample}/${sample}.cutadapt.deduped.barcoded.fasta \
	-out ${directory}${sample}/${sample}.cutadapt.deduped.barcoded.blast \
	-word_size 11 \
	-outfmt 6 \
	-num_threads 1 \
	-strand plus

#Filter reads
#length ($4) >= 14
#evalue ($11) < 0.05
#mismatch ($5) <= 2
#gapopen ($6) <= 1
#if mismatch ($5) = 2 , gapopen ($6) = 0

#Consider turning these all into variables
#minimum_length=14
#minimum_evalue=.05
#maximum_mismatch=2
#maximum_gap=1
#maximum_gap_if_maximum_mismatch=

	awk '($11 + 0) <= .05' ${directory}${sample}/${sample}.cutadapt.deduped.barcoded.blast | \
	awk '$4 >= 14' | \
	awk '$5 <= 2' | \
	awk '$6 <= 1' | \
	awk '!($5 == 2 && $6 == 1)' | \
	awk '!x[$1]++' \
	> ${directory}${sample}/${sample}.cutadapt.deduped.barcoded.blast.filtered

	rm ${directory}${sample}/${sample}.cutadapt.deduped.barcoded.blast

#Make tabular version of barcoded FASTA file (before BLAST) using FASTA_Formatter

	awk 'BEGIN{RS=">";OFS="\t"}NR>1{print $1,$2}' \
	${directory}${sample}/${sample}.cutadapt.deduped.barcoded.fasta \
	> ${directory}${sample}/${sample}.cutadapt.deduped.barcoded.tab

	rm ${directory}${sample}/${sample}.cutadapt.deduped.barcoded.fasta

#Join tabular input file with filtered BLAST results
#It is very important to include the -k1,1 tag to the sort command, I have found that sort alone sometimes produces different arrangements for different files

	join \
	<(sort -k1,1 ${directory}${sample}/${sample}.cutadapt.deduped.barcoded.blast.filtered -T ${directory}${sample}) \
	<(sort -k1,1 ${directory}${sample}/${sample}.cutadapt.deduped.barcoded.tab -T ${directory}${sample}) \
	> ${directory}${sample}/${sample}.blast.merged

	rm ${directory}${sample}/${sample}.cutadapt.deduped.barcoded.blast.filtered
	rm ${directory}${sample}/${sample}.cutadapt.deduped.barcoded.tab

##########
	echo "${sample} BLAST Aligned Reads: $(( $(wc -l ${directory}${sample}/${sample}.blast.merged | awk '{print $1}') / 1 ))" >> ${directory}${sample}/${sample}.summary.txt
##########


###Identify sncRNA-first chimeras###

#Create FASTA file with sequence of sncRNA remove and name appended to the read ID
#sncRNA start position: 1
#read length after sncRNA: >=15

	awk '$7 == 1' ${directory}${sample}/${sample}.blast.merged | \
	awk '{stop=$8;readlength=length($13);print $0,"\t", readlength-(stop)}' | \
	awk '$14 >= 15' | \
	awk '{start=$8; fulllength=length($13);print $0,"\t",substr($13, start+1, fulllength-(start));}' | \
	awk '{print">"$1"."$2"\n"$15}' \
	> ${directory}${sample}/${sample}.target.fasta

	rm ${directory}${sample}/${sample}.blast.merged

##########
	echo "${sample} sncRNA-first Chimeras: $(( $(wc -l ${directory}${sample}/${sample}.target.fasta | awk '{print $1}') / 2 ))" >> ${directory}${sample}/${sample}.summary.txt
##########


###Align reads to the mouse genome (HISAT2)###
#For further usage details: http://daehwankimlab.github.io/hisat2/manual/

#Run HISAT2

	hisat2 \
	-x /media/sf_UbuntuSharing/mouse/mouse \
	-f ${directory}${sample}/${sample}.target.fasta \
	-S ${directory}${sample}/${sample}.aligned.sam \
	--summary-file ${directory}${sample}/${sample}.hisat2summary.txt

	hisat2 \
	-x /media/sf_UbuntuSharing/mouse/mouse \
	-f ${directory}${sample}/${sample}.target.fasta \
	-S ${directory}${sample}/${sample}.aligned.max50000intron.sam \
	--summary-file ${directory}${sample}/${sample}.hisat2summary.txt \
	--max-intronlen 50000

##########
	echo "${sample} Unique Aligned Reads: $(( $(awk '$5 == 60' ${directory}${sample}/${sample}.aligned.sam | wc -l | awk '{print $1}') / 1 ))" >> ${directory}${sample}/${sample}.summary.txt
	echo "${sample} Unique Aligned Reads (max 50000 intron): $(( $(awk '$5 == 60' ${directory}${sample}/${sample}.aligned.max50000intron.sam | wc -l | awk '{print $1}') / 1 ))" >> ${directory}${sample}/${sample}.summary.txt	
##########


#Select only uniquely mapped reads

	awk '/^@/ || $5 == 60' ${directory}${sample}/${sample}.aligned.sam \
	> ${directory}${sample}/${sample}.aligned.unique.sam

#Convert to bam file

	samtools view \
	-S \
	-h \
	-b ${directory}${sample}/${sample}.aligned.unique.sam \
	> ${directory}${sample}/${sample}.aligned.unique.bam

echo "Finish:" >> ${directory}${sample}/${sample}.summary.txt
date >> ${directory}${sample}/${sample}.summary.txt

done

