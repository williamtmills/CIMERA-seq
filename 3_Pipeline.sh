#! /bin/sh

#Configure for using Miniconda environment in bash script
#If you used the installation from '1_Pipeline_Setup_Tips.sh' and installed Miniconda in your home directory, you should just have to change 'william' to your username

source /home/william/miniconda3/etc/profile.d/conda.sh

conda activate PACeR

#Structure of files
#Directory (folder) for every sample
#Each sample directory contains 2 FASTQ files, one for each end of paired-end sequencing (denoted with '_R1' and '_R2')
#Single end reads should just take the form 'sample.fastq.gz'

#    path
#    │   
#    │       
#    │
#    └───to
#        │   
#        │       
#        │
#        │───files
#        │    │   
#        │    │       
#        │    │
#        │    └───sample1
#        │    │       sample1_R1.fastq.gz
#        │    │       sample1_R2.fastq.gz
#        │    │
#        │    └───sample2
#        │    │       sample2_R1.fastq.gz
#        │    │       sample2_R2.fastq.gz
#        │    │
#        │    └───sample3
#        │            sample3_R1.fastq.gz
#        │            sample3_R2.fastq.gz
#        │  
#        │       
#        │
#        │───reference (e.g. mouse)


#Create an environment variable for the path to the directory where the sample directories are stored
#Make sure the path ends with a forward slash

#Format
#directory="/path/to/files/"

#Example
directory="/media/sf_Ubuntu_Sharing_2022/files/"

#Create an environment variable for the path to the directory where the reference files are stored

#Format
#reference_directory="/path/to/references/"

#Example
reference_directory="/media/sf_Ubuntu_Sharing_2022/mouse/"

#Create an environment variable for indexed reference genome prefix (from '2_Reference_File_Configuration.sh')

reference_genome_prefix="mouse"


#Create an environment variable containing a list of all samples
#One sample per line
#FASTQ file name suffix must be '.fastq.gz' not '.fq.gz'
#Assuming paired-end reads are in format "file1_R1.fastq" and "file1_R2.fastq"

#Format
#samples="
#sample1
#sample2
#sample3
#"

#Example
samples="
CS7
"



#Assign 5' and 3' barcode length

five_prime_barcode_length=4
three_prime_barcode_length=2


#Assign filters for BLAST results

minimum_evalue=.05
minimum_length=14
maximum_mismatch=2
maximum_gap=1
maximum_gap_if_maximum_mismatch=0

#Assign filters for sncRNA-first chimeras

maximum_sncRNA_start_position=1
minimum_length_after_sncRNA=15

#In this pipeline, intermediate files are removed to save space on the computer's memory
#If you wish to keep intermediate files (and have sufficient storage on your computer), consider removing the lines of code beginning with 'rms'

#File is generated in sample folder called "Sample".summary.txt that summarizes number of reads following each step
#Sections separated with ########## are code for counting reads

#For the most part, code modifications are not required below this line
##############################################################################

###Quality Control (FastQC)###

##for sample in $samples
##do

#Make output directory for FastQC reports

##	mkdir ${directory}FastQC_Reports

#Run FastQC
#For further usage details: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/

##for file in $(ls ${directory}${sample}/ | awk '/fastq/')
##do

##	fastqc \
##	-o ${directory}FastQC_Reports \
##	${directory}${sample}/${file}

##done

##done


###Quality Control (MultiQC)###

#Make output directory for MultiQC report

##	mkdir ${directory}MultiQC_Report

#Run MultiQC
#For further usage details: https://multiqc.info/docs/

##	multiqc ${directory}FastQC_Reports -o ${directory}MultiQC_Report


for sample in $samples
do

##########
##	[ -e ${directory}${sample}/${sample}.summary.txt ] && rm ${directory}${sample}/${sample}.summary.txt

##	echo "$sample" >> ${directory}${sample}/${sample}.summary.txt

##	echo "Start:" >> ${directory}${sample}/${sample}.summary.txt
##	date >> ${directory}${sample}/${sample}.summary.txt

##for file in $(ls ${directory}${sample}/ | awk '/fastq/')
##do

##	gunzip -c ${directory}${sample}/${file} > ${directory}${sample}/${file}.tmp
##	echo "${file} Reads: $(( $(wc -l ${directory}${sample}/${file}.tmp | awk '{print $1}') / 4 ))" >> ${directory}${sample}/${sample}.summary.txt
##	rm ${directory}${sample}/${file}.tmp

##done
##########	


###Combine Paired-End Reads (FLASh)###
#Step may be skipped if using single-end sequencing
#If step is skipped, adjust input file name in cutadapt step (currently input file ends with '.extendedFrags.fastq.gz')

#Make directory for FLASh outpout

##	mkdir ${directory}${sample}/${sample}_FLASh

#Run FLASh
#For further usage details: http://gensoft.pasteur.fr/docs/FLASH/1.2.11/flash

##	flash \
##	--allow-outies \
##	--output-directory=${directory}${sample}/${sample}_FLASh/ \
##	--output-prefix=${sample} \
##	--max-overlap=150 \
##	--min-overlap=6 \
##	--compress \
##	${directory}${sample}/${sample}_R1.fastq.gz \
##	${directory}${sample}/${sample}_R2.fastq.gz \
##	2>&1 | tee ${directory}${sample}/${sample}_FLASh/FLASh_${sample}.log

#The following line simply renames the output of FLASh to '.fastq.gz' so that the Cutadapt step doesn't need to be modified to accept single-end or FLASh'd paired-end reads

##	mv ${directory}${sample}/${sample}_FLASh/${sample}.extendedFrags.fastq.gz ${directory}${sample}/${sample}_FLASh/${sample}.fastq.gz


##########
	gunzip -k ${directory}${sample}/${sample}_FLASh/${sample}.fastq.gz
	echo "${sample} FLASh Reads: $(( $(wc -l ${directory}${sample}/${sample}_FLASh/${sample}.fastq | awk '{print $1}') / 4 ))" >> ${directory}${sample}/${sample}.summary.txt
	rm ${directory}${sample}/${sample}_FLASh/${sample}.fastq
##########


###Remove Adapters, Qualty Trimming, Minimum Length Threshold (Cutadapt)###

#Run Cutadapt
#For further usage details: https://cutadapt.readthedocs.io/en/stable/guide.html
#Sequence of adapters may need to be adjusted based on biochemical protocol used
#Single-end sequence should use unlinked adapters (-a 3'ADAPTERSEQUENCE -g 5'ADAPTERSEQUENCE)

	cutadapt \
	-g CTACAGTCCGACGATC...TGGAATTCTCGGGTGCCAAGG \
	-q 30 \
	-m 30 \
	-o ${directory}${sample}/${sample}.cutadapt.fastq.gz \
	${directory}${sample}/${sample}_FLASh/${sample}.fastq.gz

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

#	rm ${directory}${sample}/${sample}.cutadapt.fastq


##########
	echo "${sample} Deduplicated Reads: $(( $(wc -l ${directory}${sample}/${sample}.cutadapt.deduped.fasta | awk '{print $1}') / 2 ))" >> ${directory}${sample}/${sample}.summary.txt
##########


###Barcode Reads###

	awk 'BEGIN{RS=">";OFS="\t"}NR>1{print ">"$1,$2}' ${directory}${sample}/${sample}.cutadapt.deduped.fasta | \
	awk -v var1=$five_prime_barcode_length -v var2=$three_prime_barcode_length '{print $1"\t"$2"\t"substr($2,1,var1)""substr($2,length($2)-var2+1,var2)}' | \
	awk -v var1=$five_prime_barcode_length '{$2 = substr($2, var1+1); print }' | \
	awk -v var2=$three_prime_barcode_length '{$2 = substr($2, 1, length($2)-var2); print }' | \
	awk '{print $1"_"$3"\n"$2}' \
	> ${directory}${sample}/${sample}.cutadapt.deduped.barcoded.fasta

#	rm ${directory}${sample}/${sample}.cutadapt.deduped.fasta


###Identify sncRNA (BLAST)###

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

	awk -v var1=$minimum_evalue '($11 + 0) <= var1' ${directory}${sample}/${sample}.cutadapt.deduped.barcoded.blast | \
	awk -v var1=$minimum_length '$4 >= var1' | \
	awk -v var1=$maximum_mismatch '$5 <= var1' | \
	awk -v var1=$maximum_gap '$6 <= var1' | \
	awk -v var1=$maximum_mismatch -v var2=$maximum_gap_if_maximum_mismatch '!($5 == var1 && $6 > var2)' | \
	awk '!x[$1]++' \
	> ${directory}${sample}/${sample}.cutadapt.deduped.barcoded.blast.filtered

#	rm ${directory}${sample}/${sample}.cutadapt.deduped.barcoded.blast

#Make tabular version of barcoded FASTA file (from before BLAST)

	awk 'BEGIN{RS=">";OFS="\t"}NR>1{print $1,$2}' \
	${directory}${sample}/${sample}.cutadapt.deduped.barcoded.fasta \
	> ${directory}${sample}/${sample}.cutadapt.deduped.barcoded.tab

#	rm ${directory}${sample}/${sample}.cutadapt.deduped.barcoded.fasta

#Join tabular input file with filtered BLAST results
#It is very important to include the -k1,1 tag to the sort command, I have found that sort alone sometimes produces different arrangements for different files

	join \
	<(sort -k1,1 ${directory}${sample}/${sample}.cutadapt.deduped.barcoded.blast.filtered -T ${directory}${sample}) \
	<(sort -k1,1 ${directory}${sample}/${sample}.cutadapt.deduped.barcoded.tab -T ${directory}${sample}) \
	> ${directory}${sample}/${sample}.blast.merged

#	rm ${directory}${sample}/${sample}.cutadapt.deduped.barcoded.blast.filtered
#	rm ${directory}${sample}/${sample}.cutadapt.deduped.barcoded.tab


##########
	echo "${sample} BLAST Aligned Reads: $(( $(wc -l ${directory}${sample}/${sample}.blast.merged | awk '{print $1}') / 1 ))" >> ${directory}${sample}/${sample}.summary.txt
##########


###Identify sncRNA-first chimeras###

#Create FASTA file with sequence of sncRNA remove and name appended to the read ID

#Double check

	awk -v var1=$maximum_sncRNA_start_position '$7 <= var1' ${directory}${sample}/${sample}.blast.merged | \
	awk '{stop=$8;readlength=length($13);print $0,"\t", readlength-(stop)}' | \
	awk -v var1=$minimum_length_after_sncRNA '$14 >= var1' | \
	awk '{print $0,"\t",substr($13, $8+1, length($13)-($8))}' | \
	awk '{print">"$1"."$2"\n"$15}' \
	> ${directory}${sample}/${sample}.target.fasta

##	rm ${directory}${sample}/${sample}.blast.merged


##########
	echo "${sample} sncRNA-first Chimeras: $(( $(wc -l ${directory}${sample}/${sample}.target.fasta | awk '{print $1}') / 2 ))" >> ${directory}${sample}/${sample}.summary.txt
##########


###Align reads to the mouse genome (HISAT2)###
#For further usage details: http://daehwankimlab.github.io/hisat2/manual/

#Run HISAT2

	hisat2 \
	-x ${reference_directory}${reference_genome_prefix} \
	-f ${directory}${sample}/${sample}.target.fasta \
	-S ${directory}${sample}/${sample}.aligned.sam \
	--summary-file ${directory}${sample}/${sample}.hisat2summary.txt


##########
	echo "${sample} Unique Aligned Reads: $(( $(awk '$5 == 60' ${directory}${sample}/${sample}.aligned.sam | wc -l | awk '{print $1}') / 1 ))" >> ${directory}${sample}/${sample}.summary.txt
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

conda deactivate
