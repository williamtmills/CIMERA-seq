#! /bin/sh

#Configure for using Miniconda environment in bash script
#If you used the installation from '1_Pipeline_Setup_Tips.sh' and installed Miniconda in your home directory, you should just have to change 'william' to your username

source /home/william/miniconda3/etc/profile.d/conda.sh

conda activate PACeR

#Structure of files
#Directory (folder) for every sample
#Each sample directory contains the '.aligned.unique.bam' file obtained from '3_Pipeline.sh'

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
#        │    │       sample1.aligned.unique.bam
#        │    │
#        │    └───sample2
#        │    │       sample2.aligned.unique.bam
#        │    │
#        │    └───sample3
#        │            sample3.aligned.unique.bam
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


#Create an environment variable for the path to the reference genome

#Format
#input_genome_file="/path/to/genome/file"

#Example
input_genome_file="/media/sf_Ubuntu_Sharing_2022/mouse/mm10.fa"


#Create an environment variable containing a list of all samples
#One sample per line

samples="
CS3
CS7
W8FCS1
W8FCS6
"

#Create an environment variable for the minimum number of reads (combined from all samples) necessary to call a peak

minimum_reads=3

#Create an environment variable for the minimum number of libraries that must contain a read supporting a peak for it to be called

minimum_libraries=2

#Create an environment variable for the number of cycles of peak calling to perform

cycles=3


#For the most part, code modifications are not required below this line
##############################################################################

#Delete any temporary files

	[ -e ${directory}combined.rearranged.bed ] && rm ${directory}combined.rearranged.bed
	[ -e ${directory}PACeR.peaks.bed ] && rm ${directory}PACeR.peaks.bed

mkdir ${directory}PACeR_Peak_Calling

for sample in ${samples}
do

	gatk SplitNCigarReads \
	-R ${input_genome_file} \
	-I ${directory}${sample}/${sample}.aligned.unique.bam \
	-O ${directory}PACeR_Peak_Calling/${sample}.split.bam

	bedtools \
	bamtobed \
	-i ${directory}PACeR_Peak_Calling/${sample}.split.bam \
	> ${directory}PACeR_Peak_Calling/${sample}.split.bed


#Rearrange to place sncRNA name next to chromosome number to allow simultaneous peak calling

	sed 's/\./\t/' ${directory}PACeR_Peak_Calling/${sample}.split.bed | \
	awk '{print $1"."$5"\t"$2"\t"$3"\t""X""\t""X""\t"$7}' \
	> ${directory}PACeR_Peak_Calling/${sample}.split.rearranged.bed


#Combine bed files from each sample into a single file

	cat "${directory}PACeR_Peak_Calling/${sample}.split.rearranged.bed" >> ${directory}PACeR_Peak_Calling/combined.rearranged.bed

done


#Sort combined bed file for easier processing

	sort -k1,1 -k2,2n ${directory}PACeR_Peak_Calling/combined.rearranged.bed > ${directory}PACeR_Peak_Calling/combined.rearranged.sorted.bed


#Create a temporary copy of the combined bed file that will be modified during successive rounds of peak calling

	cp ${directory}PACeR_Peak_Calling/combined.rearranged.sorted.bed ${directory}PACeR_Peak_Calling/combined.rearranged.sorted.tmp.bed 


#Perform cycles of peak calling

for i in $(seq $cycles)
do

#Merge overlapping reads

	bedtools \
	merge \
	-i ${directory}PACeR_Peak_Calling/combined.rearranged.sorted.tmp.bed \
	-s \
	-c 6 \
	-o distinct \
	> ${directory}PACeR_Peak_Calling/combined.rearranged.sorted.merged


#Determine coverage across overlapping regions
#Remove regions with fewer than minimum read threshold (assigned at beginning of script)
#Remaining regions termed "windows"

	bedtools \
	coverage \
	-a <(awk '{print $1"\t"$2"\t"$3"\t""X""\t""X""\t"$4}' ${directory}PACeR_Peak_Calling/combined.rearranged.sorted.merged) \
	-b ${directory}PACeR_Peak_Calling/combined.rearranged.sorted.tmp.bed \
	-s \
	-d | \
	awk -v var=${minimum_reads} '$8 >= var' \
	> ${directory}PACeR_Peak_Calling/combined.rearranged.sorted.merged.coverage


#Find start and stop position of maximal coverage within window

	paste \
	<(sort -k1,1 -k2,2n -k8,8nr -k7,7n ${directory}PACeR_Peak_Calling/combined.rearranged.sorted.merged.coverage | awk '!window[$1, $2, $3, $6]++' | cut -f 1-3,6,7) \
	<(sort -k1,1 -k2,2n -k8,8nr -k7,7nr ${directory}PACeR_Peak_Calling/combined.rearranged.sorted.merged.coverage | awk '!window[$1, $2, $3, $6]++' | cut -f 7) | \
	awk '{print $1"\t"($2+$5-1)"\t"($2+$6)"\t""X""\t""X""\t"$4}' \
	> ${directory}PACeR_Peak_Calling/combined.rearranged.sorted.merged.coverage.filtered

	cat ${directory}PACeR_Peak_Calling/combined.rearranged.sorted.merged.coverage.filtered >> ${directory}PACeR_Peak_Calling/PACeR.peak.bed


#Remove reads that overlap with peak

	bedtools \
	subtract \
	-a ${directory}PACeR_Peak_Calling/combined.rearranged.sorted.tmp.bed \
	-b ${directory}PACeR_Peak_Calling/combined.rearranged.sorted.merged.coverage.filtered \
	-s \
	-A \
	> ${directory}PACeR_Peak_Calling/combined.rearranged.sorted.merged.coverage.filtered.outside


#Replace temporary bed file with bed file that reads were removed from

	mv ${directory}PACeR_Peak_Calling/combined.rearranged.sorted.merged.coverage.filtered.outside ${directory}PACeR_Peak_Calling/combined.rearranged.sorted.tmp.bed

done


#Sort bed file for easier processing

	sort -k1,1 -k2,2n ${directory}PACeR_Peak_Calling/PACeR.peak.bed > ${directory}PACeR_Peak_Calling/PACeR.peak.sorted.bed


#Determine coverage for each individual sample

for sample in ${samples}
do

	bedtools \
	coverage \
	-a ${directory}PACeR_Peak_Calling/PACeR.peak.sorted.bed \
	-b <(sort -k1,1 -k2,2n ${directory}PACeR_Peak_Calling/${sample}.split.rearranged.bed) \
	-s \
	> ${directory}PACeR_Peak_Calling/${sample}.split.rearranged.coverage.bed

done


#Combine coverages for each sample into a single file
#Each column represents a sample

	first=$(echo "${samples}" | sed -r '/^\s*$/d' | awk 'NR==1')

	[ -e ${directory}tmp.tmp ] && rm ${directory}tmp.tmp

	cut -f 1-6 ${directory}PACeR_Peak_Calling/${sample}.split.rearranged.coverage.bed > ${directory}tmp.tmp

for sample in ${samples}
do

	cut \
	-f 7 \
	${directory}PACeR_Peak_Calling/${sample}.split.rearranged.coverage.bed | \
	paste ${directory}tmp.tmp - \
	> ${directory}tmp2.tmp

	mv ${directory}tmp2.tmp ${directory}tmp.tmp

done

	mv ${directory}tmp.tmp ${directory}PACeR_Peak_Calling/PACeR.peak.sorted.coverage.bed

#Remove peaks supported by fewer than minimum number of libraries (assigned at beginning of script)
#Final peak file has information about parameters at the beginning of the bed file, if this is undesired, place a # infront of the following 3 lines

#	echo "#Minimum reads: $minimum_reads" >> ${directory}PACeR.peaks.bed
#	echo "#Cycles: $cycles" >> ${directory}PACeR.peaks.bed
#	echo "#Minimum libraries: $minimum_libraries" >> ${directory}PACeR.peaks.bed

	awk \
	-v var=${minimum_libraries} \
	-v var2=$((6 + $(echo "$samples" | sed -r '/^\s*$/d' | wc -l) )) \
	'{nz=0; for(i=7;i<=var2;i++) nz+=($i!=0)} nz>=var' \
	${directory}PACeR_Peak_Calling/PACeR.peak.sorted.coverage.bed \
	>> ${directory}PACeR_Peak_Calling/PACeR.peaks.bed

	bedtools \
	coverage \
	-a <(awk '!/#/' ${directory}PACeR_Peak_Calling/PACeR.peaks.bed | awk '{print $1"\t"$2"\t"$3"\t""X""\t""X""\t"$6}' | sort -k1,1 -k2,2n) \
	-b <(awk '{print $1"\t"$2"\t"$3"\t""X""\t""X""\t"$6}' ${directory}PACeR_Peak_Calling/combined.rearranged.sorted.bed | sort -k1,1 -k2,2n) \
	-s | \
	sed 's/\./\t/' | \
	awk '{print $1"\t"$3"\t"$4"\t"$2"\t"$8"\t"$7}' \
	> ${directory}PACeR.peaks.bed

#Intermediate files are stored in folder called 'PACeR_Peak_Calling' which are deleted after peaks are called. If you would like to keep the intermediate files you can place a # in front of the next line

	rm -r ${directory}PACeR_Peak_Calling

#Output bed format
#Column1: chromosome
#Column2: start
#Column3: stop
#Column4: sncRNA
#Column5: Number of reads supporting peak (combined from all sample)
#Column6: Strand

conda deactivate
