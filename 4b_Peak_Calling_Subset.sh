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


#Create an environment variable with the path to where reference files will be stored (for example, in a folder named after the model organism)
#This is where the motif files were saved in '2_Reference_File_Configuration.sh'

#Format
#reference_directory="/path/to/reference/directory/"

#Example
reference_directory="/media/sf_Ubuntu_Sharing_2022/mouse/"


#Create an environment variable for the path to the reference genome

#Format
#input_genome_file="/path/to/genome/file"

#Example
input_genome_file="/media/sf_Ubuntu_Sharing_2022/mouse/mm10.fa"


#Create an environment variable for the path to the reference genome lengths file

#Format
#genome_lengths="/path/to/genome/lengths"

#Example
genome_lengths="/media/sf_Ubuntu_Sharing_2022/mouse/mm10.chrom.sizes"


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


#Create an environment variable for how the output files should be named
#Name cannot not include spaces

name="let-7-family"


#Create an environment variable for the length of the motif to be used for enrichment analysis
#Number should be spelled out (Only 'six' or 'seven' currently supported)

xmer="six"


#Create an environment variable for how far peaks should be extended from their center for motif enrichment analysis
#Total width of peak for motif enrichment analysis will be twice the value of the environment variable

extension=10


#IF SUBSETTING FOR A SINGLE sncRNA: Create an environment variable with the exact name of the sncRNA (from miRbase (https://www.mirbase.org/ftp.shtml) or tRFdb (http://genome.bioch.virginia.edu/trfdb/index.php)
#IF SUBSETTING FOR A FAMILY OF sncRNAS: Code below must be modified directly (ignore this variable by placing a # in front)

#subset="miR-138-5p"


#For the most part, code modifications are not required below this line
##############################################################################

#Create an environment variable to the sixmer and sevenmer files obtained in '2_Reference_File_Configuration.sh'

sixmer_file="${reference_directory}sixmer.txt"
sevenmer_file="${reference_directory}sevenmer.txt"

sixmers=$(cat $sixmer_file)
sevenmers=$(cat $sevenmer_file)

	mkdir ${directory}${name}_PeakCalling


#Delete and temporary files

	[ -e ${directory}${name}_PeakCalling/combined.rearranged.bed ] && rm ${directory}${name}_PeakCalling/combined.rearranged.bed
	[ -e ${directory}${name}_PeakCalling/PACeR.${name}.peaks.bed ] && rm ${directory}${name}_PeakCalling/PACeR.${name}.peaks.bed
	[ -e ${directory}${name}_PeakCalling/PACeR.peak.bed ] && rm ${directory}${name}_PeakCalling/PACeR.peak.bed
	[ -e ${directory}${name}_PeakCalling/${name}.${xmer}mer.count ] && rm ${directory}${name}_PeakCalling/${name}.${xmer}mer.count

for sample in ${samples}
do

#	samtools view -h \
#	${directory}${sample}/${sample}.aligned.unique.bam | \
#	awk -v subset=${subset} '/^@/ || $0~subset' | \
#	samtools view -S -h -b - | \
#	samtools sort - \
#	-o ${directory}${name}_PeakCalling/${sample}.${name}.aligned.unique.bam


#Example for subsetting the -5p members of the let-7 family

	samtools view -h \
	${directory}${sample}/${sample}.aligned.unique.bam | \
	awk '/^@/ || /let-7.-5p/ || /miR-98-5p/' | \
	samtools view -S -h -b - | \
	samtools sort - \
	-o ${directory}${name}_PeakCalling/${sample}.${name}.aligned.unique.bam

#Example for subsetting the miR-29 family

#	samtools view -h \
#	${directory}${sample}/${sample}.aligned.unique.bam | \
#	awk '/^@/ || /miR-29a-/ || /miR-29b-/ || /miR-29c-/' | \
#	samtools view -S -h -b - | \
#	samtools sort - \
#	-o ${directory}${name}_PeakCalling/${sample}.${name}.aligned.unique.bam

	gatk SplitNCigarReads \
	-R ${input_genome_file} \
	-I ${directory}${name}_PeakCalling/${sample}.${name}.aligned.unique.bam \
	-O ${directory}${name}_PeakCalling/${sample}.split.bam

	bedtools \
	bamtobed \
	-i ${directory}${name}_PeakCalling/${sample}.split.bam \
	> ${directory}${name}_PeakCalling/${sample}.split.bed


#Rearrange to place sncRNA name next to chromosome number

	sed 's/\./\t/' ${directory}${name}_PeakCalling/${sample}.split.bed | \
	awk -v name="${name}" '{print $1"."name"\t"$2"\t"$3"\t""X""\t""X""\t"$7}' \
	> ${directory}${name}_PeakCalling/${sample}.split.rearranged.bed


#Combine bed files

	cat "${directory}${name}_PeakCalling/${sample}.split.rearranged.bed" >> ${directory}${name}_PeakCalling/combined.rearranged.bed

done


#Sort combined bed file for easier processing

	sort -k1,1 -k2,2n ${directory}${name}_PeakCalling/combined.rearranged.bed > ${directory}${name}_PeakCalling/combined.rearranged.sorted.bed


#Create a temporary copy of the combined bed file that will be modified during successive rounds of peak calling

	cp ${directory}${name}_PeakCalling/combined.rearranged.sorted.bed ${directory}${name}_PeakCalling/combined.rearranged.sorted.tmp.bed 


#Perform cycles of peak calling

#for cycles in $cycle_range
#do

for i in $(seq $cycles)
do

#Merge overlapping reads

	bedtools \
	merge \
	-i ${directory}${name}_PeakCalling/combined.rearranged.sorted.tmp.bed \
	-s \
	-c 6 \
	-o distinct \
	> ${directory}${name}_PeakCalling/combined.rearranged.sorted.merged


#Determine coverage across overlapping regions
#Remove regions with fewer than minimum read threshold (assigned at beginning of script)
#Remaining regions termed "windows"

	bedtools \
	coverage \
	-a <(awk '{print $1"\t"$2"\t"$3"\t""X""\t""X""\t"$4}' ${directory}${name}_PeakCalling/combined.rearranged.sorted.merged) \
	-b ${directory}${name}_PeakCalling/combined.rearranged.sorted.tmp.bed \
	-s \
	-d | \
	awk -v var=${minimum_reads} '$8 >= var' \
	> ${directory}${name}_PeakCalling/combined.rearranged.sorted.merged.coverage


#Find start and stop position of maximal coverage within window

	paste \
	<(sort -k1,1 -k2,2n -k8,8nr -k7,7n ${directory}${name}_PeakCalling/combined.rearranged.sorted.merged.coverage | awk '!window[$1, $2, $3, $6]++' | cut -f 1-3,6,7) \
	<(sort -k1,1 -k2,2n -k8,8nr -k7,7nr ${directory}${name}_PeakCalling/combined.rearranged.sorted.merged.coverage | awk '!window[$1, $2, $3, $6]++' | cut -f 7) | \
	awk '{print $1"\t"($2+$5-1)"\t"($2+$6)"\t""X""\t""X""\t"$4}' \
	> ${directory}${name}_PeakCalling/combined.rearranged.sorted.merged.coverage.filtered

	cat ${directory}${name}_PeakCalling/combined.rearranged.sorted.merged.coverage.filtered >> ${directory}${name}_PeakCalling/PACeR.peak.bed


#Remove reads that overlap with peak

	bedtools \
	subtract \
	-a ${directory}${name}_PeakCalling/combined.rearranged.sorted.tmp.bed \
	-b ${directory}${name}_PeakCalling/combined.rearranged.sorted.merged.coverage.filtered \
	-s \
	-A \
	> ${directory}${name}_PeakCalling/combined.rearranged.sorted.merged.coverage.filtered.outside


#Replace temporary bed file with bed file that reads were removed from

	mv ${directory}${name}_PeakCalling/combined.rearranged.sorted.merged.coverage.filtered.outside ${directory}${name}_PeakCalling/combined.rearranged.sorted.tmp.bed

done


#Sort bed file for easier processing

	sort -k1,1 -k2,2n ${directory}${name}_PeakCalling/PACeR.peak.bed > ${directory}${name}_PeakCalling/PACeR.peak.sorted.bed

#Determine coverage for each individual sample

for sample in ${samples}
do

	bedtools \
	coverage \
	-a ${directory}${name}_PeakCalling/PACeR.peak.sorted.bed \
	-b <(sort -k1,1 -k2,2n ${directory}${name}_PeakCalling/${sample}.split.rearranged.bed) \
	-s \
	> ${directory}${name}_PeakCalling/${sample}.split.rearranged.coverage.bed

done


#Combine coverages for each sample into a single file
#Each column represents a sample

	first=$(echo "${samples}" | sed -r '/^\s*$/d' | awk 'NR==1')

	[ -e ${directory}${name}_PeakCalling/tmp.tmp ] && rm ${directory}${name}_PeakCalling/tmp.tmp

	cut -f 1-6 ${directory}${name}_PeakCalling/${sample}.split.rearranged.coverage.bed > ${directory}${name}_PeakCalling/tmp.tmp

for sample in ${samples}
do

	cut \
	-f 7 \
	${directory}${name}_PeakCalling/${sample}.split.rearranged.coverage.bed | \
	paste ${directory}${name}_PeakCalling/tmp.tmp - \
	> ${directory}${name}_PeakCalling/tmp2.tmp

	mv ${directory}${name}_PeakCalling/tmp2.tmp ${directory}${name}_PeakCalling/tmp.tmp

done

	mv ${directory}${name}_PeakCalling/tmp.tmp ${directory}${name}_PeakCalling/PACeR.${name}.peak.sorted.coverage.bed


#Remove peaks supported by fewer than minimum number of libraries (assigned at beginning of script)
#Final peak file has information about parameters at the beginning of the bed file, if this is undesired, place a # infront of the following 4 lines

#	echo "#$name" >> ${directory}${name}_PeakCalling/PACeR.${name}.peaks.bed
#	echo "#Minimum reads: $minimum_reads" >> ${directory}${name}_PeakCalling/PACeR.${name}.peaks.bed
#	echo "#Cycles: $cycles" >> ${directory}${name}_PeakCalling/PACeR.${name}.peaks.bed
#	echo "#Minimum libraries: $minimum_libraries" >> ${directory}${name}_PeakCalling/PACeR.${name}.peaks.bed

	awk \
	-v var=${minimum_libraries} \
	-v var2=$((6 + $(echo "$samples" | sed -r '/^\s*$/d' | wc -l) )) \
	'{nz=0; for(i=7;i<=var2;i++) nz+=($i!=0)} nz>=var' \
	${directory}${name}_PeakCalling/PACeR.${name}.peak.sorted.coverage.bed \
	>> ${directory}${name}_PeakCalling/PACeR.${name}.peaks.bed

	bedtools \
	coverage \
	-a <(awk '!/#/' ${directory}${name}_PeakCalling/PACeR.${name}.peaks.bed | awk '{print $1"\t"$2"\t"$3"\t""X""\t""X""\t"$6}' | sort -k1,1 -k2,2n) \
	-b <(awk '{print $1"\t"$2"\t"$3"\t""X""\t""X""\t"$6}' ${directory}${name}_PeakCalling/combined.rearranged.sorted.bed | sort -k1,1 -k2,2n) \
	-s | \
	sed 's/\./\t/' | \
	awk '{print $1"\t"$3"\t"$4"\t"$2"\t"$8"\t"$7}' \
	> ${directory}${name}_PeakCalling/PACeR.${name}.peaks.final.bed

#Output bed format
#Column1: chromosome
#Column2: start
#Column3: stop
#Column4: Name assigned at beginning of script
#Column5: Number of reads supporting peak (combined from all sample)
#Column6: Strand

#Determine the center of each peak

	awk 'BEGIN { OFS = "\t" ; OFMT="%.0f" } {mid=(int($2)+int($3))/2 ; print($1,mid,mid,$4,$5,$6);}' \
	${directory}${name}_PeakCalling/PACeR.${name}.peaks.final.bed \
	> ${directory}${name}_PeakCalling/PACeR.${name}.peaks.final.bed.center

#Extend the peak in each direction by the number of bases assigned at the beginning of the script

	bedtools \
	slop \
	-i ${directory}${name}_PeakCalling/PACeR.${name}.peaks.final.bed.center \
	-g ${genome_lengths} \
	-b ${extension} \
	> ${directory}${name}_PeakCalling/PACeR.${name}.peaks.final.extended.center.bed

#Obtain the sequence of nucleotides within each peak
#Change 'T' to 'U' to represent RNA

	bedtools \
	getfasta \
	-fi ${input_genome_file} \
	-bed ${directory}${name}_PeakCalling/PACeR.${name}.peaks.final.extended.center.bed \
	-s | \
	awk '{print toupper($0)}' | \
	sed 's/T/U/g' | \
	awk '!/>/' \
	> ${directory}${name}_PeakCalling/PACeR.${name}.peaks.final.extended.sequences

	sort -T ${directory} ${directory}${name}_PeakCalling/PACeR.${name}.peaks.final.extended.sequences \
	> ${directory}${name}_PeakCalling/PACeR.${name}.peaks.final.extended.sequences.sorted


#Perform motif enrichment analysis
#Output:
#Column1: Motif
#Column2: Numer of peaks containing motif

for xmer_search in $(eval echo "\$${xmer}mers")
do

	echo "$(( $(awk -v motif="${xmer_search}" '$0~motif' ${directory}${name}_PeakCalling/PACeR.${name}.peaks.final.extended.sequences.sorted | wc -l) ))" >> ${directory}${name}_PeakCalling/${name}.${xmer}mer.count

done

	paste -d "\t" $(eval echo "\$${xmer}mer_file") ${directory}${name}_PeakCalling/${name}.${xmer}mer.count | sort -T ${directory} -k2,2nr > ${directory}${name}_PeakCalling/${name}.${xmer}mer.count.sorted

	mv ${directory}${name}_PeakCalling/${name}.${xmer}mer.count.sorted ${directory}PACeR.${name}.${xmer}mer.count.sorted

	mv ${directory}${name}_PeakCalling/PACeR.${name}.peaks.final.bed ${directory}PACeR.${name}.peaks.final.bed


#Intermediate files are stored in folder called '${name}_PeakCalling' which are deleted after peaks are called. If you would like to keep the intermediate files you can place a # in front of the next line

	rm -r ${directory}${name}_PeakCalling/

conda deactivate
