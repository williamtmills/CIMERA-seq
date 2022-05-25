#! /bin/sh

#Create an environment variable with the path to where miniconda is installed

miniconda_path="/home/william/miniconda3/"

#Create an environment variable with the path to where reference files will be stored (for example, in a folder named after the model organism)

reference_directory="/media/sf_Ubuntu_Sharing_2022/mouse/"

#Create an environment variable for the prefix of the reference genome

reference_genome_prefix="mouse"

#For the most part, code modifications are not required below this line
##############################################################################

conda activate PACeR

source ${miniconda_path}etc/profile.d/conda.sh

#Obtain miRNAs from miRbase
#As the web link sometimes changes, the 'mature.fa' file can be found at https://www.mirbase.org/ftp.shtml

	wget https://www.mirbase.org/ftp/CURRENT/mature.fa.gz -O ${reference_directory}miRbase.fasta.gz
	gunzip ${reference_directory}miRbase.fasta.gz


#This portion is written specifically for obtaining mouse miRNAs
#If working with a different model organism, the grep command can be modified to select the 3-letter abbreviation for that model organism
#"miRNA-" is added to the beginning of each miRNA name to differentiate it from other sncRNA sequences
#Final miRNA name format (example): miRNA-mmu-let-7c-5p

	awk '/^>/ {printf "%s%s ", pfx, $0; pfx="\n"; next} {printf "%s", $0} END {print ""}' ${reference_directory}miRbase.fasta | \
	grep "mmu-" | \
	awk '{print $1"\n"$6}' | \
	sed 's/>/>miRNA-/g' \
	> ${reference_directory}miRbase.mmu.fasta


#Obtain tRNA fragments (tRFs) from tRFdb (http://genome.bioch.virginia.edu/trfdb/index.php)
##Unfortunately these have to be downloaded manually via a web browser (http://genome.bioch.virginia.edu/trfdb/search.php)
###Download the tRF-1, tRF-3, and tRF-5

#Example
#	Select Organism: Mouse
#	Output format: HTML + Comma Separated Values

#For simplicity, the tRF sequences were accessed March 11, 2022 and uploaded to github and may be downloaded using the following code:

	wget https://raw.githubusercontent.com/williamtmills/PACeR/main/mouse_tRF-1.csv -O ${reference_directory}mouse_tRF-1.csv
	wget https://raw.githubusercontent.com/williamtmills/PACeR/main/mouse_tRF-3.csv -O ${reference_directory}mouse_tRF-3.csv
	wget https://raw.githubusercontent.com/williamtmills/PACeR/main/mouse_tRF-5.csv -O ${reference_directory}mouse_tRF-5.csv


#Deduplicate and convert .csv to fasta

	awk -F ', ' 'NR!=1 {print ">tRF-"$1"_"$9}' ${reference_directory}mouse_tRF-1.csv | \
	sort | \
	uniq | \
	sed 's/_/\n/' \
	> ${reference_directory}mouse_tRF-1.fasta

	awk -F ', ' 'NR!=1 {print ">tRF-"$1"_"$9}' ${reference_directory}mouse_tRF-3.csv | \
	sort | \
	uniq | \
	sed 's/_/\n/' \
	> ${reference_directory}mouse_tRF-3.fasta

	awk -F ', ' 'NR!=1 {print ">tRF-"$1"_"$9}' ${reference_directory}mouse_tRF-5.csv | \
	sort | \
	uniq | \
	sed 's/_/\n/' \
	> ${reference_directory}mouse_tRF-5.fasta


#Combine tRF files into a single file

	cat ${reference_directory}mouse_tRF-1.fasta ${reference_directory}mouse_tRF-3.fasta ${reference_directory}mouse_tRF-5.fasta \
	> ${reference_directory}tRFdb.mmu.fasta

#Combine miRNA and tRF files into a single file

	cat ${reference_directory}miRbase.mmu.fasta ${reference_directory}tRFdb.mmu.fasta > ${reference_directory}sncRNA.fasta


#Make BLAST database for guide RNA list

	makeblastdb \
	-in ${reference_directory}sncRNA.fasta \
	-dbtype nucl

#Download mouse genome (mm10)

	wget http://hgdownload.cse.ucsc.edu/goldenpath/mm10/bigZips/mm10.fa.gz -O ${reference_directory}mm10.fa.gz
	gunzip -c ${reference_directory}mm10.fa.gz > ${reference_directory}mm10.fa

#Build mouse genome (HISAT2)
#Expect ~2 hours

	hisat2-build ${reference_directory}mm10.fa ${reference_directory}${reference_genome_prefix}

#Index mouse genome

	gatk CreateSequenceDictionary -R ${reference_directory}mm10.fa
	samtools faidx ${reference_directory}mm10.fa

#Download reference genome chromosome sizes

	wget http://hgdownload.cse.ucsc.edu/goldenpath/mm10/bigZips/mm10.chrom.sizes -O ${reference_directory}mm10.chrom.sizes

#Download 6-mer and 7-mer motif files for motif enrichment analysis

	wget https://raw.githubusercontent.com/williamtmills/PACeR/main/sixmer.txt -O ${reference_directory}sixmer.txt
	wget https://raw.githubusercontent.com/williamtmills/PACeR/main/sevenmer.txt -O ${reference_directory}sevenmer.txt

conda deactivate
