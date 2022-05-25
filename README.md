# PACeR: a bioinformatic pipeline for the analysis of chimeric RNA-seq data

| File Name                     | Description |
| ----------- | ----------- |
| 1_Pipeline_Setup_Tips.sh      | Instructions for installing VirtualBox and Ubuntu, creating a shared folder between your computer and VirtualBox, installing Miniconda, and setting up a Miniconda environment with the requisite packages from running PACeR.       |
| environment.yml   | File for configuring the Miniconda environment with the requisite packages from running PACeR.        |
| 2_Reference_File_Configuration.sh      | Script for downloading and configuring the sncRNA reference list and the reference genome.       |
| 3_Pipeline.sh   | Script for processing compressed FASTQ files to generate a BAM file containing uniquely aligned reads annotated with the corresponding sncRNA within the read.        |
| 3\*\_Pipeline_for_CLEAR-CLIP.sh   | Script for processing compressed FASTQ files from a previously published dataset ([Moore et al. 2015](https://www.nature.com/articles/ncomms9864)) to generate a BAM file containing uniquely aligned reads annotated with the corresponding sncRNA within the read.        |
| 4a_Peak_Calling_Total.sh   | Script for calling peaks from BAM files containing uniquely aligned reads annotated with the corresponding sncRNA within the read.        |
| 4b_Peak_Calling_Subset.sh   | Script for calling peaks from BAM files containing uniquely aligned reads annotated with the corresponding sncRNA within the read. Subsets reads by individual sncRNAs or sncRNA families and includes motif enrichment analysis.        |
