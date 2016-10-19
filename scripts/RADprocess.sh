#!/bin/bash
#$ -cwd
#$ -V
#$ -N rad
#$ -pe shared 1
#$ -l highp,h_data=8G,time=100:00:00
#$ -M rachaelbay@gmail.com
#$ -m bea

. /u/local/Modules/default/init/modules.sh
module load samtools

###Uncompress the raw fastqs
gunzip WIFL-Plate1_S28_L003_R1_001.fastq.gz
gunzip WIFL-Plate1_S28_L003_R2_001.fastq.gz

###For bestRAD only: flip read pairs that are reversed and trim extra "GG"
perl ~/scripts/flip_trim_werrors.pl ~/scripts/barcodes.txt WIFL-Plate1_S28_L003_R1_001.fastq WIFL-Plate1_S28_L003_R2_001.fastq \
	WIFL_Plate1_R1.fq WIFL_Plate1_R2.fq true 2

###Recompress fastq. This isn't strictly necessary, but I did it to save space
gzip WIFL_Plate1_R1.fq
gzip WIFL_Plate1_R2.fq

###This calls the newer g++ compiler, necessary for running the newer stacks version on Hoffman2
source ~/scripts/mypath.sh

###Make directory for demultiplexed fastq files
mkdir demultiplexed

###Run stacks to demultiplex and trim adapter sequences
~/programs/stacks-1.39/process_radtags -i gzfastq -1 WIFL_Plate1_R1.fq.gz -2 WIFL_Plate1_R2.fq.gz \
	-o ./demultiplexed -b WIFL1barcodes.txt -c -q -r -e sbfI --filter_illumina \
	--adapter_1 GATCGGAAGAGCACACGTCTGAACTCCAGTC --adapter_2 CACTCTTTCCCTACACGACGCTCTTCCGATCT

###Move orphaned reads to their own directory
mkdir orphans
mv demultiplexed/*rem* orphans 

###Make directory for duplicate filtered reads
mkdir dupfiltered
cd demultiplexed

###Filter duplicate reads from demultiplexed fastq files, creating new files within the dupfiltered directory
for sample in `ls *.1.fq.gz | cut -f1 -d'.'`
do
~/programs/stacks-1.39/clone_filter -1 $sample.1.fq.gz -2 $sample.2.fq.gz -i gzfastq -o ../dupfiltered
done
