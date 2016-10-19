#!/bin/bash
#$ -cwd
#$ -V
#$ -N radMap
#$ -pe shared 8
#$ -l highp,h_data=4G,time=100:00:00
#$ -M rachaelbay@gmail.com
#$ -m bea

. /u/local/Modules/default/init/modules.sh
module load samtools
module load bowtie2

PLATE="WIFL1"

###Make directory for raw bam files
mkdir bam
cd dupfiltered


###Align each sample to genome. Note that genome reference must already be built through bowtie2-build
for sample in `ls *.1.1.fq.gz | cut -f1 -d'.'`
do
	ID="$PLATE.$sample"
	bowtie2 -x ~/nobackup-klohmuel/References/WIFL/WIFL \
	--threads 8 -1 $sample.1.1.fq.gz -2 $sample.2.2.fq.gz \
	--rg-id $ID --rg SM:$sample --rg LB:$PLATE --rg PU:$PLATE --rg PL:illumina | \
	samtools view -bhS - | \
	samtools sort - ../bam/$ID

done