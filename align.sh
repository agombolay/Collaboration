#!/usr/bin/env bash
#Author: Alli Gombolay

program=$0

function usage () {
        echo "Usage: $program [-i] '/path/to/file1.fastq etc.' [-b] 'Bowtie index' [-d] 'Ribose-seq directory' [-h]
          -i Filepaths of input FASTQ files 
          -b Basename of Bowtie index to be searched
          -d Location to save local Ribose-seq directory"
}

while getopts "i:b:d:h" opt;
do
    case $opt in
        i ) fastq=($OPTARG) ;;
	b ) index=$OPTARG ;;
	d ) directory=$OPTARG ;;
        h ) usage ;;
    esac
done

if [ "$1" == "-h" ];
then
        exit
fi

for samples in ${fastq[@]};
do
	
	filename=$(basename "$fastq")
	samples="${filename%.*}"
	
	inputDirectory=$(dirname "${fastq}")
	
	UMI=NNNNNNNN

	input=$inputDirectory/$samples.fastq

	output=$directory/ribose-seq/results/$samples/alignment

	if [[ ! -d $output ]];
	then
    		mkdir -p $output
	fi

	umiTrimmed=$output/$samples.umiTrimmed.fastq.gz

	intermediateSAM=$output/$samples.intermediate.sam
	intermediateBAM=$output/$samples.intermediate.bam

	sortedBAM=$output/$samples.sorted.bam

	finalBAM=$output/$samples.bam

	statistics=$output/$samples.statistics.txt

	BED=$output/$samples.bed.gz

	umitools.py trim $input $UMI | gzip -c > $umiTrimmed

	zcat $umiTrimmed | bowtie -m 1 --sam $index - 2> $statistics 1> $intermediateSAM

	samtools view -ShuF4 $intermediateSAM > $intermediateBAM

	samtools sort $intermediateBAM > $sortedBAM

	umitools.py rmdup $sortedBAM $finalBAM | gzip -c > $BED

done
