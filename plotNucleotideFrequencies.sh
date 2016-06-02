#! /usr/bin/env bash
#Author: Alli Gombolay

program=$0

function usage () {
        echo "Usage: $program [-i] 'sample1 etc.' [-d] 'Ribose-seq directory' [-h]
        -i Sample names of input BAM files (i.e, sample1 for sample1.bam)
        -d Location to save local Ribose-seq directory"
}

while getopts "i:d:h" opt;
do
    case $opt in
    	i ) samples=($OPTARG) ;;
	d ) directory=$OPTARG ;;
    	h ) usage ;;
    esac
done

if [ "$1" == "-h" ];
then
        exit
fi

offset_values=(100 50 15)

modes=("all" "only-mitochondria" "no-mitochondria" "only-2micron")

input=$directory/ribose-seq/results/$samples/nucleotideFrequencies/

output="$directory/ribose-seq/results/$samples/plots/nucleotideFrequencies"

if [[ ! -d $output ]];
then
	mkdir -p $output
fi

for index in ${!modes[@]};
do

	mode=${modes[$index]}

	for value in ${offset_values[@]};
	do
    
		sampleID="$samples.subset-$mode"
    
		tables="$input/$samples.$mode.nucleotideFrequencies.tab"
    
		Rscript plotNucleotideFrequencies.R -n "$sampleID" -d $output --offsetmax $value $tables
    
	done

done
