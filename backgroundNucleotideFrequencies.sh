#! /usr/bin/env bash
#Author: Alli Gombolay

program=$0

function usage () {
        echo "Usage: $program [-r] 'reference genome' [-d] 'Ribose-seq directory' [-h]
          -r Reference genome of interest (i.e., sacCer2)
          -d Location to save local Ribose-seq directory"
}

while getopts "r:d:h" opt;
do
    case $opt in
        r ) reference=$OPTARG ;;
	d ) directory=$OPTARG ;;
        h ) usage ;;
    esac
done

if [ "$1" == "-h" ];
then
        exit
fi

FASTA=$directory/ribose-seq/data/reference/$reference.fa

output="$directory/ribose-seq/results/backgroundNucleotideFrequencies"

if [[ ! -d $output ]];
then
  mkdir $output
fi

genome="$output/$reference.genome.nucleotide.frequencies.tab"
python2.7 backgroundNucleotideFrequencies.py $FASTA --region-size-minimum 1 --region-size-maximum 3 --verbose > $genome
