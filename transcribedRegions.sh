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

#VARIABLE SPECIFICATION

strands="+ -"

BED=$directory/ribose-seq/data/reference/$reference.bed

chromosomeSizes=$directory/ribose-seq/data/reference/$reference.chrom.sizes

output=$directory/ribose-seq/results/transcribedRegions

if [[ ! -d $output ]];
then
	mkdir -p $output
fi

genes="$output/$(basename $BED .bed).bed"

complement="$output/$(basename $BED .bed).complement.bed"

for strand in $strands; do

	awk -v strand=$strand '$6 == strand' < $BED |
	
	bedtools sort -i - |

	bedtools merge -i - |
		
	awk -v strand=$strand 'BEGIN {OFS="\t"} {print $0, ".", ".", strand}'

done | bedtools sort -i - > $genes

bedtools complement -i $genes -g $chromosomeSizes | bedtools sort -i - |
awk 'BEGIN {OFS="\t"} {print $0, ".", ".", "."}' > $complement

grep '^chrM' $genes > "$output/$(basename $genes .bed).mito.bed"

grep '^chrM' $complement > "$output/$(basename $complement .bed).mito.bed"

grep -v '^chrM' $genes | grep -v '^2micron' > "$output/$(basename $genes .bed).nuc.bed"

grep -v '^chrM' $complement | grep -v '^2micron' > "$output/$(basename $complement .bed).nuc.bed"
