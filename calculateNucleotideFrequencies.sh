#! /usr/bin/env bash

program=$0

function usage () {
        echo "Usage: $program [-i] 'sample1 etc.' [-r] 'reference genome' [-d] 'Ribose-seq directory' [-h]
        -i Sample names of input BAM files (i.e, sample1 for sample1.bam)
        -r File containing sizes in base pairs of chromosomes (i.e, sacCer2)
        -d Location to save local Ribose-seq directory"
}

while getopts "i:r:d:h" opt;
do
    case $opt in
        i ) names=($OPTARG) ;;
	r ) reference=$OPTARG ;;
	d ) directory=$OPTARG ;;
    	h ) usage ;;
    esac
done

if [ "$1" == "-h" ];
then
        exit
fi

# Mononucleotides, dinucleotides, and trinucleotides
sizes="1 2 3"

modes=("all")

arguments=("")

input=$directory/ribose-seq/results/$samples/alignment

output=$directory/ribose-seq/results/$samples/nucleotideFrequencies

if [[ ! -d $output ]]; then
    mkdir -p $output
fi

offset_minimum=-100
offset_maximum=100

BAM=$input/$sample.bam

FASTA=$directory/ribose-seq/data/reference/$reference.fa

for index in ${!modes[@]};
do

        mode=${modes[$index]}
        
        argument=${arguments[$index]}

        tables="$output/$sample.$mode.nucleotideFrequencies.tab"
        
        if [[ -f $tables ]];
        then
            rm -f $tables
        fi

        if [[ $mode == "all" ]];
        then
            BackgroundFrequencies="$directory/ribose-seq/results/backgroundNucleotideFrequencies/$reference.genome.nucleotide.frequencies.tab"
        fi

        #Signals need to be reverse complemented since sequence is reverse complement of the captured strand
        for size in $sizes;
        do
            python2.7 calculateNucleotideFrequencies.py $BAM $FASTA --verbose --region-size $size $arguments --revcomp-strand \
            --background-freq-table $BackgroundFrequencies --offset-min $offset_minimum --offset-max $offset_maximum >> $tables
        done
    
done
