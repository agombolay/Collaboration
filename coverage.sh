#!/usr/bin/env bash
#Author: Alli Gombolay

#COMMAND LINE OPTIONS

#Name of the program (1_alignment.sh)
program=$0

#Usage statement of the program
function usage () {
        echo "Usage: $program [-i] 'sample1 etc.' [-r] 'reference genome' [-d] 'Ribose-seq directory' [-h]
          -i Sample names of input BAM files (i.e, sample1 for sample1.bam)
          -r File containing sizes in base pairs of chromosomes (i.e, sacCer2)
          -d Location of user's local Ribose-seq directory"
}

#Use getopts function to create the command-line options ([-i], [-r], [-d], and [-h])
while getopts "i:r:d:h" opt;
do
    case $opt in
        #Specify input as arrays to allow multiple input arguments
        i ) names=($OPTARG) ;;
	#Specify input as variable to allow only one input argument
	r ) reference=$OPTARG ;;
	d ) directory=$OPTARG ;;
        #If user specifies [-h], print usage statement
        h ) usage ;;
    esac
done

#Exit program if user specifies [-h]
if [ "$1" == "-h" ];
then
        exit
fi

for samples in ${names[@]};
do

	input=$directory/ribose-seq/results/$samples/alignment/$samples.bam
	
	chromosomeSizes=$directory/ribose-seq/data/reference/$reference.chrom.sizes
	
	output=$directory/ribose-seq/results/$samples/bedgraphs

	if [[ ! -d $output ]];
	then
		mkdir -p $output
	fi

	BothStrands=$directory/ribose-seq/results/$samples/bedgraphs/$samples.bothStrands.coverage.bg
	PositiveStrands=$directory/ribose-seq/results/$samples/bedgraphs/$samples.positiveStrands.coverage.bg
	NegativeStrands=$directory/ribose-seq/results/$samples/bedgraphs/$samples.negativeStrands.coverage.bg

	bedtools genomecov -ibam $input -g $chromosomeSizes -5 -bg > $BothStrands
	bedtools genomecov -ibam $input -g $chromosomeSizes -5 -strand + -bg > $PositiveStrands
	bedtools genomecov -ibam $input -g $chromosomeSizes -5 -strand - -bg > $NegativeStrands

	cd $output
	
	LC_COLLATE=C sort -k1,1 -k2,2n $BothStrands > $(basename $BothStrands .bg).sorted.bg
	LC_COLLATE=C sort -k1,1 -k2,2n $PositiveStrands > $(basename $PositiveStrands .bg).sorted.bg
	LC_COLLATE=C sort -k1,1 -k2,2n $NegativeStrands > $(basename $NegativeStrands .bg).sorted.bg
done

for files in $(ls $output/*.sorted.bg);
do
	bigWig="$output/$(basename $files .bg).bw"
	bedGraphToBigWig $files $chromosomeSizes $bigWig
done
