#! /usr/bin/env bash
#Author: Alli Gombolay

program=$0

function usage () {
        echo "Usage: $program [-i] '/path/to/file1.bam etc.' [-r] 'reference genome' [-d] 'Ribose-seq directory' [-h]
          -i Filepaths of input BAM files 
          -d Location to save local Ribose-seq directory"
}

while getopts "i:r:d:h" opt;
do
    case $opt in
        i ) BAM=($OPTARG) ;;
        r ) reference=$OPTARG ;;
        d ) directory=$OPTARG ;;
        h ) usage ;;
    esac
done

if [ "$1" == "-h" ];
then
        exit
fi

strands=("positive" "negative")

flags=("-F 0x10" "-f 0x10")

maximumPeakLength=5000

asFile=$directory/ribose-seq/data/narrowPeak.as

chromosomeSizes=$directory/ribose-seq/data/reference/$reference.chrom.sizes

for samples in ${BAM[@]};
do

	filename=$(basename "$BAM")
	samples="${filename%.*}"
	
	inputDirectory=$(dirname "${BAM}")
	
	input=$inputDirectory/$samples.bam
	
	output=$directory/ribose-seq/results/$samples/peaks

	if [[ ! -f $output ]]; then
		mkdir -p $output
	fi
	
		for index in ${!strands[@]};
		do
	
			strandBAM=$output/$samples.${strands[index]}.strand.bam
        
			samtools view -hb ${flags[$index]} $input > $strandBAM
	
			experiment=$samples.${strands[index]}.strand
	
			narrowPeak=${experiment}_peaks.narrowPeak

			macs2 callpeak -t $strandBAM -n $experiment -s 25 --keep-dup all --nomodel --extsize 5 --call-summits
		
        		narrowPeak_temporary="$narrowPeak.tmp"

        		awk 'BEGIN {OFS="\t"} { if ($5 > 1000) $5 = 1000; print $0}' < $narrowPeak |
            		awk -v maxlen=$maximumPeakLength '$3 - $2 < maxlen' > $narrowPeak_temporary
        		mv $narrowPeak_temporary $narrowPeak
        
		done
done
