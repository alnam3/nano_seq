#!/bin/bash

#Removing gaps from alignment in fasta sequences
while getopts p:c: flag
do
    case "${flag}" in
        p) plaque=${OPTARG};;
        c) cidAorB=${OPTARG};;
    esac
done

[ ! -d "../${plaque}/output_refs/${cidAorB}/refs_no_gaps" ] && mkdir -p "../${plaque}/output_refs/${cidAorB}/refs_no_gaps"

rm ../$plaque/output_refs/$cidAorB/refs_no_gaps/*.fa* #remove previously present fai files

for fic in ../${plaque}/output_refs/${cidAorB}/sorted/*.fa
do
    sample=`echo ${fic} | cut -d '/' -f 6 | cut -d '_' -f 1-4` 
    ~/apps/seqkit replace -s -p "\-" -r '' ${fic} | awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' > ../${plaque}/output_refs/${cidAorB}/refs_no_gaps/${sample}_ng.fas 
done 

