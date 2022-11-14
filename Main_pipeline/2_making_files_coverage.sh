#!/bin/bash

#Obtaining the coverage of each reference.


while getopts p: flag
do
    case "${flag}" in
        p) plaque=${OPTARG};;
    esac
done

for cidAorB in cidA cidB
do
    [ ! -d "../${plaque}/coverage" ] && mkdir -p "../${plaque}/coverage"

    for fic in ../${plaque}/mapping_all_refs/${cidAorB}/*.bam 
    do
        echo $fic
        name=`basename ${fic} .bam | cut -d '_' -f 2-4`
        samtools coverage ${fic} > ../${plaque}/coverage/cov_${name}
    done
done

