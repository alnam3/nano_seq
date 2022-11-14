#!/bin/bash
#This script maps all cid reads on the cid references
#References have been shortened to keep variable regions only

while getopts r:p:c: flag
do
    case "${flag}" in
        r) refs=${OPTARG};;
        p) plaque=${OPTARG};;
        c) cidAorB=${OPTARG};;
    esac
done

#Does mapping and SNP calling for cidA


for fic in ../${plaque}/${cidAorB}/*.fq 
    do
    [ ! -d "../${plaque}/mapping_all_refs/${cidAorB}" ] && mkdir -p "../${plaque}/mapping_all_refs/${cidAorB}"

    name=`basename ${fic} .fq`
    ~/apps/minimap2/minimap2 -ax map-ont --secondary=no ${refs}/refs_${cidAorB}_combin_310122.fas ${fic} | samtools view -b | samtools sort - > ../${plaque}/mapping_all_refs/${cidAorB}/mapping_${name}_all_refs.bam
    echo ${fic} >> ../${plaque}/mapping_all_refs/${cidAorB}/coverage_${cidAorB}_${plaque}
    samtools coverage ../${plaque}/mapping_all_refs/${cidAorB}/mapping_${name}_all_refs.bam >> ../${plaque}/mapping_all_refs/${cidAorB}/coverage_${cidAorB}_${plaque}
done

