#!/bin/bash

#This script does the mapping + snp calling on the restricted pool of references which have been obtained through the previous steps.

module load bcftools

while getopts p:c:l: flag
do
    case "${flag}" in
        p) plaque=${OPTARG};;
        c) cidAorB=${OPTARG};;
    esac
done

[ ! -d "../${plaque}/Final_mapping/${cidAorB}" ] && mkdir -p "../${plaque}/Final_mapping/${cidAorB}"

rm ../${plaque}/Final_mapping/${cidAorB}/*.bam
rm ../${plaque}/Final_mapping/${cidAorB}/*.vcf

for fic in ../${plaque}/${cidAorB}/selected_reads/*.fq
do
    sample=`basename ${fic} .fqi | cut -d '_' -f 1-3`
    $HOME/apps/minimap2/minimap2 -ax map-ont --secondary=no ../${plaque}/output_refs/${cidAorB}/refs_no_gaps/refs_${sample}_ng.fas ${fic} | samtools view -b | samtools sort - > ../${plaque}/Final_mapping/${cidAorB}/mapping_${sample}.bam
    $HOME/apps/bcftools-1.15/bcftools mpileup -q 10 -d 75000 -B -a AD --config ont -f ../${plaque}/output_refs/${cidAorB}/refs_no_gaps/refs_${sample}_ng.fas ../${plaque}/Final_mapping/${cidAorB}/mapping_${sample}.bam | $HOME/apps/bcftools-1.15/bcftools call -m -P 0.5 > ../${plaque}/Final_mapping/${cidAorB}/${sample}_P05.vcf
done

