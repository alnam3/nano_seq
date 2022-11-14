#!/bin/bash

#This script creates separate cidA/cidB files and folders.
#Fields have to be given in input:
# - path to data folder 
# - plate specification
# - path to reference files

while getopts d:r:p: flag
do
    case "${flag}" in
        d) data=${OPTARG};;
        r) refs=${OPTARG};;
        p) plaque=${OPTARG};;
    esac
done

for fic in ${data}/${plaque}/*.fastq.gz
do 
    #checking if required folders exist and creating it if not
    [ ! -d "../${plaque}/cidA/" ] && mkdir -p "../${plaque}/cidA" 
    [ ! -d "../${plaque}/cidB/" ] && mkdir -p "../${plaque}/cidB" 

	sample=`basename ${fic} .fastq.gz`
	#cidA

	#creating the files
	$HOME/apps/minimap2/minimap2 --for-only -ax map-ont ${refs}/common_up_cidA_short.fa ${fic} | samtools fastq -n -F 4 - > ../${plaque}/cidA/cidA_for_${sample}.fq
	$HOME/apps/minimap2/minimap2 --rev-only -ax map-ont ${refs}/common_up_cidA_short.fa ${fic} | samtools fastq -n -F 4 - | $HOME/apps/seqkit seq --reverse --complement -t DNA > ../${plaque}/cidA/cidA_rev_${sample}.fq
	cat ../${plaque}/cidA/cidA_for_${sample}.fq ../${plaque}/cidA/cidA_rev_${sample}.fq > ../${plaque}/cidA/cidA_${sample}.fq
	rm ../${plaque}/cidA/cidA_for_${sample}.fq
	rm ../${plaque}/cidA/cidA_rev_${sample}.fq

	#cidB
	#creating the files
	$HOME/apps/minimap2/minimap2 --for-only -ax map-ont ${refs}/common_mid_cidB.fa ${fic} | samtools fastq -n -F 4 - > ../${plaque}/cidB/cidB_for_${sample}.fq
	$HOME/apps/minimap2/minimap2 --rev-only -ax map-ont ${refs}/common_mid_cidB.fa ${fic} | samtools fastq -n -F 4 - | $HOME/apps/seqkit seq --reverse --complement -t DNA > ../${plaque}/cidB/cidB_rev_${sample}.fq
	cat ../${plaque}/cidB/cidB_for_${sample}.fq ../${plaque}/cidB/cidB_rev_${sample}.fq > ../${plaque}/cidB/cidB_${sample}.fq
	rm ../${plaque}/cidB/cidB_for_${sample}.fq
	rm ../${plaque}/cidB/cidB_rev_${sample}.fq

done	
