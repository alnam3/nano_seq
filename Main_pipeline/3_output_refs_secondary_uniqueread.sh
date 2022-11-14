#!/bin/bash

#This script aims at making the reference file from the coverage file.
#Inputs : cidAorB, path to reference, plaque, cidB version, dividor if want to change it (5th argument, no name -- is used to divide the full coverage. Defaults to 15)
#Sorting the references and removing the highly similar references to prevent SNP calling issues.

module load R/4.0.5

#Variables given as input
while getopts r:p:v:c: flag
do
    case "${flag}" in
        r) refs=${OPTARG};;
        p) plaque=${OPTARG};;
        v) version=${OPTARG};;
        c) cidAorB=${OPTARG};;
    esac
done

#Clearing the existing files and making the required directories
rm ../$plaque/output_refs/$cidAorB/*.fa*
rm ../$plaque/$cidAorB/selected_reads/*
rm ../$plaque/$cidAorB/selected_reads/*.sam
rm ../$plaque/$cidAorB/selected_reads/*.bam


[ ! -d "../${plaque}/output_refs/$cidAorB" ] && mkdir -p "../${plaque}/output_refs/$cidAorB"
[ ! -d "../${plaque}/output_refs/$cidAorB/sorted" ] && mkdir -p "../${plaque}/output_refs/$cidAorB/sorted"
[ ! -d "../${plaque}/$cidAorB/selected_reads" ] && mkdir -p "../${plaque}/$cidAorB/selected_reads"


#Chosing the right files of refs
if [ $cidAorB == "cidB" ]
then
    ref_fic=${refs}"/refs_"${cidAorB}"_"${version}".fas"
else 
    ref_fic=${refs}"/refs_"${cidAorB}"_combin_310122.fas"
fi

echo $ref_fic
#Outputting the right references
for fic in ../${plaque}/coverage/cov_*
do
    #calculating threshold
    div=${9:-15}
    thresh=$(awk -v d=$div 'BEGIN{FS="\t";sum=0} {sum+=$7} END {print int(sum/d)}' $fic)
    echo $div $thresh
    #outputting references
    echo ${fic}
    name=`echo ${fic} | cut -d '_' -f 3-4`
    awk -F "\t" -v t=${thresh} 'BEGIN {OFS=FS} $7 > t && NR > 1 {print $1}' ${fic} > ../$plaque/name_${cidAorB}_output_refs_thresh
    $HOME/apps/seqkit grep -f ../$plaque/name_${cidAorB}_output_refs_thresh ${ref_fic} > ../${plaque}/output_refs/${cidAorB}/refs_${cidAorB}_${name}.fas
    Rscript  --vanilla ./sort_refs.R ../$plaque/output_refs/$cidAorB/refs_${cidAorB}_$name.fas ../$plaque/output_refs/$cidAorB/sorted/
    
    #Keeping only reads mapping on the selected references and making a fastq file for further steps
    samtools index ../${plaque}/mapping_all_refs/${cidAorB}/mapping_${cidAorB}_${name}_all_refs.bam
    while read -r line
    do 
        samtools view -hS -o ../${plaque}/${cidAorB}/selected_reads/${cidAorB}_${name}_$line.sam ../${plaque}/mapping_all_refs/${cidAorB}/mapping_${cidAorB}_${name}_all_refs.bam $line
        samtools view ../$plaque/${cidAorB}/selected_reads/${cidAorB}_${name}_${line}.sam | awk -F "\t" '{print $1}' >> ../$plaque/$cidAorB/mappedID_${cidAorB}_${name}
    done < ../$plaque/name_${cidAorB}_output_refs_thresh
    sort ../$plaque/$cidAorB/mappedID_${cidAorB}_${name} | uniq > ../$plaque/$cidAorB/mappedID_${cidAorB}_${name}_uniq 
    
    #Making new fastq files with selected reads only
    ~/apps/seqtk/seqtk subseq ../$plaque/$cidAorB/${cidAorB}_${name}.fq ../$plaque/$cidAorB/mappedID_${cidAorB}_${name}_uniq > ../$plaque/$cidAorB/selected_reads/${cidAorB}_${name}_sel.fq 
done

rm ../$plaque/$cidAorB/selected_reads/*.sam
