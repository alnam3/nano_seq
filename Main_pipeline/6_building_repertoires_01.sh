#!/bin/bash
#This script aims at obtaining repertoires. There are 3 cases (i) cases without SNPs (easy); (2) cases with either only alternative alleles or a single polymorphic SNP (ok); (3) cases with several polymorphic SNPs

while getopts p:c: flag
do
    case "${flag}" in
        p) plaque=${OPTARG};;
        c) cidAorB=${OPTARG};;
    esac 
done

[ ! -d "../${plaque}/repertoires/${cidAorB}/no_SNPs" ] && mkdir -p "../${plaque}/repertoires/${cidAorB}/no_SNPs"
[ ! -d "../${plaque}/repertoires/$cidAorB/SNPs" ] && mkdir -p "../$plaque/repertoires/$cidAorB/SNPs"

#Removing previous files
rm ../$plaque/repertoires/$cidAorB/no_SNPs/*
rm ../$plaque/repertoires/$cidAorB/SNPs/*
rm ../$plaque/Final_mapping/$cidAorB/Variant_calling/*

#Listing files which have or do not have no SNPs
grep -RL "GT:PL" ../${plaque}/Final_mapping/${cidAorB}/*.vcf > ../${plaque}/no_SNPs_${cidAorB}
grep -Rl "GT:PL" ../${plaque}/Final_mapping/${cidAorB}/*.vcf > ../${plaque}/SNPs_${cidAorB}

#Removing previous files
rm ../$plaque/repertoires/$cidAorB/no_SNPs/*
rm ../$plaque/repertoires/$cidAorB/SNPs/*

###For files which do not have SNPs
while read -r line; do 
    name=`basename $line .vcf | cut -d '_' -f 1-3` 
    samtools coverage ../${plaque}/Final_mapping/${cidAorB}/mapping_${name}.bam | awk -F "\t" 'BEGIN {OFS=","} {print $1,$7}' > ../${plaque}/repertoires/${cidAorB}/no_SNPs/${name}
done < ../${plaque}/no_SNPs_$cidAorB

###for files which do have SNPs
##Also making the coverage file, but names will have to be updated
while read -r line; do
    name=`basename $line .vcf | cut -d '_' -f 1-3`
    samtools coverage ../${plaque}/Final_mapping/${cidAorB}/mapping_${name}.bam | awk -F "\t" 'BEGIN {OFS=","} {print $1,$7}' > ../${plaque}/repertoires/${cidAorB}/SNPs/${name}
done < ../${plaque}/SNPs_$cidAorB


#Files which have multiple polymorphic SNPs
grep -c "0/1" ../$plaque/Final_mapping/$cidAorB/*.vcf | awk 'BEGIN{FS=":"} $2>1 {print $1}' > ../$plaque/multiple_SNPs_$cidAorB

#Other files, for which a consensus can be used (having either multiple alt SNPs or a single polymorphic SNP)
diff ../$plaque/multiple_SNPs_${cidAorB} ../$plaque/SNPs_${cidAorB} | grep '>' | sed 's/> //' > ../$plaque/cns_SNPs_$cidAorB

while read -r line; do 
    name=`basename $line .vcf | cut -d '_' -f 1-3`
    echo $name
    #call variants
    /home/anamias/apps/bcftools-1.15/bcftools mpileup -q 10 -d 75000 -f ../$plaque/output_refs/$cidAorB/refs_no_gaps/refs_${name}_ng.fas ../$plaque/Final_mapping/$cidAorB/mapping_${name}.bam | /home/anamias/apps/bcftools-1.15/bcftools call -mv -Ob -o ../$plaque/Final_mapping/$cidAorB/Variant_calling/calls_${name}.bcf
    /home/anamias/apps/bcftools-1.15/bcftools index ../$plaque/Final_mapping/$cidAorB/Variant_calling/calls_${name}.bcf
    #normalizing indels
    /home/anamias/apps/bcftools-1.15/bcftools norm -f ../$plaque/output_refs/$cidAorB/refs_no_gaps/refs_${name}_ng.fas ../$plaque/Final_mapping/$cidAorB/Variant_calling/calls_${name}.bcf -Ob -o ../$plaque/Final_mapping/$cidAorB/Variant_calling/calls_${name}.norm.bcf
    #filter adjacent indels within 5 bp
    /home/anamias/apps/bcftools-1.15/bcftools filter --IndelGap 5 ../$plaque/Final_mapping/$cidAorB/Variant_calling/calls_${name}.norm.bcf -Ob -o ../$plaque/Final_mapping/$cidAorB/Variant_calling/calls_${name}.norm.flt-indels.bcf
    /home/anamias/apps/bcftools-1.15/bcftools index ../$plaque/Final_mapping/$cidAorB/Variant_calling/calls_${name}.norm.flt-indels.bcf
    #apply variants to create consensus sequence
    #option -I is used to enable recovering pseudo-heterozyous positions -- only works with bcftools > 1.13
    cat ../$plaque/output_refs/$cidAorB/refs_no_gaps/refs_${name}_ng.fas | /home/anamias/apps/bcftools-1.15/bcftools consensus --haplotype I ../$plaque/Final_mapping/$cidAorB/Variant_calling/calls_${name}.norm.flt-indels.bcf | awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' > ../$plaque/Final_mapping/$cidAorB/Variant_calling/consensus_${name}.fa
done < ../${plaque}/cns_SNPs_$cidAorB

rm ../$plaque/Final_mapping/$cidAorB/Variant_calling/*.bcf

#(3) files having more than one polymorphic SNP -> haplotype phasing
while read -r line
do
    name=`basename $line .vcf | cut -d '_' -f '1-3'`
    samtools index ../$plaque/Final_mapping/$cidAorB/mapping_$name.bam
    samtools faidx ../$plaque/output_refs/$cidAorB/refs_no_gaps/refs_${name}_ng.fas
    whatshap phase --ignore-read-groups -o ../$plaque/Final_mapping/$cidAorB/Variant_calling/phased_$name.vcf --reference ../$plaque/output_refs/$cidAorB/refs_no_gaps/refs_${name}_ng.fas ../$plaque/Final_mapping/$cidAorB/${name}_P05.vcf ../$plaque/Final_mapping/$cidAorB/mapping_$name.bam
done < ../$plaque/multiple_SNPs_$cidAorB
