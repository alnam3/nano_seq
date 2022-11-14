#!/bin/bash

while getopts p:c:l:d flag
do
    case "${flag}" in
        p) plaque=${OPTARG};;
        c) cidAorB=${OPTARG};;
        d) data=${OPTARG};;
        r) refs=${OPTARG};;
    esac
done

./0_mapping_cidAcidB.sh -d ${data} -r ${refs} -p ${plaque}
./1_mapping_all_refs_cid.sh -r ${refs} -p ${plaque} -c ${cidAorB}
./2_making_files_coverage.sh -p ${plaque}
./3_output_refs_secondary_uniqueread.sh -r ${refs} -p ${plaque} -c ${cidAorB}
./4_remove_gaps_cid.sh -p ${plaque} -c ${cidAorB}
./5_mapping_snp_calling_cid.sh -c ${cidAorB} -p ${plaque}
./6_building_repertoires_01.sh -p ${plaque} -c ${cidAorB}
