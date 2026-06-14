#!/bin/bash

data_dir=$1
output_dir=$2
this_set=$3
patient=$4
protein=$5

BAM_dir="${output_dir}BAM/${this_set}/${patient}/${protein}/"
pat_dir="${data_dir}reads/${this_set}/MSA/${patient}/${protein}/"
ref_dir="${data_dir}reads/${this_set}/ref/${patient}_ref/${protein}/"
sv_dir="${output_dir}SV/${this_set}/${patient}/${protein}/"

mkdir -p "${sv_dir}"

cd ${BAM_dir}

for dir in */; do
    base_name="${dir%/}"
    bam_file="${BAM_dir}${base_name}/${base_name}_sorted.bam"
    mkdir "${sv_dir}${base_name}"
    if [[ -f "$bam_file" ]]; then
        cuteSV \
            "${bam_file}" \
            "${ref_dir}${patient}_${protein}_ref.fa" \
            "${sv_dir}${base_name}/${base_name}.vcf" \
            "${sv_dir}${base_name}" \
            --threads 1
    fi
done

exit 0
