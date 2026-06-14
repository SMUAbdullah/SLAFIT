#!/bin/bash

data_dir=$1
output_dir=$2
this_set=$3
patient=$4
protein=$5

BAM_maindir="${output_dir}BAM/"
pat_dir="${data_dir}reads/${this_set}/MSA/${patient}/${protein}/"
ref_dir="${data_dir}reads/${this_set}/ref/${patient}_ref/${protein}/"

mkdir -p "${BAM_maindir}${this_set}/${patient}/${protein}/"

cd ${ref_dir}
bwa index "${patient}_${protein}_ref.fa"
samtools faidx "${patient}_${protein}_ref.fa"

cd ${pat_dir}

# Function to process files
process_files() {
    local r1_file=$1    
	base_name=$(echo "$r1_file" | cut -f 1 -d '.')	
    output_dir="${BAM_maindir}${this_set}/${patient}/${protein}/${base_name}/"
    mkdir -p "$output_dir"
    
    # Execute alignment and processing
    minimap2 -a "${ref_dir}${patient}_${protein}_ref.fa" "${pat_dir}${r1_file}" > "${output_dir}${base_name}.sam"
    samtools sort "${output_dir}${base_name}.sam" -o "${output_dir}${base_name}_sorted.sam"
    samtools view -bt "${ref_dir}${patient}_${protein}_ref.fa" "${output_dir}${base_name}.sam" > "${output_dir}${base_name}.bam"
    samtools sort "${output_dir}${base_name}.bam" > "${output_dir}${base_name}-sorted.bam"
    samtools index "${output_dir}${base_name}-sorted.bam"
    rm "${output_dir}${base_name}.sam" "${output_dir}${base_name}.bam"
    samtools view -bq 1 -F 4 "${output_dir}${base_name}-sorted.bam" > "${output_dir}${base_name}_sorted.bam"
    rm "${output_dir}${base_name}-sorted.bam" "${output_dir}${base_name}-sorted.bam.bai"
    samtools index "${output_dir}${base_name}_sorted.bam"
}

for file in *; do
    if [[ -f "$file" ]]; then
        process_files "$file"
    fi
done

exit 0
