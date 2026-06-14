data_dir=$1
output_dir=$2
this_set=$3
patient=$4
protein=$5
numcores=$6

aln_maindir="${output_dir}aligned_haplotypes/"
recons_maindir="${output_dir}reconstructed_haplotypes/"
BAM_maindir="${output_dir}BAM/"
pat_dir="${recons_maindir}${this_set}/${patient}/${protein}/"
ref_dir="${data_dir}reads/${this_set}/ref/${patient}_ref/${protein}/"

mkdir -p "${aln_maindir}${this_set}/${patient}/${protein}/"

cd ${pat_dir}

for dir in */; do
    base_name="${dir%/}"
    hap_file="${pat_dir}${base_name}/${base_name}_consensus.fasta"
    hap_file_aln="${pat_dir}${base_name}/${base_name}_consensus_aln.fasta"
    if [ ! -f "${hap_file}" ]; then	
        sam_dir="${BAM_maindir}${this_set}/${patient}/${protein}/${base_name}/"
        sam_file="${sam_dir}/${base_name}_sorted.sam"
        samtools consensus -f fasta "${sam_file}" -o "${hap_file}"
        sed -i -e 's/ref_1/read0_1.0000/g' "${hap_file}"
        mafft --amino --keeplength --thread -${numcores} --addfragments "${hap_file}" "${ref_dir}${patient}_${protein}_ref.fa" > "${hap_file_aln}"
    else
        mafft --amino --keeplength --thread -${numcores} --addfragments "${hap_file}" "${ref_dir}${patient}_${protein}_ref.fa" > "${hap_file_aln}"
    fi    
done
exit 0
