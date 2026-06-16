#!/bin/sh

# Directory declarations, no need to edit
#--------------------------------------------------------------------
cd ..
main_dir=$(pwd)/
data_dir=${main_dir}"data/"
script_dir=${main_dir}"src/"
output_dir=${main_dir}"output/"
logs_dir=${main_dir}"logs/"
bash_scripts_dir=${main_dir}"bash_scripts/"
#--------------------------------------------------------------------
# Variable declarations, please edit based on data
#--------------------------------------------------------------------
numcores=8 		                 # number of CPU cores to be used
this_set="10"                    # name of the dataset. Helps in naming the files
protein="synth"                  # name of the protein. Helps in naming the files
declare -a patients=("p1" "p2" ) # names of the patients
genome_start=1                   # starting index of genome
genome_end=4500                  # ending index of genome
maximum_sequence_length=10000    # maximum length of the sequence. User is welcome to adjust based on requirement
maximum_time_points=40 	         # maximum number of time points. User is welcome to adjust based on requirement
gamma=20 		                 # regularization parameter
thresh=0.10                      # threshold below which trajectories are considered noise
rvhaplo_dir=""                   # location of the RVHaplo-main directory, such as /home/software/RVHaplo-main/
#--------------------------------------------------------------------
for patient in ${patients[@]}
do
# Read alignment
#--------------------------------------------------------------------
cd $script_dir
./0_read_alignment.sh $data_dir $output_dir $this_set $patient $protein
#--------------------------------------------------------------------
# Structural variant detection
#--------------------------------------------------------------------
cd $script_dir
./1_SV_detection.sh $data_dir $output_dir $this_set $patient $protein
#--------------------------------------------------------------------
# Generate files to perform haplotype reconstruction
#--------------------------------------------------------------------
cd $script_dir
python 2_reconstruction_files_gen.py $data_dir $output_dir $bash_scripts_dir $logs_dir $rvhaplo_dir $this_set $patient $protein $numcores
#--------------------------------------------------------------------
# Perform haplotype reconstruction
#--------------------------------------------------------------------
chmod -R 700 ${bash_scripts_dir}"reconstruction_call/"${this_set}
cd ${bash_scripts_dir}"reconstruction_call/"${this_set}
for filename in *
do ./${filename} # an & sign can be placed after this command. It is used for running the commands in parallel. It is recommended to first test the code in serial. After testing, the code can be run in parallel. Running multiple processes in parallel without testing can lead to the creation of a number of unwanted processes which are difficult to kill, hence it is a good practice to test in serial first.
done
#--------------------------------------------------------------------
# Reconstructed haplotype alignment
#--------------------------------------------------------------------
cd $script_dir
./3_recons_align.sh $data_dir $output_dir $this_set $patient $protein $numcores
#--------------------------------------------------------------------
# Generate files to perform haplotype reconstruction
#--------------------------------------------------------------------
cd $script_dir
python 4_convert_format.py $output_dir $this_set $patient $protein
#--------------------------------------------------------------------
done
# Linkage disequilibrium and fitness estimation
#--------------------------------------------------------------------
cd $script_dir
matlab -nodisplay -nojvm -nosplash -nodesktop -r "MPL_7(\"${main_dir}\",\"${data_dir}\",\"${this_set}\",\"${protein}\",\"${output_dir}\",\"${genome_start}\",\"${genome_end}\",\"${gamma}\",\"${thresh}\"); exit";
exit 0
