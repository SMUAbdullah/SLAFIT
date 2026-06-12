# Overview
This repository contains the SLAFIT pipeline and an example dataset accompanying the manuscript.

### SLAFIT: Detection of structural variants, linkage disequilibrium estimation, and fitness inference from third generation sequencing data  
Syed Muhammad Umer Abdullah<sup>1</sup>, Dehan Cai<sup>1</sup>, Jonathan Daniel Ip<sup>2</sup>, Kelvin Kai-Wang To<sup>2</sup>, and Yanni Sun<sup>1*</sup>

<sup>1</sup> Department of Electrical Engineering, City University of Hong Kong, Tat Chee Avenue, Kowloon, Hong Kong Special Administrative Region, China  
<sup>2</sup> State Key Laboratory of Emerging Infectious Diseases, Carol Yu Centre for Infection, Department of Microbiology, Li Ka Shing Faculty of Medicine, The University of Hong Kong, 102 Pokfulam Road, Hong Kong Special
Administrative Region, China

<sup>*</sup> Corresponding author: E-mail: [yannisun@cityu.edu.hk](mailto:yannisun@cityu.edu.hk)

## SLAFIT pipeline

### Dependencies
> [samtools 1.8](https://github.com/samtools/samtools/releases/tag/1.8)

> [RVHaplo](https://github.com/dhcai21/RVHaplo)

> [minimap2](https://github.com/lh3/minimap2)

> [cuteSV](https://github.com/tjiangHIT/cuteSV)

> [MAFFT](https://mafft.cbrc.jp/alignment/software/)

> [Python 3.11.5](https://www.python.org/downloads/) with [pysam](https://pypi.org/project/pysam/) and [numpy](https://pypi.org/project/numpy/) libraries

> [MATLAB R2021b](https://www.mathworks.com/products/get-matlab.html) with [Bioinformatics toolbox](https://www.mathworks.com/products/bioinfo.html) and [Parallel computing toolbox](https://ww2.mathworks.cn/en/products/parallel-computing.html)

> [snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) (Optional. Only required if running the pipeline via snakemake)
### Installation
- The dependencies can be installed via the provided links.
- The dependencies (except MATLAB and RVHaplo) can also be installed via the following commands:
```console
conda create -n slafit python==3.11.5
conda activate slafit
conda update conda
conda update --all
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict
conda install bioconda::samtools==1.18
conda install bioconda::pysam
conda install conda-forge::numpy
conda install bioconda::cutesv
conda install bioconda::minimap2
conda install bioconda::mafft
conda install bioconda::snakemake
```
### Example data
Example data is present as `data/reads.zip` and needs to be extracted into the directory `data` before running the pipeline.
### Running the pipeline
- Execution privileges can be set by `chmod -R 700 SLAFIT-master`.
- The pipeline can be run on the example data by running the file `src/SLAFIT.sh`.
- The pipeline can also be run via snakemake by typing `snakemake --cores n` in the main directory, where `n` is the number of available CPU cores.
### Output files
- Output files generated during each step of the pipeline are stored in the directory `output`.
### Selection coefficient estimates
- The selection coefficient estimates are stored in the directory `output/s_estimates`.
- Two files, `s_MPL_R_*.txt` and `s_MPL_iden_*.txt`, are generated for the selection coefficient estimates for each patient.
- `s_MPL_R_*.txt` stores the selection coefficients estimated by considering genetic linkage, and `s_MPL_iden_*.txt` stores the selection coefficients estimated by ignoring genetic linkage.
### Running on user-supplied data
- Before running the pipeline on user-supplied data, please check lines 20-34 in the file `src/SLAFIT.sh` and edit the variables based on the data in use.
### Operating system
All scripts were written and tested on a Linux based Operating System
### Known issues and troubleshooting
- Newer versions of samtools may have compatibility issues with the pipeline. It is highly recommended to use the exact version specified in the section `Required software`.
- The naming convention for the input directories is `data/reads/MSA/<patient>/<protein>`, as demonstrated in the naming of the example data. Any deviation from this might lead to errors.
- The naming convention for the input FASTQ files is `<patient>_<protein>_t<time-point>.fq`, as demonstrated in the naming of the example data.  Any deviation from this might lead to errors.
- The script `0_FASTQ_to_BAM.sh` creates a subdirectory `BAM` in the directory `output`. Further subdirectories are created for the dataset, patient, protein, and time-point. The final output are 3 files, `<patient>_<protein>_t<time-point>_sorted.bam`, `<patient>_<protein>_t<time-point>_sorted.bam.bai`, and `<patient>_<protein>_t<time-point>_sorted_depth.txt`, created in the subdirectory for the time-point, i.e. `output/BAM/<dataset>/<patient>/<protein>/<patient>_<protein>_t<time-point>`. If these directories or files are not created or there are incorrect file extensions, this could indicate an issue with running the script. As a first step, please double check that the naming convention is followed exactly as in the directory `data/reads` and all files have the same extensions as the example data. Afterwards, please check that the exact versions of samtools and BWA are used as mentioned in the section `Required software`. If both these checks have been performed and there are issues, feel free to report the issue via email or GitHub.
- The script `AnalysisMPL_10.m` gives the error `'NumPat = 0. Check initialization settings and run again.` when it cannot find the FASTA files of the reconstructed sequences for any patient. This can happen if this script is run before one or more of the previous scripts has finished running or if there is an issue with the location of the files. Before running this script, please ensure that all previous scripts have finished running without errors. Please also avoid manually changing the location or name of any directory after running previous scripts, otherwise this script will return the error mentioned before.
## License
This repository is dual licensed as [GPL-3.0](https://github.com/SMUAbdullah/paper-MPL-short-reads/blob/master/LICENSE-GPL) (source code) and [CC0 1.0](https://github.com/SMUAbdullah/paper-MPL-short-reads/blob/master/LICENSE-CC0) (figure and documentation).
## Feedback and troubleshooting
For queries or comments, please email at [umer_973@hotmail.com](mailto:umer_973@hotmail.com).
## Citation
Software tools used in the SLAFIT pipeline can be cited as

**SAMtools**
- Li H, Handsaker B, Wysoker A et al. (2009). The sequence alignment/map format and SAMtools. _Bioinform._, **25**:2078–9. https://doi.org/10.1093/bioinformatics/btp352

**Minimap2**
- Li, H. (2018). minimap2: pairwise alignment for nucleotide sequences. _Bioinform._, **34**(18), 3094–3100. https://doi.org/10.1093/bioinformatics/bty191

**cuteSV**
- Jiang, T., Liu, S., Cao, S. & Wang, Y. (2022). Structural variant detection from long-read sequencing data with cutesv. _Variant Calling: Methods and Protocols_, 137–151. https://doi.org/10.1007/978-1-0716-2293-3_9 

**RVHaplo**
- Cai D., Sun Y. (2022). Reconstructing viral haplotypes using long reads. _Bioinform._, **38** (8), 2127–2134, https://doi.org/10.1093/bioinformatics/btac089

**MAFFT**
- Katoh, K., Misawa, K., Kuma, K. I., & Miyata, T. (2002). MAFFT: a novel method for rapid multiple sequence alignment based on fast Fourier transform. _Nucleic Acids Res._, **30**(14), 3059-3066. https://doi.org/10.1093/nar/gkf436

**MPL**
- Sohail, M. S., Louie, R. H., McKay, M. R. & Barton, J. P. (2021). MPL resolves genetic linkage in fitness inference from complex evolutionary histories. _Nat. Biotechnol._ **39**, 472–479. https://doi.org/10.1038/s41587-020-0737-3
