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
### Operating system
All scripts were written and tested on a Linux based Operating System
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
