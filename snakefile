import os

rule all:
    input:
        "SLAFIT_completed.txt"

rule MPL_R:
    output:
        "SLAFIT_completed.txt"
    shell:
        """
        cd src && bash SLAFIT.sh
        
        touch ../SLAFIT_completed.txt
        """
