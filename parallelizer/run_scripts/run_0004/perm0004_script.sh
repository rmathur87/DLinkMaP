#!/bin/bash

# load R
source /opt/asn/etc/asn-bash-profiles-special/modules.sh
module load R/4.0.2

# R CMD INSTALL ~/QTL/DLinkMap/DSPRqtl_2.0-5.tar.gz
Rscript ~/QTL/DLinkMaP/mapping/MAP_general.R /home/ualcpr/QTL/DLinkMaP/parallelizer/run_scripts/run_0004/params_0004.csv

