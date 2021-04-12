#!/bin/bash

# load R
module purge
sleep 5
source /opt/asn/etc/asn-bash-profiles-special/modules.sh
sleep 5
module load R/4.0.2
sleep 5

# R CMD INSTALL ~/QTL/DLinkMap/DSPRqtl_2.0-5.tar.gz
Rscript ~/QTL/DLinkMaP/mapping/MAP_general.R /home/ualcpr/QTL/DLinkMaP/parallelizer/run_scripts/run_0009/params_0009.csv
sleep 10
# uniqing multiple runs
cp /home/ualcpr/QTL/DLinkMaP/parallelizer/run_scripts/run_0009/maleWt/p-value_Male_avgbyvial.csv /home/ualcpr/QTL/DLinkMaP/parallelizer/run_scripts/run_0009/maleWt/temp-pval.csv
cat /home/ualcpr/QTL/DLinkMaP/parallelizer/run_scripts/run_0009/maleWt/temp-pval.csv | uniq > /home/ualcpr/QTL/DLinkMaP/parallelizer/run_scripts/run_0009/maleWt/p-value_Male_avgbyvial.csv
cp /home/ualcpr/QTL/DLinkMaP/parallelizer/run_scripts/run_0009/maleWt/logLike_Male_avgbyvial.csv /home/ualcpr/QTL/DLinkMaP/parallelizer/run_scripts/run_0009/maleWt/temp-log_like.csv
cat /home/ualcpr/QTL/DLinkMaP/parallelizer/run_scripts/run_0009/maleWt/temp-log_like.csv | uniq > /home/ualcpr/QTL/DLinkMaP/parallelizer/run_scripts/run_0009/maleWt/logLike_Male_avgbyvial.csv
