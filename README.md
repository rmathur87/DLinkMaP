# DLinkMaP
Drosophila Linkage Mapping Pipeline (DLinkMaP)


This R code was developed to conduct Quantitative Trait Loci (QTL) mapping based on haplotype probabilities for the Drosophila Synthetic Population Resource (DSPR). The code performs a linear mixed model analysis at 10kB intervals throughout the Drosophila genome. Based on the round robin crossing design described in Dew-Budd et al. (2018, in preparation), random effect terms are included for month, round robin group, RIL by month, and cross (all of these are included by default in the software). Using haplotype frequencies as found in the "DSPRqtlDataA", test founder probabilities corresponding to all 8 x 8 = 64 founder by maternal/paternal combinations are calculated. Six separate non-null models corresponding to additive, dominant, and full found by parent effects as main effects and diet interactions were included in the model (all of these models are included by default in the software). Statistical inference is conducted by calculating the χ2-statistics based on the difference in model log-likelihoods in a hierarchical manner, and determining p-values based on the χ2-distribution.

The software has the capability to conduct epistasis (gene-gene interaction) modeling by calculating the cross product of the haplotype frequencies for the additive model. Dominant and full models are not considered for epistatic modeling. By default the software calculates epistatic effects for all pairwise 10kB positions in the genome. As  in  the  non-epistatic  analysis,  the  statistical  inference  is  based  onχ2-statistics  and  negative  log  p-values  are  reported.

For  each  significant  peak (by default the complete genome includes 11,768 QTLs, thus a negative log p-value above 5.37 is bonferroni adjusted significant),  the  range  of  influence  for  each  peak,  95%  confidence  intervals is  computed  using  a  LOD  drop  of  two.  For  the  epistatic  models,  a  marginal  confidence  interval  was  calculated  for  each  QTL  involved  in  a  significant  epistatic  interaction  (threshold  of  six  for  the  –log  p-value).  A  Bayesian  model  was  conducted  to  estimate  the  variance  explained  for  each  non-epistatic  peak  and  the  significant  epistatic  interactions  (more  details  are  described  in  supplemental  methods). 


Running DLinkMaP:
By running the "MAP_general.R" script with appropriate inputs (see below), the QTL mapping will be conducted (using the user defined phenotype) for all QTLs in DSPR Population A. The following R scripts are included in the pipeline:

MAP_general.R:

Conduct QTL Mapping for the user defined phenotype ('metabolites', 'triglyceride', 'male weight', 'female weight', and 'trehalose' have been tested)

  Inputs:

    commDir: directory where all QTL Mapping scripts are contained

    p: Whether permutation testing should be (p=1) or should not be (p=0) conducted

    outDir: directory where all output will be written

    phenotype: which phenotype is being analyzed - 'metabolite', 'weight', 'trehalose', and 'TG' have been tested

    weight.sex: if the 'weight' phenotype is being analyzed, which sex is analyzed - for all other phenotypes, ignore this!

    weight.type: if the 'weight' phenptype is being anayzed, whether the original data or the average by vial is being analyzed - for all other phenotypes, ignore this!

    fileName: the filepath and name of the dataset to be analyzed

    epistaticModel: True/False if an epistatic model is run. F is normal mapping with additive, dominant and full effects; T is additive with epistatic effects for all pairs of QTLs

    epistaticQTL: Integer (out of the total number of QTLs tested) for the epistatic model


MAPFun_general.R; FUN.R; gradMM.R; MM_Process.R; QC.R; QTL_Process.R
These script contain key functions which are used by MAP_general.R to conduct the design matrix setup, PCA, and inference.

Necessary R packages include:
- gplots - https://cran.r-project.org/web/packages/gplots/index.html
- lme4 - https://cran.r-project.org/web/packages/lme4/index.html
- DSPRqtl - http://wfitch.bio.uci.edu/~dspr/Tools/Tutorial/index.html
- DSPRqtlDataA - http://wfitch.bio.uci.edu/~dspr/Tools/Tutorial/index.html


