# DLinkMaP
Drosophila Linkage Mapping Pipeline


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
