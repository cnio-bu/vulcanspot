#!/bin/bash

### TITLE : Identification of Genetic Dependencies based on molecular data of Cancer Cell Lines
### AUTHOR : Perales-Paton, Javier - jperales@cnio.es
### LICENSE : GPL-v3
### DESCRIPTION : It detects genetic dependencies based on molecular profiles of 
###		genetic alterations (transcriptome and DNA alteration, source: CCLE data portal) 
###		and matched gene essentialities (genome-wide loss-of-function screening, source: DEPMAP.org)
###		Basically, it takes the GFA (Genetic Function Alterations, transformed into LoF or GoF) to test
###		if cell lines that harbor a GFA in a given Gene (GeneA) are dependent on the function of other 
###		genes (GeneB). 
###		As an exception, for GeneA with GoF, it also tests if the function of this gene
###		itself is essential for cell viability. For instance, it is the case of oncogenic Kras (KRAS with GoF),
###		which is one of the cases where the function of GeneA itself when GoF tends to be essential.
###		The statistical test for GD (Gene Dependencies) identification is a Kolmogorov-Smirnov Test.

#####################################################################
### DEFINING DATASETS OF GENE ESSENTIALITY ##########################
#####################################################################
declare -A DS;

DS=( ["CRISPR"]="./data/CCL/portal-Avana-2018-05-10.csv"
	["RNAi"]="./data/CCL/portal-RNAi_merged-2018-05-10.csv")

#####################################################################
### GENE ESSENTIALITY PROPERTIES ####################################
#####################################################################
echo -e "[INFO] Gene Essentiality analysis."
Rscript ./src/__CGD_gene_essentiality_analysis.R --DSNAMES CRISPR:RNAi --DEPMATRICES ./data/CCL/portal-Avana-2018-05-10.csv:./data/CCL/portal-RNAi_merged-2018-05-10.csv &> ./log/__CGD_gene_essentiality_analysis.log;

echo -e "[INFO] Gene Essentiality table."
Rscript ./src/__CGD_gene_essentiality_built.R

#####################################################################
### TEST FOR GENETIC CANCER DEPENDENCIES ############################
#####################################################################
for dataset in "${!DS[@]}";do 
	echo -e "[INFO] Testing for GDs in '$dataset' dataset.";
	Rscript ./src/__CGD_identify_genetic_dependencies.R --ESSENTIALITY $ds_file --GENEFUNCTION ./data/CCL/geneFunc.rds --SCALEESSENTIALLITY TRUE --MINSIZE 5 --MAXSIZE 2000 --DATASETNAME $dataset --NPROC 1 --NPERMUT 10000 --minimumSKEWNESS -0.5 &> ./log/__CGD_identify_genetic_dependencies.${dataset}.log;

	echo -e "[INFO] Adjusting for GDs in '$dataset' that are associated to Oncogenic Dependencies.";
	Rscript ./src/__CGD_colinearity_correction.R --DATASETNAME $dataset --OVERLAP 0.7 --JACKARD 0.6 &> ./log/__CGD_colinearity_correction.${dataset}.log;
done

#####################################################################
### CREATE A TABLE OF GENE DEPENDENCIES FOR DRUG PRESCRIPTION #######
#####################################################################
## Note: Indeed, these GDs are the only of interest for the database.
echo -e "[INFO] Create table for gene dependencies (FDR < 0.25)"
Rscript ./src/__OUT_table_gene_essentiality.R --INDIR ./data/GD --OUTDIR ./data/GD --FDR 0.25;

# The same, but for DB
Rscript ./src/__CGD_GD_relationships_built.R --FDR 0.25;
