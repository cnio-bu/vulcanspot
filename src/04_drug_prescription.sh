#!/bin/bash

### TITLE : Drug Prescription following two strategies: knowledge-based and drug repositioning.
### AUTHOR : Perales-Paton, Javier - jperales@cnio.es
### LICENSE : GPL-v3

##############################################################
## Define universe of drugs
Rscript ./src/__DTP_drug_universe.R

## Define sources
echo -e "id\tname" > ./data/DrugPrescription/sources.tsv
echo -e "1\tPANDRUGS" >> ./data/DrugPrescription/sources.tsv
echo -e "2\tLINCS" >> ./data/DrugPrescription/sources.tsv

##############################################################
### knowledge-based 
# Please, Note: This single-line command is going to perform a 
# very complex task: the drug repositioning using PanDrugs database.
# You always could check the ./src/PanDrugs folder to understand
# How this is performed. Moreover, you can check the original paper
# to understand the approach.

	# Piñeiro-Yáñez E et al. PanDrugs: a novel method to prioritize anticancer drug treatments 
	# according to individual genomic data. Genome Med. 2018 May 31;10(1):41

bash ./src/__DTP_KB_pandrugs.sh;


##############################################################
### drug repositioning
# Please, Note: This approach is the original approach of vulcanSpot
# which is based on looking for highly similar and specific connections
# of drugs which mimic the gene knock-down effect in different contexts.
# You can check the supplementary information from the manuscirpt to better
# understand this step.

## 1 Reconstruct the 'UNIVERSAL' network (biology) #####
echo -e "[INFO] Reconstruction of a UNIVERSAL DrugPPIN (Protein-Protein-Drug Interaction Network)."
# Rscript ./src/__NET_create_DrugProteinProtein_interaction_usingPINA.R &> ./log/__NET_create_DrugProteinProtein_interaction_usingPINA.log;

## 2 Contextualize (cancer) in different scenarios the 'UNIVERSAL' network #####
UNIVERSALNET="./data/NET/UNIVERSAL_DrugPPIN.PINACLUE.igraph.rds";

echo -e "[INFO] DrugPPIN contextualization.";
contexts=("PANCANCER" "LARGE_INTESTINE" "LUNG" "BREAST" "SKIN" "KIDNEY" "LIVER" "PROSTATE");

#for context in "${contexts[@]}";do
#	echo -e "[INFO] Contextualize universal network on '$context' .";
#	Rscript ./src/__NET_contextualize_DrugPPIN.R --UNIVERSALNET $UNIVERSALNET --CONTEXT $context --CONTEXT_TYPE tissue --UPC ./data/CCL/GEP_hgnc_symbol_feb2014_GRCh37.UPC.tsv --minimumUPC 0.01 &> ./log/__NET_contextualize_DrugPPIN.${context}.log;
#done

## 3 Calculate the KDCP score ########
for context in "${contexts[@]}";do
	echo -e "[INFO] Calculating KDCP score for $context at $(date)";
	Rscript ./src/__DTP_KDCP_score.R --CONTEXT $context --GDFILE "./data/GD/table_gene_essentiality.tsv" --NRND 1000 &> ./log/__DTP_KDCP_score.${context}.log;
done

#########################################################################
### CREATE TABLE OF DRUG-GENE ASSOCIATIONS
#Rscript ./src/__OUT_table_drug-gene_associations.R --INDIR ./data/DrugPrescription --OUTDIR ./data/DrugPrescription --DSCORE 0 --KDCP_SCORE 0
Rscript ./src/__OUT_table_drug-gene_associations.R --INDIR ./data/DrugPrescription --OUTDIR ./data/DrugPrescription
Rscript ./src/__DTP_lnk_genes_drugs_built.R
