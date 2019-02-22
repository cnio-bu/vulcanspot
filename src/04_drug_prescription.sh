#!/bin/bash

### TITLE : Drug Prescription following two strategies: knowledge-based and drug repositioning.
### AUTHOR : Perales-Paton, Javier - jperales@cnio.es
### LICENSE : GPL-v3
### DESCRIPTION : Gene Dependencies were detected in the previous step, pointing out the candidate
###	targets for cancer therapies. Herein, two complementary approaches are used to prescribe
###	drugs among these candidates.

##############################################################
## Create directory if needed
if [ ! -e ./data/DrugPrescription ];then mkdir ./data/DrugPrescription;fi

##############################################################
## Define universe of drugs
echo -e "[INFO] Obtaining the universe of drugs.";
Rscript ./src/__DTP_drug_universe.R

## Define sources
echo -e "[INFO] Obtaining the universe of sources (pandrugs and lincs).";
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
echo -e "[INFO] Prescribing drugs using PanDrugs system."
bash ./src/__DTP_KB_pandrugs.sh;


##############################################################
### drug repositioning
# Please, Note: This approach is the original approach of vulcanSpot
# which is based on looking for highly similar and specific connections
# of drugs which mimic the gene knock-down effect in different contexts.
# You can check the supplementary information from the manuscirpt to better
# understand this step.
echo -e "[INFO] Prescribing drugs using KDCP approach."

## Calculated TAU matrices for the KDCP approach
declare -A TAU;

TAU=( ["BREAST"]="https://ndownloader.figshare.com/files/14341517?private_link=b34fc8e7f4d77ca6d580"
["KIDNEY"]="https://ndownloader.figshare.com/files/14341523?private_link=b34fc8e7f4d77ca6d580"
["LARGE_INTESTINE"]="https://ndownloader.figshare.com/files/14341526?private_link=b34fc8e7f4d77ca6d580"
["LIVER"]="https://ndownloader.figshare.com/files/14341532?private_link=b34fc8e7f4d77ca6d580"
["LUNG"]="https://ndownloader.figshare.com/files/14341535?private_link=b34fc8e7f4d77ca6d580"
["PANCANCER"]="https://ndownloader.figshare.com/files/14341538?private_link=b34fc8e7f4d77ca6d580"
["PROSTATE"]="https://ndownloader.figshare.com/files/14341547?private_link=b34fc8e7f4d77ca6d580"
["SKIN"]="https://ndownloader.figshare.com/files/14341550?private_link=b34fc8e7f4d77ca6d580")

## 1 Reconstruct the 'UNIVERSAL' network (biology) #####
echo -e "[INFO] Reconstruction of a UNIVERSAL DrugPPIN (Protein-Protein-Drug Interaction Network)."
# Rscript ./src/__NET_create_DrugProteinProtein_interaction_usingPINA.R &> ./log/__NET_create_DrugProteinProtein_interaction_usingPINA.log;

## 2 Contextualize (cancer) in different scenarios the 'UNIVERSAL' network #####
UNIVERSALNET="./data/NET/UNIVERSAL_DrugPPIN.PINACLUE.igraph.rds";

echo -e "[INFO] DrugPPIN contextualization.";
contexts=("PANCANCER" "LARGE_INTESTINE" "LUNG" "BREAST" "SKIN" "KIDNEY" "LIVER" "PROSTATE");

## Precalculated in the repository
#for context in "${contexts[@]}";do
#	echo -e "[INFO] Contextualize universal network on '$context' .";
#	Rscript ./src/__NET_contextualize_DrugPPIN.R --UNIVERSALNET $UNIVERSALNET --CONTEXT $context --CONTEXT_TYPE tissue --UPC ./data/CCL/GEP_hgnc_symbol_feb2014_GRCh37.UPC.tsv --minimumUPC 0.01 &> ./log/__NET_contextualize_DrugPPIN.${context}.log;
#done

## 4 Calculate the KDCP score ########
# First download the TAU matrices for context
# See 'https://github.com/jperales/geneConf_GSMM' for 
# code on how to extract signatures.
# See 'Subramanian A et al 2017' to see how to calculate TAU
echo -e "[INFO] Downloading precalculated TAU matrices."
# Define a bash function to download files using wget
download_DS() {
	context=$1;
	url=$2;
	flname=$3;

	if [ ! -e $flname ];then
		echo -e "[INFO] The file for the TAU matrix for '$context' does _NOT_ exist. Downloading it from Figshare..."
		echo -e "[CMD]\t wget -O $flname $url";
		wget -O $flname $url;
	else echo -e "[INFO] The file for the TAU matrix for '$context' file already exists: $flname";
	fi
}

for context in "${contexts[@]}";do 
	echo -e "[INFO] Checking if the precalculated TAU matrix if available for $context";
	fl_url="${TAU[$context]}";
	flname="./data/CMap-L1000/${context}.TAU_matrix.rds";
	download_DS $context $fl_url $flname;
done

# Second calculate KDCP scores, which takes TAU and adjust this
# score for the distance between drugs and targets of interest
for context in "${contexts[@]}";do
	echo -e "[INFO] Calculating KDCP score for $context at $(date)";
	Rscript ./src/__DTP_KDCP_score.R --CONTEXT $context --GDFILE "./data/GD/table_gene_essentiality.tsv" --NRND 1000 &> ./log/__DTP_KDCP_score.${context}.log;
done

#########################################################################
### CREATE TABLE OF DRUG-GENE ASSOCIATIONS
Rscript ./src/__OUT_table_drug-gene_associations.R --INDIR ./data/DrugPrescription --OUTDIR ./data/DrugPrescription
Rscript ./src/__DTP_lnk_genes_drugs_built.R
