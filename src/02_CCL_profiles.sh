#!/bin/bash

### TITLE : Data processing for the molecular profiles (GEP,GCN,MUT) of Cancer Cell Lines (CCL)
### AUTHOR : Perales-Paton, Javier - jperales@cnio.es
### LICENSE : GPL-v3
### DESCRIPTION : 



## 1 Process Gene Expression Profiles ###
echo -e "[INFO] Processing GEPs"
Rscript ./src/__GEP_hgnc_GRCh37.R &> ./log/__GEP_hgnc_GRCh37.log;
Rscript ./src/__GEP_reads2UPC.R &> ./log/__GEP_reads2UPC.log;

## 2 Process Gene Copy-Number Profiles ###
echo -e "[INFO] Processing GCN"
Rscript ./src/__GCN_hgnc_GRCh37.R &> ./log/__GCN_hgnc_GRCh37.log;

## 3 Process Point Mutations ###
echo -e "[INFO] Processing MUT"
Rscript ./src/__MUT_makeMatrix.R &> ./log/__MUT_makeMatrix.log;

## 4 Gene Functionality binarization ###
echo -e "[INFO] Binarization of distinct molecular profiles for Alterations in Gene Functionality"
Rscript ./src/__GFA_binary_gene_alteration.R &> ./log/__GFA_binary_gene_alteration.log;

## 5 Define Genetic alterations
echo -e "[INFO] Define genetic alterations : which cell lines harbor a given genetic alterations"
Rscript ./src/__ALT_define_genetic_alterations.R &> ./log/__ALT_define_genetic_alterations.log;

echo -e "[INFO] Built genetic alterations table."
Rscript ./src/__ALT_genetic_alteration_built.R

## 6 Define and download datasets of gene dependency: CRISPR & RNAi

#####################################################################
### DEFINING DATASETS OF GENE ESSENTIALITY ##########################
#####################################################################
declare -A DS;
declare -A DS_fname;
declare -A DS_REF;

## Original version, the URLs are deprecated now
## It seems the DepMap data portal have changed the service they
## bring the data, so when I get the files these were matrices
## for the live version of the data portal. Then they created
## a 'download page' instead of direct URLs. In this page the files
## are in a different format. In the meantime I adapt the pipeline,
## a mirror of the files are the source for the data.
#DS=( ["CRISPR"]="https://depmap.org/portal/dataset/download/Avana/portal-Avana-2018-05-10.csv"
#	["RNAi"]="https://depmap.org/portal/dataset/download/RNAi_merged/portal-RNAi_merged-2018-05-10.csv")

## Mirror of the files in figshare
DS=( ["CRISPR"]="https://ndownloader.figshare.com/files/14342636?private_link=1c2f1319b7c95067fad6"
	["RNAi"]="https://ndownloader.figshare.com/files/14342639?private_link=1c2f1319b7c95067fad6")

DS_fname=( ["CRISPR"]="portal-Avana-2018-05-10.csv"
["RNAi"]="portal-RNAi_merged-2018-05-10.csv")

DS_REF=( ["CRISPR"]="AVANA_18Q2"
	["RNAi"]="Combined_RNAi")


#####################################################################
### DOWNLOAD RAW DATA FOR GENE DEPENDENCIES FROM DEPMAP #############
#####################################################################
download_DS() {
	dataset=$1;
	url=$2;
	ds_file=$3;

	if [ ! -e $ds_file ];then
		echo -e "[INFO] The dataset '$dataset' file does _NOT_ exist. Downloading it from DEPMAP.org ..."
		echo -e "[CMD]\t wget -O $ds_file $url";
		wget -O $ds_file $url;
	else echo -e "[INFO] The dataset '$dataset' file already exists."
	fi
}

#####################################################################
### DOWNLOAD DATASETS OF GENE DEPENDENCIES ##########################
#####################################################################
echo -e "id\tname\treference" > ./data/CCL/datasets.tsv;

for dataset in "${!DS[@]}";do 
	let "cnt=cnt+1";
	ds_url=${DS[$dataset]};
	ds_file="./data/CCL/${DS_fname[$dataset]}";
	echo -e "[INFO] Starting processing '$dataset' at $(date)";
	
	echo -e "[INFO] Checking if it is neccesary to download the dataset...";
	download_DS $dataset $ds_url $ds_file;

	echo -e "[INFO] ";
	echo -e "$cnt\t$dataset\t${DS_REF[$dataset]}" >> ./data/CCL/datasets.tsv;
done

#####################################################################
### GENE ESSENTIALITY PROPERTIES ####################################
#####################################################################
echo -e "[INFO] Gene Essentiality analysis."
Rscript ./src/__CGD_gene_essentiality_analysis.R --DSNAMES CRISPR:RNAi --DEPMATRICES ./data/CCL/portal-Avana-2018-05-10.csv:./data/CCL/portal-RNAi_merged-2018-05-10.csv &> ./log/__CGD_gene_essentiality_analysis.log;


## 7 Define Universe of genes and add gene annotations
echo -e "[INFO] GENES: Explore all the genes across datasets to define unique entities."
Rscript ./src/__GEN_gene_universe.R --MATRICES "./data/CCL/geneFunc.rds:./data/CCL/portal-Avana-2018-05-10.csv:./data/CCL/portal-RNAi_merged-2018-05-10.csv" &> ./log/__GEN_gene_universe.log;

echo -e "[INFO] GENES: Calculate GScore (gene scores)."
python ./src/PanDrugs/calculate_gscore.py

echo -e "[INFO] GENES: Add annotation related to Refseq description, entrez id, role in cancer, etc."
Rscript ./src/__GEN_mygene_annotation.R

echo -e "[INFO] GENES: Build table."
Rscript ./src/__GEN_gene_built.R

## 8 Define Universe of Contexts
echo -e "[INFO] CONTEXTS: Explore all the contexts across datasets to define unique entities."
Rscript ./src/__CTX_context_built.R --PROFILESMATRIX ./data/CCL/geneFunc.rds --MATRICES "./data/CCL/portal-Avana-2018-05-10.csv:./data/CCL/portal-RNAi_merged-2018-05-10.csv" &> ./log/__CTX_context_built.log;


