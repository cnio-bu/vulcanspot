### TITLE : Create tables to fetch DB
### AUTHOR : Perales-Paton, Javier - jperales@cnio.es
### DESCRIPTION : After the previous steps (0-4), different tables of
###	biological entities were generated to be related in 
###	a MySQL database.

# Context (or cell lineages)
echo -e "[INFO] Create 'context' table.";
cut -f 1,2 ./data/CCL/contexts.tsv > ./DB/contexts.tsv
# Genes
echo -e "[INFO] Create 'genes' table.";
cp ./data/CCL/genes.tsv ./DB
# Datasets of LoF screenings
echo -e "[INFO] Create 'datasets' table.";
cp ./data/CCL/datasets.tsv ./DB
# Genetic Alterations in panel of cancer cell lines
echo -e "[INFO] Create 'genetic_alterations' table.";
cut -f 1,2,3,6,7 ./data/CCL/genetic_alterations.tsv > ./DB/genetic_alterations.tsv
# Gene Dependencies (GD) relationship detected
echo -e "[INFO] Create 'relationships' table.";
cp ./data/GD/relationships.tsv ./DB
# Gene Essentiality
echo -e "[INFO] Create 'gene_essentiality' table.";
cp ./data/GD/gene_essentiality.tsv ./DB
# Sources for drug prescription
echo -e "[INFO] Create 'sources' table.";
cp ./data/DrugPrescription/sources.tsv ./DB
# Drug entities
echo -e "[INFO] Create 'drugs' table.";
cp ./data/DrugPrescription/drugs.tsv ./DB
# Link of genes to drug (prescription)
echo -e "[INFO] Create 'lnk_genes_drugs' table.";
cp ./data/DrugPrescription/lnk_genes_drugs.tsv ./DB
# Make tar.gz
echo -e "[INFO] Creater TAR archive with the set of tables."
todays_date=`date +%y%m%d`;
cd ./DB
tar czvf "vulcanSpot_${todays_date}.tar.gz" *
