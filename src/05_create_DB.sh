


cut -f 1,2 ./data/CCL/contexts.tsv > ./DB/contexts.tsv
cp ./data/CCL/genes.tsv ./DB
cp ./data/CCL/datasets.tsv ./DB
cut -f 1,2,3,6,7 ./data/CCL/genetic_alterations.tsv > ./DB/genetic_alterations.tsv
cp ./data/GD/relationships.tsv ./DB
cp ./data/GD/gene_essentiality.tsv ./DB
cp ./data/DrugPrescription/sources.tsv ./DB
cp ./data/DrugPrescription/drugs.tsv ./DB
cp ./data/DrugPrescription/lnk_genes_drugs.tsv ./DB
cd ./DB
tar czvf vulcanSpot_180912.tar.gz *
