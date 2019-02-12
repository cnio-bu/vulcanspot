### TITLE : hgnc_symbol mapping between 2012 (GRCh37) and 2014 (GRCh37)
### AUTHOR : Perales-Paton, Javier - jperales@cnio.es
### LICENSE : GPL-v3
### DESCRIPTION : One of the omics data was generated in 2012. Most of the data
###     were generated in 2014-2018. Mapping human gene symbols prompt with
###     some mismatches due to renaming issues between these two time points.
###     Here, we get the mapping table to match both.


options(stringsAsFactors = FALSE)

# Load library
library(biomaRt)


# GET Hugo Symbols for hg19 - Feb2014
human.14 = useMart(host="feb2014.archive.ensembl.org",
                biomart="ENSEMBL_MART_ENSEMBL",
                dataset="hsapiens_gene_ensembl")


# Get Hugo symbol and acc, gene biotype
ENG.14 <- getBM(attributes = c("ensembl_gene_id","hgnc_symbol","hgnc_id","entrezgene"),
             # filters="ensembl_gene_id",values=GEP.ensembl_gene_id,
             mart = human.14, uniqueRows=TRUE)
ENG.14$HGNC_id <- sapply(ENG.14$hgnc_id,function(z) {if(!is.na(z)) paste0("HGNC:",as.character(z)) else z})
colnames(ENG.14) <- paste0(colnames(ENG.14),".feb2014")

# GET Hugo Symbols for hg19 - may2012
human.12 = useMart(host="may2012.archive.ensembl.org",
                   biomart="ENSEMBL_MART_ENSEMBL",
                   dataset="hsapiens_gene_ensembl")
ENG.12 <- getBM(attributes = c("ensembl_gene_id","hgnc_symbol","hgnc_id","entrezgene"),
                # filters="ensembl_gene_id",values=GEP.ensembl_gene_id,
                mart = human.12, uniqueRows=TRUE)
ENG.12$HGNC_id <- sapply(ENG.12$hgnc_id,function(z) {if(!is.na(z)) paste0("HGNC:",as.character(z)) else z})
colnames(ENG.12) <- paste0(colnames(ENG.12),".may2012")

ENG <- merge(x=ENG.14,y=ENG.12,by.x="ensembl_gene_id.feb2014",by.y="ensembl_gene_id.may2012",all=TRUE)

# Print the unmatched gene symbols
cat("- This gene symbols do not match between the two years:\n",file=stdout())
table(ENG$hgnc_symbol.feb2014 == ENG$hgnc_symbol.may2012)

# Save it
write.table(ENG,file = paste0("./data/Annotation/","ANN_may2012-feb2014_GRCh37.tsv"),sep="\t",
            row.names = FALSE,col.names = TRUE)

cat("DONE.\n\n",file=stdout())
sessionInfo()
