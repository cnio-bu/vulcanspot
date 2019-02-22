### TITLE: Create 'gene_essentiality' table
### AUTHOR : Perales-Paton, Javier - jperales@cnio.es
### LICENSE : GPL-v3

library(reshape2)

Genes.df <- read.table("./data/CCL/genes.tsv",sep="\t",header=TRUE,quote="",fill=TRUE)
id_genes <- setNames(Genes.df$id,Genes.df$symbol)

Datasets.df <- read.table("./data/CCL/datasets.tsv",sep="\t",header=TRUE,quote="",fill=TRUE)
id_datasets <- setNames(Datasets.df$id,Datasets.df$name)


skw <- readRDS("./data/GD/gene_essentiality.skewness.rds")
skw2 <- melt(skw)

gene_essentiality.df <- data.frame("id"=1:nrow(skw2),
                                   "id_genes"=id_genes[as.character(skw2$Var1)],
                                   "id_datasets"=id_datasets[as.character(skw2$Var2)],
                                   "skewness"=skw2$value)

write.table(gene_essentiality.df, file="./data/GD/gene_essentiality.tsv",sep="\t",
            col.names = TRUE,row.names = FALSE,quote=FALSE)
