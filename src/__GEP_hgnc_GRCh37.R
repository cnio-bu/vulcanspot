#!/usr/bin/env Rscript

### TITLE : Mapping ensembl_gene_id to hgnc_symbol for read counts
### AUTHOR : Perales-Paton, Javier
### LICENSE : GPL-v3
### DESCRIPTION : read counts on gene symbols will be used for multiple tools and further
###             steps in the study.

options(stringsAsFactors = FALSE)

# 1 Load gene metadata ######
cat("[INFO] Loading HUGO symbols. It must be used to retrieve genes of interest.\n",file=stdout())
ENG <- read.table("./data/Annotation/ANN_hgnc_feb2014_GRCh37.tsv",sep="\t",header=TRUE)
# Create a dictionary of HGNC symbols by ensembl_gene_id
ENG2SYM <- setNames(ENG$hgnc_symbol,ENG$ensembl_gene_id)


# 2 Load Gene Expr data ####
# If the Gene Expression Profile (GEP) file does not exist, then download it.
if(!file.exists("./data/CCL/CCLE_DepMap_18Q2_RNAseq_reads_20180502.gct")) {
  cat("[INFO] GEP file from CCLE does not exist. Downloading it from CCLE data portal FTP repository...\n",
      file=stdout())
  download.file(destfile = "./data/CCL/CCLE_DepMap_18Q2_RNAseq_reads_20180502.gct",
                url = "https://data.broadinstitute.org/ccle/CCLE_DepMap_18Q2_RNAseq_reads_20180502.gct")
}

# Load GEP
cat("[INFO] Loading Gene Expression raw data.\n",file=stdout())
CNT <- read.table(file = "./data/CCL/CCLE_DepMap_18Q2_RNAseq_reads_20180502.gct",
                          sep="\t",skip = 2,
                          header=TRUE,check.names = FALSE)
CNT$Name <- gsub("\\.[0-9]+","",CNT$Name)

# 3 Subset to only those en HGNC
ENG <- ENG[!is.na(ENG$hgnc_id),] # Remove them
# Retrieve
cat("[INFO] Only HGNC genes were retrieved:\n", file=stdout())
table(CNT$Name %in% ENG$ensembl_gene_id)
CNT <- CNT[CNT$Name %in% ENG$ensembl_gene_id,]

# NOTE: The CCLE gene expression data is huge: the quantifiation of gene expression
#   was done over 55k features (genes). The point is that it was included snoRNAs, snRNAs,
#   MicroRNAs and so on. Conventional RNAseq is not suitable for
#   these types of transcripts. Therefore, most of these "features" are zero-inflated
#   in the matrix of gene expression.
#   To avoid the zero-inflated bias in the UPC/RPKM cutoff, we will focus only in HGNC genes.
#   In fact, the metabolic model is built using ONLY:
#         - "protein_coding"
#         - "processed_transcript",
#         - "polymorphic_pseudogene"
#         - "pseudogene"

# # Gene biotypes that are included in the metabolic model
# gene_biotype.model <- c("protein_coding","processed_transcript",
#                       "polymorphic_pseudogene","pseudogene")

## Remove RNAs not captured by conventional RNAseq
# gene_biotype.out <- c("miRNA","snRNA","snoRNA")

# ENG <- ENG[!ENG$gene_biotype %in% gene_biotype.out,]

# Rename the hgnc symbol based on ensembl_gene_id
      # table(ENG2SYM[CNT$Name] == CNT$Description)
      # FALSE  TRUE 
      # 447 32570
CNT$Description <- ENG2SYM[CNT$Name]

# There are a few duplicated gene symbols
cat("[WARN]\t There are a few duplicated gene symbols in the matrix:\n", file=stdout())
table(duplicated(CNT$Description))

# Because gene duplication, we will aggregate their values. We will use the function 'sum' because
# read quantification usually add all reads mapping on the exons of a entity (gene).
# During aggregation, we will remove the first column that is the "ensembl_gene_id".
CNT$Description <- as.factor(CNT$Description)
CNT <- aggregate(. ~ Description, data = CNT[,-1], sum)
rownames(CNT) <- CNT$Description
CNT <- CNT[,colnames(CNT)!="Description"]

# Write it
write.table(cbind("hgnc_symbol"=rownames(CNT),CNT),
            file=paste0("./data/CCL/","GEP_hgnc_symbol_feb2014_GRCh37.reads.tsv"),
            sep = "\t",quote=FALSE,row.names = FALSE,col.names = TRUE)

cat("[INFO] DONE.\n\n",file=stdout())
sessionInfo()
