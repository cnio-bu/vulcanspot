#!/usr/bin/env Rscript

### TITLE : Transform reads to Universal exPression Codes only for HGNC genes
### AUTHOR : Perales-Paton, Javier - jperales@cnio.es
### LICENSE : GPL-v3
### DESCRIPTION : Transform Gene Expression Profiles (GEP, read counts) into
###                 Universal exPression Code (Transcriptionally Active scale, 0-1).

options(stringsAsFactors = FALSE)

# 1 Load libraries #######
suppressPackageStartupMessages(require(SCAN.UPC))

# 3 Load gene metadata ######
cat(paste0("[INFO] Loading HUGO symbols, GC content and transcript length.",
           " It must be used for two purposes:","\n",
           "\t To retrieve genes of interest.","\n",
           "\t To transform into UPC values, which requires GC and transcript length","\n"),
    file=stdout())
ENG <- read.table("./data/Annotation/ANN_hgnc_feb2014_GRCh37.tsv",sep="\t",header=TRUE)

# We realized that there are a few duplicated gene symbols, with different ensembl_gene_id
#   that do not have length and GC metadata. So we retrieve only those w/ ensembl_gene_id
ENG <- ENG[grepl("ENSG",ENG$ensembl_gene_id),]

HGNC.len <- setNames(ENG$exon_length,ENG$hgnc_symbol)
HGNC.GC <- setNames(ENG$percentage_gc_content,ENG$hgnc_symbol)

# 4 Load Gene Expr data ####
cat("Loading Gene Expression raw data.\n",file=stdout())
CNT <- read.table(file=paste0("./data/CCL/","GEP_hgnc_symbol_feb2014_GRCh37",".reads.tsv"),
                  sep="\t",header=TRUE,check.names = FALSE)
rownames(CNT) <- CNT$hgnc_symbol
CNT <- CNT[,colnames(CNT)!="hgnc_symbol"]

# 5 Sanity check: check if all gene symbols in the matrix have length info.
stopifnot(all(rownames(CNT) %in% ENG$hgnc_symbol))

# 6 UPC transformation
cat("UPC transformation...\n",file=stdout())
UPC.mat <- apply(CNT,2,function(z) {
  UPC_RNASeq_Single(z,rownames(CNT),
                    # ignoreZeroes=TRUE,
                    lengths = HGNC.len[rownames(CNT)],
                    gcContent = HGNC.GC[rownames(CNT)],
                    verbose=FALSE);
})
# How bad... verbose=FALSE but it still shows verbose...

# Write it
write.table(cbind("hgnc_symbol"=rownames(UPC.mat),UPC.mat),
            file=paste0("./data/CCL/","GEP_hgnc_symbol_feb2014_GRCh37",".UPC.tsv"),
            sep = "\t",quote=FALSE,row.names = FALSE,col.names = TRUE)

cat("DONE.\n\n",file=stdout())
sessionInfo()
