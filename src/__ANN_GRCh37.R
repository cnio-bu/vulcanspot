### TITLE : Gene features for GRCh37/hg19 last release
### AUTHOR : Perales-Paton, Javier - jperales@cnio.es
### LICENSE: GPL-v3
### DESCRIPTION : Different Omics data and analysis requires metadata to map biological entities
###             between databases and resources. This script will get this annotation for future
###             implementations.
###
###           All this data will be obtained from the ensembl feb2014 (release 75):
###             ensembl_gene_id
###             hgnc_symbol
###             hgnc_id
###             gene_biotype
###             HGNC_id
###             exon.length : sum of non-overlaping exon regions for a given gene (aka. gene length)
###             percentage_gc_content
###             uniprot_sptrembl
###             uniprot_genename
###             uniprot_swissprot_accession

options(stringsAsFactors = FALSE)

# Load library
library(biomaRt)

# GET Hugo Symbols for hg19 - Feb2014
  human = useMart(host="feb2014.archive.ensembl.org",
                  biomart="ENSEMBL_MART_ENSEMBL",
                  dataset="hsapiens_gene_ensembl")

# Get Hugo symbol and acc, gene biotype
  # Metadata related to hgnc
  ENG <- getBM(attributes = c("ensembl_gene_id","hgnc_symbol","hgnc_id","gene_biotype"),
               mart = human, uniqueRows=TRUE)
  ENG$HGNC_id <- sapply(ENG$hgnc_id,function(z) {if(!is.na(z)) paste0("HGNC:",as.character(z)) else z})
  # Metadata related to uniprot
  UNI <- getBM(attributes = c("ensembl_gene_id",
                              "uniprot_sptrembl","uniprot_genename",
                              "uniprot_swissprot_accession"),
               mart = human, uniqueRows=TRUE)
  UNI2 <- merge(x=ENG[,c("ensembl_gene_id","hgnc_symbol")],y=UNI,
                by.x="ensembl_gene_id",by.y="ensembl_gene_id",
                all.x=FALSE,all.y=TRUE)
# Get GC content for each gene
  GC <- getBM(attributes = c("ensembl_gene_id","percentage_gc_content"),
                      mart = human, uniqueRows=TRUE)

# Extract Total Non-Overlapping Exon Length Per Gene With Bioconductor
  # source: https://www.biostars.org/p/83901/
  # First, import the GTF-file that you have also used as input for htseq-count
  library(GenomicFeatures)
  txdb <- makeTxDbFromGFF("./data/Annotation/Homo_sapiens.GRCh37.75.gtf",format="gtf")
  # then collect the exons per gene id
  exons.list.per.gene <- exonsBy(txdb,by="gene")
  # then for each gene, reduce all the exons to a set of non overlapping exons, calculate their lengths (widths) and sum then
  exonic.gene.sizes <- lapply(exons.list.per.gene,function(x){sum(width(reduce(x)))})
  # as.data.frame
  stopifnot(all(names(exonic.gene.sizes)==unlist(names(exonic.gene.sizes))))
  Len <- data.frame("ensembl_gene_id"=names(exonic.gene.sizes),
                    "exon_length"=unlist(exonic.gene.sizes))

# Merge the three data.frames into one
  ENG2 <- merge(x = ENG,y=Len,by.x="ensembl_gene_id",by.y="ensembl_gene_id",all=TRUE)
  ENG3 <- merge(x = ENG2, y=GC, by.x="ensembl_gene_id",by.y="ensembl_gene_id",all=TRUE)
  
# Save it
write.table(ENG3,file = "./data/Annotation/ANN_hgnc_feb2014_GRCh37.tsv",
            sep = "\t",row.names = FALSE,col.names = TRUE,quote=FALSE)
write.table(UNI2,file = "./data/Annotation/ANN_uniprot_feb2014_GRCh37.tsv",
            sep = "\t",row.names = FALSE,col.names = TRUE,quote=FALSE)


cat("DONE.\n\n",file=stdout())
sessionInfo()
