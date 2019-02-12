
suppressPackageStartupMessages(require(GenomicFeatures)) # used to create Genomic Intervals objects and get distances between genes.
suppressPackageStartupMessages(require(TxDb.Hsapiens.UCSC.hg19.knownGene)) # human genes annotation: genomic regions
suppressPackageStartupMessages(require(AnnotationDbi)) # SQL query
suppressPackageStartupMessages(require(org.Hs.eg.db)) # To map UCSC ids to Gene Symbols (HGNC)

# Create a GenomicRegions (GR) object with the human genes from UCSC
GR_genes <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
seqnames(GR_genes)
# Get a dictionary of UCSC ids to gene symbols
gene_symbol<- AnnotationDbi::select(org.Hs.eg.db, keys=GR_genes$gene_id,columns="SYMBOL", keytype="ENTREZID")
# Sanity check
all.equal(GR_genes$gene_id, gene_symbol$ENTREZID)
# Rename GR
names(GR_genes) <- gene_symbol$SYMBOL

# Create the same for EnsEMBL:
txdb <- makeTxDbFromGFF("./data/Annotation/Homo_sapiens.GRCh37.75.gtf",format="gtf",dbxrefTag=)
ENG <- read.table("./data/Annotation/ANN_hgnc_feb2014_GRCh37.tsv",sep="\t",header=TRUE)
dic2 <- setNames(ENG$hgnc_symbol,ENG$ensembl_gene_id)
GR_genes2 <- genes(txdb)
names(GR_genes2) <- dic2[names(GR_genes2)]
GR_genes2 <- GR_genes2[names(GR_genes2)!=""]

seqlevels(GR_genes2) <- sapply(seqlevels(GR_genes2),function(Gseq) ifelse(Gseq %in% c(1:23,"X","Y"),paste0("chr",Gseq),Gseq))
seqlevels(GR_genes2)[which(seqlevels(GR_genes2)=="MT")] <- "chrM"
seqnames(GR_genes2) <- sapply(seqnames(GR_genes2),function(Gseq) ifelse(Gseq %in% c(1:23,"X","Y"),paste0("chr",Gseq),Gseq))
seqnames(GR_genes2)[which(seqnames(GR_genes2)=="MT")] <- "chrM"

# Combine both
GR_genes <- c(GR_genes,GR_genes2[names(GR_genes2)[!names(GR_genes2) %in% names(GR_genes)]])

saveRDS(GR_genes,"./data/Annotation/hg19_GRCh37_genomicRegions.rds")
