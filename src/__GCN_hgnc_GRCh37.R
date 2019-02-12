#!/usr/bin/env Rscript

### TITLE : Get log2(CN ratio) mapped on 2014 hgnc symbols.
### AUTHOR : Perales-Paton, Javier - jperales@cnio.es
### LICENSE : GPL-v3
### DESCRIPTION : Create the log2(CN ratio) matrix of Genes X cancer cell lines by mapping hgnc_symbols from 2012 to 2014.

options(stringsAsFactors = FALSE)

# 1 Load gene metadata ######
cat("[INFO] Loading HUGO symbols. It must be used to retrieve genes of interest.\n",file=stdout())
ENG <- read.table("./data/Annotation/ANN_may2012-feb2014_GRCh37.tsv",sep="\t",header=TRUE)
# Create a dictionary of HGNC symbols from 2014 by HGNC symbols from 2012
SYM2012_to_SYM2014 <- setNames(ENG$hgnc_symbol.feb2014,ENG$hgnc_symbol.may2012)
SYM2012_to_SYM2014[SYM2012_to_SYM2014==""] <- NA

# 2. Load the CN table from CCLE data portal ######
if(!file.exists("./data/CCL/CCLE_copynumber_byGene_2013-12-03.txt")) {
  cat("[INFO] GCN file from CCLE does not exist. Downloading it from CCLE data portal FTP repository...\n",
      file=stdout())
  download.file(destfile = "./data/CCL/CCLE_copynumber_byGene_2013-12-03.txt",
                url="https://data.broadinstitute.org/ccle_legacy_data/dna_copy_number/CCLE_copynumber_byGene_2013-12-03.txt")
}

cat("[INFO] Load the CN table from CCLE data portal.\n",file=stdout())
GCN <- read.table(paste0("./data/CCL/CCLE_copynumber_byGene_2013-12-03.txt"),
                  sep="\t",header=TRUE,check.names = FALSE)

GCN$SYMBOL14 <- GCN$SYMBOL
GCN$SYMBOL14[!is.na(SYM2012_to_SYM2014[GCN$SYMBOL])] <- SYM2012_to_SYM2014[GCN$SYMBOL][!is.na(SYM2012_to_SYM2014[GCN$SYMBOL])]

# There are just 26 duplicated gene symbols, but to be really careful
for(dup.gene in GCN$SYMBOL14[duplicated(GCN$SYMBOL14)]) {
  old.gene.ids <- which(GCN$SYMBOL14==dup.gene)
  old.uniq.id <- old.gene.ids[which(GCN$SYMBOL[old.gene.ids]!=dup.gene)]
  GCN$SYMBOL14[old.uniq.id] <- GCN$SYMBOL[old.uniq.id]
}
# MiRNAs come up and down with the same identifier several times during the years
for(dup.gene in GCN$SYMBOL14[duplicated(GCN$SYMBOL14)]) {
  old.gene.ids <- which(GCN$SYMBOL14==dup.gene)
  old.uniq.id <- old.gene.ids[which(GCN$SYMBOL[old.gene.ids]!=dup.gene)]
  GCN$SYMBOL14[old.uniq.id] <- GCN$SYMBOL[old.uniq.id]
}

if(any(duplicated(GCN$SYMBOL14))) {
  stop("[ERROR] 1 : Mapping came up with duplications.")
} else {
  rownames(GCN) <- GCN$SYMBOL14
  # REMOVE GENES FROM THE chrY
  GCN <- GCN[which(GCN$CHR!="Y"),]
  # Dropped-out the columns we are not interested in
  GCN <- GCN[!colnames(GCN) %in% c("EGID","SYMBOL","SYMBOL14","CHR","CHRLOC","CHRLOCEND"),]
}

MAP.res <- GCN[,c("SYMBOL","SYMBOL14")]
colnames(MAP.res) <- c("original.hgnc_symbol","actual.hgnc_symbol")
write.table(MAP.res,file="./data/CCL/GCN_2012-2014mapping.log",
            sep="\t",row.names = FALSE,col.names=FALSE,quote=FALSE)

# 3. write it ######
GCN <- GCN[,!colnames(GCN) %in% c("EGID","SYMBOL","CHR","CHRLOC","CHRLOCEND","SYMBOL14")]
write.table(cbind("hgnc_symbol"=rownames(GCN),GCN),
            file="./data/CCL/GCN_hgnc_symbol_feb2014_GRCh37.log2CNratio.tsv",
            sep="\t",quote = FALSE,row.names = FALSE,col.names = TRUE)

cat("[INFO] DONE.\n\n",file=stdout())
sessionInfo()
