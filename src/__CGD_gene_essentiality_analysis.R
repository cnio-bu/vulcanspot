#!/usr/bin/env Rscript

### TITLE : Analysis on Gene Essentiality
### AUTHOR : Perales-Paton, Javier - jperales@cnio.es
### LICENSE : GPL-v3
### DESCRIPTION : This script

options(stringsAsFactors = FALSE)


suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(e1071)) # For skewness measurement calculation


option_list = list(
  make_option(c("--DSNAMES"), action="store",
              default="CRISPR:RNAi",
              type='character',
              help="DataSet Names for the corresponding matrices of gene dependencies"),
  make_option(c("--DEPMATRICES"), action="store",
              default="./data/CCL/portal-Avana-2018-05-10.csv:./data/CCL/portal-RNAi_merged-2018-05-10.csv",
              type='character',
              help="colon-separated string with the paths to matrices files with all the contexts of the database.")
)

opt = parse_args(OptionParser(option_list=option_list))

# Print Parameters
cat("PARAMETERS:\n",file=stdout())
for(opt.idx in 1:length(opt)) {
  cat(paste0(names(opt)[opt.idx],"\t",opt[[opt.idx]],"\n"),file = stdout())
}
cat("---\n",file=stdout())


DS_names <- as.vector(sapply(opt$DSNAMES,function(z) strsplit(z,split=":")[[1]]))
fls <- as.vector(sapply(opt$DEPMATRICES,function(z) strsplit(z,split=":")[[1]]))


# 2 Define functions ######
# Project Achilles proposed GCT format to store the Gene Essentiality Scores as a matrix.
readGCT <- function(file,col.symbols=1) {
  DF <- read.table(file = file,header=TRUE,skip=2,sep="\t",check.names = FALSE)
  rownames(DF) <- DF[,col.symbols]
  DF <- DF[,-c(1,2)]
  return(DF)
}

# Project DEPMAP proposed CSV format to store the Gene Essentiality Scores as a matrix.
readCSV <- function(file,col.symbols=1) {
  DF <- read.table(file=file,header=TRUE,sep=",",check.names = FALSE)  
  rownames(DF) <- DF[,col.symbols]
  DF <- DF[,-col.symbols]
  return(DF)
}

# 3 Load matrices #########
GF.list <- setNames(vector("list",length=length(fls)),DS_names)

for(idx in 1:length(fls)) {
  GS.fl <- fls[idx]
  cat("[INFO] Reading Gene Solutions matrix...\n",file=stdout())
  if(grepl("\\.gct",basename(GS.fl))) {
    GS <- readGCT(GS.fl)
  } else if (grepl("\\.csv",basename(GS.fl))) {
    GS <- readCSV(GS.fl)
  } else {
    stop("[ERROR] #2 : GS.fl could not be read.")
  }
  
  GF.list[[DS_names[idx]]] <- GS
}

# 4 Universe of genes #####################
genes <- sort(unique(unlist(lapply(GF.list,function(mat) rownames(mat)))))

# 5 Analysis on Gene Solution distributions #########
## 5.1 Skewness #####
skw.mat <- matrix(NA,ncol=length(GF.list),nrow=length(genes),
                  dimnames=list(genes,names(GF.list)))

for(DS in colnames(skw.mat)) {
  skw <- apply(GF.list[[DS]],1,function(z) skewness(z,na.rm=TRUE))
  skw.mat[names(skw),DS] <- skw
}

## 5.2 Top 10 quantile lethal genes in 90% cell lines
lethal.mat <- matrix(NA,ncol=length(GF.list),nrow=length(genes),
                     dimnames=list(genes,names(GF.list)))

for(DS in colnames(lethal.mat)) {
  depl01 <- apply(GF.list[[DS]],2,function(z) quantile(z,0.1,na.rm=TRUE)) # thresholds for individual top 10 quantile lethals
  depl01 <- apply(GF.list[[DS]],2,function(z) quantile(z,0.1,na.rm=TRUE)) # thresholds for individual top 10 quantile lethals
  potentialPAN <- t(apply(GF.list[[DS]],1,function(dep) dep <= depl01))   # Logical. Wether the gene is lethal in the individual CCL
  Ntested <- apply(GF.list[[DS]],1,function(dep) sum(!is.na(dep)))        # How many CCL were tested for that gene
  # rowSums(potentialPAN, na.rm=TRUE) >= ( ncol(potentialPAN)*0.9 )
  propLethal <- rowSums(potentialPAN, na.rm=TRUE) / Ntested               # Proportion of CCL where the gene is lethal
  lethal.mat[names(propLethal),DS] <- propLethal                          # store it
}

# cat(na.omit(rownames(lethal.mat)[lethal.mat[,1] >= 0.9]),sep="\n",file="./TMP/mypan.txt")

## 5.3 median
med.mat <- matrix(NA,ncol=length(GF.list),nrow=length(genes),
                  dimnames=list(genes,names(GF.list)))

for(DS in colnames(med.mat)) {
  med <- apply(GF.list[[DS]],1,function(z) median(z,na.rm=TRUE))
  med.mat[names(med),DS] <- med
}

# 6 Save it ####
saveRDS(skw.mat,"./data/GD/gene_essentiality.skewness.rds")
saveRDS(lethal.mat,"./data/GD/gene_essentiality.commonlethal.rds")
saveRDS(med.mat,"./data/GD/gene_essentiality.median.rds")
