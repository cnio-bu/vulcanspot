#!/usr/bin/env Rscript

### TITLE : Define Genetic Alterations
### AUTHOR : Perales-Paton, Javier - jperales@cnio.es
### LICENSE : GPL-v3

options(stringsAsFactors = FALSE)

# Functions
readCSV <- function(file,col.symbols=1,nrows=-1) {
  DF <- read.table(file=file,header=TRUE,sep=",",check.names = FALSE,nrows=nrows) 
  rownames(DF) <- DF[,col.symbols]
  DF <- DF[,-col.symbols]
  return(DF)
}

# Libraries
suppressPackageStartupMessages(require(optparse))

# 1 Get parameters ############
option_list = list(
  make_option(c("--GENEFUNCTION"), action="store", type='character',
              default = "./data/CCL/geneFunc.rds",
              help="GoF/LoF matrix in a RDS file for ccl in the screening."),
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

GF.fl <- opt$GENEFUNCTION # GF.fl <- "./data/CCL/geneFunc.rds"
fls <- as.vector(sapply(opt$DEPMATRICES,function(z) strsplit(z,split=":")[[1]]))

### 2 Read the Gene Funcionality Data (by GoF/LoF) #####
cat("[INFO] Reading Gene Function status matrix (GoF/LoF)...\n",file=stdout())
GF <- readRDS(GF.fl)

### 3 Create a library of Cell Sets based on Genetic alterations ##### 
cat("[INFO] Creating a library of Cell Sets that harbor alterations in Gene A...\n",file=stdout())

# First, remove all genes without any alterations.
# This is a easy step, since we start with a signed binary matrix,
# We just have to remove those genes with all 0s in a row
  # > table(rowSums(abs(GF)) == 0)
  # 
  # FALSE  TRUE 
  # 19692 15725
  # > table(rowSums(GF==0) != ncol(GF))
  # 
  # FALSE  TRUE 
  # 15725 19692 
GF <- GF[rowSums(GF==0) != ncol(GF),]

# Let start with the definition of PanCancer cell sets
GoF <- apply(GF,1,function(z) colnames(GF)[z == +1])
names(GoF) <- paste0("PANCANCER","::","GoF","::",names(GoF))

LoF <- apply(GF,1,function(z) colnames(GF)[z == -1])
names(LoF) <- paste0("PANCANCER","::","LoF","::",names(LoF))

CellSets <- c(GoF,LoF)
# Now, we will create subsets of these PanCancer cell sets based on cancer type
i <- 0;
CellSets.cancertype <- vector(mode = "list")
for(idx in names(CellSets)) {
  i <- i+1;
  if( (i %% 1000)==0 | (i==length(CellSets)) ) cat(paste0("#",i,"/",length(CellSets),"\n"),file=stdout())
  
  CellSet <- CellSets[[idx]]
  if(length(CellSet)==0) next;
  
  # Unique tissues across cells
  tissues <- sort(unique(sapply(CellSet,function(z) paste(strsplit(z,split="_")[[1]][-1],collapse = "_"))))
  for(tissue in tissues) {
    newidx <- gsub("PANCANCER::",paste0(tissue,"::"),idx)
    tmp <- list(grep(paste0("_",tissue,"$"),CellSet,value=TRUE))
    names(tmp)[1] <- newidx
    CellSets.cancertype <- c(CellSets.cancertype,tmp)
    rm(tmp)
  }
  rm(CellSet,tissues)
}

CellSets <- c(CellSets,CellSets.cancertype)
rm(CellSets.cancertype)

### 4 Subset these alterations based on the cell lines which were screenned in LoF assays ######
cat("[INFO]: Reading matrices from DEPMAP\n",file=stdout())
ccls <- vector()
for(fl in fls) {
  cat(paste0("[INFO]: \tReading File:",fl,"\n"),file=stdout())
  if(grepl("\\.rds$",fl)) {
    mat <- readRDS(fl)
  } else if(grepl("\\.csv",fl)) {
    mat <- readCSV(fl,nrows = 1)
  }
  stopifnot(is.matrix(mat) | is.data.frame(mat))
  
  # CellSets2 <- lapply(CellSets2,function(ccl) intersect(ccl,colnames(mat)))
  ccls <- unique(c(ccls,colnames(mat)))
}

CellSets2 <- lapply(CellSets,function(ccl) intersect(ccl,ccls))

### 5 Save files ###############
cat("[INFO]: Saving Cell Sets of Genetic Alterations in './data/CCL/'\n",file=stdout())
saveRDS(CellSets, file="./data/CCL/cellsets.genetics.rds")
saveRDS(CellSets2, file="./data/CCL/cellsets.anydep.rds")

