#!/usr/bin/env Rscript

### TITLE : Cancer Gene dependency by Kolmogorov-Smirnov test
### AUTHOR : Perales-Paton, Javier - jperales@cnio.es
### LICENSE : GPL-v3
### DESCRIPTION : This script

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(data.table))
suppressPackageStartupMessages(require(fgsea)) # For kolgomorov-smirnov testing
suppressPackageStartupMessages(require(e1071)) # For skewness measurement calculation


# 1 Get parameters ############
option_list = list(
  make_option(c("--ESSENTIALITY"), action="store", type='character',
              help="File containing scores for decrease cell fitness by LoF-gene screening."),
  make_option(c("--GENEFUNCTION"), action="store", type='character',
              help="GoF/LoF matrix in a RDS file for ccl in the screening."),
  make_option(c("--SCALEESSENTIALLITY"), action="store",default = TRUE, type='logical',
              help="Logical whether scaling the Gene Essentiality scores."),
  make_option(c("--MINSIZE"), action="store", default=5, type='numeric',
              help="Cell line from the core of CMap-L1000 core"),
  make_option(c("--MAXSIZE"), action="store", default=150, type='numeric',
              help="Cell line from the core of CMap-L1000 core"),
  make_option(c("--DATASETNAME"), action="store", default="TS2017", type='character',
              help="Dataset name for the Gene Essentiality dataset screening."),
  make_option(c("--NPERMUT"), action="store", default=10000, type='numeric',
              help="Number of permutations for fgsea."),
  make_option(c("--NPROC"), action="store", default=10, type='numeric',
              help="Number of processes to be used for fgsea."),
  make_option(c("--minimumSKEWNESS"), action="store", default=NULL, type='numeric',
              help="Minimum SKEWNESS measure of the distribution to be considered for statistical testing"),
  make_option(c("--minimumAVGscore"), action="store", default=NULL, type='numeric',
              help="Minimum Average Essentiality Score to be considered for statistical testing"),
  make_option(c("--minimumGDscore"), action="store", default=NULL, type='numeric',
              help="Minimum GD Essentiality Score to be considered for statistical testing"),
  make_option(c("--OUTDIR"), action="store", default="./data/GD", type='character',
              help="Output directory for the results")
)

opt = parse_args(OptionParser(option_list=option_list))

# Print Parameters
cat("PARAMETERS:\n",file=stdout())
for(opt.idx in 1:length(opt)) {
  cat(paste0(names(opt)[opt.idx],"\t",opt[[opt.idx]],"\n"),file = stdout())
}
cat("---\n",file=stdout())

GS.fl <- opt$ESSENTIALITY # GS.fl <- "./data/CCL/portal-RNAi_merged-2018-05-10.csv"
GF.fl <- opt$GENEFUNCTION # GF.fl <- "./data/CCL/geneFunc.rds"

min.size <- opt$MINSIZE # min.size <- 5
max.size <- opt$MAXSIZE # max.size <- 150

DS.name <- opt$DATASETNAME

RNKscale <- opt$SCALEESSENTIALLITY

npermut <- opt$NPERMUT
nproc <- opt$NPROC

SKEWNESS_min <- opt$minimumSKEWNESS
AVGscore_min <- opt$minimumAVGscore
GDscore_min <- opt$minimumGDscore
outdir <- opt$OUTDIR

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

# This is a the core function of this script.
# It will get a GS matrix (Gene Essentiality Score)
performFGSEA <- function(CellSets,RNK,geneB,RNKscale=TRUE,
                         min.size,max.size,nperm=10000,nproc=1,
                         AVGscore_min=NULL,GDscore_min=NULL) {
  
  # Average Gene Essentiality Score from the raw data for each Cell Set
  #   NOTE: We take the avg BEFORE scaling the data. This is important because
  #   it allows us to rank the significant GDs based on the observation of 
  #   how much decrease cell viability in the context
  avgGES.cellsets <- unlist(lapply(CellSets,function(cells) mean(RNK[intersect(names(RNK),cells)])))
  
  # Calculate the GSscore as the difference between the median Gene Essentiality Score between the two groups
  #   NOTE: We take the avg BEFORE scaling the data. This is important because
  #   it allows us to rank the significant GDs based on the observation of 
  #   how much decrease cell viability in the context
  GDscore.cellsets <- unlist(lapply(CellSets,function(cells) {
    -1*(median(RNK[intersect(names(RNK),cells)],na.rm = TRUE) - median(RNK[setdiff(names(RNK),cells)],na.rm = TRUE))
    }))
  
  # Calculate skewness of the gene B dependency
  #   NOTE: We take the avg BEFORE scaling the data. This is important because
  #   it allows us to rank the significant GDs based on the observation of 
  #   how much decrease cell viability in the context
  skewness.gene <- skewness(RNK,na.rm = TRUE)
  
  # Scale ranking using a z-score
  if(RNKscale) RNK <- scale(RNK)[,1];
  
  # Remove CellSets that are not relevant to avoid multiple testing issues
  if(!is.null(AVGscore_min)) {
    CellSets <- CellSets[avgGES.cellsets <= AVGscore_min & !is.na(avgGES.cellsets)]
    avgGES.cellsets <- avgGES.cellsets[avgGES.cellsets <= AVGscore_min & !is.na(avgGES.cellsets)]
  }
  
  if(!is.null(GDscore_min)) {
    CellSets <- CellSets[GDscore.cellsets >= GDscore_min & !is.na(GDscore.cellsets)]
    GDscore.cellsets <- GDscore.cellsets[GDscore.cellsets >= GDscore_min & !is.na(GDscore.cellsets)]
  }
 
  # Perform the KS test along all the CellSets
  res <- fgsea(pathways = CellSets,stats = RNK,
               minSize = min.size,maxSize = max.size,
               nperm = nperm,nproc=nproc)
  
  # Unfold the metadata from the cell sets for annotation
  CS.metainfo <- sapply(res$pathway,function(z) strsplit(z,split="::"))
  
  # Add metadata for the records
  res$GeneA <- unlist(lapply(CS.metainfo,function(z) z[[3]]))
  res$GeneB <- geneB
  res$Context <- unlist(lapply(CS.metainfo,function(z) z[[1]]))
  res$GeneAFunction <- unlist(lapply(CS.metainfo,function(z) z[[2]]))
  res$avgGES <- avgGES.cellsets[res$pathway]
  res$GDscore <- GDscore.cellsets[res$pathway]
  res$freq <- unlist(res$size)/length(RNK)
  res$skewness <- skewness.gene
  
  # Calculate the GScore based on LeadingEdge
  # NOTE: Discarded bc some Cell Sets are significant just with one cell line in the leading edge which
  # inflate so far the outliers with a very high GDscore
  # res$GDscore2 <- sapply(res$pathway,function(cs_name,RNK,CellSets) {
  #   leadingEdge_CCLs <- unlist(res[which(res$pathway==cs_name),"leadingEdge"])
  #   all_CCL <- CellSets[[cs_name]]
  #   
  #   -1*(median(RNK[leadingEdge_CCLs]) - median(RNK[setdiff(names(RNK),all_CCL)],na.rm = TRUE))
  #   
  # },RNK=RNK,CellSets=CellSets)
  
  return(res)
}



### 3 Read the Gene Essentiality data #####
# GS stands for Gene Solution
if(!file.exists(GS.fl)) {
  stop("[ERROR] #1 GS.fl does not exist.")
}

cat("[INFO] Reading Gene Solutions matrix...\n",file=stdout())
if(grepl("\\.gct",basename(GS.fl))) {
  GS <- readGCT(GS.fl)
} else if (grepl("\\.csv",basename(GS.fl))) {
  GS <- readCSV(GS.fl)
} else {
  stop("[ERROR] #2 : GS.fl could not be read.")
}

### 4 Read the Gene Funcionality Data (by GoF/LoF) #####
cat("[INFO] Reading Gene Function status matrix (GoF/LoF)...\n",file=stdout())
GF <- readRDS(GF.fl)

### 5 Read the Gene Roles in cancer based on Cancer Gene Census ####
CGC <- read.table("./data/Annotation/Census_allWed_Oct_4_09-53-19_2017.tsv",sep="\t",header = TRUE,fill=TRUE)
CGC.Role <- setNames(sapply(gsub(" +","",CGC$Role.in.Cancer),function(z) strsplit(z,split=",")[[1]]),CGC$Gene.Symbol)
ONC <- names(CGC.Role)[unlist(lapply(CGC.Role,function(z) any(z=="oncogene") & !any(z=="TSG")))]
TSG <- names(CGC.Role)[unlist(lapply(CGC.Role,function(z) !any(z=="oncogene") & any(z=="TSG")))]
AMB <- names(CGC.Role)[unlist(lapply(CGC.Role,function(z) any(z=="oncogene") & any(z=="TSG")))]

### 6 Create library of Cell Sets ##### 
# Process the GoF/LoF data to create a library of Cell Sets where those cell lines
#   harbor an alteration (GoF/LoF) on a Gene A. This library is going to be used
#   to perform Kolmogorov-Smirnov tests.

cat("[INFO] Creating a library of Cell Sets that harbor alterations in Gene A...\n",file=stdout())

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

### 7 Avoid Massive Multiple Testing (order of magniture of millions of tests) #####
# NOTE: We decided to use an emperical approach to avoid multiple testing. The source
# is the fact that we test for every GeneA-GeneB pair, so as long as the datasets grows
# (both, molecular profiling or gene essentiality), more a more massive multiple testing
# are added. Instead of empowering the statistical test association, we have an exponential
# increase of the number of test and we loss associations when adjusting for multiple testing.
# The addiction of hundres of cell lines in the new releases of the datasets have promote this
# event in the analysis before.
# To avoid this, we follow two criteria:
#   7.1 We only consider Gene Dependencies on recurrently altered genes in cancer considering
#     a minimum and maximum sample size.
#   7.2 [optional] We only consider Gene Dependencies for genes that present a skewed distribution towards
#     dependency.
#   7.3 [optional] We only consider Gene Dependencies for genes that clearly show a functional essentiality
#     in the data. That is to say, the cell lines that harbor the alteration on gene A must
#     present a tendency to essentiality to geneB. This is not the statistical test per-se, but
#     it is for removing genes with flat patterns on essentiality (that are the vast majority of
#     the cases).
# Later, we are going to statistical test whether cell lines that harbor alteration in Gene A
# are the cell lines with the STRONGEST essentiality to the function of the Gene B.

### 7.1 Remove Non-recurrently altered genes from the statistical testing ######
# Retrieve only those Cell lines that are present in the actual screening
CellSets <- lapply(CellSets, function(CCL) CCL[CCL%in%colnames(GS)])
# Remove those CellSets that are empty (length==0)
CellSets <- CellSets[unlist(lapply(CellSets,length))!=0]

# Save all the CellSets defined for the analysis of this dataset
saveRDS(CellSets,
        paste0(outdir,"/contexts.",DS.name,".All_CellSets.rds"))

# # Remove CellSets that do not fit the criteria for the test
CS.pass <- unlist(lapply(CellSets,function(z) length(z) >= min.size & length(z) <= max.size))
CellSets <- CellSets[CS.pass]

# Remove Cell Sets where GeneA is confounding context
# uniq_cell_lines <- unique(sort(unlist(CellSets)))
# cnt_contexts <- table(sapply(uniq_cell_lines,function(z) paste(strsplit(z,split="_")[[1]][-1],collapse = "_")))
# 
# CS.pass <- sapply(names(CellSets),function(z) {
#   
# })

# Save the CellSets object used for the analysis of this dataset
saveRDS(CellSets,
        paste0(outdir,"/contexts.",DS.name,".CellSets.rds"))

### 7.2 Retrieve only those genes that present a skewed distribution on gene dependency scores ######

    # # Maybe for future implementations, instead of an absolute value of skewness, it would be interesting to set a
    # # relative threshold based on the within-dataset distribution of skewness.
    # if(!is.null(SKEWNESS_belowQT)) {
    #   cat("[INFO] Calculating skewness for a set quantile...\n",file=stdout())
    #   # cat(paste0("\t Skewness Quantile =",SKEWNESS_belowQT,"\n"),file=stdout())
    #   skw_tmp <- apply(GS,1,function(z) skewness(z,na.rm=TRUE))
    #   SKEWNESS_min <- quantile(skw_tmp,SKEWNESS_belowQT)
    #   cat(paste0("\t Skewness Quantile =",SKEWNESS_belowQT," => minimum Skewness measurement",SKEWNESS_min,"\n"),file=stdout())
    # }

if(!is.null(SKEWNESS_min)) {
  cat("[INFO] Retrieving only genes with a minimum skewness from their gene dependency score distribution...\n",file=stdout())
  cat(paste0("\t minimum Skewness measurement =",SKEWNESS_min,"\n"),file=stdout())
  
  skw <- apply(GS,1,function(z) skewness(z,na.rm=TRUE))
  
  cat(paste0("\t- Retrieved: ",sum(skw <= SKEWNESS_min),"\n"),file=stdout())
  cat(paste0("\t- Discarded: ",sum(skw > SKEWNESS_min),"\n"),file=stdout())
  GS <- GS[skw <= SKEWNESS_min,]
}

### 7.2 Remove those genes from the LoF screening that have flat patterns #####
if(!is.null(AVGscore_min)) {
  cat("[INFO] Retrieving only genes with a minimum AVG gene dependency on the cell sets...\n",file=stdout())
  GS_rel <- apply(GS,1, function(essentiality,avgGES) {
    sum(unlist(lapply(CellSets,function(cells,arg,avgGES) {
      RNK <- arg[cells]
      if(!all(is.na(RNK))) {
        mean(na.omit(RNK),na.rm = TRUE) <= avgGES;
      } else {
        FALSE
      }
    }, arg=essentiality,avgGES=avgGES)))
  }, avgGES=AVGscore_min)
  
  GS <- GS[GS_rel!=0,]  
}

### 8 Statistical test for GeneA-GeneB dependencies #####
cat("[INFO] Performing systematic pre-ranked GSEAs for all GeneA-GeneB pairs. This will take a long time...\n",file=stdout())

# To concatenate different data.tables, first we need to create the first one. 
# This first one is going to be empty
res <- data.table("pathway"=vector("character",length=0),
           "pval"=vector("numeric",length=0),
           "padj"=vector("numeric",length=0),
           "ES"=vector("numeric",length=0),
           "NES"=vector("numeric",length=0),
           "nMoreExtreme"=vector("numeric",length=0),
           "size"=vector("numeric",length=0),
           "leadingEdge"=vector("list",length=0),
           "GeneA"=vector("character",length=0),
           "GeneB"=vector("character",length=0),
           "Context"=vector("character",length=0),
           "GeneAFunction"=vector("character",length=0),
           "avgGES"=vector("numeric",length=0),
           "GDscore"=vector("numeric",length=0),
	   "freq"=vector("numeric",length=0),
	   "skewness"=vector("numeric",length=0))

i <- 0;
for(gene in rownames(GS)) {
  i <- i+1;
  if( (i %% 100)==0 | (i==nrow(GS)) ) cat(paste0("#",i,"/",nrow(GS),"\n"),file=stdout())
  # Sys.sleep(0.5)
  
  RNK <- sort(unlist(GS[gene,]),decreasing = TRUE)
  # Remove NAs values from the ranking. There are only NAs values in the case of shRNA screening due
  # to bad quality metrics or sparse data when combining different shRNA screening.
  RNK <- na.omit(RNK) 
  
  tmp <- performFGSEA(CellSets,               # List with all Cell Sets (stratified by Tissue + GeneA altered + type of Alteration)
                      RNK = RNK,geneB = gene, # Ranking and Gene B name for the Kolmogorov-smirnov test
                      RNKscale=RNKscale,      # Scale the Ranking of scores
                      min.size,max.size,      # Restrictions in the cell sets
                      npermut,                # No. permutations to get a p-value
                      nproc,                  # No. cores
                      AVGscore_min=AVGscore_min,
                      GDscore_min=GDscore_min)
  
  # If there is nrow==0, the test couldn't be performed becase of a lack of cell lines
  # in the screening (only in the case of datasets inflated with NA values, like McDonald2017)
  if(nrow(tmp)>0) {
    res <- rbind(res,tmp) 
  }
  
  rm(tmp,RNK)
}

### 9 Correct for multiple testing #####
cat(paste0("[INFO] Correcting for Multiple testing (method: fdr)\n"),file=stdout())
for(context in unique(res$Context)) {
  cat(paste0("[INFO] \tStarting with '",context,":'\n"),file=stdout())
  GeneA_byContext <- unique(res[which(res$Context==context),]$pathway)
  i <- 0
  for(GeneA in GeneA_byContext) {
    i <- i+1
    if(i%%100 == 0) {
      cat(paste0("[INFO]\t\t",i," out of ",
                 length(GeneA_byContext)," GeneA alterations\n"),
                 file=stdout())
    }
    
    idx <- which(res$pathway==GeneA & res$Context==context)
    res[idx,]$padj <- p.adjust(res[idx,]$pval, method="fdr")
    
  }
}

### 10 Add cancer roles to Gene A and Gene B #####
cat(paste0("[INFO] Add annotation to the table related to cancer gene roles\n"),file=stdout() );

res$GeneA.Role <- "unknown"
res$GeneA.Role[res$GeneA %in% ONC] <- "oncogene"
res$GeneA.Role[res$GeneA %in% TSG] <- "TSG"
res$GeneA.Role[res$GeneA %in% AMB] <- "oncogene/TSG"

res$GeneB.Role <- "unknown"
res$GeneB.Role[res$GeneB %in% ONC] <- "oncogene"
res$GeneB.Role[res$GeneB %in% TSG] <- "TSG"
res$GeneB.Role[res$GeneB %in% AMB] <- "oncogene/TSG"


### 11 Save results #####
cat(paste0("[INFO] Save results\n"),file=stdout())
saveRDS(res,
        paste0(outdir,"/GD.",DS.name,".rds"))

cat("[INFO] Finished\n",file=stdout())

sessionInfo()
