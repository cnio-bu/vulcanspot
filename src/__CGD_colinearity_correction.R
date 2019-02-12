### TITLE : Adjustment for colinearity in cancer GDs
### AUTHOR : Perales-Paton, Javier - jperales@cnio.es
### LICENSE : GPL-v3
### DESCRIPTION : Exploratory analysis on Genetic Dependencies reveals colinearity of significant GDs 
###         due to co-occurrence of genetic alterations associated to the dependency of a Gene B.
###         The source of this colinearity is:
###           1. Gene Copy-number Alterations on Genomic regions (GCN on GR): Copy-Number alterations usually
###           spans large segment of DNA, including several neighbouring genes. When GDs are tested by pairs,
###           neighbour genes take similar statistics.
###           2. Co-occurrent Passenger alterations (passenger): Cancer typically show dependency to some driver
###           events, such as Gain of Function of BRAF in melanoma or GoF of ERBB2 en breast cancer. However,
###           other co-ocurrent events with these other events are going to be statistically significant with 
###           dependency of the driver event.
### Adjustment of these colinearity events are going to follow the assumption of a parsimonius solution.


options(stringsAsFactors = FALSE)

### 1 Load libraries #######
suppressPackageStartupMessages(require(optparse)) # Parsing arguments
suppressPackageStartupMessages(require(data.table)) # load data.tables (GD object)
suppressPackageStartupMessages(require(GenomicFeatures)) # used to create Genomic Intervals objects and get distances between genes.
# suppressPackageStartupMessages(require(TxDb.Hsapiens.UCSC.hg19.knownGene)) # human genes annotation: genomic regions
# suppressPackageStartupMessages(require(AnnotationDbi)) # SQL query
# suppressPackageStartupMessages(require(org.Hs.eg.db)) # To map UCSC ids to Gene Symbols (HGNC)

### 2 Get Parameters #######
option_list = list(
  make_option(c("--DATASETNAME"), action="store", default="RNAi", type='character',
              help="Dataset name for the Gene Essentiality dataset screening."),
  make_option(c("--OVERLAP"), action="store", default=0.7, type='numeric',
              help="Grade of Overlap (0,1) for Passenger alterations of driver events."),
  make_option(c("--JACKARD"), action="store", default=0.7, type='numeric',
              help="Jackard Index threshold for Genomic Regions partnerships."),
  make_option(c("--INDIR"), action="store", default="./data/GD", type='character',
              help="Input directory for the results"),
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

dataset_name <- opt$DATASETNAME # dataset_name <- "RNAi"
in_dir <- opt$INDIR         # in_dir <- "./data/GD"
out_dir <- opt$OUTDIR       # out_dir <- "./data/GD"
jIdx_cuttoff <- opt$JACKARD # jIdx_cuttoff <- 0.65
oIdx_cutoff <- opt$OVERLAP  # oIdx_cutoff <- 0.7

GD_fl <- paste0(in_dir,"/GD.",dataset_name,".rds")               # GD_fl <- "./data/GD/contexts.CRI.CellSets.rds"
CS_fl <- paste0(in_dir,"/contexts.",dataset_name,".CellSets.rds") # CS_fl <- "./data/GD/contexts.CRI.CellSets.rds"

# Sanity check
if(!file.exists(GD_fl)) {
  stop("[ERROR]: GD file does not exist.")
}
if(!file.exists(CS_fl)) {
  stop("[ERROR]: Cell Sets file does not exist.")
}

### 3 Prepare GenomicRegions object for human genes #####
cat(paste0("[INFO] :","Prepare a GenomicRegions object combining the UCSC and EnsEMBL datasets.","\n"),file=stdout())
GR_genes <- readRDS("./data/Annotation/hg19_GRCh37_genomicRegions.rds")

# Load cancer genes
cat(paste0("[INFO] :","Loading cancer genes","\n"),file=stdout())
CGC <- read.table("./data/Annotation/Census_allWed_Oct_4_09-53-19_2017.tsv",sep="\t",header=TRUE,quote="")
cancer_genes <- CGC$Gene.Symbol

Driver_data <- read.table("./data/Annotation/srep02650-s3.csv",sep=",",header=TRUE,quote="")
HighConf_drivers <- Driver_data[Driver_data$Putative.Driver.Category=="High Confidence Driver","Gene.Symbol"]

cancer_genes <- unique(c(cancer_genes,HighConf_drivers,"KLHL9","YAP1","MTAP"))

### 4 Function definition #######

# 'jackard_index' calculates the jarckard index for two vectors, where each vector is a set
# of Cancer Cell lines harboring an alteration in Gene A.
jackard_index <- function(set1,set2) {
  jidx <- length(intersect(set1,set2))/length(unique(c(set1,set2)))
  return(jidx)
}

# DESCRIPTION HERE
Overlap_index <- function(set1,set2) {
  IN_idx <- sum(set1 %in% set2)/length(set1) 
  return(IN_idx)
}

# 'pairwise_dist' provides a matrix of Gene Pair-wise distances along the genome.
# 'NA' values are obtained when two genes are in different chromosomes.
# The diagnoal of the matrix is always going to be 0 because it is the distance to gene itself.
# The numbers are base-pairs along the genome.
pairwise_dist <- function(cogenes,GR=GR_genes,verbose=FALSE,include.unk=FALSE) {
  cogenes_iset <-   intersect(cogenes,names(GR))
  cogenes_left <- setdiff(cogenes,names(GR))
  
  if(verbose) {
    if(length(cogenes_left)>0) {
      cat("[WARN] : Not all the genes are in the GR object. See verbose:\n",file=stderr())
      cat(paste0("\t",setdiff(cogenes,cogenes_iset),"\n"),file=stderr())
    } 
  }
  
  PW_dist <- matrix(NA,ncol=length(cogenes),nrow=length(cogenes),dimnames=list(cogenes,cogenes))
  
  if(length(cogenes_iset)>0) {
    subGR <- GR[cogenes_iset]
    cogenes_iset_dist <- sapply(names(subGR),function(gene_name) distance(subGR[gene_name],subGR,ignore.strand=TRUE))
    rownames(cogenes_iset_dist) <- colnames(cogenes_iset_dist)
    PW_dist[rownames(cogenes_iset_dist),colnames(cogenes_iset_dist)] <- cogenes_iset_dist
    
  }
  
  if(include.unk) {
    PW_dist[cogenes_left,] <- 0
    PW_dist[,cogenes_left] <- 0
  }

  return(PW_dist)
}


# DESCRIPTION HERE
stratifyByIdx <- function(CS_Idx,threshold,n,selectRows=NULL,avoid_itself=FALSE,avoid_triangle=NULL) {
  
  if(avoid_itself) diag(CS_Idx) <- -666
  
  if(!is.null(avoid_triangle)) {
    if(avoid_triangle=="upper") {
      CS_Idx[upper.tri(CS_Idx, diag = avoid_itself)] <- -666
    } else if(avoid_triangle=="lower") {
      CS_Idx[lower.tri(CS_Idx, diag = avoid_itself)] <- -666
    }
  }
  
  if(!is.null(selectRows)) CS_Idx <- CS_Idx[selectRows,]
    
  coOccur <- apply(CS_Idx,1, function(jIdx,jIdx_cuttoff) {
    jIdx <- sort(jIdx,decreasing = TRUE)
    jIdx[jIdx > jIdx_cuttoff]
  },jIdx_cuttoff=threshold)
  
  coOccur <- coOccur[unlist(lapply(coOccur,length))>n]
  
  
  return(coOccur)
}

# DESCRITION HERE
retrieveCS <- function(CS_list,feature_sel) {
  stopifnot(is.list(CS_list))
  stopifnot(is.vector(feature_sel))
  if(any(!feature_sel %in% c("Context","Alteration","GeneA"))) {
    stop("[ERROR]: feature selection only works with 'Context,Alteration,GeneA'");
  }
  
  for(cs_name in names(CS_list)) {
    features <- strsplit(cs_name,split="::")[[1]]
    coCSs <- CS_list[[cs_name]]
    coFeatures <- sapply(names(coCSs),function(z) strsplit(z,split="::"))
    
    if("Context" %in% feature_sel) {
      coCSs_valid <- unlist(lapply(coFeatures,function(z) z[1])) == features[1]
      coCSs <- coCSs[coCSs_valid]
      coFeatures <- coFeatures[coCSs_valid]
    }
    if("Alteration" %in% feature_sel) {
      coCSs_valid <- unlist(lapply(coFeatures,function(z) z[2])) == features[2]
      coCSs <- coCSs[coCSs_valid]
      coFeatures <- coFeatures[coCSs_valid]
    }
    if("GeneA" %in% feature_sel) {
      coCSs_valid <- unlist(lapply(coFeatures,function(z) z[3])) == features[3]
      coCSs <- coCSs[coCSs_valid]
      coFeatures <- coFeatures[coCSs_valid]
    }
    CS_list[[cs_name]] <- coCSs
  }
  
  return(CS_list)
}

# 'pairwise_minDist' get the minimum distance across genes to create clusters
pairwise_minDist <- function(PW_dist,knowngenes=NULL){
  
  if(is.matrix(PW_dist)) {
    diag(PW_dist) <- Inf
    PW_dist[is.na(PW_dist)] <- Inf
    PW_minDist <- apply(PW_dist,2,function(g2g_dist) min(g2g_dist[!is.infinite(g2g_dist)],na.rm = TRUE))
  } else {
    PW_minDist <- min(PW_dist)
  }
  
  if(!is.null(knowngenes)) {
    PW_minDist[!names(PW_minDist)%in%knowngenes] <- 0
  }
  
  if(length(PW_minDist)==1) {
    if(is.infinite(PW_minDist)) PW_minDist[1] <- 0
  }
  
  
  return(PW_minDist)
}


### 5 Load Cell Sets and calculate jackard indexes #####
cat(paste0("[INFO] :","Loading Gene Dependencies and Cell Sets for adjustment","\n"),file=stdout())
CS <- readRDS(CS_fl)
GD <- readRDS(GD_fl)

# Calculate Jackard index by pairs
cat(paste0("[INFO] : Calculating jackard indexes between Cell sets\n"),file=stdout())
CS_jIdx <- matrix(NA,ncol=length(CS),nrow=length(CS),dimnames=list(names(CS),names(CS)))
diag(CS_jIdx) <- 1

ridx_cnt <- 0
for(ridx in rownames(CS_jIdx)) {
  ridx_cnt <- ridx_cnt+1
  if(ridx_cnt%%100 == 0 ) cat(paste0("#",ridx_cnt," out of ",nrow(CS_jIdx),"\n"),file=stdout())
  for(cidx in colnames(CS_jIdx)) {
    if(ridx==cidx) break;
    
    if(is.na(CS_jIdx[ridx,cidx])) {
      jidx <- jackard_index(CS[[ridx]],CS[[cidx]])
      CS_jIdx[ridx,cidx] <- jidx
      CS_jIdx[cidx,ridx] <- jidx
    }
  }
}

# Calculate overlap between pairs (it takes a long time)
cat(paste0("[INFO] : Calculating Overlap between Cell sets\n"),file=stdout())
CS_oIdx <- matrix(NA,ncol=length(CS),nrow=length(CS),dimnames=list(names(CS),names(CS)))
diag(CS_oIdx) <- 1
ridx_cnt <- 0
for(ridx in rownames(CS_oIdx)) {
  ridx_cnt <- ridx_cnt+1
  if(ridx_cnt%%100 == 0 ) cat(paste0("#",ridx_cnt," out of ",nrow(CS_oIdx),"\n"),file=stdout())
  for(cidx in colnames(CS_oIdx)) {
    if(ridx==cidx) next;
      CS_oIdx[ridx,cidx] <- Overlap_index(CS[[ridx]],CS[[cidx]])
  }
}

# Temporary code for testing
# CS_jIdx <- readRDS("./TMP/contexts.RNAi.CellSets.jackardIdx.rds")
# CS_oIdx <- readRDS("./TMP/contexts.RNAi.CellSets.overlap.rds")

### 6 Remove passenger events which co-occur with oncogenic dependencies #####
cat(paste0("[INFO] :","Removing passenger alterations together with oncogenic dependencies (actual driver events)","\n"),file=stdout())
# Get Oncogenic Dependency by GoF (significant ones)
OncoGD <- unlist(GD[which(GD$GeneA == GD$GeneB & GD$GeneA.Role=="oncogene" & GD$pval < 0.05),
                                      "pathway"])
names(OncoGD) <- NULL

# Get which Cell Sets overlap with these Oncogenic GoFs
coONCs <- stratifyByIdx(t(CS_oIdx), 0.7, 1,
                        selectRows=OncoGD,
                        avoid_itself=TRUE,avoid_triangle = NULL)

# Discard them
discardCS <- vector()
for(coONC in names(coONCs)) {
  features <- strsplit(coONC,split="::")[[1]]
  idx <- which(GD$pathway %in% names(coONCs[[coONC]]) & GD$Context == features[1] & GD$GeneB == features[3])
  discardCS <- unique(c(discardCS,idx))
}

GD <- GD[-discardCS,]


# Check if dependencies over BRAF GoF -> BRAF have dissapeared
# SKIN_BRAF <- GD[GD$GeneB=="BRAF" & GD$Context=="SKIN",]
# SKIN_BRAF[SKIN_BRAF$padj < 0.05,]

### 7 Discard passenger due to GenomicRegions GCN  #####
cat(paste0("[INFO] :","Removing collateral events in genomic clusters due to deletions/amplifications of cancer genes.","\n"),file=stdout())
coCS <- stratifyByIdx(CS_jIdx, jIdx_cuttoff, 1)

coCS <- retrieveCS(coCS,feature_sel = c("Context","Alteration"))
coCS <- coCS[unlist(lapply(coCS,length))>1]

# Calculate Pair-wise gene distances along the genome
coCS_dist <- lapply(coCS,function(CSs,GR) {
  cogenes <- sapply(names(CSs),function(j) strsplit(j,split="::")[[1]][3])
  pairwise_dist(unique(cogenes),GR)
},GR=GR_genes)

coCS_minDist <- lapply(coCS_dist,function(PW_dist, knowngenes) pairwise_minDist(PW_dist,knowngenes), knowngenes=names(GR_genes))

coCS_GR <- lapply(coCS_minDist,function(Dists) names(Dists[Dists<2e6]))
coCS_GR <- coCS_GR[unlist(lapply(coCS_GR,length))>0]

  # Recover Context for the GR genes?
  coCS_GR <- sapply(names(coCS_GR),function(cs_name) {
    # Remove GeneA
    Context_GR <- paste(strsplit(cs_name,split="::")[[1]][-3],collapse="::")
    
    # Get members
    coCS_GRmembers <- coCS_GR[[cs_name]]
    # Complete with context
    coCS_GRmembers <- paste0(Context_GR,"::",coCS_GRmembers)
    # Return
    coCS_GRmembers
  })
  
# First, we are going to treat the problem selecting Cancer Genes
coCS_GRbyCancerGene <- lapply(coCS_GR,function(CS_GRmember,cancer_genes) {
  
  GRmembers <- unlist(sapply(CS_GRmember,function(z) strsplit(z,split="::")[[1]][3]))
  CS_GRmember <- CS_GRmember[GRmembers %in% cancer_genes]
  CS_GRmember
  },cancer_genes=cancer_genes)
coCS_GRbyCancerGene <- unique(unlist(coCS_GRbyCancerGene))

discarded <- setdiff(unlist(coCS_GR[names(coCS_GR) %in% coCS_GRbyCancerGene]),coCS_GRbyCancerGene)

discarded_idx <- which(GD$pathway %in% discarded)

GD <- GD[-discarded_idx,]

### 8 Save the adjusted object ######
cat(paste0("[INFO] :","Saving the adjusted GDs (data.table of R, RDS format).","\n"),file=stdout())
saveRDS(GD,paste0(out_dir,"/","GD_adj.",dataset_name,".rds"))