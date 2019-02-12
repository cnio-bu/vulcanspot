#!/usr/bin/env Rscript

### TITLE : Calculate KDCP score
### AUTHOR : Perales-Paton, Javier - jperales@cnio.es
### LICENSE : GPL-v3


options(stringsAsFactors = FALSE)

### Load libraries #######
suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(igraph))
suppressPackageStartupMessages(require(reshape2))


# 1 Get parameters ############
option_list = list(
  make_option(c("--CONTEXT"), action="store", default="VCAP", type='character',
              help="PANCANCER OR (Tissue or cell line from the Cancer Cell Line Encyclopedia)"),
  make_option(c("--GDFILE"), action="store", default=NULL, type='character',
              help="TSV file with the GeneB that are significant in 'gene dependencies' relationships."),
  make_option(c("--NRND"), action="store", default=500, type='numeric',
              help="Number of Randomizations to calculate Probability of finding a shortest path in randomized models."),
  make_option(c("--OUTDIR"), action="store", default="./data/DrugPrescription/", type='character',
              help="Output directory.")
  )


opt = parse_args(OptionParser(option_list=option_list))

# Print Parameters
cat("PARAMETERS:\n",file=stdout())
for(opt.idx in 1:length(opt)) {
  cat(paste0(names(opt)[opt.idx],"\t",opt[[opt.idx]],"\n"),file = stdout())
}
cat("---\n",file=stdout())

context <- opt$CONTEXT
GD_fl <- opt$GDFILE
outdir <- opt$OUTDIR
N_rnds <- opt$NRND

# Define a function to get molten matrices
melt_genesXdrugs <- function(mat,cols=c("Gene","Drug","Value")) {
  stopifnot(is.matrix(mat))
  stopifnot(all(is.numeric(unlist(mat))))
  stopifnot(length(cols)==3)
  
  df <- setNames(reshape2::melt(mat,variable.factor = FALSE),cols)
  df[,cols[1]] <- as.character(df[,cols[1]])
  df[,cols[2]] <- as.character(df[,cols[2]])
  
  return(df)
}

### 1 Load the TAU of KD over CP ####
cat(paste0("[INFO] Loading tau matrix from '",context,"'\n"),file=stdout())
tau.fl <- paste0("./data/CMap-L1000/",context,".TAU_matrix.rds")
if(file.exists(tau.fl)) {
  tau <- readRDS(tau.fl)
  tau <- tau[,!is.na(tau[1,])] # Remove KDs with no data
} else {
  stop("ERROR #2 : The cell line has no KD-CP perturbation data in CMap-L1000 dataset.")
}

### 2 Load the table of CP signatures #######
# The drug info from the CMap-L1000 dataset
cat(paste0("[INFO] Reading the metadata of gene expression signatures.\n"),file=stdout())
tab <- read.table("./data/Annotation/sig_id_table_LINCS.tsv",
                  sep="\t",header=TRUE,quote="",fill=TRUE,colClasses = "character")
# short_list <- read.table("./Tidy/sig_ids_short.txt",sep="\t",header=FALSE,colClasses = "character")[,1]
# tab <- tab[tab$sig_id %in% short_list,]
tab$sig_id <- paste0("sig_",tab$sig_id)

sig_id2common_name <- setNames(tab$common_name,tab$sig_id)

# We collect those drug signatures that are redundant on the perturbation.
# The source of these perturbations is the different vendors that the chemical compounds
# were obtained from. Anyway, we have created a consensus signature for these cases that are
# those we are going to retrieve for further analysis.
common_names.duplicated <- tab$common_name[which(tab$consensus=="1")]
pert.replicate <- tab$sig_id[which(tab$common_name%in%common_names.duplicated & 
                                              tab$consensus=="0")]

### 3 Remove these replicate signatures from the KD-CP matrices ######
tau <- tau[,which(!colnames(tau) %in% pert.replicate)]

### 4 Load the Contextualized network of drug-targets and protein-protein interactions ####
cat(paste0("[INFO] Loading the Drug-PPIN from '",context,"'\n"),file=stdout())
net <- readRDS(paste0("./data/NET/",context,"_DrugPPIN.PINACLUE.igraph.rds"))
# From the network
drugs.net <- names(V(net))[which(vertex_attr(net)$entity=="drug")]
genes.net <- names(V(net))[which(vertex_attr(net)$entity=="gene")]

# From the CMap-L1000 dataset, which includes some drugs without known targets
drugs.pert <- sig_id2common_name[colnames(tau)]
genes.pert <- rownames(tau)

# The universe of Drugs and Genes
all.drugs <- sort(unique(c(drugs.net,drugs.pert)))
all.genes <- sort(unique(c(genes.net,genes.pert)))
all.genes <- all.genes[all.genes!=""]

### 5 Calculate distances between drug-targets all around the network #####
calculate_distances <- function(g,all.genes, all.drugs) {
  drugs.net <- names(V(g))[which(vertex_attr(g)$entity=="drug")]
  genes.net <- names(V(g))[which(vertex_attr(g)$entity=="gene")]
  
  # Create the matrix which will contain all the distances for
  #   - Genes and Drugs in the graph
  #   - Genes and Drugs out of the graph (by avg)
  Genes2Drug_dist <- matrix(NA,
                            nrow=length(all.genes),
                            ncol=length(all.drugs),
                            dimnames=list(all.genes,
                                          all.drugs))
  
  # Calculate the actual distances through the network
  g_dist <- distances(graph = g,v = genes.net,to = drugs.net)
  # Those unconnected drugs (i.e. wo known target) has Inf values.
  # We will set them to NA for calculation purposes
  g_dist[is.infinite(g_dist)] <- NA
  # It could happen that there is one row that is "", so we remove it
  g_dist <- g_dist[rownames(g_dist)!="",]
  
  # Fill the matrix of distances with those calculated in the network
  Genes2Drug_dist[rownames(g_dist),colnames(g_dist)] <- g_dist
  
  # For now, this is so far we can do with the ppi+dt network...
  # Besides this, we have unconnected nodes with not available properties that we should infer.
  # So we decided to replace those unconnected nodes with averages tendencies following:
  #   1.- Drugs without known target : the average of each gene through the network
  #   2.- Genes with no drugs targeting them : the average of all genes through the network
  
  # Identify them
  missing_genes <- rownames(Genes2Drug_dist)[rowSums(is.na(Genes2Drug_dist)) == ncol(Genes2Drug_dist)]
  missing_drugs <- colnames(Genes2Drug_dist)[colSums(is.na(Genes2Drug_dist)) == nrow(Genes2Drug_dist)]
  
  # 1. Do similarly with the drugs
  drugs_AVGdist <- setNames(rep(NA,length(all.drugs)),all.drugs)
  drugs_AVGdist[colnames(g_dist)] <- colMeans(g_dist,na.rm=TRUE)
  
  # Many of these drugs are unconnected in the network
  drugs_AVGdist[is.na(drugs_AVGdist)] <- mean(drugs_AVGdist,na.rm=TRUE)
  
  for(missing_drug in missing_drugs) {
    Genes2Drug_dist[,missing_drug] <- drugs_AVGdist[missing_drug]
  }
  
  # 2. We start with the genes
  genes_AVGdist <- setNames(rep(NA,length(all.genes)),all.genes)
  genes_AVGdist[rownames(g_dist)] <- rowMeans(g_dist,na.rm = TRUE)
  
  # Many of these Genes are not in the network
  genes_AVGdist[is.na(genes_AVGdist)] <- mean(genes_AVGdist,na.rm = TRUE)
  
  for(missing_gene in missing_genes) {
    Genes2Drug_dist[missing_gene,] <- genes_AVGdist[missing_gene]
  }
  
  # Fill remaining NA values, which belong to drugs that have a few gene targets and
  # those targets are totally unconnected from the rest of the PPIN. We fill the NAs
  # with the average distance from any drug to that gene with the NA.
  NAs.idx <- which(is.na(Genes2Drug_dist),arr.ind = TRUE)
  if(nrow(NAs.idx)>0) {
    for(idx in 1:nrow(NAs.idx)) {
      Genes2Drug_dist[NAs.idx[idx,"row"],NAs.idx[idx,"col"]] <- genes_AVGdist[rownames(NAs.idx)[idx]]
    }
  }

  return(Genes2Drug_dist)
}

# Get distances
cat(paste0("[INFO] Calculating Drug-Target gene distances thorugh Drug-PPIN.\n"),file=stdout())
G2D.dist <- calculate_distances(g=net,all.genes = all.genes, all.drugs=all.drugs)

# This block of code was reprogrammed to avoid overload of RAM memory with
#   temporary random networks allocation.
    # rnd_nets <- vector(mode = "list",length=N_rnds)
    # for(iter in 1:length(rnd_nets)) {
    #   rnd_nets[[iter]] <- rewire(net,with = keeping_degseq(loops = TRUE,
    #                                                        niter = 10*ecount(net)))
    # }
    # 
    # rnd_nets.dist <- lapply(rnd_nets,function(rnd_net) {
    #   calculate_distances(g=rnd_net,all.genes = all.genes, all.drugs=all.drugs)
    # })
    # 
    # rnd_count <- rnd_nets.dist[[1]] < G2D.dist
    # for(rnd_net.idx in 2:length(rnd_nets.dist)) {
    #   rnd_count <- rnd_count + (rnd_nets.dist[[rnd_net.idx]] < G2D.dist )
    # }
    # rnd_count <- rnd_count/N_rnds

cat(paste0("[INFO] Calculating Probability of finding a shortest distance",
           " between two Drug-target pairs in randomized networks preserving",
           " the original degree.",
           " Randomizations=",N_rnds,
           "\n"),file=stdout())

rnd_count <- matrix(0,nrow=nrow(G2D.dist),ncol=ncol(G2D.dist),dimnames=dimnames(G2D.dist))
for(iter in 1:N_rnds) {
  if(iter%%10 == 0) cat(paste0("#",iter," out of ",N_rnds," randomizations.\n"),file=stdout());
  
  rnd_net <- rewire(net,with = keeping_degseq(loops = TRUE,
                                              niter = 10*ecount(net)))
  rnd_net.dist <- calculate_distances(g=rnd_net,all.genes = all.genes, all.drugs=all.drugs)
  rnd_count <- rnd_count + (rnd_net.dist < G2D.dist )
  
  rm(rnd_net, rnd_net.dist)
}
rnd_count <- rnd_count/N_rnds


# We melt the matrix into a data.frame
D2G.features <- melt_genesXdrugs(G2D.dist,cols=c("Gene","Drug","Dist"))
rnd_count2 <- melt_genesXdrugs(rnd_count,cols=c("Gene","Drug","LowerDist_prob"))
stopifnot(all(D2G.features$Gene == rnd_count2$Gene))
stopifnot(all(D2G.features$Drug == rnd_count2$Drug))
D2G.features$LowerDist_prob <- rnd_count2$LowerDist_prob

### 6 Calculate degree for the Drugs #####
cat(paste0("[INFO] Obtain the Connectivity degree of all Drugs.\n"),file=stdout())
drugs.degree <- setNames(rep(NA,length(all.drugs)),all.drugs)
drugs.degree.net <- igraph::degree(net)[drugs.net]

drugs.degree[names(drugs.degree.net)] <- drugs.degree.net
drugs.degree[is.na(drugs.degree)] <- mean(drugs.degree,na.rm=TRUE)
drugs.degree[which(drugs.degree==0)] <- mean(drugs.degree,na.rm=TRUE)


D2G.features$Degree <- drugs.degree[D2G.features$Drug]

D2G.features$key <- paste0(D2G.features$Gene,"::",D2G.features$Drug)

### 7 Melt tau #####
cat(paste0("[INFO] Create table of drug-target gene features.\n"),file=stdout())
tau <- melt_genesXdrugs(tau,cols=c("Gene","id_drug","tau"))
tau$Drug <- sig_id2common_name[tau$id_drug]
tau$key <- paste0(tau$Gene,"::",tau$Drug)

tab <- merge(tau,D2G.features[,c("Dist","Degree","LowerDist_prob","key")],
             by.x="key",by.y="key",all.x = TRUE)

cat(paste0("[INFO] ... adding KDCP score.\n"),file=stdout())

# The second version of KDCP score, described in the manuscript
tab$KDCP_score <- tab$tau / sqrt(1+tab$LowerDist_prob)

### 8 [optional] Retrieve only Drug Prescription on a white list of genes,
### related to Gene Dependencies

if(!is.null(GD_fl)) {
  if(file.exists(GD_fl)) {
    GDs <- read.table(GD_fl, sep="\t",header=TRUE,quote = "")
    GenesGD <- as.character(GDs$GeneB[which(GDs$Context==context)])
    # Perform the filter-out
    tab <- tab[which(tab$Gene %in% GenesGD),]
  } else {
    stop("[ERROR] File with Gene Dependencies (GDs) does not exist.")
  }
}


### 9 Save it ####
cat(paste0("[INFO] Save table\n"),file=stdout())
write.table(tab[,colnames(tab)!="key"],
            paste0("./data/DrugPrescription/repositioning.",context,".features.tsv"),sep="\t",
            col.names=TRUE,row.names=FALSE,quote=FALSE)

cat(paste0("[INFO] Finished\n"),file=stdout())