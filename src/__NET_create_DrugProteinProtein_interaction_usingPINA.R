### TITLE : Create the Protein-Protein-Drug interaction network based on PINA
### AUTHOR: Perales-Paton, Javier - jperales@cnio.es
### LICENSE : GPL-v3

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(require(igraph))

# We initialize the directory where the graph is going to be saved
if(!dir.exists("./data/NET/")) dir.create("./data/NET/")

### 1 Create a dictionary of UniprotKB to HUGNC gene symbols #####
### PINA network is described using PPI between UniProtKB. We are going to transform these entities
### into human gene symbols.
cat("[INFO] Create a dictionary to map UniprotID to HGNC gene symbol in homo sapiens (source: ensEMBL 2014)\n",file=stdout())
# Load Annotation
human.BM <- read.table("./data/Annotation/ANN_uniprot_feb2014_GRCh37.tsv",sep="\t",header=TRUE,quote="")
# We are only interested in Uniprot and HGNC
human.BM <- unique(human.BM[,c("hgnc_symbol","uniprot_swissprot_accession")])
human.BM <- human.BM[which(human.BM$uniprot_swissprot_accession!=""),]

# Create a dictionary to match uniprot accs to gene symbols in human
uniprot2symbol <- setNames(vector("list",length = length(human.BM$uniprot_swissprot_accession)),
                           human.BM$uniprot_swissprot_accession)
for(uniprot.id in names(uniprot2symbol)) {
  uniprot2symbol[[uniprot.id]] <- human.BM[which(human.BM$uniprot_swissprot_accession==uniprot.id),"hgnc_symbol"]
}

# The dictionary looks good because the vast majority of the entries from UniProt has a unique gene symbol,
# and almost all uniprot ids have at least one gene symbol:
  # table(unlist(lapply(uniprot2symbol,length)))
  # 0     1     2     3     4     5     7    10    12    14 
  # 171 18984    90    12     2     3     1     1     1     1

### 2 build the PPI network based on PINA database #######
cat("[INFO] Build Protein-Protein Interaction network based on supported evidences (source:PINA)\n",file=stdout())

# Check if it needs to download 
if(!file.exists("./data/NET//PINA_Homosapiens-20140521.tsv")) {
  cat("[INFO] Downloading PINA 2014 (interactome) in MITAB format.\n",file=stdout())
  download.file("./data/NET/PINA_Homosapiens-20140521.tsv",
                url = "http://cbg.garvan.unsw.edu.au/pina/download/Homo%20sapiens-20140521.tsv")
}

# load PINA network
cat("[INFO] Reading PPI network (source: PINA)\n",file=stdout())
PINA <- read.table("data/NET/PINA_Homosapiens-20140521.tsv",sep="\t",header=TRUE,quote = "",fill=TRUE)
NoDatabases <- sapply(PINA$X.Source.database.s..,function(z) length(unique(strsplit(z,split="\\|")[[1]])))
      # > table(NoDatabases)
      # NoDatabases
      # 1      2      3      4      5 
      # 111230  39198   9333   6604    411
MI.interest <- c("MI:0218", # physical interaction
                 "MI:0407") # direct interaction
MI.bol <- grepl(paste(MI.interest,collapse = "|"),PINA$X.Interaction.type.s..)

# Subset PINA network based on Criteria:
#   - at least two databases supporting the assoc.
#   - Evidence of PPI : Physical and direct interaction of proteins

cat("[INFO] Subsetting the whole PPIN network (source: PINA) following this criteria:\n",file=stdout())
cat("[INFO] \t - At least PPI supported by two different original databases (IntAct, BioGRID, MINT, DIP, HPRD, MIPS/MPact)\n",file=stdout())
cat("[INFO] \t - Evidence of PPI supported: MI:0218 (physical interaction) OR MI:0407 (direct interaction)\n",file=stdout())
NET <- PINA[NoDatabases >= 2 & MI.bol,]
rownames(NET) <- NULL

# Rename identifiers to map in the dictionary
cat("[INFO] Mapping proteins (Uniprot) to genes (gene symbol)...\n",file=stdout())
NET$X.ID.s..interactor.A. <- gsub("^uniprotkb:","",NET$X.ID.s..interactor.A.)
NET$X.ID.s..interactor.B. <- gsub("^uniprotkb:","",NET$X.ID.s..interactor.B.)

NET$Symbol.A <- sapply(NET$X.ID.s..interactor.A.,function(z) if(length(uniprot2symbol[[z]])==0) NA else uniprot2symbol[[z]][1])
NET$Symbol.B <- sapply(NET$X.ID.s..interactor.B.,function(z) if(length(uniprot2symbol[[z]])==0) NA else uniprot2symbol[[z]][1])

# We remove those Uniprot Ids which couldn't map in hgnc
cat("[INFO] Remove Vertices which couldn't be mapped...\n",file=stdout())
NET <- NET[!is.na(NET$Symbol.A) & !is.na(NET$Symbol.B),]

# Finally, build the igraph object with the network
cat("[INFO] Build igraph object and save it.\n",file=stdout())
g <- graph_from_data_frame(NET[,c("Symbol.A","Symbol.B")],
                           directed = FALSE)
# PPI are mutual
g <- as.directed(g,mode = "mutual")

# For the records, we are going to set some graph attributes (vertex and edges)
vertex_attr(g) <- c(vertex_attr(g),list("entity"=rep("gene",length(V(g)))))
edge_attr(g) <- list("link"=rep("ppi",length(E(g))))
  # Take a look to the ppi 
  # plot(induced_subgraph(g,vids = c("CALR",names(neighbors(g,v = "CALR")))))

saveRDS(g,"./data/NET//UNIVERSAL_PPIN.PINA.igraph.rds")

### 3 CREATE DRUG -> gene target network #####
cat("[INFO] Build Drug-Gene target Interaction network (source:CLUE)\n",file=stdout())

source("./src/CLUE/api_clue.R")

# Read the user_key password. You can get your own password in www.clue.io for the use
# of the API RestFull service. We are going to use this API to get the information of
# Drug -> target interactions
user_key <- scan(file="./src/CLUE/user_key.txt",what = "character")

cat("[INFO] Performing several Queries to the API RestFullservice at CLUE to get the drug targets relationships.\n",file=stdout())
sink(file="/dev/null") # The user_key (private) is shown as part of the log, so we try to get it hiden
tmp <- GETpert_info(pert_iname = NULL,pert_ids = NULL,where = list("pert_type"="trt_cp"),
                    fields = c("pert_iname","target"),
                    user_key = user_key,
                    service="perts")
sink() # disable sink

# Renaming element names of the list: avoid 1:1000 by each query as element names of the list
names(tmp) <- NULL

# Reshape the drug-target information as a data.frame
DRT.df <- data.frame(A=unlist(sapply(1:length(tmp),function(idx) rep(tmp[[idx]]$pert_iname,length(tmp[[idx]]$target)))),
                  B=unlist(lapply(tmp,function(z) z$target)))
DRT.df <- unique(DRT.df)

# Build the igraph object with drug-target info
cat("[INFO] Build igraph and save it\n",file=stdout())
dt <- graph_from_data_frame(DRT.df,directed = TRUE)

# Add graph attributes
vertex_attr(dt) <- c(vertex_attr(dt),list("entity"=ifelse(names(V(dt)) %in% DRT.df$A,"drug","gene")))
edge_attr(dt) <- list("link"=rep("dt",length(E(dt))))

write.table(DRT.df,"./data/NET/CLUE_drugtargets.tsv",sep="\t",row.names=FALSE,col.names = TRUE,quote=FALSE)
saveRDS(dt,"./data/NET//UNIVERSAL_DrugTIN.CLUE.igraph.rds")

### 4 Append a unique graph with the two elements: Protein-Protein Interactions and Drug-Target intractions ######
cat("[INFO] Build Appended graph of Protein-Protein & Drug target Interaction network (combined PINA+CLUE)\n",file=stdout())
net <- igraph::union(g,dt)

# Vertex attributes
vattr_sub <- grepl("entity_[0-9]$",names(vertex_attr(net)))
vertex_attr(net) <- c(vertex_attr(net)[!vattr_sub],
                      setNames(vertex_attr(net)[vattr_sub],
                               paste0(gsub("[0-9]$","",names(vertex_attr(net)[vattr_sub])),c("ppi","dt"))))
entity.mat.tmp <- matrix(c(vertex_attr(net)$entity_ppi,
                           vertex_attr(net)$entity_dt),ncol=2,byrow = FALSE)
# unique(entity.mat.tmp) # All posibilities. I have checked it out that is correct
entity.tmp <- apply(entity.mat.tmp,1,function(z) {
  if (z[2]=="drug" & !is.na(z[2])) {
    "drug"
  } else {
    "gene"
  }
})

vertex_attr(net) <- c(vertex_attr(net),list("entity"=entity.tmp))
V(net)$color <- c("green","gold", "tomato")[as.integer(as.factor(vertex_attr(net)$entity))]
# vertex_attr(net) <- c(vertex_attr(net),"entity"=unlist(lapply(vertex_attr(net)[vattr_sub],function(z) z[!is.na(z)])))
V(net)$size <- as.vector(sqrt(igraph::degree(net))+1)

# edges attributes
net <- set.edge.attribute(net,"color",value=ifelse(is.na(edge_attr(net)$link_2),"grey60","red4"))

# Save it
cat("[INFO] Save it\n",file=stdout())
saveRDS(net,file = "./data/NET//UNIVERSAL_DrugPPIN.PINACLUE.igraph.rds")

### 5 Just try some examples to show how the network looks #######
cat("[INFO] Saving an example of BIRC5 and its neighbors from the built PPDI network.\n",file=stdout())
png("./log/DrugPPIN_BIRC5_neighbors.png",height = 800*3,width = 800*3,res=280)
plot(induced_subgraph(net,vids = c("BIRC5",names(neighbors(as.undirected(net),v = "BIRC5")))))
dev.off()

cat("[INFO] Saving another example of AURKA and its neighbors from the built PPDI network.\n",file=stdout())
png("./log/DrugPPIN_AURKA_neighbors.png",height = 800*3,width = 800*3,res=280)
plot(induced_subgraph(net,vids = c("AURKA",names(neighbors(as.undirected(net),v = "AURKA")))))
dev.off()

cat("[INFO] Saving an example of BIRC5 and AURKA and its neighborhod from the built PPDI network.\n",file=stdout())
png("./log/DrugPPIN_BIRC5_AURKA_neighbors.png",height = 800*3,width = 800*3,res=280)
plot(induced_subgraph(net,
                      vids = c("BIRC5",names(neighbors(as.undirected(net),v = "BIRC5")),
                               "AURKA",names(neighbors(as.undirected(net),v = "AURKA")))))
dev.off()

