### TITLE : Generate the table of drug-gene associations 
### AUTHOR: Perales-Paton, Javier - jperales@cnio.es
### LICENSE : GPL-v3


options(stringsAsFactors = FALSE)

### Load libraries #######
suppressPackageStartupMessages(require(optparse))

# 1 Get parameters ############
option_list = list(
  make_option(c("--INDIR"), action="store", default="./data/DrugPrescription/", type='character',
              help="Input directory which contains drug prescription files"),
  make_option(c("--OUTDIR"), action="store", default="./data/DrugPrescription/", type='character',
              help="Output directory"),
  make_option(c("--DSCORE"), action="store", default=NULL, type='numeric',
              help="Minimum threshold for Drug score (PanDrugs)"),
  make_option(c("--KDCP_SCORE"), action="store", default=NULL, type='numeric',
              help="Minimum threshold for KDCP score (LINCS)")
)


opt = parse_args(OptionParser(option_list=option_list))

# Print Parameters
cat("PARAMETERS:\n",file=stdout())
for(opt.idx in 1:length(opt)) {
  cat(paste0(names(opt)[opt.idx],"\t",opt[[opt.idx]],"\n"),file = stdout())
}
cat("---\n",file=stdout())

input.dir <- opt$INDIR
output.dir <- opt$OUTDIR
DScore.cutoff <- opt$DSCORE
KDCP_score.cutoff <- opt$KDCP_SCORE

### 2 Load cross-reference database entities ######
  # Drugs
  # drugs <- read.table("./data/DrugPrescription/drugs.tsv",sep="\t",quote="",header=TRUE)
  drugs <- readRDS("./data/DrugPrescription/drugs.rds")
  drugsID <- setNames(drugs$id,drugs$internal_id)
  # Genes
  genes <- read.table("./data/CCL/genes.tsv",sep="\t",header=TRUE,quote="")
  genesID <- setNames(genes$id,genes$symbol)
  # Contexts
  Contexts <- read.table("./data/CCL/contexts.tsv",sep="\t",header=TRUE, quote="")
  ContextsID <- setNames(Contexts$id, Contexts$name)
  # Sources
  sources <- read.table("./data/DrugPrescription/sources.tsv",sep="\t",header=TRUE,quote="")
  sourcesID <- setNames(sources$id,sources$name)
  
# dga.colnames <- c("Gene","Drug","Score","Pval","Context","Source","id_drug")
dga.colnames <- c("id_drugs"="id_drugs","id_genes"="id_genes","id_contexts"="id_contexts","id_sources"="id_sources",
                  "score"="KDCP_score")


### 2 Collect all Drug Repositioning files (individual for each lineage) #######
# The context (lineage or cell line) is automatically detected by the pattern of the filename
repurposing.fls <- list.files(input.dir,pattern="repositioning",full.names = TRUE)
names(repurposing.fls) <- gsub("repositioning\\.(.*)\\.features\\.tsv","\\1",basename(repurposing.fls))


KDCP <- vector()
for(context in names(repurposing.fls)) {
  fl <- repurposing.fls[context]
  tmp <- read.table(fl,sep="\t",header=TRUE,quote="")
  if(nrow(tmp) > 0) {
    tmp$Pval <- NA
    tmp$Context <- context
    tmp$Source <- "LINCS"
    if(!is.null(KDCP_score.cutoff)) {
      tmp <- tmp[tmp$KDCP_score >= KDCP_score.cutoff,
                 c("Gene","Drug","KDCP_score","Context","Source","id_drug")]
    } else {
      tmp <- tmp[,c("Gene","Drug","KDCP_score","Context","Source","id_drug")]
    }
    
    KDCP <- rbind(KDCP,tmp) 
  }
}

KDCP$id_drugs <- drugsID[KDCP$id_drug]
KDCP$id_genes <- genesID[KDCP$Gene]
KDCP$id_contexts <- ContextsID[KDCP$Context]
KDCP$id_sources <- sourcesID[KDCP$Source]


KDCP <- KDCP[order(KDCP$Gene,KDCP$Context,-KDCP$KDCP_score), dga.colnames]
colnames(KDCP) <- names(dga.colnames)
rownames(KDCP) <- NULL

### 3 Add the PanDrugs prescription ##########
# Note: In this invidiual case, file contains all the precription amongst context-specific GDs

# PanDrugs developer told me the header of the file
PD_header <- c("tumor_type","donor","affected_gene","gscore","mut_branch","cnv_branch",
               "exp_branch","cons_branch","alert","gene_symbol","sources","show_drug_name",
               "status","pathology","cancer","extra","extra2","pathways","target_marker",
               "resistance","alteration","ind_pathway","relation","ind-gene","dscore")

# Read the PanDrugs prescription
PD <- read.table(paste0(input.dir,"/knowledgebased_pandrugs.tsv"),sep="\t",
                 header=FALSE,quote="",fill=TRUE,col.names = PD_header)
# Remove those Gene Bs that are not present in PanDrugs database
PD <- PD[PD$gene_symbol!="",]

# Remove those grug-gene assoc. without a drug (?)
PD <- PD[PD$show_drug_name!="",]

# Retrieve only drugs with a direct target, no experimental
PD <- PD[PD$target_marker=="target" & 
         PD$relation=="DIRECT" & 
         PD$resistance=="sensitivity",]

# [Optional] Filter-out drug-gene associations based on DScore
if(!is.null(DScore.cutoff)) {
  PD <- PD[PD$dscore >= DScore.cutoff,]
}

# Load the metadata that match both drug prescription approaches: LINCS and PanDrugs
pd_id <- read.table("./data/Annotation/PD_table.tsv",sep="\t",header=TRUE,quote="")

# Add this metadata to the original table
PD <- merge(x=PD, y=pd_id, by.x="show_drug_name", by.y="source_id",all.x=TRUE, all.y=FALSE)

# Next, we are going to reformat the PanDrugs table
# Rename some columns
colnames(PD)[which(colnames(PD)=="pd_id")] <- "id_drug"
colnames(PD)[which(colnames(PD)=="tumor_type")] <- "Context"

# Add some dummy features per instance
PD$Source <- "PANDRUGS"
PD$id_drugs <- drugsID[PD$id_drug]
PD$id_genes <- genesID[PD$gene_symbol]
PD$id_contexts <- ContextsID[PD$Context]
PD$id_sources <- sourcesID[PD$Source]

# Subset the table to the information that is relevant
PD <- PD[,c("id_drugs","id_genes","id_contexts","id_sources","dscore")]
# Unique table
PD <- unique(PD)

#Rename
colnames(PD) <- names(dga.colnames)
rownames(PD) <- NULL

### 4 Write down the drug-gene associations
# KDCP approach
write.table(KDCP,paste0(output.dir,"/lnk_genes_drugs.tsv"),
            sep="\t",row.names = FALSE,col.names = TRUE, quote=FALSE)

# PanDrugs
first_bool <- !file.exists(paste0(output.dir,"/lnk_genes_drugs.tsv"))
write.table(PD,paste0(output.dir,"/lnk_genes_drugs.tsv"),
            sep="\t",row.names = FALSE,col.names = first_bool,quote=FALSE,
            append = !first_bool)