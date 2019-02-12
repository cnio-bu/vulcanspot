



cols.interest <- c("symbol","name","role","description","entrez_id","driver","druggable","score")

all_genes <- scan("./data/CCL/genes.txt",what = "character")

mygene.info <- read.table("./data/CCL/mygene.tsv",sep="\t",header=TRUE,quote="",fill=TRUE)


Gscore <- read.table("./data/CCL/gscores.tsv",sep="\t",header=FALSE,quote="",col.names = c("symbol","Gscore"))
Gscore <- setNames(Gscore$Gscore,Gscore$symbol)

# Load PanDrugs
PD <- read.table("./src/PanDrugs/PanDrugsFiles/Pandrugs_Feb2018.tsv",sep="\t",header=TRUE,
                 check.names = FALSE,quote="",fill=TRUE)
# Subset for only target
PD <- PD[PD$target_marker=="target" & 
           PD$resistance=="sensitivity",]
Gdruggable <- unique(PD$gene_symbol)

# Load the Role
CGC <- read.table("./data/Annotation/Census_allWed_Oct_4_09-53-19_2017.tsv",sep="\t",header = TRUE,fill=TRUE)
CGC.Role <- setNames(sapply(gsub(" +","",CGC$Role.in.Cancer),function(z) strsplit(z,split=",")[[1]]),CGC$Gene.Symbol)
CGC.ONC <- names(CGC.Role)[unlist(lapply(CGC.Role,function(z) any(z=="oncogene") & !any(z=="TSG")))]
CGC.TSG <- names(CGC.Role)[unlist(lapply(CGC.Role,function(z) !any(z=="oncogene") & any(z=="TSG")))]
CGC.AMB <- names(CGC.Role)[unlist(lapply(CGC.Role,function(z) any(z=="oncogene") & any(z=="TSG")))]

genes.role <- sapply(all_genes,function(gene) {
  if(gene %in% CGC.ONC) {
    "oncogene"
  } else if(gene %in% CGC.TSG) {
    "TSG"
  } else if(gene %in% CGC.AMB) {
   "oncogene/TSG"
  } else {
    "unknown"
  }
})

# Load Drivers
Gdrivers <- read.table("./data/Annotation/srep02650-s3.csv",sep=",",header=TRUE)
Gdrivers <- setNames(Gdrivers$Putative.Driver.Category, Gdrivers$Gene.Symbol)

Drivers <- sapply(all_genes,function(gene) {
  if(gene %in% names(Gdrivers)) {
    Gdrivers[gene]
  } else {
    "ND"
  }
})


stopifnot(all(all_genes==mygene.info$gene_symbol))
# stopifnot(all(all_genes==Gscore$symbol))
stopifnot(all(all_genes==names(genes.role)))

gene.tab <- data.frame("id"=1:length(all_genes),
                      "symbol"=all_genes,
                       "name"=mygene.info$gene_name,
                       "role"=genes.role,
                       "description"=mygene.info$gene_description,
                       "entrez_id"=mygene.info$entrez_id,
                       "driver"=Drivers,
                       "druggable"=all_genes %in% Gdruggable,
                       "score"=Gscore[all_genes])

write.table(gene.tab,file="./data/CCL/genes.tsv",sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)
