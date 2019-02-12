### TITLE : Generate the table of genes 
### AUTHOR: Perales-Paton, Javier - jperales@cnio.es
### LICENSE : GPL-v3

options(stringsAsFactors = FALSE)

### Load libraries #######
suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(mygene))

### 1 Get parameters ############
option_list = list(
  make_option(c("--INPUT"), action="store", default="./data/CCL/genes.txt", type='character',
              help="TXT file with the genes to be annotated by mygene."))

opt = parse_args(OptionParser(option_list=option_list))

# Print Parameters
cat("PARAMETERS:\n",file=stdout())
for(opt.idx in 1:length(opt)) {
  cat(paste0(names(opt)[opt.idx],"\t",opt[[opt.idx]],"\n"),file = stdout())
}
cat("---\n",file=stdout())

fl <- opt$INPUT

GeneX <- scan(fl,what = "character")

### 2  #####################
table_columns <- c("gene_symbol","gene_name","gene_description","entrez_id")
table_GeneX <- data.frame(matrix(NA,ncol=length(table_columns),nrow=length(GeneX),
                                dimnames=list(GeneX,table_columns)))
table_GeneX$gene_symbol <- GeneX

# Reannotate
i <- 0
for(Symbol in GeneX) {
  i <- i+1;
  if(i%%100==0) cat(paste0("#",i,"/",length(GeneX),"\n"),file=stdout())
  tmp <- query(q=Symbol,scopes="symbol",species="human",fields="entrezgene")$hits
  if(length(tmp)!=0) {
    # Get the "EntrezGeneID" by REST service
    # gene.id <- tmp[1,"entrezgene"]
    if("entrezgene"%in%colnames(tmp)) {
      gene.id <- na.omit(tmp[,"entrezgene"])[1]
    } else {
      gene.id <- tmp[1,"_id"]
    }
    rm(tmp)
    # Get some annotations
    status <- tryCatch(
      data.mygene <- mygene::getGene(gene.id, fields = c("summary","name")),
      error = function(e) e
    )
    
    if(inherits(status,  "error")) {
      cat(paste0("[ERROR]: GeneSymbol=",Symbol," failed in getting metadata.\n"),file=stdout())
      data.mygene <- list(c("_id"=gene.id,"name"=gene.id,"summary"="N/A"))
    }
    # data.mygene <- mygene::getGene(gene.id, fields = c("summary","name"))
    
    table_GeneX[Symbol,"entrez_id"] <- data.mygene[[1]][["_id"]]
    if(!is.null(data.mygene[[1]][["name"]])) {
      table_GeneX[Symbol,"gene_name"] <- data.mygene[[1]][["name"]]
    }
    
    if(!is.null(data.mygene[[1]][["summary"]])) {
      table_GeneX[Symbol,"gene_description"] <- data.mygene[[1]][["summary"]]
    } 
  }
}


write.table(table_GeneX, file="./data/CCL/mygene.tsv",sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
