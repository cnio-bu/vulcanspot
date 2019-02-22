### TITLE: Contextualized Universal DPPIN by cell lineage (context)
### AUTHOR : Perales-Paton, Javier - jperales@cnio.es
### LICENSE : GPL-v3

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(igraph))
suppressPackageStartupMessages(require(ComplexHeatmap))

### 1 Get parameters ########
option_list = list(
  make_option(c("--UNIVERSALNET"), action="store",
              default="./data/NET/UNIVERSAL_DrugPPIN.PINACLUE.igraph.rds",
              type='character',
              help="RDS file with the igrpah object for the network"),
  make_option(c("--CONTEXT"), action="store", default="PANCREAS", type='character',
              help="PANCANCER OR (Tissue or cell line from the Cancer Cell Line Encyclopedia)"),
  make_option(c("--CONTEXT_TYPE"), action="store", default = "tissue",type='character',
              help="Two options: tissue or cell"),
  make_option(c("--UPC"), action="store",
              default = "./data/CCL/GEP_hgnc_symbol_feb2014_GRCh37.UPC.tsv", type='character',
              help="TSV file with UPC values"),
  make_option(c("--minimumUPC"), action="store",
              default = 0.01, type='numeric',
              help="Minimum value"),
  make_option(c("--OVERWRITE"), action="store",
              default = FALSE, type='logical',
              help="Whether downstream files must be over-written.")
)


opt = parse_args(OptionParser(option_list=option_list))

# Print Parameters
cat("PARAMETERS:\n",file=stdout())
for(opt.idx in 1:length(opt)) {
  cat(paste0(names(opt)[opt.idx],"\t",opt[[opt.idx]],"\n"),file = stdout())
}
cat("---\n",file=stdout())


igraph_fl <- opt$UNIVERSALNET
context <- opt$CONTEXT
ctype <- opt$CONTEXT_TYPE
upc_fl <- opt$UPC
upc_min <- opt$minimumUPC
context_fl <- gsub("UNIVERSAL",context,igraph_fl)

if(!file.exists(context_fl) | opt$OVERWRITE) {

### 2 Load UPC matrix ######
UPC <- read.table(upc_fl,sep="\t",header=TRUE,quote="",check.names = FALSE)
rownames(UPC) <- UPC$hgnc_symbol
UPC <- UPC[,which(colnames(UPC)!="hgnc_symbol")]

### 3 Extract upc values for that context #######
if(ctype=="tissue") {
  if(context=="PANCANCER") {
    upc <- apply(UPC,1,median)
  } else {
    upc <- apply(UPC[,grep(paste0("_",context,"$"),colnames(UPC))],1,median) 
  }
} else if (ctype=="cell") {
  upc <- setNames(UPC[,grep(paste0("^",context,"_"),colnames(UPC))],rownames(UPC))
} else {
  stop("[ERROR] : CONTEXT_TYPE must be 'tissue' or 'cell'\n")
}

### 4 Load igraph object ######
g <- readRDS(igraph_fl)
vertex.names <- vertex_attr(g)$name

### 5 Filter out those vertex that are not expressed #####
table(vertex.names %in% names(upc))
commongenes <- intersect(vertex.names,names(upc))
filteredout <- commongenes[upc[commongenes] < upc_min]

### 6 Remove nodes that were filtered out #####
g <- g - filteredout

### 7 Save it ######
saveRDS(g,file=context_fl)

} else {
  cat("[WARN] There is a igraph object file already for this contextualization.",file=stdout())
  cat(" You could use '--OVERWRITE TRUE' to re-do this process. Skip it...\n",file=stdout())
}
