### TITLE : Get Genes A 
### AUTHOR: Perales-Paton, Javier - jperales@cnio.es
### LICENSE : GPL-v3


options(stringsAsFactors = FALSE)


# Functions
readCSV <- function(file,col.symbols=1) {
  DF <- read.table(file=file,header=TRUE,sep=",",check.names = FALSE)  
  rownames(DF) <- DF[,col.symbols]
  DF <- DF[,-col.symbols]
  return(DF)
}
### Load libraries #######
suppressPackageStartupMessages(require(optparse))


option_list = list(
  make_option(c("--MATRICES"), action="store",
              default="./data/CCL/geneFunc.rds:./data/CCL/portal-Avana-2018-05-10.csv:./data/CCL/portal-RNAi_merged-2018-05-10.csv",
              type='character',
              help="colon-separated string with the paths to matrices files with all the genes of the database.")
  )

opt = parse_args(OptionParser(option_list=option_list))

# Print Parameters
cat("PARAMETERS:\n",file=stdout())
for(opt.idx in 1:length(opt)) {
  cat(paste0(names(opt)[opt.idx],"\t",opt[[opt.idx]],"\n"),file = stdout())
}
cat("---\n",file=stdout())

fls <- as.vector(sapply(opt$MATRICES,function(z) strsplit(z,split=":")[[1]]))


### Load the data to extract the genes ####
cat("[INFO]: Reading matrices with the original genes use to build the database\n",file=stdout())
genes <- vector()
for(fl in fls) {
  cat(paste0("[INFO]: \tReading File:",fl,"\n"),file=stdout())
  
  if(grepl("\\.rds$",fl)) {
    mat <- readRDS(fl)
    stopifnot(is.matrix(mat))
    genes <- sort(unique(c(genes,rownames(mat))))
    
  } else if(grepl("\\.csv",fl)) {
    mat <- readCSV(fl)
    genes <- sort(unique(c(genes,rownames(mat))))
  }
}

### Save it ####
cat("[INFO]: Save the total list of genes:\n",file=stdout())
cat(paste0("[INFO]: \t-",length(genes),"\n"),file=stdout())
cat(genes,sep="\n",file="./data/CCL/genes.txt")

sessionInfo()