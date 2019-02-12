### TITLE : Get Contexts 
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
  make_option(c("--PROFILESMATRIX"), action="store",
              default="./data/CCL/geneFunc.rds",
              type='character',
              help="Path to RDS file with the molecular profiles from the Cancer Cell Lines"),
  make_option(c("--DEPMATRICES"), action="store",
              default="./data/CCL/portal-Avana-2018-05-10.csv:./data/CCL/portal-RNAi_merged-2018-05-10.csv",
              type='character',
              help="colon-separated string with the paths to matrices files with all the contexts of the database."),
  make_option(c("--CMAPCONTEXTS"), action="store",
              default="./data/CMap-L1000/CMap_contexts.txt",
              type='character',
              help="TXT file with the contexts included in CMap for drug repositioning (limited).")
  )

opt = parse_args(OptionParser(option_list=option_list))

# Print Parameters
cat("PARAMETERS:\n",file=stdout())
for(opt.idx in 1:length(opt)) {
  cat(paste0(names(opt)[opt.idx],"\t",opt[[opt.idx]],"\n"),file = stdout())
}
cat("---\n",file=stdout())

mol <- opt$PROFILESMATRIX
fls <- as.vector(sapply(opt$DEPMATRICES,function(z) strsplit(z,split=":")[[1]]))


### Initilize the vector which will contain the different contexts ####
contexts <- vector()

### Load the molecular profile matrix to get contexts and size by context ######
cat("[INFO]: Reading Matrix with Molecular profile of Cancer Cell Lines\n",file=stdout())
if(!file.exists(mol)) {
  stop("ERROR: Molecular profiles file does not exit.")
}
mat <- readRDS(mol)
stopifnot(is.matrix(mat) | is.data.frame(mat))
cntxts <- sapply(colnames(mat),function(z) paste(strsplit(z,split="_")[[1]][-1],collapse = "_"))

# Contexts for the vulcanspot database
contexts <- sort(unique(c(contexts,cntxts)))

# Get sizes for each context
cntxts_size <- table(cntxts)
# Get a list of cell lines
cntxts_list <- sapply(unique(cntxts),function(cntx) names(cntxts)[cntxts==cntx])

### Load the data to extract the contexts ####
cat("[INFO]: Reading matrices with the original contexts use to build the database\n",file=stdout())
cntxts_commonDEP_wPROFILE <- cntxts_list
for(fl in fls) {
  cat(paste0("[INFO]: \tReading File:",fl,"\n"),file=stdout())
  
  if(grepl("\\.rds$",fl)) {
    mat <- readRDS(fl)
  } else if(grepl("\\.csv",fl)) {
    mat <- readCSV(fl)
  }
    stopifnot(is.matrix(mat) | is.data.frame(mat))
    # Contexts in dataset (depmap)
    cntxts <- sapply(colnames(mat),function(z) paste(strsplit(z,split="_")[[1]][-1],collapse = "_"))
    
    # Contexts for the vulcanspot database
    contexts <- sort(unique(c(contexts,cntxts)))
    
    # Cell lines profiled in both genetic alterations and gene dependency in all DEPMAP datasets
    cntxts_commonDEP_wPROFILE <- lapply(cntxts_commonDEP_wPROFILE,function(ccl) intersect(ccl,colnames(mat)))
}

size_cntxts_commonDEP_wPROFILE <- unlist(lapply(cntxts_commonDEP_wPROFILE,length))


#### Reading contexts form CMAP ############
cat("[INFO]: Reading contexts form CMAP:\n",file=stdout())
cntxts_cmap <- scan(opt$CMAPCONTEXTS,what = "character")

contexts <- sort(unique(c(contexts,cntxts_cmap)))

### Fill with zeros between datasets
cntxts_size <- c(cntxts_size,
                 setNames(rep(0,sum(!contexts%in%names(cntxts_size))),setdiff(contexts,names(cntxts_size)))
                 )
size_cntxts_commonDEP_wPROFILE <- c(size_cntxts_commonDEP_wPROFILE,
                                    setNames(rep(0,sum(!contexts%in%names(size_cntxts_commonDEP_wPROFILE))),
                                             setdiff(contexts,names(size_cntxts_commonDEP_wPROFILE))))

boolean_cntxts_cmap <- setNames(contexts %in% cntxts_cmap,contexts)

# Summarize PANCANCER context, a special case where everything add
cntxts_size["PANCANCER"] <- sum(cntxts_size[names(cntxts_size)!="PANCANCER"])
size_cntxts_commonDEP_wPROFILE["PANCANCER"] <- sum(size_cntxts_commonDEP_wPROFILE[names(size_cntxts_commonDEP_wPROFILE)!="PANCANCER"])
  
### Create data.frame to save it #######
contexts.df <- data.frame(id=1:length(contexts),
                           name=contexts,
                           size_genetics=cntxts_size[contexts],
                           size_bothdep=size_cntxts_commonDEP_wPROFILE[contexts],
                            in_CMap=boolean_cntxts_cmap[contexts]
                           )


# Minimum size to be included in the testing relationships
contexts.df$included <- contexts.df$size_bothdep >= 5

### Save it ####
cat("[INFO]: Save the total list of contexts:\n",file=stdout())
cat(paste0("[INFO]: \t",length(contexts),"\n"),file=stdout())
cat(contexts,sep="\n",file="./data/CCL/contexts.txt")
write.table(contexts.df,file="./data/CCL/contexts.tsv",sep="\t",col.names = TRUE,row.names = FALSE,quote=FALSE)

sessionInfo()