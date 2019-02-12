### TITLE : Generate table_gene_essentiality.tsv for the DB
### AUTHOR : Perales-Paton, Javier - jperales@cnio.es
### LICENSE : GPL-v3


options(stringsAsFactors = FALSE)

### Load libraries #######
suppressPackageStartupMessages(require(optparse))

# 1 Get parameters ############
option_list = list(
  make_option(c("--INDIR"), action="store", default="./data/GD", type='character',
              help="Input directory where RDS files with GDs are stored."),
  make_option(c("--OUTDIR"), action="store", default="./DB", type='character',
              help="Output directory"),
  make_option(c("--FDR"), action="store", default=NULL, type='numeric',
              help="FDR cutoff to retrieve GDs."),
  make_option(c("--PVAL"), action="store", default=NULL, type='numeric',
              help="P-value cutoff to retrieve GDs.")
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
FDR.cutoff <- opt$FDR
PVAL.cutoff <- opt$PVAL


# 2 Collect all the RDS file which contain GDs (different datasets, prob only two) ####
fls <- list.files(input.dir,pattern = "GD_adj\\..*\\.rds",full.names = TRUE)
# It takes the own filename to get the name of the dataset
names(fls) <- gsub("GD_adj\\.(.*)\\.rds","\\1",basename(fls))

# 3 Define the cols of interest and how to rename these cols for the final table #####
cols_interest <- c("GeneA"="GeneA","Alteration"="GeneAFunction","GeneB"="GeneB",
                   "NES"="NES","Pval"="pval","FDR"="padj",
                   "Context"="Context",
                   "GDscore"="avgGES", # Average of the Gene Essentiality Score for the set of ccl with GeneA altered
                   "GeneA.Role"="GeneA.Role","GeneB.Role"="GeneB.Role",
                   "Dataset"="Dataset")

# 4 Load data #######
for(dataset in names(fls)) {
  fl <- fls[dataset]
  GD <- readRDS(fl)
  
  # Add the dataset of origin
  GD$Dataset <- dataset
  
  # Transform into a data.frame
  GD <- as.data.frame(GD)
  
  # Order by significance
  GD <- GD[order(GD$padj,decreasing = FALSE),]
  
  # Filter-out the NES >= 0 : those cases mean that the Alteration on GeneA promotes
  # resistance to the functionality of GeneB. This is biological interesting, but not
  # the goal of the tool
  GD <- GD[which(GD$NES <= 0),]
  
  # Filter-out by significance
  # by FDR cutoff
  if(!is.null(FDR.cutoff)) {
    GD <- GD[GD$padj < FDR.cutoff,]
  }
  
  # by P-value cutoff
  if(!is.null(PVAL.cutoff)) {
    GD <- GD[GD$pval < PVAL.cutoff,]
  }
  
  # Reformat the original table
  GD <- GD[,cols_interest]
  colnames(GD) <- names(cols_interest)
  
  # The first dataset must start and write the header of the table
  first_bool <- ifelse(which(names(fls)==dataset)==1, TRUE, FALSE)
  
  # Write the table
  write.table(GD,file=paste0(output.dir,"/table_gene_essentiality.tsv"),sep="\t",
              col.names=first_bool,row.names = FALSE,quote=FALSE,
              append = !first_bool)
  
  # Remove it
  rm(GD)
}
