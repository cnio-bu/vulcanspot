#!/usr/bin/env Rscript

### TITLE : Make a LoF/GoF Matrix based on multi-omics profiles from a CCL (gene expr, point mut, gene CN).
### AUTHOR : Perales-Paton, Javier - jperales@cnio.es
### LICENSE : GPL-v3
### DESCRIPTION : To define Cellular Context (tissue + GeneA altered) for the identification of Cancer Genetic Dependencies (CGD),
#               we will create a matrix of Gene Function Alteration (i.e. Loss-of-Function and Gain-of-Function) across cell lines.
#               To this aim, we have to process all the multi-omics layers:
#                 - GEP : Gene Expression Profiles in UPC values.
#                 - MUT : Point Mutations as the consequence of SNVs and indels.
#                 - GCN : Gene Copy-Number as log2(Copy-number ratio)= log2(Gene copies / 2). 
#                         We assume diploidy of the human genome.
#             The output will be a R matrix of genes X CCL (cancer cell lines), with the following values:
#                 +1 = GoF : Gain-of-Function
#                 -1 = LoF : Loss-of-Function
#                  0 = wt : Wildtype (gene does not harbor any alteration)
#             This script takes into account certain thresholds to define if a genetic alteration lead to a GoF or a LoF:
#                 - GEP : Currently, this omic layer is not included because author considered that the transcriptional active state
#                       of a given gene in a cell line, if it is wildtype for GEP and GCN; probably could be modulated to avoid
#                       cell death due to a gene dependency.
#                 - MUT : There is a list of mutation consequences that lead to a LoF (see below):
#                         "De_novo_Start_OutOfFrame","Frame_Shift_Del","Frame_Shift_Ins",
#                         "In_Frame_Del","In_Frame_Ins","Nonsense_Mutation",
#                         "Nonstop_Mutation","Splice_Site","Start_Codon_Del","Start_Codon_Ins",
#                         "Stop_Codon_Del","Stop_Codon_Ins".
#                         So if any of these mutations are in a VAF (Variant Allele Frequency) of certain threshold,
#                         we consider the mutation consequence as a LoF:
#                           - if mut consequence is LoF & VAF > 0.7 = LoF.
#                         However, missense_mutation could have the two consequences (GoF/LoF) based on the protein position 
#                         and amino-acid change affected. We use an 'Annotation list' of cancer gene role,
#                           - the Cancer Gene Census (COSMIC)
#                         to decipher the most likely consequence of missense_mutation:
#                           - if a missense_mutation affects an 'oncogene' & VAF > 0.2 = GoF.
#                           - if a missense_mutation affects a 'TSG (tumour suppressor gene)' & VAF > 0.7 = LoF.
#                 - GCN : In brief, we consider a GoF/LoF depending on the type of GNV affecting the gene:
#                           - LoF= Deep deletion (log2(CNratio)=log2(gene copies/2) <= log2(0.5/2))
#                           - GoF= High amplification (log2(CNratio)=log2(gene copies/2) >= log2(8/2))

options(stringsAsFactors = FALSE)

# What mutations are considered lof
# NOTE: For the particular case of TSG (Tumour Suppressor Genes), we also consider 'missense_mutation'
#       at a very high VAF (Variant Allele Frequency) as a LoF mutation consequence (loss-of-function).
#       For these particular cases, we add 'missense_mutation' to the vector below in a temporary variable.
LoF <- c("De_novo_Start_OutOfFrame","Frame_Shift_Del","Frame_Shift_Ins",
         "In_Frame_Del","In_Frame_Ins","Nonsense_Mutation",
         "Nonstop_Mutation","Splice_Site","Start_Codon_Del","Start_Codon_Ins",
         "Stop_Codon_Del","Stop_Codon_Ins")

# multi-omics data for all cancer cell lines
GEP.UPC <- read.table(file = "./data/CCL/GEP_hgnc_symbol_feb2014_GRCh37.UPC.tsv",
                      sep="\t",check.names = FALSE,header=TRUE)
rownames(GEP.UPC) <- GEP.UPC$hgnc_symbol
GEP.UPC <- GEP.UPC[,colnames(GEP.UPC)!="hgnc_symbol"]

GCN.log2R <- read.table(file = "./data/CCL/GCN_hgnc_symbol_feb2014_GRCh37.log2CNratio.tsv",
                        sep="\t",check.names = FALSE,header=TRUE)
rownames(GCN.log2R) <- GCN.log2R$hgnc_symbol
GCN.log2R <- GCN.log2R[,colnames(GCN.log2R)!="hgnc_symbol"]

MUT.CSQ <- read.table(file = "./data/CCL/MUT_hgnc_GRCh37.CSQ.tsv",sep="\t",
                      check.names = FALSE,header=TRUE)
rownames(MUT.CSQ) <- MUT.CSQ$hgnc_symbol
MUT.CSQ <- MUT.CSQ[,colnames(MUT.CSQ)!="hgnc_symbol"]

MUT.VAF <- read.table(file = "./data/CCL/MUT_hgnc_GRCh37.VAF.tsv",sep="\t",
                      check.names = FALSE,header=TRUE)
rownames(MUT.VAF) <- MUT.VAF$hgnc_symbol
MUT.VAF <- MUT.VAF[,colnames(MUT.VAF)!="hgnc_symbol"]

# Universe of features : all the Genes described in the different multi-omics
all.genes <- sort(unique(c(rownames(GEP.UPC),rownames(GCN.log2R),row.names(MUT.CSQ))))

# Please note that we set TRUE by default in 'MUT' because they only reported those
# genes that are altered in at least one cell line
STAT.all.genes <- data.frame("GEP"=all.genes %in% row.names(GEP.UPC),
                             "GCN"=all.genes %in% row.names(GCN.log2R),
                             "MUT"=rep(TRUE,length(all.genes)),
                             row.names=all.genes)

# Universe of features : cell lines
all.ccl <- sort(unique(c(colnames(GEP.UPC),colnames(GCN.log2R),colnames(MUT.CSQ))))

# OK, it seems that CNV data introduces contexts that musnt be included because these are replicates
individual_artifacts <- c("COLO201_COLO205",
                          "KELLY_p_CCLE_AffySNP_Oct2012_01_GenomeWideSNP_6_B01_1217624_KMS18_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE",
                          "KELLY_p_CCLE_AffySNP_Oct2012_01_GenomeWideSNP_6_B03_1217476_TALL1_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE",
                          "MEL202_EYE","MOLT3_MOLT4","NHA_HT_DD")
# Moreover, both GEP and MUT includes non-carcinogenic cell lines such as fibroblasts and matched normal samples
individual_artifacts <- unique(c(individual_artifacts,
                                grep("_(MATCHED_NORMAL_TISSUE|FIBROBLAST)$",colnames(MUT.CSQ),value=TRUE)))
individual_artifacts <- unique(c(individual_artifacts,
                                 grep("_(MATCHED_NORMAL_TISSUE|FIBROBLAST)$",colnames(GEP.UPC),value=TRUE)))


# We are going to remove these cell lines from the set of all cancer cell lines
all.ccl <- setdiff(all.ccl,individual_artifacts)

# Please note that we expect that all cell lines have at least one point mutation (MUT),
# so we considered are present those that are described in the original MAF file.
STAT.all.ccl <- data.frame("GEP"=all.ccl %in% colnames(GEP.UPC),
                           "GCN"=all.ccl %in% colnames(GCN.log2R),
                           "MUT"=all.ccl %in% colnames(MUT.CSQ),
                           row.names = all.ccl)

# Get the CGC list for annotation
CGC <- read.table("./data/Annotation/Census_allWed_Oct_4_09-53-19_2017.tsv",sep="\t",header = TRUE,fill=TRUE)
CGC.Role <- setNames(sapply(gsub(" +","",CGC$Role.in.Cancer),function(z) strsplit(z,split=",")[[1]]),CGC$Gene.Symbol)
CGC.ONC <- unlist(lapply(CGC.Role,function(z) any(z=="oncogene") & !any(z=="TSG")))
CGC.TSG <- unlist(lapply(CGC.Role,function(z) !any(z=="oncogene") & any(z=="TSG")))
CGC.AMB <- unlist(lapply(CGC.Role,function(z) any(z=="oncogene") & any(z=="TSG")))

# We are going to join Ambigous (oncogene+TSG) with TSG because we realized that many
#   genes in AMB set, such as TP53 are really more TSG rather than oncogenes... 
CGC.TSG <- sort(c(CGC.TSG,CGC.AMB))


# Create matrix of LoF/GoF
FUNC <- matrix(0,nrow=length(all.genes),ncol=length(all.ccl),dimnames=list(all.genes,all.ccl))


# Define frequencies for point mutations for both classes
LoF_VAF.mut <- 0.7
GoF_VAF.mut <- 0.2

# Initialize loop
i <- 0
for(gene in rownames(FUNC)) {
  i <- i+1;
  if( (i %% 100)==0 | (i==nrow(FUNC)) ) cat(paste0("#",i,"/",nrow(FUNC),"\n"),file=stdout())
  
  # First define which set of point mutations would lead to LoF depending on gene role
  if( gene %in% names(CGC.ONC[CGC.ONC]) ) {
    LoF.mut <- LoF
    GoF.mut <- "Missense_Mutation"
  } else if (gene %in% names(CGC.TSG[CGC.TSG]) ) {
    LoF.mut <- c("Missense_Mutation",LoF) # We add missense_mutation only for TSGs, so it will always satisfy the next if for TSGs.
    GoF.mut <- "NOTfeasible"
  } else {
    LoF.mut <- c("Missense_Mutation",LoF) # We add missense_mutation
    GoF.mut <- "NOTfeasible"
  }
  
  # We start for each omic layer
  # Then, for each omic layer, define the LoF/GoF of the gene
  ind.geneFunc <- setNames(rep(0,length=ncol(FUNC)),colnames(FUNC))
  
  # Point Mutations
  if( (gene %in% rownames(MUT.CSQ)) ) {
    MUT_ccl.common <- intersect(colnames(FUNC),colnames(MUT.CSQ))
    
    # Get mutational status for all cell lines
    ind.muts <- unlist(MUT.CSQ[gene,MUT_ccl.common])
    ind.mutsVAF <- unlist(MUT.VAF[gene,MUT_ccl.common])
    if(is.character(ind.mutsVAF)) ind.mutsVAF <- sapply(ind.mutsVAF,function(z) max(as.numeric(strsplit(z,split=",")[[1]])))
    
    # Then, define the LoF/GoF of the gene based on point mutations
    # NOTE: Why to split by ';'? The reason is that for cell lines that harbor SEVERAL alterations,
    #       we have concatenate them by semicolon chars. Otherwise there is only one and the split does
    #       not make any effect.
    MUT_ccl.LoF <- unlist(lapply(strsplit(ind.muts,split=";"),function(csq) any(csq %in% LoF.mut))) & ind.mutsVAF >= LoF_VAF.mut
    MUT_ccl.GoF <- unlist(lapply(strsplit(ind.muts,split=";"),function(csq) any(csq %in% GoF.mut))) & ind.mutsVAF >= GoF_VAF.mut
    MUT_ccl.GoF <- MUT_ccl.GoF==TRUE & MUT_ccl.LoF==FALSE # Necessary step, because multiple mutations on a gene per sample. Truncating is always gonna trucates
    
    # From boolean vector to which particular cell lines have the alteration
    # that is a vector with the names of cell lines that harbor the alteration
    MUT_ccl.LoF <- names(MUT_ccl.LoF)[MUT_ccl.LoF]
    MUT_ccl.GoF <- names(MUT_ccl.GoF)[MUT_ccl.GoF]
    
    # Assignation
    ind.geneFunc[names(ind.geneFunc) %in% MUT_ccl.LoF] <- -1
    ind.geneFunc[names(ind.geneFunc) %in% MUT_ccl.GoF] <- +1
    
    # Remove temporary variables
    rm(MUT_ccl.common,ind.muts,ind.mutsVAF,MUT_ccl.LoF,MUT_ccl.GoF)
  }
  
  # Gene Copy-Number
  if(gene %in% rownames(GCN.log2R)) {
    GCN_ccl.common <- intersect(colnames(FUNC),colnames(GCN.log2R))
    
    # Get mutational status for all cell lines
    ind.gcns <- unlist(GCN.log2R[gene,GCN_ccl.common])
    
    # Then, define the LoF/GoF
    GCN_ccl.GoF <- ind.gcns >= log2(8/2) & ind.geneFunc[GCN_ccl.common] == 0
    GCN_ccl.LoF <- ind.gcns <= log2(0.5/2) & ind.geneFunc[GCN_ccl.common] == 0
    
    # From boolean vector to which particular cell lines have the alteration,
    # that is a vector with the names of cell lines that harbor the alteration
    GCN_ccl.LoF <- names(GCN_ccl.LoF)[GCN_ccl.LoF]
    GCN_ccl.GoF <- names(GCN_ccl.GoF)[GCN_ccl.GoF]
    
    # Assignation
    ind.geneFunc[names(ind.geneFunc) %in% GCN_ccl.LoF] <- -1
    ind.geneFunc[names(ind.geneFunc) %in% GCN_ccl.GoF] <- +1
    
    # Remove temporary variables
    rm(GCN_ccl.common,ind.gcns,GCN_ccl.GoF,GCN_ccl.LoF)
    
  }
  
  # Write it in the matrix
  FUNC[gene,names(ind.geneFunc)] <- ind.geneFunc
  # DONE
}

# For aesthatic purposes, we are going to reorder the columns (biological samples) by their type
sample_type <- unlist(sapply(colnames(FUNC), function(sample_name) paste(strsplit(sample_name,split="_")[[1]][-1],collapse = "_")))
FUNC <- FUNC[,order(sample_type)]

# Recalculate again the contexts
sample_type <- unlist(sapply(colnames(FUNC), function(sample_name) paste(strsplit(sample_name,split="_")[[1]][-1],collapse = "_")))

# Remove those cases that are context=="", that are duplicated form MUT data or GCN (I checked manually this point)
# These are only 10 ccl
    # > table(sample_type=="")
    # 
    # FALSE  TRUE 
    # 1525    10
FUNC <- FUNC[,sample_type!=""]

# Save the matrix
saveRDS(FUNC,file = "./data/CCL/geneFunc.rds")

