#!/usr/bin/env Rscript

### TITLE : Make a Matrices of point mutation 
### AUTHOR : Perales-Paton, Javier - jperales@cnio.es
### LICENSE : GPL-v3
### DESCRIPTION : It makes several matrices (genes x ccl) with information related to the 
###               point mutation using the MAF file from CCLE data portal. These matrices
###               will present the same dimension, but contain different information:
###                   - comma separated point mutations per gene and per sample.
###                   - comma separated allele frequency for each mutation above.
###                   - whether it is considered deleterious (1) or not (0) by the MAF annotation

options(stringsAsFactors = FALSE)

# 1 Set of mutations whose consequence is loss-of-function (LoF) ######
LoF <- c("De_novo_Start_OutOfFrame","Frame_Shift_Del","Frame_Shift_Ins",
         "In_Frame_Del","In_Frame_Ins","Missense_Mutation","Nonsense_Mutation",
         "Nonstop_Mutation","Splice_Site","Start_Codon_Del","Start_Codon_Ins",
         "Stop_Codon_Del","Stop_Codon_Ins")


# 2. Load the Mutations table from CCLE data portal ######
if(!file.exists("./data/CCL/CCLE_DepMap_18Q2_maf_20180502.txt")) {
  cat("[INFO] GCN file from CCLE does not exist. Downloading it from CCLE data portal FTP repository...\n",
      file=stdout())
  download.file(destfile = "./data/CCL/CCLE_DepMap_18Q2_maf_20180502.txt",
                url="https://data.broadinstitute.org/ccle/CCLE_DepMap_18Q2_maf_20180502.txt")
}

cat("[INFO] Load the Mutations table from CCLE data portal.\n",file=stdout())
MUT <- read.table("./data/CCL/CCLE_DepMap_18Q2_maf_20180502.txt",sep="\t",header=TRUE)

# 3 Remove 'non-truncating' mutations ######
cat(paste0("[INFO] Retrieve mutations by their consequences (",paste(LoF,collapse = ";"),")\n"),
    file=stdout())
MUT <- MUT[MUT$Variant_Classification%in%LoF,]
# MUT <- MUT[!grepl("Silent|Start_Codon_SNP",MUT$Variant_Classification),]

MUT <- MUT[order(MUT$Hugo_Symbol,MUT$Start_position),]
rownames(MUT) <- NULL

# 4 Transform it into a matrix of genes x ccl ###########
cat("[INFO] Transform it into a matrix of genes x ccl.\n",file=stdout())
# Mutation Consequences
MUT.CSQ <- matrix(NA,nrow=length(unique(MUT$Hugo_Symbol)),ncol=length(unique(MUT$Tumor_Sample_Barcode)),
                  dimnames=list(unique(MUT$Hugo_Symbol),unique(MUT$Tumor_Sample_Barcode)))
MUT.CSQ <- as.data.frame(MUT.CSQ)

# Variant Allele Frequency
MUT.VAF <- matrix(0,nrow=length(unique(MUT$Hugo_Symbol)),ncol=length(unique(MUT$Tumor_Sample_Barcode)),
                  dimnames=list(unique(MUT$Hugo_Symbol),unique(MUT$Tumor_Sample_Barcode)))
MUT.VAF <- as.data.frame(MUT.VAF)

# Is deleterious
MUT.isDLT <- matrix(FALSE,nrow=length(unique(MUT$Hugo_Symbol)),ncol=length(unique(MUT$Tumor_Sample_Barcode)),
                    dimnames=list(unique(MUT$Hugo_Symbol),unique(MUT$Tumor_Sample_Barcode)))
MUT.isDLT <- as.data.frame(MUT.isDLT)


cat("[INFO] Processing point mutations: from tidy observation to matrix. This will take a while...\n",file=stdout())
i <- 0
for(mut.id in rownames(MUT)) {
  i <- i+1;
  if( (i %% 100)==0 | (i==nrow(MUT)) ) cat(paste0("#",i,"/",nrow(MUT),"\n"),file=stdout())
  mut.info <- MUT[mut.id,]
  mut.info.VAF <- mut.info[grep("_AC$",names(mut.info))]
  mut.info.VAF[mut.info.VAF==""] <- NA
  VAF <- median(sapply(mut.info.VAF[!is.na(mut.info.VAF)],function(z) {locus.reads <- as.numeric(strsplit(z,split = ":")[[1]]); locus.reads[1] / sum(locus.reads)}))
  CSQ <- mut.info$Variant_Classification
  isDLT <- mut.info[["isDeleterious"]]
  
  # Record the consenquence
  if( is.na(MUT.CSQ[mut.info$Hugo_Symbol,mut.info$Tumor_Sample_Barcode]) ){
    MUT.CSQ[mut.info$Hugo_Symbol,mut.info$Tumor_Sample_Barcode] <- CSQ;
  } else {
    MUT.CSQ[mut.info$Hugo_Symbol,
            mut.info$Tumor_Sample_Barcode] <- paste0(MUT.CSQ[mut.info$Hugo_Symbol,mut.info$Tumor_Sample_Barcode],";",CSQ)
  }
  
  # Record the variant allele frequency
  if( MUT.VAF[mut.info$Hugo_Symbol,mut.info$Tumor_Sample_Barcode] == 0 ){
    MUT.VAF[mut.info$Hugo_Symbol,mut.info$Tumor_Sample_Barcode] <- VAF;
  } else {
    MUT.VAF[mut.info$Hugo_Symbol,
            mut.info$Tumor_Sample_Barcode] <- mean(MUT.VAF[mut.info$Hugo_Symbol,mut.info$Tumor_Sample_Barcode],VAF)
  }
  
  # Record if the gene harbors ANY deleterious mutation
  if ( is.logical(isDLT) ) {
    if( isDLT ) {
      MUT.isDLT[mut.info$Hugo_Symbol,mut.info$Tumor_Sample_Barcode] <- isDLT
    }
  }
  
  rm(VAF,CSQ,mut.info,mut.info.VAF,isDLT)
}


MUT.CSQ[is.na(MUT.CSQ)] <- "" # All that is NA still, it is wildtype

### 5 Write these tables in TSV files #######
cat("[INFO] Saving the results in TSV files at './data/CCL/MUT_hgnc_GRCh37*'.\n",file=stdout())
write.table(cbind("hgnc_symbol"=rownames(MUT.CSQ),MUT.CSQ),
            file = paste0("./data/CCL/MUT_hgnc_GRCh37",".CSQ.tsv"),
            sep = "\t",row.names = FALSE,col.names=TRUE,quote=FALSE)

write.table(cbind("hgnc_symbol"=rownames(MUT.VAF),MUT.VAF),
            file = paste0("./data/CCL/MUT_hgnc_GRCh37",".VAF.tsv"),
            sep = "\t",row.names = FALSE,col.names=TRUE,quote=FALSE)

write.table(cbind("hgnc_symbol"=rownames(MUT.isDLT),MUT.isDLT),
            file = paste0("./data/CCL/MUT_hgnc_GRCh37",".isDLT.tsv"),
            sep = "\t",row.names = FALSE,col.names=TRUE,quote=FALSE)
