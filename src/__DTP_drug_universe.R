options(stringsAsFactors = FALSE)

PANDRUGS_cols <- c("name"="common_name", "internal_id"="pd_id")
PanDrugs <- read.table("./data/Annotation/PD_table.tsv",sep="\t",header=TRUE,quote="")
PanDrugs <- PanDrugs[,PANDRUGS_cols]
colnames(PanDrugs) <- names(PANDRUGS_cols)



LINCS_cols <- c("name"="common_name", "internal_id"="sig_id")
LINCS <- read.table("./data/Annotation/sig_id_table_LINCS.tsv",sep="\t",header=TRUE,quote="")
LINCS$sig_id <- paste0("sig_",LINCS$sig_id)
LINCS <- LINCS[,LINCS_cols]
colnames(LINCS) <- names(LINCS_cols)

# Merge both
drugs <- rbind(PanDrugs,LINCS)
drugs_unique <- sort(unique(c(PanDrugs$name,LINCS$name)))

drugsID <- setNames(1:length(drugs_unique),drugs_unique)
drugs$id <- drugsID[drugs$name]

drugs.df <- data.frame(id=drugsID,
                       name=names(drugsID))

# Write
saveRDS(drugs,file="./data/DrugPrescription/drugs.rds")
write.table(drugs.df,file="./data/DrugPrescription/drugs.tsv",sep="\t",row.names = FALSE,col.names = TRUE,quote=FALSE)
