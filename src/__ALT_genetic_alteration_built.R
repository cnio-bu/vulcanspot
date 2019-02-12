

# UNDER CONSTRUCTION

CellSets <- readRDS(file="./data/CCL/cellsets.genetics.rds")
CellSets2 <- readRDS(file="./data/CCL/cellsets.anydep.rds")

Genes.df <- read.table("./data/Genes/genes.tsv",sep="\t",header=TRUE,quote="",fill=TRUE)
id_genes <- setNames(Genes.df$id,Genes.df$symbol)

Contexts.df <- read.table("./data/Contexts/contexts.tsv",sep="\t",header=TRUE)
id_contexts <- setNames(Contexts.df$id, Contexts.df$name)
Contexts.size_genetics <- setNames(Contexts.df$size_genetics,Contexts.df$name)

CS_context <- sapply(names(CellSets),function(CS_name) strsplit(CS_name,split="::")[[1]][1])
CS_alt <- sapply(names(CellSets),function(CS_name) strsplit(CS_name,split="::")[[1]][2])
CS_genes <- sapply(names(CellSets),function(CS_name) strsplit(CS_name,split="::")[[1]][3])
CS_freq <- unlist(lapply(CellSets,length)) / Contexts.size_genetics[CS_context]


stopifnot(all(CS_context %in% names(id_contexts)))
stopifnot(all(CS_genes %in% names(id_genes)))


CellSets.df <- data.frame(id=1:length(CellSets),
                          id_genes=id_genes[CS_genes],
                          id_contexts=id_contexts[CS_context],
                          name=names(CellSets),
                          context = CS_context,
                          alteration=CS_alt,
                          freq=CS_freq,
                          size_genetics=unlist(lapply(CellSets,length)),
                          size_anydep=unlist(lapply(CellSets2,length)))

write.table(CellSets.df, file="./data/Genetic_alterations/genetic_alterations.tsv",sep="\t",
            col.names = TRUE,row.names = FALSE,quote=FALSE)


