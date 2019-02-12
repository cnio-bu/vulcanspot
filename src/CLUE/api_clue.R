source("/drives/slave/TBU/LINCS_RawData/Functions_l1ktools.R")

# GETquery.api.CLUE (v.1.0) ################
## Description ::: QUERY to the https://api.clue.io/ via RESTful service.
###   By default, it will try to perform the query and retrieve if there 
###     was ANY response from the RESTful service. (status of the query).
###     If there was no response, it could be any of the following situations:
###       a) an Authentification failed because the user_key from the query is not correct.
###       b) an Authentification failed because the server banned the user due to massive 
###       queries in a short time (>2500 per hour). For instance, this happens if you loop 
###       too many small queries without any delay.
###       c) The Service is temporaly down.
## Parameters
### query.url : (character) The complete URL for the query RESTful (including key_user).
### tryGETresponse : (boolean) Choose if it try to get any response from the server to 
###                     avoid missing results due to beeing banned
### tryNtimesIfFailed : (positive integer) How many times it must try to get any response if it fail.
### tryAgainAfter_Xmins : (positive numeric) How long it must wait until the next try.
## Dependencies:
###   library(rjson) # Read json into a list in R
###   require(RCurl) # To get text from URL as json object
GETquery.api.CLUE <- function(query.url,tryGETresponse=TRUE,tryNtimesIfFailed=7,tryAgainAfter_Xmins=30) {
  require(rjson)
  require(RCurl)
  # Parameters sanity check:
  if(!is.character(query.url) | length(query.url)!=1)
    stop("ERROR #1 : The 'query.url' parameter must be a character vector of length 1");
  if(!is.logical(tryGETresponse) | length(tryGETresponse)!=1)
    stop("ERROR #2 : The 'tryGETresponse' parameter must be a logical of length 1");
  if(tryNtimesIfFailed%%1!=0 | length(tryNtimesIfFailed)!=1 | tryNtimesIfFailed<0)
    stop("ERROR #3 : The 'tryNtimesIfFailed' parameter must be an positive integer of length 1");
  if(!is.numeric(tryAgainAfter_Xmins) | length(tryAgainAfter_Xmins)!=1)
    stop("ERROR #4 : The 'tryAgainAfter_Xmins' parameter must be a positive number of length 1");
  
  
  # Using tryGETresponse, it avoids missing results because of massive queries to the API.
  # So, if the query fails, it will wait 'tryAgainAfter_Xmins' seconds until tryNtimesIfFailed gets 0.
  info.got <- vector()
  if(tryGETresponse) {
    tryAgain <- TRUE # First time
    while(tryAgain) {
      #It launch the query
      status <- tryCatch(
        info.got <- fromJSON(getURLContent(query.url, customrequest = "GET")),
        error = function(e) e
      )
      if(inherits(status,  "error")) {
        warning("api.clue.io's API did NOT response");
        if(tryNtimesIfFailed!=0) {
          tryNtimesIfFailed <- tryNtimesIfFailed-1;
          print(paste0("Go to sleep and try Again after ",tryAgainAfter_Xmins,"minutes"));
          Sys.sleep(60*tryAgainAfter_Xmins);
        } else {
          tryAgain <- FALSE; # Get out, it failed but the user do not want to try again
        }
      } else { # It works out well
        tryAgain <- FALSE; # Get out
      }
    } # END of while
  } else { # If tryGETresponse is disabled, just do the query:
    info.got <- fromJSON(getURLContent(query.url, customrequest = "GET"));
  }
  
  # it returns a list, empty if the query fails or got 0 results, or with something
  return(info.got); 
}


# BUILDquery.filter ############# 
## Description: Build a Query for GETquery.api.CLUE
###   Generate a URL to be used by GETquery.api.CLUE (i.e. query.url) using
###   the normal feature (w/ or wo/ filtering). It is useful for listing documents.
BUILDquery.filter <- function(service="profiles",
                              where=c("pert_type"="trt_cp"),
                              fields=NULL,
                              features=c("limit"=1000,"skip"=0),
                              user_key="9332b818ed68dd90d5915a83d5cf753a") {
  # Internal variables
  web.url <- "https://api.clue.io/api/";
  
  # Process params
  
  if (is.list(where)) {
    # Remove NULL items
    where <- where[!unlist(lapply(where,is.null))]
    
    where.url_items <- sapply(1:length(where), function(z) {
      if(length(where[[z]])==1) {
        paste0("%22",names(where)[z],"%22:%22",where[[z]],"%22")    
      } else if (length(where[[z]])>1) {
        paste0("%22",names(where)[z],"%22:{%22inq%22:[",paste0("%22",where[[z]],"%22",collapse=","),"]}")
      }
    })
  } else if(is.vector(where)) {
    where.url_items <- sapply(1:length(where),function(z) paste0("%22",names(where)[z],"%22:%22",where[z],"%22"))
  }
  where.url <- paste(where.url_items,collapse=",")
  
  # Process fields
  if(!is.null(fields)) {
    fields.url_items <- sapply(1:length(fields),function(z) paste0("%22",fields[z],"%22:",1))
    fields.url <- paste(fields.url_items,collapse=",")
  }
  
  # Process features
  features.url_items <- sapply(1:length(features),function(z) paste0("%22",names(features)[z],"%22:",features[z]))
  features.url <- paste(features.url_items,collapse=",")
  
  # Build the query
  query.url <- paste0(web.url,service,"?filter={%22where%22:","{",where.url,"}")
  if(!is.null(fields)) query.url <- paste0(query.url,",%22fields%22:{",fields.url,"}");
  query.url <- paste0(query.url,",",features.url)
  query.url <- paste0(query.url,"}&user_key=",user_key)
  
  
  # Return the URL for the query
  return(query.url)
}

# BUILDquery.count ############
## Description: Build a Query for GETquery.api.CLUE
###   Generate a URL to be used by GETquery.api.CLUE (i.e. query.url) using
###   the normal feature (w/ or wo/ filtering). It is useful for counting documents.
BUILDquery.count <- function(service="profiles",
                              where=c("pert_type"="trt_cp"),
                              user_key="9332b818ed68dd90d5915a83d5cf753a") {
  # Internal variables
  web.url <- "https://api.clue.io/api/";
  
  # Process params
  if (is.list(where)) {
    # Remove NULL items
    where <- where[!unlist(lapply(where,is.null))]
    
    where.url_items <- sapply(1:length(where), function(z) {
      if(length(where[[z]])==1) {
        paste0("%22",names(where)[z],"%22:%22",where[[z]],"%22")    
      } else if (length(where[[z]])>1) {
        paste0("%22",names(where)[z],"%22:{%22inq%22:[",paste0("%22",where[[z]],"%22",collapse=","),"]}")
      }
    })
  } else if(is.vector(where)) {
    where.url_items <- sapply(1:length(where),function(z) paste0("%22",names(where)[z],"%22:%22",where[z],"%22"))
  }
  where.url <- paste(where.url_items,collapse=",")

  # Build the query :
  query.url <- paste0(web.url,service,"/count?where","={",where.url,"}")
  query.url <- paste0(query.url,"&user_key=",user_key)
  
  
  # Return the URL for the query
  return(query.url)
}

# GETpert_info #############
## Description :: Get metadata from all entries for a given pert_iname/pert_type
### This function get all metadata available for a given drug. For instance, all distil_ids for DMSO.
GETpert_info <- function(pert_iname=NULL,pert_ids=NULL,where=list("pert_type"="trt_cp"),
                         fields=NULL,metadata=NULL,
                         user_key="lincsdemo",waitXsecs=0,service="profiles") {
  
  # Internal variables
  # service <- "profiles"
  
  if(!is.null(pert_iname) & !is.null(pert_ids)) {
    stop("# ERROR #1: Conflict of input parameters")
  }
  
  if(!is.null(pert_iname)) {
    q1 <- BUILDquery.filter(service = "perts",
                            where = c(where,c("pert_iname"=pert_iname)),
                            fields = c("pert_id"),
                            features = c("limit"=1000,"skip"=0),
                            user_key=user_key)
    g1 <- GETquery.api.CLUE(q1)
    pert_ids <- unlist(lapply(g1,function(z) z$pert_id))
    pert_iname <- NULL
    
    if( is.null(pert_iname) & is.null(pert_ids)) {
      stop("# ERROR #2 : The user (or sometimes the database) has specified that both pert_iname and pert_id are null. We cannot generate a signature with that.\n")
    }
    
  }
  
  # The CLUE api does not include pert_iname in the service profiles, so we need to get it from the other
  # service called perts
  
  # From old LINCScloud documentation
  block_size <- 1000;
  pert_type.opt <- data.frame(Class=c("Treatment","Treatment","Treatment","Treatment","Treatment","Treatment",
                                      "Control","Control","Control","Control","Control","Control","Control"),
                              Perturbation=c("Chemical compound","Gene knockdown","Consensus gene knockdown","Gene over expression",
                                             "Mutant gene over expression","Ligand treatment","Untreated","Consensus untreated",
                                             "Control vector","Consensus control vector","Consensus seed signature","Vehicle control",
                                             "Consensus vehicle control"),
                              pert_type=c("trt_cp","trt_sh","trt_sh.cgs","trt_oe","trt_oe.mut","trt_lig","ctl_untrt","ctl_untrt.cns",
                                          "ctl_vector","ctl_vector.cns","trt_sh.css","ctl_vehicle","ctl_vehicle.cns"))
  
  q0 <- BUILDquery.count(service = service,where = c(where,list("pert_id"=pert_ids)),user_key = user_key)
  # cat(paste0("#QUERY counts:",q0,"\n"))
  g0 <- GETquery.api.CLUE(q0)
  NumDocs <- g0$count
  
  if(NumDocs==0) stop("ERROR #1 : It seems that the perturbation does not exist.\n")
  
  cat(paste0("#","pert_id(s):",paste(pert_ids,collapse = ";")," (Count=",NumDocs,")\n"),file=stdout())

    num_blocks <- ceiling(NumDocs/block_size)
    info.got <- vector()
    if(num_blocks==1) {
      q1 <- BUILDquery.filter(service = service,
                              where = c(where,list("pert_id"=pert_ids)),
                              fields = c(fields,metadata),
                              features = c("limit"=1000,"skip"=0),
                              user_key=user_key)
      g1 <- GETquery.api.CLUE(q1)
      
      info.got <- g1
    } else {
      for(i in 1:num_blocks) {
        qi <- BUILDquery.filter(service = service,
                                where = c(where,list("pert_id"=pert_ids)),
                                fields = c(fields,metadata),
                                features = c("limit"=block_size,"skip"=(i-1)*block_size),
                                user_key=user_key)
        cat(paste0("   #",i,"/",num_blocks,"::QUERY::",qi,"\n"),file=stdout())
        gi <- GETquery.api.CLUE(qi)

        # Store the output
        info.got <- c(info.got,gi)
        
        
          # Avoid BAN before getting blocked
        if(waitXsecs!=0)
          Sys.sleep(waitXsecs) # Go to sleep some seconds
      } 
    }
  
    # If only is one metadata of interest
    if(!is.null(metadata)) info.got <- unlist(lapply(info.got,function(z) z[[metadata]]))
    
    
    return(info.got)
}

# checkIFperturbationEXISTS()
# Description
#   It will check if a pert_iname exists or not.
###   pert_iname : he pert_iname for the drug/gene symbol of interest or NULL for control experiments
###   user_key : your user key
checkIFperturbationEXISTS <- function(pert_iname=NULL,pert_ids=NULL,pert_type=NULL,cell_id=NULL,user_key) {
  # Sanity Check
  if(any(length(pert_iname)>1) )
    stop("ERROR #1 : just 1 argument for parameter")
  if(is.null(pert_iname) & is.null(pert_ids)) {
    if(!(pert_type%in%c("ctl_vehicle","ctl_vector","ctl_untrt")))
      stop("ERROR #2 : If pert_iname is null. It only allows all distil_ids which were control experiments.")
  }

  # Internal variables
  service <- "profiles"
  
  if(!is.null(pert_iname) & !is.null(pert_ids)) {
    stop("Conflict of input parameters")
  }
  
  if(!is.null(pert_iname)) {
    q1 <- BUILDquery.filter(service = "perts",
                            where = c("pert_iname"=pert_iname),
                            fields = c("pert_id"),
                            features = c("limit"=1000,"skip"=0),
                            user_key=user_key)
    g1 <- GETquery.api.CLUE(q1)
    pert_ids <- unlist(lapply(g1,function(z) z$pert_id))
    pert_iname <- NULL
  }
  
  
  q0 <- BUILDquery.count(service = service,
                         where = list("pert_type"=pert_type,"pert_id"=pert_ids,"cell_id"=cell_id),
                        user_key = user_key)
  # cat(paste0("#QUERY counts:",q0,"\n"))
  g0 <- GETquery.api.CLUE(q0)
  NumDocs <- g0$count
  
  if(NumDocs==0) {
    warning(paste0("pert_iname : ",pert_iname," does not exist for this cell line."))
    res <- FALSE
  } else {
    res <- TRUE
  }
  
  return(res)
}



# GEM() ##############
# Description:
#   This function get the Gene Expression Matrix (GEM) from the raw GCTX file for a given set of distil_ids.
#   Optionally, it can transform a GEM into Gene Set Matrix. This is done only when GeneSet.file is not NULL.
GEM <- function(GCTX="/drives/slave/TBU/LINCS_RawData/q2norm_n1328098x22268.gctx",
                distil_ids,
                probe2symbol=TRUE,
                GeneSet.file=NULL,
                ssgsea.method="gsva") {
  
  
  if(!file.exists(GCTX)) {
    stop("ERROR #1 : GCTX does not exist.\n");
  }
  
  #--- Get the gene expression matrix
  require("rhdf5")
  # Parse the colnames
  id_distil <- read.gctx.ids(GCTX, "col")
  # Get the ids for those colnames
  col_IDs <- which(id_distil %in% distil_ids)
  # Get the final matrix
  MAT.GCT <- parse.gctx(GCTX,cid = col_IDs)
  # The actual Matrix
  MAT <- MAT.GCT@mat
  # Close h5
  H5close()
  
  # CollapseRows by maxmean?
  if(probe2symbol) {
    #--- collapseRows by MaxMean
    require("hgu133a.db")
    require("WGCNA")
    gene.symbols <- as.character(mget(as.character(rownames(MAT)), hgu133aSYMBOL))
    
    if(grepl("^q2norm",basename(GCTX))) {
      chooseProbe.method <- "MaxMean"
    } else if(grepl("^zspc",basename(GCTX))) {
      chooseProbe.method <- "absMaxMean"
    } else {
      chooseProbe.method <- "MaxMean"
    }
    
    cat(paste0("#Using collapseRows (method=",chooseProbe.method,")\n"),file=stdout())
    cr_obj.MaxMean <- collapseRows(MAT,
                                   rowGroup = gene.symbols,
                                   rowID = rownames(MAT),
                                   method=chooseProbe.method)
    datETcollapsed <- cr_obj.MaxMean$datETcollapsed
    MAT <- datETcollapsed[which(rownames(datETcollapsed)!="NA"),]
    
  }
  
  # Transform Genes into Gene Sets
  if(!is.null(GeneSet.file)){
    require(GSVA)
    require(GSEABase)
    cat("Transforming Gene-level expression values into Gene-Set expression values.\n",file=stdout())
    cat(paste0("The dimension of the matrix is ncol=",ncol(MAT),". "),file=stdout())
    cat("This will take a lot of time... be patient.\n",file=stdout())
    GMT <- getGmt(GeneSet.file)
    MAT <- gsva(MAT,GMT,method=ssgsea.method,rnaseq=FALSE,parallel.sz=1)
    if(ssgsea.method=="gsva")
      MAT <- MAT$es.obs
  }
  
  return(MAT)
}

# GEM2eSet() ############
# Description
#   This function generate an entire Expression-Set with the Gene Expression Matrix (GEM) as expression and 
#   metada in phenotype data. You can acess to the data as following:
#     exprs(eSet) : Gene Expression Matrix
#     pData(eSet) : Phenotype data, usually 
#                     - treatment : trt=perturbated group; cntl=control group
#                     - det_palte : Original det_plate
# Parameters
#
#
GEM2eSet <- function(GCTX="/drives/slave/TBU/LINCS_RawData/q2norm_n1328098x22268.gctx",
                     trt.distil_ids,cntl.distil_ids,probe2symbol=TRUE,
                     GeneSet.file=NULL,ssgsea.method="gsva") {
  library(Biobase)
  
  #--- Build the Expression Set
  #-- Expression Matrix
  # Get the Gene expression matrix,
  #     if the user used a GeneSet.file!=NULL, then
  #     the Gene Expression matrix will be actually a
  #     Gene-Set Expression Matrix. This allows to use limma too. See vignette(gsva)
  mat <- GEM(GCTX = GCTX,
             distil_ids = c(trt.distil_ids,cntl.distil_ids),
             probe2symbol=probe2symbol,
             GeneSet.file=GeneSet.file,
             ssgsea.method=ssgsea.method)
  # To rename future data, the original colnames will be stored.
  originalRowNames <- rownames(mat)
  
  # # Sanity check
  # if(!any(trt.distil_ids%in%colnames(mat))) {
  #   stop("There is no samples in the treatment group (trt)");
  # }
  
  #-- Phenotype Data
  # Parse some metadata
  #NOTE: Since distil_ids selfcontains useful information of plate_id, we could use it without any query
  # Treatment covariate
  treatment <- as.factor(sapply(colnames(mat),function(z) ifelse(z%in%trt.distil_ids,"trt","cntl")))
  treatment <- relevel(treatment,ref="cntl")
  # plate identifier covariate
  det_plate <- as.factor(unlist(sapply(colnames(mat),function(z) strsplit(z,split=":")[[1]][1])))
  # det_plate <- unlist(sapply(colnames(mat),function(z) strsplit(z,split=":")[[1]][1]))
  
  # Cell line
  cell_id <- as.factor(unlist(sapply(colnames(mat),function(z) strsplit(z,split="_")[[1]][2])))
  
  # Sanity check
  stopifnot(all(colnames(mat)%in% c(trt.distil_ids,cntl.distil_ids)))
  stopifnot(!any(trt.distil_ids%in%cntl.distil_ids))
  
  # Phenotype data
  pheno <- data.frame("treatment"=treatment,
                      "cell_id"=cell_id,
                      "det_plate"=det_plate,
                      row.names=colnames(mat),
                      stringsAsFactors = TRUE)
  
  #--- Create the eSet
  colnames(mat) <-  rownames(pheno)
  rownames(mat) <- make.names(rownames(mat))
  
  mat.es =  ExpressionSet(assayData=mat, phenoData = new("AnnotatedDataFrame", data=pheno))
  
  featureNames(mat.es) <- originalRowNames
  
  return(mat.es)
}

# eSet_woBatchEffect ######
eSet_woBatchEffect <- function(eSet,method="limma") {
  if(is(eSet)[1]!="ExpressionSet") {
    stop("ERROR #1 : The input must be a expression set");
  }
  
  MAT <- exprs(eSet)
  
  if(method=="limma") {
    # Design
    if(length(levels(pData(eSet)$treatment))>1) {
      treatment <- pData(eSet)$treatment
      des <- model.matrix(~treatment)
    } else {
      des <- model.matrix(~1,data=pData(eSet))
    }
    batch <- pData(eSet)$det_plate
    
    MAT2 <- removeBatchEffect(MAT,design = des,batch = batch)
  } else if (method=="combat") {
    require(sva)
    modcombat <- model.matrix(~1,data=pData(eSet))
    batch <- pData(eSet)$det_plate
    
    MAT2 = ComBat(dat=MAT, batch=batch, mod=modcombat)
  }

  eSet2 <- eSet
  exprs(eSet2) <- MAT2
  return(eSet2)
}

# limma_LVL3() ##################
# Description :
#   This function peforms a limma comparison (see `limmaUsersGuide()`) for a trt Vs. cntl experiment from l1k.
#   It works with small molecule compounds, knock-down, and over-expression perturbations. The only point
#     is that you need to specify which distil_ids are for trt condition, and which ones are for cntl (CONTROL).
#   This function only works for the level 3 data from LINCS L1K
# Arguments :
#   GCTX : GCTX file with the level 4 data.
#   trt.distil_ids : TREATED CONDITION - vector of distil_ids which will be matched in the GCTX expression matrix.
#   cntl.distil_ids : CONTROL CONDITION - vector of distil_ids which will be matched in the GCTX expression matrix.
#   GeneSet.file : The GMT file which contains the gene set to be transformed into from the Gene Expression matrix
#   ssgsea.method : choose a ssgsea method between gsva or ssgsea
limma_LVL3 <- function(GCTX="/drives/slave/TBU/LINCS_RawData/q2norm_n1328098x22268.gctx",
                       trt.distil_ids,cntl.distil_ids,probe2symbol=TRUE,
                       GeneSet.file=NULL,ssgsea.method="gsva") # Only if it must be transformed into gene sets
{
  require(limma)
  
  # Get the Gene expression matrix,
  #     if the user used a GeneSet.file!=NULL, then
  #     the Gene Expression matrix will be actually a 
  #     Gene-Set Expression Matrix. This allows to use limma too. See vignette(gsva)
  mat <- GEM(GCTX = GCTX,distil_ids = c(trt.distil_ids,cntl.distil_ids),
             probe2symbol=probe2symbol,
             GeneSet.file=GeneSet.file,
             ssgsea.method=ssgsea.method)
  
  # Sanity check
  if(!any(trt.distil_ids%in%colnames(mat))) {
    stop("There are no samples in the treatment group (trt)");
  }
  
  if(!any(cntl.distil_ids%in%colnames(mat))) {
    stop("There are no samples in the control group (cntl)")
  }
  
  #--- Parse some metadata
  #NOTE: Since distil_ids selfcontains useful information of plate_id, we could use it without any query
  # Treatment covariate
  treatment <- as.factor(sapply(colnames(mat),function(z) ifelse(z%in%trt.distil_ids,"trt","cntl")))
  treatment <- relevel(treatment,ref="cntl")
  # plate identifier covariate
  det_plate <- as.factor(unlist(sapply(colnames(mat),function(z) strsplit(z,split=":")[[1]][1])))
  # det_plate <- unlist(sapply(colnames(mat),function(z) strsplit(z,split=":")[[1]][1]))
  
  # Sanity check
  stopifnot(all(colnames(mat)%in% c(trt.distil_ids,cntl.distil_ids)))
  stopifnot(!any(trt.distil_ids%in%cntl.distil_ids))
  
  #--- Build the design model
  # Generate the formula
  form.str <- "~ treatment";
  if(length(levels(det_plate))>1) 
    form.str <- paste(form.str,"det_plate",sep=" + ");
  # [... you can add some covariates into the model using the code above]
  des <- model.matrix(as.formula(form.str))
  
  # VERBOSE
  cat(paste0("#Genome-wide Differential Gene Expression Analysis verbose","\n"),file=stdout())
  cat(paste0("\t Group A (#trt=",sum(treatment=="trt"),")\n"),file=stdout())
  cat(paste0("\t Group B (#cntl=",sum(treatment=="cntl"),")\n"),file=stdout())
  cat(paste0("\t Block Factor 'det_plate' (#unique det_plate=",length(levels(det_plate)),")\n"),file=stdout())
  cat(paste0("\t Formula for Design Matrix: ",form.str,"\n"),file=stdout())
  cat(paste0("\t Coefficient of Contrast: ",colnames(des)[2],"\n"),file=stdout())
  
  # Perform the gene-wide limma test
  fit <- lmFit(mat,design = des)
  fit.cont <- contrasts.fit(fit,coefficients = 2) # NOTE: Coef=2 is the 'treatment' covariate, 
  #   because the 1st is the intercept
  eBay <- eBayes(fit.cont)
  
  return(eBay)
}

# GETeSet()
# Description : pending
# Parameters : pending
GETeSet <- function(pert_iname=NULL,user_key,pert_id=NULL,
                    trt.pert_type="trt_cp",
                    cntl.pert_type="ctl_vehicle",
                    cell_ids=NULL,
                    probe2symbol=TRUE,
                    GeneSet.file=NULL,ssgsea.method="ssgsea") {
  
  require("rhdf5")
  require("WGCNA")
  
  #--- Generate the signature
  cat("GENERATING Expression-Set(s) for a given perturbation:\n",file=stdout())
  cat(paste0("\tpert_iname:\t",pert_iname,"\n"),file=stdout())
  
  # There is a common block in this function which consists in getting the distil_ids
  if(!is.null(trt.pert_type)) {
    #   for the Treatment Group (TRT).
    TRT.where <- list("pert_type"=trt.pert_type)
    if(!is.null(cell_ids)) TRT.where <- c(TRT.where,list("cell_id"=cell_ids))
    
    TRT.info.got<- GETpert_info(pert_iname = pert_iname,pert_ids = pert_id,
                                where = TRT.where,fields = c("det_plate","distil_id"),
                                metadata = NULL,user_key = user_key)
    TRT.distil_ids <- unlist(lapply(TRT.info.got,function(z) z$distil_id))
    TRT.det_plates <- unlist(lapply(TRT.info.got,function(z) z$det_plate))
  } else {
    TRT.distil_ids <- NULL
  }
  
  # GET CONTROL distil_ids
  CNTL.where <- list("pert_type"=cntl.pert_type)
  if(!is.null(trt.pert_type)) CNTL.where <- c(CNTL.where,list("det_plate"=TRT.det_plates));
  if(!is.null(cell_ids)) CNTL.where <- c(CNTL.where,list("cell_id"=cell_ids));
  
  CNTL.distil_ids<- GETpert_info(pert_iname = NULL,pert_ids = NULL,
                                 where = CNTL.where,fields = NULL,
                                 metadata = "distil_id",user_key = user_key)
  
  #--- Getting the Expression Set
    # Just one eSet
    eSet <- GEM2eSet(trt.distil_ids = TRT.distil_ids,
                     cntl.distil_ids = CNTL.distil_ids,
                     probe2symbol=probe2symbol,
                     GeneSet.file=GeneSet.file,ssgsea.method=ssgsea.method)
  
  return(eSet);
} # END f(x)

# GETsignature.LIMMA_LVL3 ###############
# This is a wrapper function for getting signatures by LIMMA_LVL3.
# Please, see the following functions:
#     - GET_SMCsignature()
#     - GET_KDsignature()
#     - GET_OEsignature()
GETsignature.LIMMA_LVL3 <- function(pert_iname=NULL,pert_id=NULL,
                                    trt.pert_type="trt_cp",
                                    cntl.pert_type="ctl_vehicle",
                                    cell_ids=NULL,
                                    probe2symbol=TRUE,
                                    user_key,
                                    GeneSet.file=NULL,
                                    ssgsea.method="ssgsea") {
  
  rawsignatures.type <- "LIMMA_LVL3"
  
  #--- Generate the signature
  cat("GENERATING signatures by LIMMA\n",file=stdout())
  require(limma)
  require(WGCNA)
  
  # There is a common block in this function which consists in getting the distil_ids
  #   for the Treatment Group (TRT).
  TRT.where <- list("pert_type"=trt.pert_type)
  if(!is.null(cell_ids)) TRT.where <- c(TRT.where,list("cell_id"=cell_ids))
  
  TRT.info.got<- GETpert_info(pert_iname = pert_iname,pert_ids = pert_id,
                              where = TRT.where,fields = c("det_plate","distil_id"),
                              metadata = NULL,user_key = user_key)
  TRT.distil_ids <- unlist(lapply(TRT.info.got,function(z) z$distil_id))
  TRT.det_plates <- unlist(lapply(TRT.info.got,function(z) z$det_plate))
  
  # GET CONTROL distil_ids
  CNTL.where <- list("pert_type"=cntl.pert_type,"det_plate"=TRT.det_plates)
  if(!is.null(cell_ids)) CNTL.where <- c(CNTL.where,list("cell_id"=cell_ids))
  
  CNTL.distil_ids<- GETpert_info(pert_iname = NULL,pert_ids = NULL,
                                 where = CNTL.where,fields = NULL,
                                 metadata = "distil_id",user_key = user_key)
  
  #--- Perform eBayes
  # Just a limma including all pert_ids
  eBay <- limma_LVL3(trt.distil_ids = TRT.distil_ids,
                     cntl.distil_ids = CNTL.distil_ids,
                     probe2symbol = probe2symbol,
                     GeneSet.file=GeneSet.file,ssgsea.method=ssgsea.method)
  # # Sanity check
  # if(length(list.eBays)==0) {
  #   stop(paste0("There are no eBayes generated for this pert_iname:",pert_iname,pert_id))
  # }
  
  return(eBay);
} # END f(x)


# GETsignature.default ::
# This is a wrapper function implemented with the aim of adding more rawsignatures.types in the future.
# Please, see the following functions:
#     - GET_SMCsignature()
#     - GET_KDsignature()
#     - GET_OEsignature()
GETsignature.default <- function(pert_iname=NULL,pert_id=NULL,
                                 rawsignatures.type=c("LIMMA_LVL3","PRL_LVL4"),
                                 user_key,
                                 trt.pert_type="trt_cp",
                                 cntl.pert_type="ctl_vehicle",
                                 cell_ids=NULL,
                                 probe2symbol=TRUE,
                                 GeneSet.file=NULL,ssgsea.method="gsva") {
  # Load libraries
  require(GeneExpressionSignature)
  require(GSEABase)
  require(Biobase)
  require(RCurl)
  require(rjson)
  require(rhdf5) # Access to GCTX files
  require("hgu133a.db")
  require(AnnotationDbi)
  require(annotate)
  #
  
  if(rawsignatures.type=="LIMMA_LVL3") {
    eBay <- GETsignature.LIMMA_LVL3(pert_iname=pert_iname,
                                         pert_id=pert_id,
                                         user_key=user_key,
                                         trt.pert_type=trt.pert_type,
                                         cntl.pert_type=cntl.pert_type,
                                         cell_ids=cell_ids,
                                         probe2symbol=probe2symbol,
                                         GeneSet.file=GeneSet.file,ssgsea.method=ssgsea.method)
    
  } else if (awsignatures.type=="PRL_LVL4"){
    stop("DEPRECATED PRL FUNCTION")
    if(!is.null(GeneSet.file))
      stop("ERROR #2 : GeneSet transformation is only available when LIMMA_LVL3");
    # OBJ.final <- GETsignature.PRL_LVL4(pert_iname=pert_iname,
    #                                    Ngenes=Ngenes,user_key=user_key,
    #                                    trt.pert_type=trt.pert_type,
    #                                    cell_ids=cell_ids,
    #                                    stratifyBy.pert_id=stratifyBy.pert_id,
    #                                    return.GSC=return.GSC,
    #                                    probe2symbol=probe2symbol)
  }
  
  return(eBay)
}


# GET_SMCsignature() ::
# Description :
#   This function executes a battery of primary functions in order to get a gene expression signature for a
#   Small Molecule Compound
#   It can be used by following two approaches: limma using lvl3, and PRL using lvl4.
# Arguments :
#   drug : The pert_iname for the knocked-down gene. It must be a gene symbol.
#   rawsignatures.type : "LIMMA_LVL3" or "PRL_LVL4", depending on the method to get the signature.
#   Ngenes : Number of genes in the final rank aggregated collapsed gene sets.
#   stratifyBy.pert_id : logical. If different pert_id for a given pert_iname must be gathered individually.
#   cntl_cache.distil_ids : file. The RDS file of control distil_ids.
#   probe2symbol : logical. If probes must be collapse to symbols.
#   user_key : Your user key for API LINCS L1k

GET_SMCsignature <- function(drug=NULL,pert_id=NULL,
                             rawsignatures.type=c("LIMMA_LVL3","PRL_LVL4"),
                             user_key,
                             cell_ids=NULL,
                             probe2symbol=TRUE,
                             GeneSet.file=NULL,ssgsea.method="gsva") {
  
  eBay <- GETsignature.default(pert_iname=drug,
                                    pert_id=pert_id,
                                    rawsignatures.type=rawsignatures.type,
                                    user_key=user_key,
                                    trt.pert_type="trt_cp",cntl.pert_type="ctl_vehicle",
                                    cell_ids=cell_ids,
                                    probe2symbol=probe2symbol,
                                    GeneSet.file=GeneSet.file,ssgsea.method=ssgsea.method)
  return(eBay)  
}


# GET_KDsignature() ::
# Description :
#   This function executes a battery of primary functions in order to get a gene expression signature for a
#   knock-down gene.
#   It can be used by following two approaches: limma using lvl3, and PRL using lvl4.
# Arguments :
#   symbol : The pert_iname for the knocked-down gene. It must be a gene symbol.
#   rawsignatures.type : "LIMMA_LVL3" or "PRL_LVL4", depending on the method to get the signature.
#   Ngenes : Number of genes in the final rank aggregated collapsed gene sets.
#   stratifyBy.pert_id : logical. If different pert_id for a given pert_iname must be gathered individually.
#   cntl_cache.distil_ids : file. The RDS file of control distil_ids.
#   probe2symbol : logical. If probes must be collapse to symbols.
#   user_key : Your user key for API LINCS L1k

GET_KDsignature <- function(symbol,
                            pert_id=NULL,
                            rawsignatures.type=c("LIMMA_LVL3","PRL_LVL4"),
                            user_key,
                            cell_ids=NULL,
                            probe2symbol=TRUE,
                            GeneSet.file=NULL,ssgsea.method="gsva") {
  
  OBJ.final <- GETsignature.default(pert_iname=symbol,
                                    pert_id=pert_id,
                                    rawsignatures.type=rawsignatures.type,
                                    user_key=user_key,
                                    trt.pert_type="trt_sh",cntl.pert_type="ctl_vector",
                                    cell_ids=cell_ids,
                                    probe2symbol=probe2symbol,
                                    GeneSet.file=GeneSet.file,ssgsea.method=ssgsea.method)
  return(OBJ.final)  
}

# GET_OEsignature() ::
# Description :
#   This function executes a battery of primary functions in order to get a gene expression signature for a
#   over-expression gene.
#   It can be used by following two approaches: limma using lvl3, and PRL using lvl4.
# Arguments :
#   symbol : The pert_iname for the over-expressed gene. It must be a gene symbol.
#   rawsignatures.type : "LIMMA_LVL3" or "PRL_LVL4", depending on the method to get the signature.
#   Ngenes : Number of genes in the final rank aggregated collapsed gene sets.
#   stratifyBy.pert_id : logical. If different pert_id for a given pert_iname must be gathered individually.
#   cntl_cache.distil_ids : file. The RDS file of control distil_ids.
#   probe2symbol : logical. If probes must be collapse to symbols.
#   user_key : Your user key for API LINCS L1k

GET_OEsignature <- function(symbol,
                            pert_id=NULL,
                            rawsignatures.type=c("LIMMA_LVL3","PRL_LVL4"),
                            user_key,
                            cell_ids=NULL,
                            probe2symbol=TRUE,
                            GeneSet.file=NULL,ssgsea.method="gsva") {
  
  OBJ.final <- GETsignature.default(pert_iname=symbol,
                                    pert_id=pert_id,
                                    rawsignatures.type=rawsignatures.type,
                                    user_key=user_key,
                                    trt.pert_type="trt_oe",cntl.pert_type="ctl_untrt",
                                    cell_ids=cell_ids,
                                    probe2symbol=probe2symbol,
                                    GeneSet.file=GeneSet.file,ssgsea.method=ssgsea.method)
  return(OBJ.final)  
}

#### PERMUTATIONS -----#######
# limma_LVL3.Npermut
# Description :
#   This function peforms a permutation of N limma comparison (see `limmaUsersGuide()`) for a trt Vs. cntl experiment from l1k.
#   It works with small molecule compounds, knock-down, and over-expression perturbations. The only point
#     is that you need to specify which distil_ids are for trt condition, and which ones are for cntl (CONTROL).
#   This function only works for the level 3 data from LINCS L1K
# Arguments :
#   GCTX : GCTX file with the level 4 data.
#   trt.distil_ids : TREATED CONDITION - vector of distil_ids which will be matched in the GCTX expression matrix.
#   cntl.distil_ids : CONTROL CONDITION - vector of distil_ids which will be matched in the GCTX expression matrix.
#   no.permutations : Number of permutations to perform
limma_LVL3.Npermut <- function(GCTX="/drives/slave/TBU/LINCS_RawData/q2norm_n1328098x22268.gctx",
                               trt.distil_ids,cntl.distil_ids,probe2symbol=TRUE,
                               GeneSet.file=NULL,ssgsea.method="gsva", # Only if it must be transformed into gene sets
                               no.permutations=1000,statistic.name="t") 
{
  require(limma)
  
  # Get the Gene expression matrix,
  #     if the user used a GeneSet.file!=NULL, then
  #     the Gene Expression matrix will be actually a 
  #     Gene-Set Expression Matrix. This allows to use limma too. See vignette(gsva)
  mat <- GEM(GCTX = GCTX,distil_ids = c(trt.distil_ids,cntl.distil_ids),
             probe2symbol=probe2symbol,
             GeneSet.file=GeneSet.file,
             ssgsea.method=ssgsea.method)
  
  # Sanity check
  if(!any(trt.distil_ids%in%colnames(mat))) {
    stop("There is no samples in the treatment group (trt)");
  }
  
  #--- Parse some metadata
  #NOTE: Since distil_ids selfcontains useful information of plate_id, we could use it without any query
  # Treatment covariate
  treatment.original <- as.factor(sapply(colnames(mat),function(z) ifelse(z%in%trt.distil_ids,"trt","cntl")))
  treatment.original <- relevel(treatment.original,ref="cntl")
  # plate identifier covariate
  det_plate <- as.factor(unlist(sapply(colnames(mat),function(z) strsplit(z,split=":")[[1]][1])))
  # det_plate <- unlist(sapply(colnames(mat),function(z) strsplit(z,split=":")[[1]][1]))
  
  # Sanity check
  stopifnot(all(colnames(mat)%in% c(trt.distil_ids,cntl.distil_ids)))
  stopifnot(!any(trt.distil_ids%in%cntl.distil_ids))
  
  #--- Build the design model
  # Generate the formula
  form.str <- "~ treatment";
  if(length(levels(det_plate))>1) 
    form.str <- paste(form.str,"det_plate",sep=" + ");
  # [... you can add some covariates into the model using the code above]
  
  
  # Create the matrix which will contain the statistics from the permutation
  PERMmat <- matrix(NA,nrow=nrow(mat),ncol=no.permutations,dimnames=list(rownames(mat)))
  
  cat(paste0("Performing n=",no.permutations," permutations. Statistic:",statistic.name),
      sep="\n",file=stdout())
  for(permut.idx in 1:no.permutations) {
    # For each loop, sampling the trt/cntl tags. Then generate the design matrix
    treatment <- sample(treatment.original)
    names(treatment) <- names(treatment.original)
    
    des <- model.matrix(as.formula(form.str))
    
    # Perform the gene-wide limma test
    fit <- lmFit(mat,design = des)
    fit.cont <- contrasts.fit(fit,coefficients = 2) # NOTE: Coef=2 is the 'treatment' covariate, 
    #   because the 1st is the intercept
    # eBay <- eBayes(fit.cont)
    
    PERMmat[,permut.idx] <- eBayes(fit.cont)[[statistic.name]]
    
  }
  
  return(PERMmat)
}




# GETpermutSignature.LIMMA_LVL3 ::
# This is a wrapper function for getting signatures by LIMMA_LVL3.Npermut
# Please, see the following functions:
#     - GET_SMCpermutSignature()
#     - GET_KDpermutSignature()
#     - GET_OEspermutSignature()
GETpermutSignature.LIMMA_LVL3 <- function(pert_iname=NULL,pert_id=NULL,
                                          user_key,
                                          trt.pert_type="trt_cp",
                                          cntl.pert_type="ctl_vehicle",
                                          cell_ids=NULL,
                                          probe2symbol=TRUE,ssgsea.method="ssgsea",GeneSet.file=NULL,
                                          no.permutations=1000) {
  
  rawsignatures.type <- "LIMMA_LVL3"
  
  #--- Generate the signature
  cat("GENERATING signatures by LIMMA\n",file=stdout())
  require(limma)
  require(WGCNA)
  
  # There is a common block in this function which consists in getting the distil_ids
  #   for the Treatment Group (TRT).
  TRT.where <- list("pert_type"=trt.pert_type)
  if(!is.null(cell_ids)) TRT.where <- c(TRT.where,list("cell_id"=cell_ids))
  
  TRT.info.got<- GETpert_info(pert_iname = pert_iname,pert_ids = pert_id,
                              where = TRT.where,fields = c("det_plate","distil_id"),
                              metadata = NULL,user_key = user_key)
  TRT.distil_ids <- unlist(lapply(TRT.info.got,function(z) z$distil_id))
  TRT.det_plates <- unlist(lapply(TRT.info.got,function(z) z$det_plate))
  
  # GET CONTROL distil_ids
  CNTL.where <- list("pert_type"=cntl.pert_type,"det_plate"=TRT.det_plates)
  if(!is.null(cell_ids)) CNTL.where <- c(CNTL.where,list("cell_id"=cell_ids))
  
  CNTL.distil_ids<- GETpert_info(pert_iname = NULL,pert_ids = NULL,
                                 where = CNTL.where,fields = NULL,
                                 metadata = "distil_id",user_key = user_key)
  
  #--- Perform eBayes
    # Just a limma including all pert_ids
    PERMmat <- limma_LVL3.Npermut(trt.distil_ids = TRT.distil_ids,
                                  cntl.distil_ids = CNTL.distil_ids,
                                  probe2symbol = probe2symbol,
                                  GeneSet.file=GeneSet.file,ssgsea.method=ssgsea.method,
                                  no.permutations=no.permutations)


  return(PERMmat);
} # END f(x)

# GETpermutSignature.default ::
# This is a wrapper function implemented with the aim of adding more rawsignatures.types in the future.
# Please, see the following functions:
#     - GET_SMCsignature()
#     - GET_KDsignature()
#     - GET_OEsignature()
GETpermutSignature.default <- function(pert_iname=NULL,pert_id=NULL,
                                       rawsignatures.type=c("LIMMA_LVL3","PRL_LVL4"),
                                       user_key,
                                       trt.pert_type="trt_cp",
                                       cntl.pert_type="ctl_vehicle",
                                       cell_ids=NULL,
                                       probe2symbol=TRUE,
                                       GeneSet.file=NULL,ssgsea.method="gsva",
                                       no.permutations=1000) {
  # Load libraries
  require(GeneExpressionSignature)
  require(GSEABase)
  require(Biobase)
  require(RCurl)
  require(rjson)
  require(rhdf5) # Access to GCTX files
  require("hgu133a.db")
  require(AnnotationDbi)
  require(annotate)
  #
  
  
  if(!is.null(GeneSet.file))
    stop("ERROR #1 : 'return.GSC'=TRUE is NOT available when the gene expression is transformed into Gene Sets.")
  
  if(rawsignatures.type=="LIMMA_LVL3") {
    MAT.final <- GETpermutSignature.LIMMA_LVL3(pert_iname=pert_iname,
                                               pert_id=pert_id,
                                               user_key=user_key,
                                               trt.pert_type=trt.pert_type,
                                               cntl.pert_type=cntl.pert_type,
                                               probe2symbol=probe2symbol,
                                               cell_ids=cell_ids,
                                               GeneSet.file=GeneSet.file,ssgsea.method=ssgsea.method,
                                               no.permutations=no.permutations)
  } else {
    # if(!is.null(GeneSet.file))
    #   stop("ERROR #2 : GeneSet transformation is only available when LIMMA_LVL3");
    # OBJ.final <- GETsignature.PRL_LVL4(pert_iname=pert_iname,
    #                                    Ngenes=Ngenes,user_key=user_key,
    #                                    trt.pert_type=trt.pert_type,
    #                                    stratifyBy.pert_id=stratifyBy.pert_id,
    #                                    return.GSC=return.GSC,
    #                                    probe2symbol=probe2symbol)
    stop("ERROR #1 under construction")
  }
  
  return(MAT.final)
}

# GET_SMCpermutSignature() ::
# Description :
#   This function executes a battery of primary functions in order to get a gene expression signature for a
#   Small Molecule Compound
#   It can be used by following two approaches: limma using lvl3, and PRL using lvl4.
# Arguments :
#   drug : The pert_iname for the knocked-down gene. It must be a gene symbol.
#   rawsignatures.type : "LIMMA_LVL3" or "PRL_LVL4", depending on the method to get the signature.
#   Ngenes : Number of genes in the final rank aggregated collapsed gene sets.
#   stratifyBy.pert_id : logical. If different pert_id for a given pert_iname must be gathered individually.
#   cntl_cache.distil_ids : file. The RDS file of control distil_ids.
#   probe2symbol : logical. If probes must be collapse to symbols.
#   user_key : Your user key for API LINCS L1k

GET_SMCpermutSignature <- function(drug=NULL,pert_id=NULL,
                                   rawsignatures.type=c("LIMMA_LVL3","PRL_LVL4"),
                                   user_key,
                                   probe2symbol=TRUE,
                                   GeneSet.file=NULL,ssgsea.method="gsva",
                                   cell_ids=NULL,
                                   no.permutations=1000) {
  
  MAT.final <- GETpermutSignature.default(pert_iname=drug,
                                          pert_id=pert_id,
                                          rawsignatures.type=rawsignatures.type,
                                          user_key=user_key,
                                          trt.pert_type="trt_cp",cntl.pert_type="ctl_vehicle",
                                          cell_ids=cell_ids,
                                          probe2symbol=probe2symbol,
                                          GeneSet.file=GeneSet.file,ssgsea.method=ssgsea.method,
                                          no.permutations=no.permutations)
  return(MAT.final)  
}

#### END of PERMUTATIONS -----#######


#### Gene Regulatory Networks -----################

### GET_msVIPER ###############
### Description : Get the Virtual Inferente of Protein Activities for a perturbation given a GRN model.
### Parameters:
### pert_iname :: pert_iname
### pert_id :: pert_ids
### trt.pert_type :: pert_type for the treatment group
### cntl.pert_type :: pert_type for the control group
### user_key :: user key for the API (api.clue.io)
### cell_ids :: cell_id(s) of interest. Null is the consensus (incl. all cell ids tested in the dataset)
### regulon :: regulon R object or file.rds with it.
### batchRemoval.method :: Method for the batch effect removal. See eSet_woBatchEffect()
### Comments: Make sure you use the ARACNe GRN appropiate for the query.

GET_msVIPER <- function(pert_iname=NULL,pert_id=NULL,
                     trt.pert_type="trt_cp",cntl.pert_type = "ctl_vehicle",
                     user_key,cell_ids=NULL,probe2symbol,regulon,batchRemoval.method="combat") {
  require(viper)
  
  # Sanity checks
  if(is(regulon)[1] != "regulon") {
    if(file.exists(regulon) & grep("\\.rds$",regulon)) {
      regulon <- readRDS(regulon)
      if(is(regulon)[1] != "regulon") {
        stop("ERROR #2 : The .rds file does not contain a regulon object\n")
      }
    } else {
      stop("ERROR #1 : The regulon object must be a .RDS file or a 'regulon' R object\n")
    }
  }
  
  if(is.null(trt.pert_type) | is.null(cntl.pert_type)) {
    stop("ERROR #3 : The msVIPER algorithm requires 2 treatment conditions.")
  }
  
  cat(paste0("#Obtaining the Expression-set object","\n"),file=stdout())
  eset <- GETeSet(pert_iname = pert_iname,user_key = user_key,pert_id = pert_id,
                  trt.pert_type = trt.pert_type,cntl.pert_type = cntl.pert_type,
                  cell_ids = cell_ids,probe2symbol = TRUE)
  
  cat(paste0("#Batch effect removal","\n"),file=stdout())
  eset_woBatch <- eSet_woBatchEffect(eSet = eset,method=batchRemoval.method)
  
  if(length(levels(pData(eset_woBatch)$treatment))!=2) {
    stop("ERROR #4 : The msVIPER algorithm requires 2 treatment conditions.")
  }

  #--- Following the vignette from VIPER paackage
  # Gerating the gene expression signature
  cat(paste0("#Generating the signature","\n"),file=stdout())
  signature <- rowTtest(eset_woBatch, "treatment", "trt", "cntl")
  signature <- (qnorm(signature$p.value/2, lower.tail = FALSE) * sign(signature$statistic))[, 1]
  
  # NULL model by sample permutations
  cat(paste0("#Generating the null model","\n"),file=stdout())
  nullmodel <- ttestNull(eset_woBatch, "treatment", "trt", "cntl", per = 1000,repos = TRUE, verbose = FALSE)
  
  # msVIPER
  cat(paste0("#Calculating VIPER","\n"),file=stdout())
  mrs <- msviper(signature,regulon,nullmodel,verbose=FALSE)
  
  # Return the object
  return(mrs)
}

###


test_TES <- function(A=list("pert_iname"="palbociclib",
                            "pert_id"=NULL,
                            "trt.pert_type"="trt_cp",
                            "cntl.pert_type"="ctl_vehicle",
                            "cell_ids"=NULL),
                     B=list("pert_iname"="CDK4",
                            "pert_id"=NULL,
                            "trt.pert_type"="trt_sh",
                            "cntl.pert_type"="ctl_vector",
                            "cell_ids"=NULL),
                     Ngenes=250,
                     # rawsignatures.type="LIMMA_LVL3",
                     probe2symbol = TRUE,
                     user_key="lincsdemo",
                     no.permutations=100,
                     show.plot=TRUE) {
  require(fgsea)
  
  # Sanity check
  if(!all(c("pert_iname","pert_id","trt.pert_type","cntl.pert_type","cell_ids") %in% names(A))) {
    stop("ERROR #1: Not all params are included in A.")
  }
  if(!all(c("pert_iname","pert_id","trt.pert_type","cntl.pert_type","cell_ids") %in% names(B))) {
    stop("ERROR #1: Not all params are included in A.")
  }
  
  # First get the signatures to be compared each other
  cat(paste0("# Generating the 'A signature' and 1000permutations by sampling samples\n"),file=stdout())
  A.eBay <- GETsignature.LIMMA_LVL3(pert_iname=A$pert_iname,
                                    pert_id=A$pert_id,
                                    user_key=user_key,
                                    trt.pert_type=A$trt.pert_type,
                                    cntl.pert_type=A$cntl.pert_type,
                                    cell_ids=A$cell_ids,
                                    probe2symbol=probe2symbol,
                                    GeneSet.file=NULL,ssgsea.method=NULL)
  
  A.p1k <- GETpermutSignature.LIMMA_LVL3(pert_iname=A$pert_iname,
                                         pert_id=A$pert_id,
                                         user_key=user_key,
                                         trt.pert_type=A$trt.pert_type,
                                         cntl.pert_type=A$cntl.pert_type,
                                         probe2symbol=probe2symbol,
                                         cell_ids=A$cell_ids,
                                         GeneSet.file=NULL,ssgsea.method=NULL,
                                         no.permutations=no.permutations)
  
  cat(paste0("# Generating the 'B signature' and 1000permutations by sampling samples\n"),file=stdout())
  B.eBay <- GETsignature.LIMMA_LVL3(pert_iname=B$pert_iname,
                                    pert_id=B$pert_id,
                                    user_key=user_key,
                                    trt.pert_type=B$trt.pert_type,
                                    cntl.pert_type=B$cntl.pert_type,
                                    cell_ids=B$cell_ids,
                                    probe2symbol=probe2symbol,
                                    GeneSet.file=NULL,ssgsea.method=NULL)
  
  B.p1k <- GETpermutSignature.LIMMA_LVL3(pert_iname=B$pert_iname,
                                         pert_id=B$pert_id,
                                         user_key=user_key,
                                         trt.pert_type=B$trt.pert_type,
                                         cntl.pert_type=B$cntl.pert_type,
                                         probe2symbol=probe2symbol,
                                         cell_ids=B$cell_ids,
                                         GeneSet.file=NULL,ssgsea.method=NULL,
                                         no.permutations=no.permutations)
  
  # Prepare the data to perform a pre-ranked GSEA
  A.p1k_NgenesUP <- sapply(1:ncol(A.p1k),function(j) {
    rownames(A.p1k)[order(A.p1k[,j],decreasing = TRUE)][1:Ngenes]},
    simplify = FALSE)
  names(A.p1k_NgenesUP) <- paste("A.Npermut",1:no.permutations,"UP",sep="_")
  A.NgenesUP <- c(list("A.actual_UP"=names(sort(A.eBay$t[,1],decreasing = TRUE))[1:Ngenes]),A.p1k_NgenesUP) 
  
  A.p1k_NgenesDN <- sapply(1:ncol(A.p1k),function(j) {
    rownames(A.p1k)[order(A.p1k[,j],decreasing = FALSE)][1:Ngenes]},
    simplify = FALSE)
  names(A.p1k_NgenesUP) <- paste("A.Npermut",1:no.permutations,"DN",sep="_")
  A.NgenesDN <- c(list("A.actual_DN"=names(sort(A.eBay$t[,1],decreasing = FALSE))[1:Ngenes]),A.p1k_NgenesDN) 
  
  
  
  
  B.p1k_NgenesUP <- sapply(1:ncol(B.p1k),function(j) {
    rownames(B.p1k)[order(B.p1k[,j],decreasing = TRUE)][1:Ngenes]},
    simplify = FALSE)
  names(B.p1k_NgenesUP) <- paste("B.Npermut",1:no.permutations,"UP",sep="_")
  B.NgenesUP <- c(list("B.actual_UP"=names(sort(B.eBay$t[,1],decreasing = TRUE))[1:Ngenes]),B.p1k_NgenesUP) 
  
  
  B.p1k_NgenesDN <- sapply(1:ncol(B.p1k),function(j) {
    rownames(B.p1k)[order(B.p1k[,j],decreasing = FALSE)][1:Ngenes]},
    simplify = FALSE)
  names(B.p1k_NgenesDN) <- paste("B.Npermut",1:no.permutations,"DN",sep="_")
  B.NgenesDN <- c(list("B.actual_DN"=names(sort(B.eBay$t[,1],decreasing = FALSE))[1:Ngenes]),B.p1k_NgenesDN) 
  
  
  AonB.NES_UP <- fgsea(pathways = B.NgenesUP,stats = A.eBay$t[,1],nperm = 1000)$NES
  AonB.NES_DN <- fgsea(pathways = B.NgenesDN,stats = A.eBay$t[,1],nperm = 1000)$NES
  AonB.cNES <- setNames((AonB.NES_UP - AonB.NES_DN),c("Actual",paste("Npermut",1:no.permutations,sep = "_")))
  
  BonA.NES_UP <- fgsea(pathways = A.NgenesUP,stats = B.eBay$t[,1],nperm = 1000)$NES
  BonA.NES_DN <- fgsea(pathways = A.NgenesDN,stats = B.eBay$t[,1],nperm = 1000)$NES
  BonA.cNES <- setNames( (BonA.NES_UP - BonA.NES_DN),c("Actual",paste("Npermut",1:no.permutations,sep = "_")))
  
  TES <- (AonB.cNES + BonA.cNES)/2
  
  if(show.plot) {
    hist(TES[-1],xlim=c(min(TES),max(TES)),
         main=paste0(A$pert_iname,A$pert_id," (",A$trt.pert_type,")"," ~ ",
                     B$pert_iname,B$pert_id," (",B$trt.pert_type,")"),
         xlab="Null Distribution of Total Enrichment Score (TES)",las=1)
    abline(v=TES[1],col="red")
  }
  
  pval = ifelse(TES[1] <0, (sum(TES[-1] < TES[1]) +1)/length(TES[-1]), (sum(TES[-1] > TES[1]) +1)/length(TES[-1]))
  names(pval) <- "pval"
  
  
  return(c("TES"=TES[1],pval))
}


# test_TES(A = list("pert_iname"="palbociclib",
#                   "pert_id"=NULL,
#                   "trt.pert_type"="trt_cp",
#                   "cntl.pert_type"="ctl_vehicle",
#                   "cell_ids"=NULL),
#          B=list("pert_iname"="CDK4",
#                 "pert_id"=NULL,
#                 "trt.pert_type"="trt_sh",
#                 "cntl.pert_type"="ctl_vector",
#                 "cell_ids"=NULL),
#          Ngenes=250,
#          # rawsignatures.type="LIMMA_LVL3",
#          probe2symbol = TRUE,
#          user_key="9ca8871996ddd3f41d4062a446504d7b",
#          no.permutations=100,
#          show.plot=TRUE)