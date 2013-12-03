
genLabel <- function(dbTmp) {
  seqNm <- paste(dbTmp$genusorhigher, dbTmp$modifier, dbTmp$species, dbTmp$Loc,
                 paste(dbTmp$"Collection.Code", dbTmp$Catalog_number, sep=""), sep="_")
  seqNm <- gsub("_{2,}", "_", seqNm)
  seqNm <- gsub("_$", "", seqNm)
  seqNm <- gsub("\\s+$", "", seqNm)
  seqNm
}
genSp <- function(dbTmp) {
  seqNm <- paste(dbTmp$genusorhigher, dbTmp$modifier, dbTmp$species, sep="_")
  seqNm <- gsub("_{2,}", "_", seqNm)
  seqNm
}

genFasta <- function(db, out=file.path("/tmp", paste(format(Sys.time(), "%Y%m%d-%H%M%S"), "seq.fas", sep="_"))) {
### db -- database in which the data is stored
### out -- file name of the fasta that will be generated
  ##if (file.exists(out)) stop("file ", out, "already exists")
  for (i in 1:nrow(db)) {
    dbTmp <- db[i, ]
    seqNm <- genLabel(dbTmp)
    seqNm <- paste(">", seqNm, sep="")
    seqNm <- gsub("\\s+", "", seqNm)
    seqTmp <- dbTmp$Sequence
    cat(seqNm, "\n", dbTmp$Sequence, "\n", file=out,
        append=ifelse(i == 1, FALSE, TRUE), sep="")
  }
  TRUE
}

getUniqueSeq <- function(s, ...) {
### returns the number of unique sequences in an alignment (and the names of the
###   duplicated seq as an attribute "whichDup")
### s -- is an DNA alignements
### ... -- further arguments to be passed to dist.dna

  d <- dist.dna(s, ..., as.matrix=TRUE)
  d[upper.tri(d, diag=TRUE)] <- NA
  lbl <- cbind(dimnames(d)[[1]][row(d)], dimnames(d)[[1]][col(d)])
  whichDup <- lbl[which(d == 0), ]
  res <- nrow(s) - length(unique(whichDup[,1]))
  attr(res, "whichDup") <- whichDup
  return(res)
}

getIntraInterDist <- function(tr, d, check.boot=NULL) {
### tr -- phylogenetic tree that gives the groups (output of findGroups)
### d -- distance *matrix* (from dist.dna(..., as.matrix=TRUE))
### check.boot -- should interspecific distances be returned only for nodes
###   above a certain bootstrap value? If NULL all nodes are considered,
###   otherwise a numerical value indicating bootstrap values to be considered
###   e.g., 90: only nodes with bootstrap values >= 90 will be included
  
  if(!inherits(tr, "phylo4"))
    stop("tr must be a phylo4(d) object")
  if(!inherits(d, "matrix"))
    stop("d must be a matrix")
  if(is.null(tipData(tr)$Group))
    stop("The groups must be specified in a column named ", sQuote("Group"),
         " in the ", sQuote("tr"), " object.")
  if(!is.null(check.boot) && !inherits(tr, "phylo4d") &&
     !is.null(tdata(tr)$labelValues))
    stop(sQuote("tr"), "is not properly formatted if you want to use",
         "the bootstrap values. Make sure it's a phylo4d object and that",
         "the node labels are stored as *data* with the name",
         sQuote("labelValues"))

  iDimNm <- match(dimnames(d)[[1]], tipLabels(tr))
  stopifnot(all(!is.na(iDimNm)))
  tipLabels(tr) <- paste(tipData(tr)$Group, tipLabels(tr), sep="_")
  dimnames(d)[[1]][iDimNm] <- tipLabels(tr)
  dimnames(d)[[2]][iDimNm] <- tipLabels(tr)
  
  alrdy <- numeric(nTips(tr))
  resIntra <- resInter <- vector("list", nTips(tr))
  lAlrdy <- 1
  for (i in 1:nTips(tr)) {
    if (tipLabels(tr)[i] %in% alrdy) next
    ancI <- ancestors(tr, i)
    for (j in 1:length(ancI)) {
      descI <- descendants(tr, ancI[j], "tips")
      nmDescI <- names(descI)
      prefix <- gsub("_.+", "", nmDescI)
      if (length(unique(prefix)) == 2) {
        allNd <- names(getNode(tr, nmDescI))
        intraNd1 <- grep(paste("^", unique(prefix)[1], "_", sep=""), allNd, value=T)
        intraNd2 <- grep(paste("^", unique(prefix)[2], "_", sep=""), allNd, value=T)
        alrdy[lAlrdy:(lAlrdy+length(allNd)-1)] <- allNd
        lAlrdy <- lAlrdy+length(allNd)+1        
        if (length(intraNd1) < 2 || length(intraNd2) < 2)
          next
        else if (!is.null(check.boot)) {
          mr <- MRCA(tr, descI)
          bs <- tdata(tr)$labelValues[MRCA(tr, descI)]
          if (bs < check.boot)
            next
          else {
            iToSub <- match(allNd, dimnames(d)[[1]])
            tmpD <- d[iToSub, iToSub]
            tmpD[upper.tri(tmpD, diag=TRUE)] <- NA
            interD <- tmpD[prefix[col(tmpD)] != prefix[row(tmpD)]]
            intraD <- tmpD[prefix[col(tmpD)] == prefix[row(tmpD)]]
            resInter[[i]] <- na.omit(interD)
            names(resInter)[i] <- paste(unique(prefix), collapse="-")
            resIntra[[i]] <- na.omit(intraD)
            names(resIntra)[i] <- paste(unique(prefix), collapse="-")
            break
          }
        }
      }
      else if (length(unique(prefix)) > 2) {
        break
      }
      else next    
    }
  }
  toRmIntra <- sapply(resIntra, is.null)
  toRmInter <- sapply(resInter, is.null)
  tmpIntra <- resIntra[!toRmIntra]
  tmpInter <- resInter[!toRmInter]
  toRepIntra <- lapply(tmpIntra, length)
  toRepInter <- lapply(tmpInter, length)
  nmIntra <- rep(names(tmpIntra), toRepIntra)
  nmInter <- rep(names(tmpInter), toRepInter)
  tmpIntra <- unlist(tmpIntra)
  tmpInter <- unlist(tmpInter)
  data.frame(typeDist=c(rep("intra", length(tmpIntra)),
                       rep("inter", length(tmpInter))),
             namePair=c(nmIntra, nmInter),
             dist=c(tmpIntra, tmpInter))

}


#### ----- End of functions
