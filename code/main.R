genLabel <- function(dbTmp) {
  seqNm <- paste(dbTmp$genusorhigher, dbTmp$modifier, dbTmp$species, dbTmp$Loc,
                 paste(dbTmp$"Collection.Code", dbTmp$Catalog_number, sep=""), sep="_")
  seqNm <- gsub("_{2,}", "_", seqNm)
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
  if (file.exists(out)) stop("file ", out, "already exists")
  for (i in 1:nrow(db)) {
    dbTmp <- db[i, ]
    seqNm <- genLabel(dbTmp)
    seqNm <- paste(">", seqNm, sep="")
    cat(seqNm, "\n", dbTmp$Sequence, "\n", file=out, append=TRUE) 
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
library(doMC)
registerDoMC()
library(ape)
library(phylobase)
source("write.dna.R")
source("findGroups.R")
assignInNamespace("write.dna", write.dna, "ape")
source("dist.topo.R")
assignInNamespace("boot.phylo", boot.phylo, "ape")
set.seed(10101)

db <- read.csv(file="MARBoL_Echinos_VI_13.csv", header=T, stringsAsFactors=F)

## only pass == "yes" and Seq_length > 550
db <- db[db$pass == "yes" & db$Seq_length > 550, ]
dbH <- subset(db, class_ == "Holothuroidea")
dbA <- subset(db, class_ == "Asteroidea")
dbE <- subset(db, class_ == "Echinoidea")
dbO <- subset(db, class_ == "Ophiuroidea")
dbC <- subset(db, class_ == "Crinoidea")

#### ----- Identify problematic sequences

genFasta(dbHol, out="/tmp/allEchino.fas")
system("mafft --auto --op 10 --thread -1 /tmp/allEchino.fas > /tmp/allEchino.afa")
allEch <- read.dna(file="/tmp/allEchino.afa", format="fasta")

### Remove sequences with stop codons
### Remove sequences with more than 1 ambiguous (incl. N) base pair
## stop codons
tranE <- foreach (i = 1:nrow(allEch)) %dopar% {
  translate(as.character(allEch[i, ]), frame=1, numcode=9)
}
seqWithStop <- dimnames(allEch)[[1]][grep("\\*", tranE)]

## ambiguous nucleotides
funnyLetters <- numeric(nrow(allEch))
for (i in 1:nrow(allEch)) {
  tmpTbl <- table(as.character(allEch[i, ]))
  if (length(tmpTbl) == 4 && all(names(tmpTbl) %in% c("a", "c", "g", "t")))
    funnyLetters[i] <- 0
  else {
    tmpRes <- try(sum(tmpTbl[-match(c("-", "a", "c", "g", "t"), names(tmpTbl))]))
    funnyLetters[i] <- ifelse(inherits(tmpRes, "try-error"), browser(), tmpRes)
  }
}
seqWithAmb <- dimnames(allEch)[[1]][funnyLetters > 1]
toRm <- union(seqWithStop, seqWithAmb)

toRmInd <- match(toRm, dimnames(allEch)[[1]])
allEch <- allEch[-toRmInd, ]

write.dna(allEch, file="~/Documents/echinoBarcode/20120621.allEchino.fas", format="fasta", colsep="")



#########
######### ---- Can start from here
#########
allEch <- read.dna(file="20120621.allEchino.fas", format="fasta")

## All Echinos
## seqAllEch <- allEch
## dimnames(seqAllEch)[[1]] <- make.unique(dimnames(seqAllEch)[[1]])
## treEch <- nj(dist.dna(seqAllEch))
## treEch$edge.length[treEch$edge.length < 0] <- 0
## bootEch <- boot.phylo(treEch, seqAllEch, function(xx) nj(dist.dna(xx)), B=100)
## treEch$node.label <- bootEch
## treEchr <- root(treEch, grep("Balanometra", treEch$tip.label), resolve.root=T)

## treEch4 <- phylo4d(treEchr, check.node.labels="asdata")

## grpEch1 <- findGroups(treEch4, threshold=.005)
## grpEch2 <- findGroups(treEch4, threshold=.010)
## grpEch3 <- findGroups(treEch4, threshold=.015)
## grpEch4 <- findGroups(treEch4, threshold=.020)
## grpEch5 <- findGroups(treEch4, threshold=.025)
## grpEch6 <- findGroups(treEch4, threshold=.030)

## Crinoids
lblC <- sapply(1:nrow(dbC), function(i) genLabel(dbC[i, ]))
keepC <- dimnames(allEch)[[1]] %in% lblC
seqC <- allEch[keepC, ]
dimnames(seqC)[[1]] <- make.unique(dimnames(seqC)[[1]])
treC <- nj(dist.dna(seqC))
treC$edge.length[treC$edge.length < 0] <- 0
bootC <- boot.phylo(treC, seqC, function(xx) nj(dist.dna(xx)), B=100)
treC$node.label <- bootC
treCr <- root(treC, grep("Balanometra", treC$tip.label), resolve.root=T)

treC4 <- phylo4d(treCr, check.node.label="asdata")
grpC3 <- findGroups(treC4, threshold=0.015)

nrow(seqC)
max(tipData(grpC3)$Group)
c(getUniqueSeq(seqC))

## Ophiuroids
lblO <- sapply(1:nrow(dbO), function(i) genLabel(dbO[i, ]))
keepO <- dimnames(allEch)[[1]] %in% lblO
seqO <- allEch[keepO, ]
dimnames(seqO)[[1]] <- make.unique(dimnames(seqO)[[1]])
treO <- nj(dist.dna(seqO))
treO$edge.length[treO$edge.length < 0] <- 0
bootO <- boot.phylo(treO, seqO, function(xx) nj(dist.dna(xx)), B=100)
treO$node.label <- bootO
treOr <- root(treO, grep("Gorgonocephalidae", treO$tip.label), resolve.root=T)

treO4 <- phylo4d(treOr, check.node.label="asdata")
grpO3 <- findGroups(treO4, threshold=0.015)

nrow(seqO)
max(tipData(grpO3)$Group)
c(getUniqueSeq(seqO))
(table(table(tipData(grpO3)$Group)))

save.image(file="20120703.CriOphBootGroup.RData")

## Holothuroids
lblH <- sapply(1:nrow(dbH), function(i) genLabel(dbH[i, ]))
keepH <- dimnames(allEch)[[1]] %in% lblH
seqH <- allEch[keepH, ]
dimnames(seqH)[[1]] <- make.unique(dimnames(seqH)[[1]])
treH <- nj(dist.dna(seqH))
treH$edge.length[treH$edge.length < 0] <- 0
bootH <- boot.phylo(treH, seqH, function(xx) nj(dist.dna(xx)), B=100)
treH$node.label <- bootH
treHr <- root(treH, grep("^Pseudostichopus", treH$tip.label), resolve.root=TRUE)

treH4 <- phylo4d(treHr, check.node.label="asdata")

grpH1 <- findGroups(treH4, threshold=.005)
grpH2 <- findGroups(treH4, threshold=.010)
grpH3 <- findGroups(treH4, threshold=.015)
grpH4 <- findGroups(treH4, threshold=.020)
grpH5 <- findGroups(treH4, threshold=.025)
grpH6 <- findGroups(treH4, threshold=.030)


## Dendro
dbDen <- subset(dbH, order_ == "Dendrochirotida")
lblDen <- sapply(1:nrow(dbDen), function(i) genLabel(dbDen[i, ]))
seqDen <- allEch[dimnames(allEch)[[1]] %in% lblDen, ]
dimnames(seqDen)[[1]] <- make.unique(dimnames(seqDen)[[1]])
treDen <- nj(dist.dna(seqDen))
treDen$edge.length[treDen$edge.length < 0] <- 0
bootDen <- boot.phylo(treDen, seqDen, function(xx) nj(dist.dna(xx)), B=100)
treDen$node.label <- bootDen
treDenr <- root(treDen, grep("^Lissothuria", treDen$tip.label), resolve.root=TRUE)
treDen4 <- phylo4d(treDenr, check.node.label="asdata")

grpDen1 <- findGroups(treDen4, threshold=.005)
grpDen2 <- findGroups(treDen4, threshold=.010)
grpDen3 <- findGroups(treDen4, threshold=.015)
grpDen4 <- findGroups(treDen4, threshold=.020)
grpDen5 <- findGroups(treDen4, threshold=.025)
grpDen6 <- findGroups(treDen4, threshold=.030)

## Apodid
dbApo <- subset(dbH, order_ == "Apodida")
lblApo <- sapply(1:nrow(dbApo), function(i) genLabel(dbApo[i, ]))
seqApo <- allEch[dimnames(allEch)[[1]] %in% lblApo, ]
dimnames(seqApo)[[1]] <- make.unique(dimnames(seqApo)[[1]])
treApo <- nj(dist.dna(seqApo))
treApo$edge.length[treApo$edge.length < 0] <- 0
bootApo <- boot.phylo(treApo, seqApo, function(xx) nj(dist.dna(xx)), B=100)
treApo$node.label <- bootApo
treApor <- root(treApo, grep("^Polycheira", treApo$tip.label), resolve.root=TRUE)
treApo4 <- phylo4d(treApor, check.node.label="asdata")

grpApo1 <- findGroups(treApo4, threshold=.005)
grpApo2 <- findGroups(treApo4, threshold=.010)
grpApo3 <- findGroups(treApo4, threshold=.015)
grpApo4 <- findGroups(treApo4, threshold=.020)
grpApo5 <- findGroups(treApo4, threshold=.025)
grpApo6 <- findGroups(treApo4, threshold=.030)

## Aspido
dbAsp <- subset(dbH, order_ == "Aspidochirotida")
lblAsp <- sapply(1:nrow(dbAsp), function(i) genLabel(dbAsp[i, ]))
seqAsp <- allEch[dimnames(allEch)[[1]] %in% lblAsp, ]
dimnames(seqAsp)[[1]] <- make.unique(dimnames(seqAsp)[[1]])
treAsp <- nj(dist.dna(seqAsp))
treAsp$edge.length[treAsp$edge.length < 0] <- 0
bootAsp <- boot.phylo(treAsp, seqAsp, function(xx) nj(dist.dna(xx)), B=100)
treAsp$node.label <- bootAsp
treAspr <- root(treAsp, grep("^Pseudostichopus", treAsp$tip.label), resolve.root=TRUE)
treAsp4 <- phylo4d(treAspr, check.node.label="asdata")

grpAsp1 <- findGroups(treAsp4, threshold=.005)
grpAsp2 <- findGroups(treAsp4, threshold=.010)
grpAsp3 <- findGroups(treAsp4, threshold=.015)
grpAsp4 <- findGroups(treAsp4, threshold=.020)
grpAsp5 <- findGroups(treAsp4, threshold=.025)
grpAsp6 <- findGroups(treAsp4, threshold=.030)

## Holoth
dbHol <- subset(dbH, family == "Holothuriidae")
lblHol <- sapply(1:nrow(dbHol), function(i) genLabel(dbHol[i, ]))
seqHol <- allEch[dimnames(allEch)[[1]] %in% lblHol, ]
dimnames(seqHol)[[1]] <- make.unique(dimnames(seqHol)[[1]])
treHol <- nj(dist.dna(seqHol))
treHol$edge.length[treHol$edge.length < 0] <- 0
bootHol <- boot.phylo(treHol, seqHol, function(xx) nj(dist.dna(xx)), B=100)
treHol$node.label <- bootHol
treHolr <- root(treHol, grep("platei", treHol$tip.label), resolve.root=TRUE)
treHol4 <- phylo4d(treHolr, check.node.label="asdata")

grpHol1 <- findGroups(treHol4, threshold=.005)
grpHol2 <- findGroups(treHol4, threshold=.010)
grpHol3 <- findGroups(treHol4, threshold=.015)
grpHol4 <- findGroups(treHol4, threshold=.020)
grpHol5 <- findGroups(treHol4, threshold=.025)
grpHol6 <- findGroups(treHol4, threshold=.030)


nrow(seqH)
c(getUniqueSeq(seqH))
max(tipData(grpH3))
table(table(tipData(grpH3)$Group))

## holnm <- sapply(1:nrow(dbH), function(i) genSp(dbH[i, ]))
## holnm <- cbind(lblH, holnm)
## holnm <- holnm[lblH %in% dimnames(seqH)[[1]], ]
## write.csv(unique(holnm[,2]), file="~/Documents/echinoBarcode/uniqueHolNm.csv")


## Echinoids
lblE <- sapply(1:nrow(dbE), function(i) genLabel(dbE[i, ]))
keepE <- dimnames(allEch)[[1]] %in% lblE
seqE <- allEch[keepE, ]
dimnames(seqE)[[1]] <- make.unique(dimnames(seqE)[[1]])
treE <- nj(dist.dna(seqE))
treE$edge.length[treE$edge.length < 0] <- 0
bootE <- boot.phylo(treE, seqE, function(xx) nj(dist.dna(xx)), B=100)
treE$node.label <- bootE
treEr <- root(treE, grep("Echinothrix|m_Echinoidea_Eparses", treE$tip.label), resolve.root=T)

treE4 <- phylo4d(treEr, check.node.label="asdata")

grpE1 <- findGroups(treE4, threshold=.005)
grpE2 <- findGroups(treE4, threshold=.010)
grpE3 <- findGroups(treE4, threshold=.015)
grpE4 <- findGroups(treE4, threshold=.020)
grpE5 <- findGroups(treE4, threshold=.025)
grpE6 <- findGroups(treE4, threshold=.030)

nrow(seqE)
max(tipData(grpE3))
c(getUniqueSeq(seqE))
table(table(tipData(grpE3)))

## echnm <- sapply(1:nrow(dbE), function(i) genSp(dbE[i, ]))
## echnm <- cbind(lblE, echnm)
## echnm <- echnm[lblE %in% dimnames(seqE)[[1]], ]
## write.csv(unique(echnm[,2]), file="~/Documents/echinoBarcode/uniqueEchNm.csv")


## Asteroids
lblA <- sapply(1:nrow(dbA), function(i) genLabel(dbA[i, ]))
keepA <- dimnames(allEch)[[1]] %in% lblA
seqA <- allEch[keepA, ]
dimnames(seqA)[[1]] <- make.unique(dimnames(seqA)[[1]])
treA <- nj(dist.dna(seqA))
treA$edge.length[treA$edge.length < 0] <- 0
bootA <- boot.phylo(treA, seqA, function(xx) nj(dist.dna(xx)), B=100)
treA$node.label <- bootA
treAr <- root(treA, grep("UF226$|UF3915$", treA$tip.label), resolve.root=T)

treA4 <- phylo4d(treAr, check.node.labels="asdata")

grpA1 <- findGroups(treA4, threshold=.005)
grpA2 <- findGroups(treA4, threshold=.010)
grpA3 <- findGroups(treA4, threshold=.015)
grpA4 <- findGroups(treA4, threshold=.020)
grpA5 <- findGroups(treA4, threshold=.025)
grpA6 <- findGroups(treA4, threshold=.030)

nrow(seqA)
c(getUniqueSeq(seqA))
max(tipData(grpA3))
table(table(tipData(grpA3)))

## astnm <- sapply(1:nrow(dbA), function(i) genSp(dbA[i, ]))
## astnm <- cbind(lblA, astnm)
## astnm <- astnm[lblA %in% dimnames(seqA)[[1]], ]
## astnm <- gsub("^m_", "", astnm)
## astnm <- gsub("^X_", "", astnm)
## astnm <- gsub("_$", "", astnm)
## write.csv(unique(astnm[,2]), file="~/Documents/echinoBarcode/uniqueAstNm.csv")

distA1 <- getIntraInterDist(grpA1, dist.dna(seqA, as.matrix=T), check.boot=80)
distA2 <- getIntraInterDist(grpA2, dist.dna(seqA, as.matrix=T), check.boot=80)
distA3 <- getIntraInterDist(grpA3, dist.dna(seqA, as.matrix=T), check.boot=80)
distA4 <- getIntraInterDist(grpA4, dist.dna(seqA, as.matrix=T), check.boot=80)
distA5 <- getIntraInterDist(grpA5, dist.dna(seqA, as.matrix=T), check.boot=80)
distA6 <- getIntraInterDist(grpA6, dist.dna(seqA, as.matrix=T), check.boot=80)

distE1 <- getIntraInterDist(grpE1, dist.dna(seqE, as.matrix=T), check.boot=80)
distE2 <- getIntraInterDist(grpE2, dist.dna(seqE, as.matrix=T), check.boot=80)
distE3 <- getIntraInterDist(grpE3, dist.dna(seqE, as.matrix=T), check.boot=80)
distE4 <- getIntraInterDist(grpE4, dist.dna(seqE, as.matrix=T), check.boot=80)
distE5 <- getIntraInterDist(grpE5, dist.dna(seqE, as.matrix=T), check.boot=80)
distE6 <- getIntraInterDist(grpE6, dist.dna(seqE, as.matrix=T), check.boot=80)

distH1 <- getIntraInterDist(grpH1, dist.dna(seqH, as.matrix=T), check.boot=80)
distH2 <- getIntraInterDist(grpH2, dist.dna(seqH, as.matrix=T), check.boot=80)
distH3 <- getIntraInterDist(grpH3, dist.dna(seqH, as.matrix=T), check.boot=80)
distH4 <- getIntraInterDist(grpH4, dist.dna(seqH, as.matrix=T), check.boot=80)
distH5 <- getIntraInterDist(grpH5, dist.dna(seqH, as.matrix=T), check.boot=80)
distH6 <- getIntraInterDist(grpH6, dist.dna(seqH, as.matrix=T), check.boot=80)

distDen1 <- getIntraInterDist(grpDen1, dist.dna(seqDen, as.matrix=T), check.boot=80)
distDen2 <- getIntraInterDist(grpDen2, dist.dna(seqDen, as.matrix=T), check.boot=80)
distDen3 <- getIntraInterDist(grpDen3, dist.dna(seqDen, as.matrix=T), check.boot=80)
distDen4 <- getIntraInterDist(grpDen4, dist.dna(seqDen, as.matrix=T), check.boot=80)
distDen5 <- getIntraInterDist(grpDen5, dist.dna(seqDen, as.matrix=T), check.boot=80)
distDen6 <- getIntraInterDist(grpDen6, dist.dna(seqDen, as.matrix=T), check.boot=80)

distApo1 <- getIntraInterDist(grpApo1, dist.dna(seqApo, as.matrix=T), check.boot=80)
distApo2 <- getIntraInterDist(grpApo2, dist.dna(seqApo, as.matrix=T), check.boot=80)
distApo3 <- getIntraInterDist(grpApo3, dist.dna(seqApo, as.matrix=T), check.boot=80)
distApo4 <- getIntraInterDist(grpApo4, dist.dna(seqApo, as.matrix=T), check.boot=80)
distApo5 <- getIntraInterDist(grpApo5, dist.dna(seqApo, as.matrix=T), check.boot=80)
distApo6 <- getIntraInterDist(grpApo6, dist.dna(seqApo, as.matrix=T), check.boot=80)

distAsp1 <- getIntraInterDist(grpAsp1, dist.dna(seqAsp, as.matrix=T), check.boot=80)
distAsp2 <- getIntraInterDist(grpAsp2, dist.dna(seqAsp, as.matrix=T), check.boot=80)
distAsp3 <- getIntraInterDist(grpAsp3, dist.dna(seqAsp, as.matrix=T), check.boot=80)
distAsp4 <- getIntraInterDist(grpAsp4, dist.dna(seqAsp, as.matrix=T), check.boot=80)
distAsp5 <- getIntraInterDist(grpAsp5, dist.dna(seqAsp, as.matrix=T), check.boot=80)
distAsp6 <- getIntraInterDist(grpAsp6, dist.dna(seqAsp, as.matrix=T), check.boot=80)

distHol1 <- getIntraInterDist(grpHol1, dist.dna(seqHol, as.matrix=T), check.boot=80)
distHol2 <- getIntraInterDist(grpHol2, dist.dna(seqHol, as.matrix=T), check.boot=80)
distHol3 <- getIntraInterDist(grpHol3, dist.dna(seqHol, as.matrix=T), check.boot=80)
distHol4 <- getIntraInterDist(grpHol4, dist.dna(seqHol, as.matrix=T), check.boot=80)
distHol5 <- getIntraInterDist(grpHol5, dist.dna(seqHol, as.matrix=T), check.boot=80)
distHol6 <- getIntraInterDist(grpHol6, dist.dna(seqHol, as.matrix=T), check.boot=80)


save.image(file="20120630.echinoBarcode.RData")

### number of ESUs vs threshold
nbEsu <- data.frame(Threshold=rep(1:6, 7),
                    Class=c(rep("Asteroidea", 6),
                      rep("Echinoidea", 6),
                      rep("Holothuroidea", 6),
                      rep("Dendrochirotida", 6),
                      rep("Apodida", 6),
                      rep("Aspidochirotida", 6),
                      rep("Holothuriidae", 6)),                    
                    NbESU=c(
                      max(tipData(grpA1)$Groups),
                      max(tipData(grpA2)$Groups),
                      max(tipData(grpA3)$Groups),
                      max(tipData(grpA4)$Groups),
                      max(tipData(grpA5)$Groups),
                      max(tipData(grpA6)$Groups),
                      max(tipData(grpE1)$Groups),
                      max(tipData(grpE2)$Groups),
                      max(tipData(grpE3)$Groups),
                      max(tipData(grpE4)$Groups),
                      max(tipData(grpE5)$Groups),
                      max(tipData(grpE6)$Groups),
                      max(tipData(grpH1)$Groups),
                      max(tipData(grpH2)$Groups),
                      max(tipData(grpH3)$Groups),
                      max(tipData(grpH4)$Groups),
                      max(tipData(grpH5)$Groups),
                      max(tipData(grpH6)$Groups),
                      max(tipData(grpDen1)$Groups),
                      max(tipData(grpDen2)$Groups),
                      max(tipData(grpDen3)$Groups),
                      max(tipData(grpDen4)$Groups),
                      max(tipData(grpDen5)$Groups),
                      max(tipData(grpDen6)$Groups),
                      max(tipData(grpApo1)$Groups),
                      max(tipData(grpApo2)$Groups),
                      max(tipData(grpApo3)$Groups),
                      max(tipData(grpApo4)$Groups),
                      max(tipData(grpApo5)$Groups),
                      max(tipData(grpApo6)$Groups),
                      max(tipData(grpAsp1)$Groups),
                      max(tipData(grpAsp2)$Groups),
                      max(tipData(grpAsp3)$Groups),
                      max(tipData(grpAsp4)$Groups),
                      max(tipData(grpAsp5)$Groups),
                      max(tipData(grpAsp6)$Groups),
                      max(tipData(grpHol1)$Groups),
                      max(tipData(grpHol2)$Groups),
                      max(tipData(grpHol3)$Groups),
                      max(tipData(grpHol4)$Groups),
                      max(tipData(grpHol5)$Groups),
                      max(tipData(grpHol6)$Groups)
                      ))

pdf(file="nbESU.pdf", width=8, height=3)
ggplot(subset(nbEsu, Class %in% c("Echinoidea", "Asteroidea", "Holothuroidea")),
              aes(x=Threshold, y=NbESU, group=Class, colour=Class)) + geom_line()
dev.off()

pdf(file="nbESU_hol.pdf", width=8, height=3)
ggplot(subset(nbEsu, Class %in% c("Aspidochirotida", "Dendrochirotida", "Apodida", "Holothuriidae")),
              aes(x=Threshold, y=NbESU, group=Class, colour=Class)) + geom_line()
dev.off()


## ### visualize sequence differences
## varSA <- apply(seqA, 2, function(x) max(base.freq(x)))
## varSH <- apply(seqH, 2, function(x) max(base.freq(x)))
## varSE <- apply(seqE, 2, function(x) max(base.freq(x)))

## par(mfrow=c(3,1))
## plot(1-varSA, col="blue", type="h", pch=15, cex=.5)
## plot(1-varSH, col="red", type="h", pch=15, cex=.5)
## plot(1-varSE, col="green", type="h", pch=15, cex=.5)

## matVar <- array(, dim=c(3,max(c(length(varSA), length(varSE), length(varSH)))))
## matVar[1,] <- 1-varSA
## matVar[2,] <- 1-varSE
## matVar[3,] <- 1-varSH
## barplot(matVar[,1:100], col=c("blue", "red", "green"), beside=T)

allDist <- rbind(cbind(Order=rep("Asteroidea", nrow(distA1)), Threshold=rep(1, nrow(distA1)), distA1),
        cbind(Order=rep("Asteroidea", nrow(distA2)), Threshold=rep(2, nrow(distA2)), distA2),
        cbind(Order=rep("Asteroidea", nrow(distA3)), Threshold=rep(3, nrow(distA3)), distA3),
        cbind(Order=rep("Asteroidea", nrow(distA4)), Threshold=rep(4, nrow(distA4)), distA4),
        cbind(Order=rep("Asteroidea", nrow(distA5)), Threshold=rep(5, nrow(distA5)), distA5),
        cbind(Order=rep("Asteroidea", nrow(distA6)), Threshold=rep(6, nrow(distA6)), distA6),        
        cbind(Order=rep("Echinoidea", nrow(distE1)), Threshold=rep(1, nrow(distE1)), distE1),
        cbind(Order=rep("Echinoidea", nrow(distE2)), Threshold=rep(2, nrow(distE2)), distE2),
        cbind(Order=rep("Echinoidea", nrow(distE3)), Threshold=rep(3, nrow(distE3)), distE3),
        cbind(Order=rep("Echinoidea", nrow(distE4)), Threshold=rep(4, nrow(distE4)), distE4),
        cbind(Order=rep("Echinoidea", nrow(distE5)), Threshold=rep(5, nrow(distE5)), distE5),
        cbind(Order=rep("Echinoidea", nrow(distE6)), Threshold=rep(6, nrow(distE6)), distE6),
        cbind(Order=rep("Holothuroidea", nrow(distH1)), Threshold=rep(1, nrow(distH1)), distH1),
        cbind(Order=rep("Holothuroidea", nrow(distH2)), Threshold=rep(2, nrow(distH2)), distH2),
        cbind(Order=rep("Holothuroidea", nrow(distH3)), Threshold=rep(3, nrow(distH3)), distH3),
        cbind(Order=rep("Holothuroidea", nrow(distH4)), Threshold=rep(4, nrow(distH4)), distH4),
        cbind(Order=rep("Holothuroidea", nrow(distH5)), Threshold=rep(5, nrow(distH5)), distH5),
        cbind(Order=rep("Holothuroidea", nrow(distH6)), Threshold=rep(6, nrow(distH6)), distH6),
        cbind(Order=rep("Dendrochirotida", nrow(distDen1)), Threshold=rep(1, nrow(distDen1)), distDen1),
        cbind(Order=rep("Dendrochirotida", nrow(distDen2)), Threshold=rep(2, nrow(distDen2)), distDen2),
        cbind(Order=rep("Dendrochirotida", nrow(distDen3)), Threshold=rep(3, nrow(distDen3)), distDen3),
        cbind(Order=rep("Dendrochirotida", nrow(distDen4)), Threshold=rep(4, nrow(distDen4)), distDen4),
        cbind(Order=rep("Dendrochirotida", nrow(distDen5)), Threshold=rep(5, nrow(distDen5)), distDen5),
        cbind(Order=rep("Dendrochirotida", nrow(distDen6)), Threshold=rep(6, nrow(distDen6)), distDen6),
        cbind(Order=rep("Apodida", nrow(distApo1)), Threshold=rep(1, nrow(distApo1)), distApo1),
        cbind(Order=rep("Apodida", nrow(distApo2)), Threshold=rep(2, nrow(distApo2)), distApo2),
        cbind(Order=rep("Apodida", nrow(distApo3)), Threshold=rep(3, nrow(distApo3)), distApo3),
        cbind(Order=rep("Apodida", nrow(distApo4)), Threshold=rep(4, nrow(distApo4)), distApo4),
        cbind(Order=rep("Apodida", nrow(distApo5)), Threshold=rep(5, nrow(distApo5)), distApo5),
        cbind(Order=rep("Apodida", nrow(distApo6)), Threshold=rep(6, nrow(distApo6)), distApo6),
        cbind(Order=rep("Aspidochirotida", nrow(distAsp1)), Threshold=rep(1, nrow(distAsp1)), distAsp1),
        cbind(Order=rep("Aspidochirotida", nrow(distAsp2)), Threshold=rep(2, nrow(distAsp2)), distAsp2),
        cbind(Order=rep("Aspidochirotida", nrow(distAsp3)), Threshold=rep(3, nrow(distAsp3)), distAsp3),
        cbind(Order=rep("Aspidochirotida", nrow(distAsp4)), Threshold=rep(4, nrow(distAsp4)), distAsp4),
        cbind(Order=rep("Aspidochirotida", nrow(distAsp5)), Threshold=rep(5, nrow(distAsp5)), distAsp5),
        cbind(Order=rep("Aspidochirotida", nrow(distAsp6)), Threshold=rep(6, nrow(distAsp6)), distAsp6),
        cbind(Order=rep("Holothuriidae", nrow(distHol1)), Threshold=rep(1, nrow(distHol1)), distHol1),
        cbind(Order=rep("Holothuriidae", nrow(distHol2)), Threshold=rep(2, nrow(distHol2)), distHol2),
        cbind(Order=rep("Holothuriidae", nrow(distHol3)), Threshold=rep(3, nrow(distHol3)), distHol3),
        cbind(Order=rep("Holothuriidae", nrow(distHol4)), Threshold=rep(4, nrow(distHol4)), distHol4),
        cbind(Order=rep("Holothuriidae", nrow(distHol5)), Threshold=rep(5, nrow(distHol5)), distHol5),
        cbind(Order=rep("Holothuriidae", nrow(distHol6)), Threshold=rep(6, nrow(distHol6)), distHol6))


svg(file="compareIntraInter.svg", width=9, height=5)
ggplot(subset(allDist, Order %in% c("Asteroidea", "Echinoidea", "Holothuroidea") &
              Threshold > 1),
       aes(x=dist, fill=typeDist, alpha=.5)) + geom_density() +
  facet_grid(Order ~ Threshold) + scale_fill_manual(values=c("yellow", "red"))
dev.off()

svg(file="compareIntraInterHol.svg", width=9, height=5)
ggplot(subset(allDist, Order %in% c("Dendrochirotida", "Apodida", "Aspidochirotida", "Holothuriidae") &
              Threshold > 1),
       aes(x=dist, fill=typeDist, alpha=.5)) + geom_density() +
  facet_grid(Order ~ Threshold) + scale_fill_manual(values=c("yellow", "red"))
dev.off()

### Geographic barriers
treE3 <- grpE3
tipLabels(treE3) <- paste(tipData(grpE3)$Group, tipLabels(treE3), sep="_")
bs <- nodeData(treE3)[,1, drop=F]
bsToShow <- rownames(bs)[bs > 80]
pdf(file="echinoids-015.pdf", height=100)
plot(as(treE3, "phylo"))
nodelabels(rep("", length(bsToShow)), as.integer(bsToShow), frame="circ", cex=.5)
dev.off()

treA3 <- grpA3
tipLabels(treA3) <- paste(tipData(grpA3)$Group, tipLabels(treA3), sep="_")
bs <- nodeData(treA3)[,1, drop=F]
bsToShow <- rownames(bs)[bs > 80]
pdf(file="asteroids-015.pdf", height=100)
plot(as(treA3, "phylo"))
nodelabels(rep("", length(bsToShow)), as.integer(bsToShow), frame="circ", cex=.5)
dev.off()

treH3 <- grpH3
tipLabels(treH3) <- paste(tipData(grpH3)$Group, tipLabels(treH3), sep="_")
bs <- nodeData(treH3)[,1, drop=F]
bsToShow <- rownames(bs)[bs > 80]
pdf(file="holothuroids-015.pdf", height=100)
plot(as(treH3, "phylo"), cex=.25)
nodelabels(rep("", length(bsToShow)), as.integer(bsToShow), frame="circ", cex=.125)
dev.off()

geo <- read.csv(file="geographicBarriers.csv")
propGeo <- with(geo, tapply(Allopatric, Class, function(x) {x <- x[!is.na(x)]; sum(x)/length(x)}))

svg(file="rateGeo.svg")
barplot(propGeo, ylim=c(0,1))
dev.off()

#############

stdSeqLength <- function(alg) {

}




## ## if only ufid

## ## yy <- sapply(xx, function(x) {
## ##   xt <- unlist(strsplit(x, "_"))
## ##   xt <- xt[length(xt)]
## ##   if (length(grep("^[0-9]+$", xt)) == 1)
## ##     TRUE
## ##   else FALSE
## ## })


## ## zz <- sapply(xx[yy], function(x) {
## ##   t <- unlist(strsplit(x, "_"))
## ##   t[length(t)]
## ##   })
