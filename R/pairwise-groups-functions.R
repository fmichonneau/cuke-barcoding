source("R/load.R")

pairwiseGrps <- function(d, threshold) {
    iGrp <- lapply(1:ncol(d), function(j) dimnames(d)[[1]][d[, j] < threshold])
    sngl <- sapply(iGrp, function(x) length(x) == 1)
    edg <- do.call("rbind", lapply(iGrp, function(x) {
        if (length(x) > 1) cbind(head(x, -1), tail(x, -1)) else NULL
    }))
    g <- graph.data.frame(edg)
    c(split(V(g)$name, clusters(g)$membership), iGrp[sngl])
}

getPairwiseGrp <- function(alg=cukeAlg, model=c("raw", "K80"), db=cukeDB,
                           taxonomySubset, threshold,
                           taxonomyData=taxonomyDf) {
    model <- match.arg(model)
    if (missing(taxonomySubset) || taxonomySubset == "all")
        distAlg <- ape::dist.dna(alg, model=model, as.matrix=TRUE)
    else {
        taxLvl <- taxonomyData[taxonomyData$taxa == taxonomySubset, "rank"]
        if (length(taxLvl) > 1)
            stop("Something is wrong with this taxonomic name.")
        if (taxLvl == "Order") {
            distLbl <- subset(db, order == taxonomySubset)$Labels_withAmb

        } else if (taxLvl == "Family") {
            distLbl <- subset(db, family == taxonomySubset)$Labels_withAmb

        }
        else stop("Houston, we have a problem.")
        distLbl <- gsub("\\\"", "", distLbl)
        indexLbl <- match(distLbl, dimnames(alg)[[1]])
        stopifnot(! any(is.na(indexLbl)))
        stopifnot(length(indexLbl) > 3)
        distAlg <- dist.dna(alg[indexLbl, ], model=model, as.matrix=TRUE)
    }
    lapply(threshold, function(thres) pairwiseGrps(distAlg, thres))
}

getPropSngl <- function(x) {
    tt <- sapply(x, length)
    sum(tt == 1)/length(tt)
}

build_pairwiseGrpRes <- function(pairwiseGrpMdl, pairwiseGrpTax) {
    pairwiseGrpMdl <- c("raw", "K80")
    taxonomyDf <- load_taxonomyDf()
    pairwiseGrpTax <- taxonomyDf$taxa
    pairwiseGrpRes <- vector("list", length(pairwiseGrpMdl) *
                             length(pairwiseGrpTax))
    nmGrpRes <- character(length(pairwiseGrpMdl) * length(pairwiseGrpTax))
    i <- 1
    cukeAlg <- load_cukeAlg()
    cukeDB <- load_cukeDB()
    thresVec <- load_thresholdPairwise()
    
    for (eachMdl in 1:length(pairwiseGrpMdl)) {
        for (eachTax in 1:length(pairwiseGrpTax)) {
            pairwiseGrpRes[[i]] <- getPairwiseGrp(alg=cukeAlg,
                                                  model=pairwiseGrpMdl[eachMdl],
                                                  db=cukeDB, threshold=thresVec,
                                                  taxonomySubset=pairwiseGrpTax[eachTax])
            nmGrpRes[i] <- paste(pairwiseGrpMdl[eachMdl],
                                 pairwiseGrpTax[eachTax], sep="-")
            i <- i + 1
        }
    }
    names(pairwiseGrpRes) <- nmGrpRes
    pairwiseGrpRes
}

load_pairwiseGrpRes <- function(overwrite=FALSE) {
    fnm <- "data/pairwiseGrpRes.rds"
    if (file.exists(fnm) && !overwrite) {
        pairwiseGrpRes <- readRDS(file=fnm)
    }
    else {
        pairwiseGrpRes <- build_pairwiseGrpRes()
    }
    pairwiseGrpRes
}
