source("R/packages.R")
source("R/build.R")

load_echinoDB <- function(overwrite=FALSE) {
    fnm <- "data/raw/cukeBarcodes.csv.rds"
    if (file.exists(fnm) && !overwrite) {
        echinoDB <- readRDS(fnm)
    }
    else {
        echinoDB <- read.csv(file="data/raw/MARBoL_Echinos_VIII_2013.csv",
                             stringsAsFactors=FALSE)
        saveRDS(echinoDB, file=fnm)
    }
    echinoDB
}

load_cukeDB <- function(overwrite=FALSE) {
    source("R/genFasta.R")
    
    fnm <- "data/cukeDB_withLabels.rds"
    if (file.exists(fnm) && !overwrite)
        cukeDB_lbls <- readRDS(file=fnm)
    else {
        treeH <- load_cukeTree_raw()
        echinoDB <- load_echinoDB()
        cukeDB <- subset(echinoDB, class_ == "Holothuroidea")
        dataLbls <- character(nrow(cukeDB))
        for (i in 1:nrow(cukeDB)) {
            dataLbls[i] <- genLabel(cukeDB[i, ])
        }
        
        ## using the same pattern as in seqManagement::cleanSeqLabels
        cukeDB$Labels <- gsub(":|,|\\(|\\)|;|\\[|\\]|\\'|\\s|\t", "", dataLbls)
        
        treeTips <- data.frame(Labels_withAmb = tipLabels(treeH),
                               Labels = gsub("_\\d+amb$", "", tipLabels(treeH)),
                               stringsAsFactors=FALSE)
        
        stopifnot(all(treeTips$treeLabels %in% cukeDB$Labels))
        
        cukeDB_lbls <- merge(cukeDB, treeTips, by="Labels")
        saveRDS(cukeDB_lbls, fnm)
    }
    cukeDB_lbls
}

load_cukeTree_raw <- function(overwrite=FALSE, ...) {
    fnm <- "data/cukeTree-raw.rds"
    if (file.exists(fnm) && !overwrite) {
        cukeTree <- readRDS(file=fnm)
    }
    else {
        cukeAlg <- load_cukeAlg()
        cukeTree <- build_cukeTree(cukeAlg, model="raw", Nrep=200)
        saveRDS(cukeTree, file=fnm)
    }
    cukeTree
}

load_cukeTree_k2p <- function(overwrite=FALSE, ...) {
    fnm <- "data/cukeTree-k2p.rds"
    if (file.exists(fnm) && !overwrite) {
        cukeTree <- readRDS(file=fnm)
    }
    else {
        cukeAlg <- load_cukeAlg()
        cukeTree <- build_cukeTree(cukeAlg, model="K80", Nrep=200)
        saveRDS(cukeTree, file=fnm)
    }
    cukeTree
}

load_cukeAlg <- function(overwrite=FALSE) {
    fnm <- "data/cukeBarcodes-flagAmb.rds"
    if (file.exists(fnm) && !overwrite) {
        cukeAlg <- readRDS(file=fnm)
    }
    else {
        if (overwrite)  message("Re-creating", fnm)
        algFile <- "data/cukeBarcodes-cleaned.fas"

        ## identify sequences with ambiguities and rename them
        ambSeq <- checkAmbiguity(file=algFile)
        oldNm <- names(ambSeq)
        newNm <- paste(oldNm, "_", sapply(ambSeq, length), "amb", sep="")

        cukeAlg <- read.dna(file=algFile, format="fasta")
        dimnames(cukeAlg)[[1]][match(oldNm, dimnames(cukeAlg)[[1]])] <- newNm
        cukeAlg <- cleanSeqLabels(cukeAlg, software="RAxML")

        saveRDS(cukeAlg, file="data/cukeBarcodes-flagAmb.rds")
    }
    cukeAlg
}

load_cukeTree_raw_phylo4 <- function(overwrite=FALSE) {
    fnm <- "data/cukeTree-raw-phylo4.rds"
    if (file.exists(fnm) && !overwrite) {
        cukeTree4 <- readRDS(fnm)
    }
    else {
        tmpTr <- load_cukeTree_raw()
        cukeTree4 <- build_cukeTree_phylo4(tmpTr)
        saveRDS(cukeTree4, file=fnm)
    }
    cukeTree4
}

load_cukeTree_k2p_phylo4 <- function(overwrite=FALSE) {
    fnm <- "data/cukeTree-k2p-phylo4.rds"
    if (file.exists(fnm) && !overwrite) {
        cukeTree4 <- readRDS(fnm)
    }
    else {
        tmpTr <- load_cukeTree_k2p()
        cukeTree4 <- build_cukeTree_phylo4(tmpTr)
        saveRDS(cukeTree4, file=fnm)
    }
    cukeTree4
}

load_taxonomyDf <- function(overwrite=FALSE) {
    fnm <- "data/taxonomyDf.rds"
    if (file.exists(fnm) && !overwrite) {
        taxDf <- readRDS(fnm)
    }
    else {
        cukeDB <- load_cukeDB()

        uniqOrder <- unique(cukeDB$order)
        uniqFamily <- unique(cukeDB$family)
        uniqFamily <- uniqFamily[! uniqFamily %in%
                                 c("Dactylochirotida", "?", "Uncertain")]
        uniqTaxa <- c(uniqOrder, uniqFamily)
        stopifnot(! any(duplicated(uniqTaxa)))
        taxonomyDf <- data.frame(rank = c(rep("Order", length(uniqOrder) + 1),
                                     rep("Family", length(uniqFamily))),
                                 taxa = c("all", uniqOrder, uniqFamily))
    }
}
