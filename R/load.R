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
        dimnames(cukeAlg)[[1]] <- gsub("\\\"", "", dimnames(cukeAlg)[[1]])
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
        taxonomyDf <- readRDS(fnm)
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
        stopifnot(! any(duplicated(taxonomyDf$taxa)))
        saveRDS(taxonomyDf, file=fnm)
    }
    taxonomyDf
}

load_labelsFromTaxa <- function(taxa="all") {
    taxonomyDf <- load_taxonomyDf()
    taxa <- match.arg(taxa, taxonomyDf$taxa)
    cukeDB <- load_cukeDB()

    if (taxa == "all") {
        lbls <- cukeDB[, "Labels_withAmb"]
    } else if (taxonomyDf[taxonomyDf$taxa == taxa, "rank"] == "Order") {
        lbls <- cukeDB[cukeDB$order == taxa, "Labels_withAmb"]
    } else if (taxonomyDf[taxonomyDf$taxa == taxa, "rank"] == "Family") {
        lbls <- cukeDB[cukeDB$family == taxa, "Labels_withAmb"]
    } else {
        stop("something is wrong with ", taxa)
    }
    lbls
}

load_thresholdPairwise <- function() {
    c(seq(1, 5, by=.5), 6:8)/100
}

load_thresholdClusters <- function() {
    load_thresholdPairwise()/2
}

load_tree_phylo4 <- function(distance="raw", taxa="all") {
    distance <- match.arg(distance, c("raw", "K80"))
    if(identical(distance, "raw")) {
        tree <- load_cukeTree_raw_phylo4()
    }
    else {
        tree <- load_cukeTree_k2p_phylo4()
    }

    if (taxa == "all") {
        tree
    }
    else {
        toKeep <- load_labelsFromTaxa(taxa)
        tree <- subset(tree, tips.include=toKeep)
        stopifnot(all(toKeep %in% tipLabels(tree)) ||
                  all(!is.na(toKeep)))
        tree
    }
}

load_species_pairwiseGrps <- function(distance="raw", taxa="all",
                                      threshold=0.03) {

    taxonomyDf <- load_taxonomyDf()
    thres <- load_thresholdPairwise()

    taxa <- match.arg(as.character(taxa), taxonomyDf$taxa)
    distance <- match.arg(distance, c("raw", "K80"))
    stopifnot(length(threshold) == 1 && threshold %in% thres)

    pairwseGrpRes <- load_pairwiseGrpRes()
    nmRes <- paste(distance, taxa, sep="-")
    res <- pairwiseGrpRes[[match(nmRes, names(pairwiseGrpRes))]][which(thres == threshold)]
    stopifnot(! is.null(res[[1]]))
    res[[1]]
}

load_tree_pairwiseGrps <- function(distance="raw", taxa="all",
                                   threshold=0.03) {
    spp <- load_species_pairwiseGrps(distance=distance, taxa=taxa,
                                     threshold=threshold)
    lSpp <- sapply(spp[[1]], length)
    tmpGrps <- mapply(rep, 1:length(lSpp), lSpp)
    tmpGrps <- data.frame(Groups=unlist(tmpGrps))
    rownames(tmpGrps) <- unlist(spp)

    tree <- load_tree_phylo4(distance=distance, taxa=taxa)
    addData(tree, tip.data=tmpGrps)
}

load_tree_clusterGrps <- function(distance="raw", taxa="all",
                                  threshold=0.015) {

    taxonomyDf <- load_taxonomyDf()
    thres <- load_thresholdClusters()
    taxa <- match.arg(as.character(taxa), taxonomyDf$taxa)
    distance <- match.arg(distance, c("raw", "K80"))
    stopifnot(length(threshold) == 1 && threshold %in% thres)

    distance <- ifelse(identical(distance, "raw"), "raw", "k2p")

    lstFiles <- list.files(path="data", pattern="cukeTree-.+-\\d+.rds$",
                           full.names=TRUE)

    nmRes <- paste("cukeTree", distance, taxa, gsub("\\.", "", threshold), sep="-")
    nmRes <- file.path("data", paste0(nmRes, ".rds"))

    if (file.exists(nmRes))
        readRDS(file=nmRes)
    else
        stop("can't find ", nmRes)
}

load_species_clusterGrps <- function(distance="raw", taxa="all",
                                     threshold=0.015) {
    tr <- load_tree_clusterGrps(distance, taxa, threshold)
    grps <- tdata(tr, "tip")[, "Groups", drop=FALSE]
    gg <- factor(grps$Groups)
    split(rownames(grps), gg)
}
