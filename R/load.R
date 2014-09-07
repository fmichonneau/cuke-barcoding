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
        echinoDB <- load_echinoDB(overwrite)
        cukeDB <- subset(echinoDB, class_ == "Holothuroidea")

        ## ## Select sequences TODO -- need to reincorporate this!
        ## holDB <- subset(allDB, class_ == "Holothuroidea")  # nrow = 4385
        ## holDB <- subset(holDB, pass.seq != "GenBank")      # nrow = 4360
        ## holDB <- subset(holDB, pass.seq != "fix")          # nrow = 4358
        ## holDB <- subset(holDB, pass.seq != "no_seq_yet")   # nrow = 3466
        ## holDB <- subset(holDB, pass.seq != "no")           # nrow = 3402
        ## holDB <- subset(holDB, Notes != "MH sequence")     # nrow = 3379
        ## holDB <- subset(holDB, pass.seq != "duplicate")    #

        ## lSeq <- sapply(holDB$Sequence, function(x) length(gregexpr("[actgACTG]", x)[[1]])) # only non-ambiguous bp

        ## lAmb <- sapply(holDB$Sequence, function(x) length(gregexpr("[^-]", x)[[1]]))       # all bp

        ## holDB <- holDB[lAmb > 500, ] # nrow = 2894 -- this also takes care of empty sequences (only -)
        ## check for duplicated samples

        ## dup <- holDB[duplicated(holDB$Sample), "Sample"]
        ## stopifnot(length(dup) == 0)

        ## Remove sequences with internal gaps and stop codons
        ## toRm <- union(seqWithStop, seqWithIntGap)

        ## dimnames(seqHol)[[1]][match(toRm, dimnames(seqHol)[[1]])] <- paste("stop-intgap", toRm, sep="_")
        ## toRmInd <- match(toRm, dimnames(seqHol)[[1]])
        ## seqHol <- seqHol[-toRmInd, ]

        ## These 3 sequences are not represented by other representative

        ##  it might be worth trying to figure out if we can clean up the
        ##  sequences to deal with the issues
        ##  - FRM-194
        ##  - NMV F112128
        ##  - NIWA 38032

        ## Write working copy of fasta file
        write.dna(seqHol, file="data/workingAlg.fas", format="fasta", colsep="")

        ## fix GPS coordinates
        cukeDB[cukeDB$decimalLatitude==19.95 & cukeDB$Loc == "MexicoPac",
               c("decimalLatitude", "decimalLongitude")] <- data.frame(20.3, -105.5)
        cukeDB[cukeDB$Loc == "Tanzania", "decimalLatitude"] <- -cukeDB[cukeDB$Loc == "Tanzania", "decimalLatitude"]
        cukeDB[cukeDB$Loc == "Eparses", "decimalLatitude"] <- -cukeDB[cukeDB$Loc == "Eparses", "decimalLatitude"]
        cukeDB[cukeDB$Sample == "MOLAF_0139",
               c("decimalLatitude", "decimalLongitude")] <- data.frame(-9.536, 147.289)
        cukeDB[cukeDB$Sample == "8928",
               c("decimalLatitude", "decimalLongitude")] <- data.frame(-6.14623, 39.13786)
        cukeDB[cukeDB$Sample == "FRM069",
               c("decimalLatitude", "decimalLongitude")] <- data.frame(NA, NA) ## prob cont.
        cukeDB[cukeDB$Sample == "8919",
               c("decimalLatitude", "decimalLongitude")] <- data.frame(-6.125, 39.191)
        cukeDB[cukeDB$Sample == "RUMF-ZE-00072",  ## not from okinawa but Xmas Island
               c("decimalLatitude", "decimalLongitude")] <- data.frame(-10.501823, 105.685488)
        cukeDB[cukeDB$Sample == "8858F",
               c("decimalLatitude", "decimalLongitude")] <- data.frame(NA, NA) ## prob cont.
         cukeDB[cukeDB$Sample == "9166",
                c("decimalLatitude", "decimalLongitude")] <- data.frame(-22.33917, 40.3388) ## prob cont.
        cukeDB[cukeDB$Sample == "MOLAF_0108",
               c("decimalLatitude", "decimalLongitude")] <- data.frame(-9.536, 147.289)
        cukeDB[cukeDB$Sample == "8931",
               c("decimalLatitude", "decimalLongitude")] <- data.frame(-6.14623, 39.13786)
        cukeDB[cukeDB$Sample == "8932",
               c("decimalLatitude", "decimalLongitude")] <- data.frame(-6.14623, 39.13786)
         cukeDB[cukeDB$Sample == "9190",
               c("decimalLatitude", "decimalLongitude")] <- data.frame(-22.34657, 40.33203)
         cukeDB[cukeDB$Sample == "8937",
               c("decimalLatitude", "decimalLongitude")] <- data.frame(-6.38733, 39.28727)
        cukeDB[cukeDB$Sample == "Hickman_needed3",
               c("decimalLatitude", "decimalLongitude")] <- data.frame(-0.41, -91.48)
         cukeDB[cukeDB$Sample == "8923",
               c("decimalLatitude", "decimalLongitude")] <- data.frame(-6.146230, 39.13786)
         cukeDB[cukeDB$Sample == "8933",
               c("decimalLatitude", "decimalLongitude")] <- data.frame(-6.38733, 39.28727)
         cukeDB[cukeDB$Sample == "8930",
               c("decimalLatitude", "decimalLongitude")] <- data.frame(-6.146230, 39.13786)
         cukeDB[cukeDB$Sample == "8938",
               c("decimalLatitude", "decimalLongitude")] <- data.frame(-6.38733, 39.28727)
        cukeDB[cukeDB$Sample == "6355",
               c("decimalLatitude", "decimalLongitude")] <- data.frame(-21.1008, 55.2437)
        cukeDB[cukeDB$Sample == "6322",
               c("decimalLatitude", "decimalLongitude")] <- data.frame(-21.0, 55.2437)

        ## check taxonomy
        testGenera <-  as.matrix(xtabs(~ genusorhigher + family, data=cukeDB, subset=family != "Uncertain"))
        resGenera <- apply(testGenera, 1, function(x) sum(x != 0))
        stopifnot(all(resGenera == 1))
        testFamily <- as.matrix(xtabs(~ family + order, data=cukeDB, subset=family != "Uncertain"))
        resFamily <- apply(testFamily, 1, function(x) sum(x != 0))
        stopifnot(all(resFamily == 1))

        ## add labels
        treeH <- load_cukeTree_raw_phylo4()

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
    invisible(cukeDB_lbls)
}

load_cukeDist_raw <- function(overwrite=FALSE, ...) {
    fnm <- "data/cukeDist-raw.rds"
    if (file.exists(fnm) && !overwrite) {
        cukeDist <- readRDS(file=fnm)
    } else {
        cukeDist <- ape::dist.dna(load_cukeAlg(), as.matrix=TRUE, model="raw")
        saveRDS(cukeDist, file=fnm)
    }
    invisible(cukeDist)
}

load_cukeTree_raw <- function(overwrite=FALSE, ...) {
    fnm <- "data/cukeTree-raw.rds"
    if (file.exists(fnm) && !overwrite) {
        cukeTree <- readRDS(file=fnm)
        cukeTree$tip.label <- gsub("\\\"", "", cukeTree$tip.label)
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
        cukeTree$tip.label <- gsub("\\\"", "", cukeTree$tip.label)
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
        testFamily <- as.matrix(xtabs(~ family + order, data=cukeDB,
                                      subset = !family %in% c("Dactylochirotida", "?", "Uncertain")))
        whichOrder <- apply(testFamily, 1, function(x) which(x != 0))
        tmpFamily <- data.frame(taxa = names(whichOrder),
                                higher=dimnames(testFamily)[[2]][whichOrder])
        taxonomyDf <- merge(taxonomyDf, tmpFamily, all.x=TRUE)
        saveRDS(taxonomyDf, file=fnm)
    }
    taxonomyDf
}

load_labelsFromTaxa <- function(taxa="all") {
    taxonomyDf <- load_taxonomyDf()
    taxa <- match.arg(as.character(taxa), taxonomyDf$taxa)
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

load_tree_raxml <- function(overwrite=FALSE) {
    fnm <- "data/cukeTree-raxml.rds"
    origTreeNm <- "data/raxml/RAxML_bipartitions.cukeBarcodes"
    if (file.exists(fnm) && !overwrite) {
        tree <- readRDS(file=fnm)
    }
    else {
        if (file.exists(origTreeNm) && !overwrite) {
            tree <- ape::read.tree(file=origTreeNm)
            saveRDS(tree, fnm)
        }
        else {
            build_raxml_tree()
            tree <- load_tree_raxml(FALSE)
        }
    }
    tree
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

load_tree_raxml_phylo4 <- function(taxa="all", overwrite=FALSE) {
    fnm <- "data/cukeTree-raxml-phylo4.rds"
    if (file.exists(fnm) && !overwrite) {
        tr <- readRDS(fnm)
    }
    else {
        tr <- load_tree_raxml()
        ## TODO -- double check that using 1 blindly doesn't cause
        ## any issues
        ## TODO -- use mid-point rooting, or rooting on longest branch
        tr <- ape::root(tr, 1, resolve.root=TRUE)
        tr <- as(tr, "phylo4")
        bs <- nodeLabels(tr)
        bs <- data.frame(bs, stringsAsFactors=FALSE)
        bs$bs <- as.numeric(bs$bs)
        tr <- phylo4d(tr, node.data=bs)
        tr <- removeNodeLabels(tr)
        saveRDS(tr, file=fnm)
    }

    if (taxa == "all") {
        invisible(tr)
    }
    else {
        toKeep <- load_labelsFromTaxa(taxa)
        tr <- subset(tr, tips.include=toKeep)
        stopifnot(all(toKeep %in% tipLabels(tr)) ||
                  all(!is.na(toKeep)))
        invisible(tr)
    }
}


load_species_pairwiseGrps <- function(distance="raw", taxa="all",
                                      threshold=0.03) {

    taxonomyDf <- load_taxonomyDf()
    thres <- load_thresholdPairwise()

    taxa <- match.arg(as.character(taxa), taxonomyDf$taxa)
    distance <- match.arg(distance, c("raw", "K80"))
    stopifnot(length(threshold) == 1 && threshold %in% thres)

    pairwiseGrpRes <- load_pairwiseGrpRes()
    nmRes <- paste(distance, taxa, sep="-")
    res <- pairwiseGrpRes[[match(nmRes, names(pairwiseGrpRes))]][[which(thres == threshold)]]
    stopifnot(! is.null(res))
    res
}

load_tree_pairwiseGrps <- function(distance="raw", taxa="all",
                                   threshold=0.03) {
    spp <- load_species_pairwiseGrps(distance=distance, taxa=taxa,
                                     threshold=threshold)
    lSpp <- sapply(spp, length)
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

load_tree_manualGrps <- function(taxa="Holothuriidae",
                                 overwrite=FALSE) {
    fnm <- "data/cukeTree-manualESUs.rds"
    ##distance <- match.arg(distance, c("raw", "K80"))
    taxa <- match.arg(taxa) ## only Holothuriidae for now
    if (file.exists(fnm) && !overwrite) {
        manESU <- readRDS(file=fnm)
    }
    else {
        ## holTree <- load_tree_clusterGrps("raw", "Holothuriidae", threshold=0.02)
        ## write.csv(tdata(holTree, "tip"), file="data/raw_manualESUs.csv")
        ## ESU coding
        ## - for ESU_genetic: name of the species, (genus_species), followed
        ##   by alphanum to distinguish ESUs, (_1, _1a, _2, _ESU1), followed
        ##   by geography as needed (_IO, _PO):
        ##   - If same complex but different ESUs (e.g., the ind from a region
        ##     form a rec. monophyletic clade): hol_ver_1a_IO and hol_ver_1b_PO
        ##   - If same complex and same ESUs: hol_ver_1_IO and hol_ver_1_PO
        ##     (IO or PO not rec. monophyletic, or sing. ind. w/ low div.)
        ##   - If different ESUs but not sisters: hol_ver_1 and hol_ver_2
        manESU <- read.csv(file="data/raw/manualESUs.csv", stringsAsFactors=FALSE)
        manESU$ESU_noGeo <- gsub("_[A-Z]{2}$", "", manESU$ESU_genetic)
        tree <- load_tree_raxml_phylo4(taxa)
        tDat <- data.frame(as.numeric(factor(manESU$ESU_noGeo)))
        names(tDat) <- "Groups"
        rownames(tDat) <- manESU$Labels
        manESU <- addData(tree, tip.data=tDat)
        saveRDS(manESU, file=fnm)
    }
    invisible(manESU)
}
