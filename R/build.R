source("R/packages.R")

build_cukeTree <- function(alg, model, Nrep) {
    treH <- ape::nj(ape::dist.dna(alg, model=model))
    bootH <- ape::boot.phylo(treH, alg, function(xx) {
        ape::nj(ape::dist.dna(xx, model=model))
    }, B=Nrep)
    treH$node.label <- bootH
    treH
}

build_cukeTree_phylo4 <- function(tree) {
    tree$edge.length[tree$edge.length < 0] <- 1e-6
    ## TODO -- double check that using 1 blindly doesn't cause
    ## any issues
    ## TODO -- use mid-point rooting, or rooting on longest branch
    treeRooted <- ape::root(tree, 1, resolve.root=TRUE)
    treeP4 <- as(treeRooted, "phylo4")
    bs <- nodeLabels(treeP4)
    is.na(bs[as.character(rootNode(treeP4))]) <- TRUE
    bs <- data.frame(bs, stringsAsFactors=FALSE)
    bs$bs <- as.numeric(bs$bs)
    treeP4 <- phylo4d(treeP4, node.data=bs)
    treeP4 <- removeNodeLabels(treeP4)
    treeP4
}

build_PTP_alg <- function() {
    alg <- ape::read.dna(file="data/cukeBarcodes-flagAmb.phy.reduced",
                         format="sequential")
    nmAlg <- dimnames(alg)[[1]]
    mtchNm <- data.frame(origNm=nmAlg, newNm=1:length(nmAlg),
                         stringsAsFactors=FALSE)
    saveRDS(mtchNm, file="data/match-cukeBarcodes-labels.rds")
    dimnames(alg)[[1]] <- mtchNm$newNm
    ape::write.dna(alg, file="data/cukeBarcodes-flagAmb-reduced-numbered.phy",
                   format="sequential", colsep="", colw=dim(alg)[2])
}

build_PTP_tree <- function() {
    ptpDir <- "data/raxml_ptp"
    if (! file.exists(ptpDir)) dir.create(ptpDir)
    raxmlCmd <- paste("raxmlHPC-PTHREADS-SSE3",
                      "-s data/cukeBarcodes-flagAmb-reduced-numbered.phy",
                      "-m GTRGAMMA -q data/cukeBarcodes-partition",
                      "-f a -p 10101 -x 10101 -# 500 -n cukeBarcodesPTP",
                      "-T8 -w", file.path(getwd(), ptpDir))

    system(raxmlCmd)
}
