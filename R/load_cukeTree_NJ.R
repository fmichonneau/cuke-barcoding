## This should be transfered to phylobase, not sure where it should go
## though...
removeNodeLabels <- function(phy) {
    intNd <- nodeId(phy, "internal")
    is.na(phy@label[intNd]) <- TRUE
    phy
}

build_cukeTree <- function(alg, model, Nrep) {
    model <- match.arg(model, c("raw", "K80"))
    if (identical(model, "raw")) {
        dMat <- load_cukeDist_raw()
    }
    else if (identical(model, "K80")) {
        dMat <- load_cukeDist_k2p()
    } else stop("Houston, we have a problem")

    treH <- ape::nj(dMat)
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

load_cukeTree_raw <- function(cukeAlg, Nrep, ...) {
    if(missing(Nrep)) Nrep <- 200
    cukeTree <- build_cukeTree(cukeAlg, model="raw", Nrep=Nrep, ...)
    invisible(cukeTree)
}

load_cukeTree_k2p <- function(cukeAlg, Nrep, ...) {
    if(missing(Nrep)) Nrep <- 200
    cukeTree <- build_cukeTree(cukeAlg, model="K80", Nrep=Nrep, ...)
    invisible(cukeTree)
}

load_cukeTree_raw_phylo4 <- function(cuke_tree_raw) {
    cukeTree4 <- build_cukeTree_phylo4(cuke_tree_raw)
    invisible(cukeTree4)
}

load_cukeTree_k2p_phylo4 <- function(cuke_tree_k2p) {
    cukeTree4 <- build_cukeTree_phylo4(cuke_tree_k2p)
    invisible(cukeTree4)
}
