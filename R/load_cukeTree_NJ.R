## This should be transfered to phylobase, not sure where it should go
## though...
removeNodeLabels <- function(phy) {
    intNd <- nodeId(phy, "internal")
    is.na(phy@label[intNd]) <- TRUE
    phy
}

build_cuke_tree <- function(alg, dist_mat, Nrep) {
    treH <- ape::nj(dist_mat)
    bootH <- ape::boot.phylo(treH, alg, function(xx) {
        ape::nj(ape::dist.dna(xx, model=model))
    }, B=Nrep)
    treH$node.label <- bootH
    treH
}

build_cuke_tree_phylo4 <- function(tree) {
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

load_cuke_tree <- function(cuke_alg, dist_mat, Nrep, ...) {
    if(missing(Nrep)) Nrep <- 200
    cukeTree <- build_cuke_tree(cuke_alg, dist_mat, Nrep=Nrep, ...)
    invisible(cukeTree)
}

load_cuke_tree_raw_phylo4 <- function(cuke_tree_raw) {
    cukeTree4 <- build_cuke_tree_phylo4(cuke_tree_raw)
    invisible(cukeTree4)
}

load_cuke_tree_k2p_phylo4 <- function(cuke_tree_k2p) {
    cukeTree4 <- build_cuke_tree_phylo4(cuke_tree_k2p)
    invisible(cukeTree4)
}
