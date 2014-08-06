
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
