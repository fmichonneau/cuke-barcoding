
library(ape)
library(phylobase)
library(doMC)
registerDoMC()

source("R/findGroups.R")

cuke_findGroups <- function(tree, threshold=0.015) {
    if (! require(ape) && require(phylobase) &&
        require(doMC)) {
        stop("problem")
    }

    tree$edge.length[tree$edge.length < 0] <- 1e-6
    ## TODO -- double check that using 1 blinding doesn't cause
    ## any issues
    treeRooted <- ape::root(tree, 1, resolve.root=TRUE)
    treeP4 <- as(treeRooted, "phylo4")
    bs <- nodeLabels(treeP4)
    is.na(bs[as.character(rootNode(treeP4))]) <- TRUE
    bs <- data.frame(bs, stringsAsFactors=FALSE)
    bs$bs <- as.numeric(bs$bs)
    treeP4 <- phylo4d(treeP4, node.data=bs)
    treeP4 <- removeNodeLabels(treeP4)
    
    treegr <- findGroups(treeP4, threshold=threshold,
                         experimental=FALSE, parallel=TRUE)
    treegr
}

inputFiles <- c("data/cukeTree-raw.rds",
                "data/cukeTree-k2p.rds")

## computer number of species for threshold at
##  1%, 1.5%, 2%, 2.5%, 3%, 3.5%, 4%, 4.5%, 5%, 6%, 7%, 8%
thresVec <- c(seq(1, 5, by=.5), 6:8)/200

for (j in 1:length(inputFiles)) {
    treeTmp <- readRDS(file=inputFiles[j])

    outputFiles <- paste(gsub(".rds$", "-", inputFiles[j]),
                         gsub("\\.", "", thresVec), ".rds", sep="")
    stopifnot(length(outputFiles) == length(thresVec))

    for (i in 1:length(thresVec)) {
        treeRawTmp <- cuke_findGroups(tree=treeTmp, threshold=thresVec[i])
        saveRDS(treeRawTmp, file=outputFiles[i])
    }
}
