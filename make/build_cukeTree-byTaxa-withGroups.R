
library(methods)
library(ape)
library(phylobase)
library(doMC)
registerDoMC()

source("R/findGroups.R")
source("R/removeNodeLabels.R")

phylobase.options(allow.duplicated.labels="ok")


cuke_convertTree <- function(tree) {
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
    treeP4
}

inputFiles <- c("data/cukeTree-raw-phylo4.rds",
                "data/cukeTree-k2p-phylo4.rds")

if (! all(file.exists(inputFiles))) {
    tmpTr <- readRDS(file="data/cukeTree-raw.rds")
    tmpTr <- cuke_convertTree(tmpTr)
    saveRDS(tmpTr, file="data/cukeTree-raw-phylo4.rds")
    tmpTr <- readRDS(file="data/cukeTree-k2p.rds")
    tmpTr <- cuke_convertTree(tmpTr)
    saveRDS(tmpTr, file="data/cukeTree-k2p-phylo4.rds")
    stopifnot(all(file.exists(inputFiles)))
}

## computer number of species for threshold at
##  1%, 1.5%, 2%, 2.5%, 3%, 3.5%, 4%, 4.5%, 5%, 6%, 7%, 8%
thresVec <- c(seq(1, 5, by=.5), 6:8)/200

## To find groups based on the partial dataset
cukeDB <- readRDS(file="data/cukeDB_withLabels.rds")
uniqOrder <- unique(cukeDB$order)
uniqFamily <- unique(cukeDB$family)
uniqFamily <- uniqFamily[! uniqFamily %in% c("Dactylochirotida", "?", "Uncertain")]
uniqTaxa <- c(uniqOrder, uniqFamily)
stopifnot(! any(duplicated(uniqTaxa)))
overwrite <- FALSE

taxonomyDf <- data.frame(rank = c(rep("Order", length(uniqOrder) + 1),
                             rep("Family", length(uniqFamily))),
                         taxa = c("all", uniqOrder, uniqFamily))
saveRDS(taxonomyDf, file="data/taxonomyDf.rds")

for (j in 1:length(inputFiles)) {
    treeTmp <- readRDS(file=inputFiles[j])
    
    for (k in 1:length(uniqTaxa)) {
        outputFiles <- paste(gsub("phylo4.rds$", "", inputFiles[j]),
                            uniqTaxa[k], "-",
                            gsub("\\.", "", thresVec), ".rds", sep="")
        if (k <= length(uniqOrder)) {
            toKeep <- cukeDB[cukeDB$order == uniqTaxa[k], "Labels_withAmb"]
        } else {
            toKeep <- cukeDB[cukeDB$family == uniqTaxa[k], "Labels_withAmb"]
        }
        
        for (i in 1:length(thresVec)) {
            if (!file.exists(outputFiles[i]) || overwrite) {
                stopifnot(all(toKeep %in% tipLabels(treeTmp)) || all(!is.na(toKeep)))
                treeTmpSub <- subset(treeTmp, tips.include=toKeep)
                treeTmpGrp <- findGroups(treeTmpSub, threshold=thresVec[i],
                                         experimental=FALSE, parallel=TRUE)
                saveRDS(treeTmpGrp, file=outputFiles[i])
            }  else {
                message(outputFiles[i], " already exists.")
            }
        }
    }
}
