build_cukeTree_clusterGrps <- function(overwrite=FALSE) {
    ## To find groups based on the partial dataset
    cukeDB <- load_cukeDB()
    taxonomyDf <- load_taxonomyDf()

    uniqTaxa <- taxonomyDf$taxa
    stopifnot(! any(duplicated(uniqTaxa)))

    inputFiles <- c("data/cukeTree-raw-phylo4.rds",
                    "data/cukeTree-k2p-phylo4.rds")

    ## computer number of species for threshold at
    ##  1%, 1.5%, 2%, 2.5%, 3%, 3.5%, 4%, 4.5%, 5%, 6%, 7%, 8%
    thresVec <- load_thresholdClusters()

    for (j in 1:length(inputFiles)) {

        for (k in 1:length(uniqTaxa)) {

            outputFiles <- paste(gsub("phylo4.rds$", "", inputFiles[j]),
                                 uniqTaxa[k], "-",
                                 gsub("\\.", "", thresVec), ".rds", sep="")

            toKeep <- load_labelsFromTaxa(uniqTaxa[k])

            for (i in 1:length(thresVec)) {
                if (!file.exists(outputFiles[i]) || overwrite) {
                    treeTmp <- readRDS(file=inputFiles[j])
                    stopifnot(all(toKeep %in% tipLabels(treeTmp)) ||
                              all(!is.na(toKeep)))
                    if (length(toKeep) == nTips(treeTmp)) {
                        treeTmpSub <- treeTmp
                    }
                    else {
                        treeTmpSub <- subset(treeTmp, tips.include=toKeep)
                    }
                    treeTmpGrp <- findGroups(treeTmpSub, threshold=thresVec[i],
                                             experimental=FALSE, parallel=TRUE)
                    saveRDS(treeTmpGrp, file=outputFiles[i])
                }  else {
                    message(outputFiles[i], " already exists.")
                }
            }
        }
    }

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
