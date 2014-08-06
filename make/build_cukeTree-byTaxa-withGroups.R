
source("R/load.R")
source("R/findGroups.R")
source("R/removeNodeLabels.R")

phylobase.options(allow.duplicated.labels="ok")

inputFiles <- c("data/cukeTree-raw-phylo4.rds",
                "data/cukeTree-k2p-phylo4.rds")

## computer number of species for threshold at
##  1%, 1.5%, 2%, 2.5%, 3%, 3.5%, 4%, 4.5%, 5%, 6%, 7%, 8%
thresVec <- c(seq(1, 5, by=.5), 6:8)/200

## To find groups based on the partial dataset
cukeDB <- load_cukeDB()

uniqOrder <- unique(cukeDB$order)
uniqFamily <- unique(cukeDB$family)
uniqFamily <- uniqFamily[! uniqFamily %in% c("Dactylochirotida", "?", "Uncertain")]
uniqTaxa <- c(uniqOrder, uniqFamily)
stopifnot(! any(duplicated(uniqTaxa)))

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
