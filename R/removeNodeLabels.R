## This should be transfered to phylobase, not sure where it should go
## though...

removeNodeLabels <- function(phy) {
    intNd <- nodeId(phy, "internal")
    is.na(phy@label[intNd]) <- TRUE
    phy
}
