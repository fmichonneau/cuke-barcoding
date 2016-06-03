load_taxonomy_df <- function(cukeDB) {
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
    taxonomyDf
}
