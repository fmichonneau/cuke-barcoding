cuke_tree_cluster_groups <- function(taxonomy, cuke_tree_phylo4, cuke_db) {

    st <- storr(driver_environment())
    uniq_taxa <- taxonomy$taxa
    thresVec <- load_thresholdClusters()

    for (j in seq_along(uniq_taxa)) {
        for (i in seq_along(thresVec)) {
            key <- paste(uniq_taxa[j], gsub("\\.", "", thresVec[i]),
                         sep = "-")
            message("Finding groups for ", sQuote(uniq_taxa[j]),
                    " with threshold of ", sQuote(thresVec[i]),
                    " ....", appendLF =  FALSE)
            tmp_grp <- build_cluster_group(tr = cuke_tree_phylo4,
                                           taxa = uniq_taxa[j],
                                           threshold = thresVec[i],
                                           taxonomy = taxonomy,
                                           cuke_db = cuke_db)
            message("DONE.")
            st$set(key, tmp_grp)
        }
    }
    st
}

build_cluster_group <- function(tr, taxa, threshold, taxonomy, cuke_db) {
    lbl_to_keep <- fetch_labels_from_taxa(taxonomy, cuke_db, taxa)
    stopifnot(all(lbl_to_keep %in% tipLabels(tr)))
    sub_tr <- subset(tr, tips.include = lbl_to_keep)
    grp_tr <- findGroups(sub_tr, threshold = threshold,
                         experimental = FALSE, parallel = TRUE)
    grp_tr
}

load_tree_cluster_groups <- function(cluster_store, taxa="all",
                                  threshold=0.015, taxonomy) {

    thres <- load_thresholdClusters()
    taxa <- match.arg(as.character(taxa), taxonomy$taxa)
    stopifnot(length(threshold) == 1 && threshold %in% thres)

    key <- paste(taxa, threshold, sep = "-")

    cluster_store$get(key)
}

load_species_cluster_groups <- function(cluster_store, taxa="all",
                                        threshold=0.015, taxonomy) {

    tr <- load_tree_cluster_groups(cluster_store, taxa, threshold, taxonomy)
    grps <- tdata(tr, "tip")[, "Groups", drop=FALSE]
    gg <- factor(grps$Groups)
    split(rownames(grps), gg)
}
