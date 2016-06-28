load_threshold_pairwise <- function() {
    c(seq(1, 5, by=.5), 6:8)/100
}

load_threshold_clusters <- function() {
    load_threshold_pairwise()/2
}


## the ... arguments is used to pass an arbitrary number of trees on
## which the clustering algorithm is used. The number of trees and the
## methods used to build them are specified by the 'methods' argument.
## the phylo4 objects
load_cuke_tree_clusters <- function(taxonomy, cuke_db, methods, ...) {

    trs <- list(...)
    if (length(methods) != length(trs)) {
        stop("The number of trees provided should match the number of methods.")
    }
    names(trs) <- methods

    st <- storr::storr_rds("data/storr_clusters")
    uniq_taxa <- taxonomy$taxa
    thres_vec <- load_threshold_clusters()

    for (k in seq_along(trs)) {
        for (j in seq_along(uniq_taxa)) {
            for (i in seq_along(thres_vec)) {
                key <- paste(uniq_taxa[j],
                             gsub("\\.", "", thres_vec[i]),
                             methods[k],
                             sep = "-")
                message("Finding groups for ", sQuote(uniq_taxa[j]),
                        " with threshold of ", sQuote(thres_vec[i]),
                        " .... ", appendLF =  FALSE)
                tmp_grp <- build_cluster_group(tr = trs[[k]],
                                               taxa = uniq_taxa[j],
                                               threshold = thres_vec[i],
                                               taxonomy = taxonomy,
                                               cuke_db = cuke_db)
                message("DONE.")
                st$set(key, tmp_grp)
            }
        }
    }
    st
}

build_cluster_group <- function(tr, taxa, threshold, taxonomy, cuke_db) {
    lbl_to_keep <- fetch_labels_from_taxa(taxonomy, cuke_db, taxa)
    lbl_to_keep <- intersect(tipLabels(tr), lbl_to_keep)
    sub_tr <- subset(tr, tips.include = lbl_to_keep)
    grp_tr <- findGroups(sub_tr, threshold = threshold,
                         experimental = FALSE, parallel = TRUE)
    grp_tr
}

load_tree_cluster_groups <- function(cluster_store, taxa="all",
                                     threshold=0.015, taxonomy) {

    thres <- load_threshold_clusters()
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
