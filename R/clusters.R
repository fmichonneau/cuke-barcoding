load_threshold_pairwise <- function() {
    c(seq(1, 5, by=.5), 6:8)/100
}

load_threshold_clusters <- function() {
    load_threshold_pairwise()/2
}


build_cluster_key <- function(taxa, threshold, distance) {
    ## the keys for the store holding the cluster info looks like
    ## Phyllophoridae-0015-k2p
    paste(taxa, gsub("\\.", "", threshold), distance, sep = "-")
}

## the ... arguments is used to pass an arbitrary number of trees on
## which the clustering algorithm is used. The number of trees and the
## methods used to build them are specified by the 'distances' argument.
## the phylo4 objects
make_cuke_tree_clusters <- function(taxonomy, cuke_db,
                                    distances = c("raw", "k2p"), ...) {

    ## extract the trees from ...
    trs <- list(...)

    distances <- match.arg(distances, several.ok = TRUE)

    ## make sure that the number of distances specified matches the
    ## number of trees (built using the distances) specified
    if (length(distances) != length(trs)) {
        stop("The number of trees provided should match the number of distances.")
    }
    names(trs) <- distances

    st <- storr::storr_rds("data/storr_clusters")
    uniq_taxa <- taxonomy$taxa
    thres_vec <- load_threshold_clusters()

    for (k in seq_along(trs)) {
        for (j in seq_along(uniq_taxa)) {
            for (i in seq_along(thres_vec)) {
                key <- build_cluster_key(uniq_taxa[j], thres_vec[i],
                                         distances[k])
                message("Finding groups for ", sQuote(uniq_taxa[j]),
                        " with threshold of ", sQuote(thres_vec[i]),
                        " on tree built using ", sQuote(distances[k]), " distance",
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
    if (hasDuplicatedLabels(tr)) {
        stop("The tree can't have duplicated labels")
    }
    lbl_to_keep <- fetch_guids_from_taxa(taxonomy, cuke_db, taxa)
    lbl_to_keep <- intersect(tipLabels(tr), lbl_to_keep)
    sub_tr <- subset(tr, tips.include = lbl_to_keep)
    grp_tr <- findGroups(sub_tr, threshold = threshold)
    grp_tr
}

get_tree_cluster_group <- function(cluster_store, taxa="all",
                                   distance = c("raw", "k2p"),
                                   threshold=0.015, taxonomy) {

    threshold <- match.arg(as.character(threshold), load_threshold_clusters())
    taxa <- match.arg(as.character(taxa), taxonomy$taxa)
    distance <- match.arg(distance)

    key <- build_cluster_key(taxa, threshold, distance)

    cluster_store$get(key)
}

load_species_cluster_groups <- function(cluster_store, taxa="all", distance,
                                        threshold=0.015, taxonomy) {
    tr <- get_tree_cluster_group(cluster_store, taxa, distance, threshold, taxonomy)
    grps <- tdata(tr, "tip")[, "Groups", drop=FALSE]
    grps$Groups <- factor(grps$Groups)
    split(rownames(grps), grps)
}
