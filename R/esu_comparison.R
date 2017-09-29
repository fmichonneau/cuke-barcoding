## compare ESU assignments to manual ESUs. Returns a list with number of ESUs
## that are split, and the number of ESUs that are lumped.
accuracy_grps <- function(tree, manual_esu) {
    clstr_grps <- tdata(tree, "tip")[, "Groups", drop=FALSE]
    clstr_grps <- cbind(guids=rownames(clstr_grps), clstr_grps=clstr_grps)
    man_grps <- manual_esu[, c("guids", "ESU_noGeo")]
    man_grps$man_grps <- as.numeric(factor(manual_esu$ESU_noGeo))
    comp_grps <- merge(clstr_grps, man_grps, by="guids")

    uniq_grps <- unique(comp_grps$man_grps)
    comp_res <- vector("list", length(uniq_grps))

    for (i in seq_along(uniq_grps)) {

        tmp_dt <- subset(comp_grps, man_grps == uniq_grps[i])
        if (length(unique(tmp_dt$Groups)) > 1)
            splits <- unique(tmp_dt$Groups)
        else splits <- NA

        tmp_dt2 <- subset(comp_grps, Groups %in% unique(tmp_dt$Groups))
        if( length(unique(tmp_dt2$man_grps)) > 1)
            lumps <- unique(tmp_dt2$man_grps)
        else lumps <- NA

        comp_res[[i]] <- list(splits=splits, lumps=lumps)
    }
    n_lumps <- sapply(comp_res, function(x) x$lumps)
    n_lumps <- unlist(n_lumps)
    n_lumps <- length(unique(n_lumps[!is.na(n_lumps)]))

    n_splits <- sapply(comp_res, function(x) x$splits)
    n_splits <- unlist(n_splits)
    n_splits <- length(unique(n_splits[!is.na(n_splits)]))
    list(n_lumps=n_lumps, n_splits=n_splits)
}


calculate_accuracy <- function(taxonomy, man_esu) {

    st_clusters <- storr::storr_rds("data/storr_clusters")
    st_pairwise <- storr::storr_rds("data/storr_pairwise")

    setup_pairwise <- expand.grid(
        method = "pairwise",
        distance = c("raw", "k2p"),
        threshold = load_threshold_pairwise(),
        stringsAsFactors = FALSE
    )
    setup_clusters <- expand.grid(
        method = "cluster",
        distance = c("raw", "k2p"),
        threshold = load_threshold_clusters(),
        stringsAsFactors = FALSE
    )
    setup <- dplyr::bind_rows(setup_pairwise, setup_clusters)

    split_lumps <- setup %>%
        purrr::transpose() %>%
        purrr::map_df(function(x) {
                   method <- match.arg(x$method, c("pairwise", "cluster"))
                   if (method == "pairwise")
                       st <- st_pairwise
                   else st <- st_clusters
                   tr <- get_tree_cluster_group(st, taxa = "Holothuriidae",
                                                distance = x$distance,
                                                threshold = x$threshold,
                                                taxonomy)
                   .res <- accuracy_grps(tr, man_esu)
                   list(n_lumps = .res$n_lumps,
                        n_splits = .res$n_splits,
                        n_groups = n_distinct(tdata(tr, "tip")[, "Groups"]))
               })

    dplyr::bind_cols(setup, split_lumps) %>%
        dplyr::mutate(all_error = n_lumps + n_splits,
                      p_lumps = n_splits/n_groups,
                      p_lumps = n_lumps/n_groups)

}


if (FALSE) {

pwiseTreeRaw <- lapply(load_thresholdPairwise(), function(x) load_tree_pairwiseGrps("raw", taxa="Holothuriidae", x))
pwiseTreeK80 <- lapply(load_thresholdPairwise(), function(x) load_tree_pairwiseGrps("K80", taxa="Holothuriidae", x))
clstrTreeRaw <- lapply(load_thresholdClusters(), function(x) load_tree_clusterGrps("raw", taxa="Holothuriidae", x))
clstrTreeK80 <- lapply(load_thresholdClusters(), function(x) load_tree_clusterGrps("K80", taxa="Holothuriidae", x))

pwiseGrpsRaw <- lapply(pwiseTreeRaw, accuracy_grps)
pwiseGrpsK80 <- lapply(pwiseTreeK80, accuracy_grps)
clstr_grpsRaw <- lapply(clstrTreeRaw, accuracy_grps)
clstr_grpsK80 <- lapply(clstrTreeK80, accuracy_grps)

pwiseDatRaw <- data.frame(method="pairwise", distance="raw", do.call("rbind", lapply(pwiseGrpsRaw, function(x) c(x[1], x[2]))))
pwiseDatK80 <- data.frame(method="pairwise", distance="K2P", do.call("rbind", lapply(pwiseGrpsK80, function(x) c(x[1], x[2]))))
clstrDatRaw <- data.frame(method="clustering", distance="raw", do.call("rbind", lapply(clstr_grpsRaw, function(x) c(x[1], x[2]))))
clstrDatK80 <- data.frame(method="clustering", distance="K2P", do.call("rbind", lapply(clstr_grpsK80, function(x) c(x[1], x[2]))))

compareManESUs <- rbind(pwiseDatRaw, pwiseDatK80, clstrDatRaw, clstrDatK80)

nGrps <- sapply(c(pwiseTreeRaw, pwiseTreeK80, clstrTreeRaw, clstrTreeK80),
                function(x) max(tdata(x, "tip")[, "Groups"]))

compareManESUs <- data.frame(threshold=rep(load_thresholdPairwise(), 4), compareManESUs, nGrps=nGrps)
compareManESUs$n_splits <- as.numeric(compareManESUs$n_splits)
compareManESUs$n_lumps <- as.numeric(compareManESUs$n_lumps)
allErrors <- compareManESUs$n_splits + compareManESUs$n_lumps
compareManESUs <- cbind(compareManESUs, allErrors=allErrors)

compareManESUs$pLumps <- compareManESUs$n_lumps/compareManESUs$nGrps
compareManESUs$pSplits <- compareManESUs$n_splits/compareManESUs$nGrps

saveRDS(compareManESUs, file="tmp/compareManESUs.rds")

minError <- compareManESUs[which(allErrors <= min(allErrors)), ]

if (nrow(minError) > 1) {
    minError <- minError[which.min(nESUs - minError$nGrps), ]
}
}
