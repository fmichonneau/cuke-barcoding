make_all_spp <- function(cuke_db, cuke_alg, taxonomy) {
    cuke_db_tmp <- cuke_db[match(dimnames(cuke_alg)[[1]],
                                 cuke_db$guid), ]
    genera <- unique(cuke_db_tmp$genusorhigher)
    genera <- genera[-grep("\\?", genera)]
    all_spp <- cuke_db_tmp[cuke_db_tmp$genusorhigher %in% genera,
                        c("family", "genusorhigher", "species")]
    all_spp <- all_spp[!is.na(all_spp$family), ]
    all_spp
}


make_uniq_spp <- function(all_spp, taxonomy) {
    uniq_spp <- unique(paste(all_spp$family, all_spp$genusorhigher, all_spp$species, sep="_"))
    uniq_spp <- uniq_spp[-grep("\\d", uniq_spp)]
    uniq_spp <- uniq_spp[-grep("n_?sp|new\\s?species|unique\\s?species", uniq_spp)]
    uniq_spp <- uniq_spp[-grep("(cf|aff)\\.?", uniq_spp)]
    uniq_spp <- uniq_spp[-grep("_.{1,3}$", uniq_spp)] # too generous, removes H. bo
    uniq_spp <- uniq_spp[-grep("_$", uniq_spp)]

    fams <- sapply(strsplit(uniq_spp, "_"), function(x) x[1])
    uniq_spp <- data.frame(family=fams, species=uniq_spp)
    uniq_spp <- merge(uniq_spp, taxonomy, by.x="family", by.y="taxa", all.x=TRUE)
    uniq_spp
}


calc_nspp <- function(uniq_spp) {
    list(
        all = nrow(uniq_spp),
        asp = nrow(subset(uniq_spp, higher == "Aspidochirotida")),
        hol = nrow(subset(uniq_spp, family=="Holothuriidae")),
        apo = nrow(subset(uniq_spp, higher=="Apodida")),
        den = nrow(subset(uniq_spp, higher=="Dendrochirotida")),
        ela = nrow(subset(uniq_spp, higher=="Elasipodida"))
    )
}

make_undesc_spp <- function(all_spp) {
    undesc_spp <- unique(paste(all_spp$family, all_spp$genusorhigher, all_spp$species, sep="_"))
    undesc_spp <- grep("(n)?_sp(_|\\.|\\s|nov)?\\d?", undesc_spp, value=TRUE)
    undesc_spp <- undesc_spp[-grep("spic|spin|spect", undesc_spp)]
    undesc_spp
}

calc_prop_hol <- function(cuke_db, cuke_alg) {
    db <- cuke_db[cuke_db$guid %in% dimnames(cuke_alg)[[1]], ]
    n_hol <- nrow(cuke_db[cuke_db$family == "Holothuriidae", ])
    100 * n_hol/nrow(db)
}
