is_monophyletic <- function(lbl, tree) {
    if (length(lbl) < 2) {
        NA
    } else {
        desc <- descendants(tree, MRCA(tree, lbl), "tips")
        if (length(desc) == length(lbl))
            TRUE
        else FALSE
    }
}


### Code to create the manual ESU file
if (FALSE) {
    ## To create the backbone of the spreadsheet used to encode the ESUs manually:
    st <- storr::storr_rds("data/storr_clusters/")
    hol_tree <- get_tree_cluster_group(st, "Holothuriidae", "raw", 0.02, fetch("taxonomy"))
    esus <- tdata(hol_tree, "tip")
    esus$labels <- make_labels_from_guids(cuke_db, rownames(esus))
    esus$guids <- rownames(esus)
    rownames(esus) <- NULL
    esus <- esus[order(esus$Groups), c("guids", "labels", "Groups")]
    write.csv(esus, file="data/raw/backbone_manualESUs.csv", row.names=FALSE)
    ## And then to plot the tree to make informed decision about the ESUs:
    tipLabels(hol_tree) <- paste(make_labels_from_guids(fetch("cuke_db"), tipLabels(hol_tree)), tdata(hol_tree, "tip")$Groups, sep = "_")
    hol_tr <- as(hol_tree, "phylo")
    pdf(file = "tmp/holothuriidae_raw_002_with_groups.pdf", height = 50)
    ggtree(hol_tree) + geom_tiplab()
    dev.off()


    ## ESU coding
    ## - for ESU_genetic: name of the species, (genus_species), followed
    ##   by alphanum to distinguish ESUs, (_1, _1a, _2, _ESU1), followed
    ##   by geography as needed (_IO, _PO):
    ##   - If same complex but different ESUs (e.g., the ind from a region
    ##     form a rec. monophyletic clade): hol_ver_1a_IO and hol_ver_1b_PO
    ##   - If same complex and same ESUs: hol_ver_1_IO and hol_ver_1_PO
    ##     (IO or PO not rec. monophyletic, or sing. ind. w/ low div.)
    ##   - If different ESUs but not sisters: hol_ver_1 and hol_ver_2
}


load_man_esu <- function(file) {
    manESU <- read.csv(file=file, stringsAsFactors=FALSE)
    manESU$ESU_noGeo <- gsub("_[A-Z]{2}$", "", manESU$ESU_genetic)
    manESU <- manESU[-grep("\\d+amb$", manESU$Labels), ]
    manESU
}

load_manual_esu_grps <- function(man_esu) {
    split(man_esu$Labels, manESU$ESU_noGeo)
}


load_local_gap <- function(taxa="Holothuriidae", cuke_dist_raw,
                           cuke_db, cuke_alg, species_manual_grp) {

    taxa <- match.arg(taxa) # only Holothuriidae for now

    summary_inter_dist <- function(list_species, cuke_alg) {
        lapply(list_species, function(x) {
            lapply(list_species, function(y) {
                inter_esu_dist(x, y, cuke_dist_raw)
            })
        })
    }

    summ_inter_dist <- summary_inter_dist(species_manual_grp, cuke_alg)

    min_inter <- mapply(function(summ, nm) {
        min_dist <- sapply(summ, function(x) x$min)
        min_dist <- min_dist[-match(nm, names(min_dist))]
        min_dist[which.min(min_dist)]
    }, summ_inter_dist, names(summ_inter_dist))

    intra_dist <- lapply(species_manual_grps, intra_esu_dist, cuke_dist_raw)
    max_intra <- vapply(intra_dist, function(x) x$max, numeric(1))

    esuSpatial <- spatial_from_species(species_manual_grps, cuke_db)
    nm_esu_spatial <- names(esu_spatial[[1]])

    rg_type <- lapply(names(min_inter), function(x) {
        spp <- unlist(strsplit(x, "\\."))
        i <- grep(paste0("^", spp[1], "-"), nm_esu_spatial)
        j <- grep(paste0("^", spp[2], "-"), nm_esu_spatial)
        rangeType(i, j, esu_spatial[[2]])
    })

    local_gap <- data.frame(min_inter, max_intra, row.names=names(min_inter))

    local_gap$rangeType <- sapply(rg_type, function(x) x$rangeType)
    local_gap$rangeType[is.na(local_gap$rangeType)] <- "unknown"
    local_gap <- local_gap[complete.cases(local_gap), ]

    local_gap$species <- NA
    local_gap$species[local_gap$max_intra > local_gap$min_inter] <-
        rownames(local_gap)[local_gap$max_intra > local_gap$min_inter]

    local_gap
}

calc_n_esus <- function(man_esu) {
    length(unique(man_esu$ESU_noGeo))
}

calc_n_cryptic <- function(man_esu) {
    esus <- strsplit(unique(man_esu$ESU_noGeo), "_")
    has_cryptic <- sapply(esus, function(x) length(x) > 2 & length(grep("nsp", x)) < 1)
    esu_cryptic <- esus[has_cryptic]
    esu_cryptic <- sapply(esu_cryptic, function(x) paste0(x[1:2], collapse="_"))
    nCryptic <- length(unique(esu_cryptic))
    nCryptic
}

calc_n_new_spp <- function(man_esu) {
    esus <- strsplit(unique(man_esu$ESU_noGeo), "_")
    sum(sapply(esus, function(x) length(grep("nsp", x)) > 0))
}


get_hol_tree <- function(cuke_db, man_esu, taxonomy) {
    to_keep <- intersect(man_esu$Labels, fetch_guids_from_taxa(taxonomy, cuke_db, taxa = "Holothuriidae"))
    hol_tree <- subset()
}
