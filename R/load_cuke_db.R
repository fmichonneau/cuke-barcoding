load_echino_db <- function(cuke_file) {
    echino_db <- read.csv(file = cuke_file, stringsAsFactors = FALSE)
    echino_db
}

load_cuke_db <- function(echino_db) {

    cuke_db <- subset(echino_db, class == "Holothuroidea")
    cuke_db <- subset(cuke_db, ! pass.seq %in%  c("GenBank",
                                                "fix",
                                                "no_seq_yet",
                                                "no",
                                                "duplicate"))
    cuke_db <- subset(cuke_db, Notes != "MH sequence")
    cuke_db <- cuke_db[nzchar(cuke_db$Sample), ]

    ## Taxonomic check
    test_genera <- as.matrix(xtabs(~ genusorhigher + family, data = cuke_db,
                                   subset = family != "Uncertain"))
    res_genera <- apply(test_genera, 1, function(x) sum(x != 0))
    stopifnot(all(res_genera == 1))

    test_family <- as.matrix(xtabs(~ family + order, data = cuke_db,
                                   subset = family != "Uncertain"))
    res_family <- apply(test_family, 1, function(x) sum(x != 0))
    stopifnot(all(res_family == 1))

    ## only non-ambiguous bp
    #l_seq <- sapply(cuke_db$Sequence, function(x)
    #    length(gregexpr("[actgACTG]", x)[[1]]))

    ## all bp
    l_amb <- sapply(cuke_db$Sequence, function(x)
        length(gregexpr("[^-]", x)[[1]]))

    ## this also takes care of empty sequences (only -)
    cuke_db <- cuke_db[l_amb >= 500, ]

    ## check for duplicated samples
    dup <- cuke_db[duplicated(cuke_db$Sample), "Sample"]
    stopifnot(length(dup) == 0)

    ## These 3 sequences are not represented by other representative

    ##  it might be worth trying to figure out if we can clean up the
    ##  sequences to deal with the issues
    ##  - FRM-194
    ##  - NMV F112128
    ##  - NIWA 38032

    ## fix GPS coordinates
    cuke_db[cuke_db$decimalLatitude==19.95 & cuke_db$Loc == "MexicoPac",
           c("decimalLatitude", "decimalLongitude")] <- data.frame(20.3, -105.5)
    cuke_db[cuke_db$Loc == "Tanzania", "decimalLatitude"] <- -cuke_db[cuke_db$Loc == "Tanzania", "decimalLatitude"]
    cuke_db[cuke_db$Loc == "Eparses", "decimalLatitude"] <- -cuke_db[cuke_db$Loc == "Eparses", "decimalLatitude"]
    cuke_db[cuke_db$Sample == "MOLAF_0139",
           c("decimalLatitude", "decimalLongitude")] <- data.frame(-9.536, 147.289)
    cuke_db[cuke_db$Sample == "8928",
           c("decimalLatitude", "decimalLongitude")] <- data.frame(-6.14623, 39.13786)
    cuke_db[cuke_db$Sample == "FRM069",
           c("decimalLatitude", "decimalLongitude")] <- data.frame(NA, NA) ## prob cont.
    cuke_db[cuke_db$Sample == "8919",
           c("decimalLatitude", "decimalLongitude")] <- data.frame(-6.125, 39.191)
    cuke_db[cuke_db$Sample == "RUMF-ZE-00072",  ## not from okinawa but Xmas Island
           c("decimalLatitude", "decimalLongitude")] <- data.frame(-10.501823, 105.685488)
    cuke_db[cuke_db$Sample == "8858F",
           c("decimalLatitude", "decimalLongitude")] <- data.frame(NA, NA) ## prob cont.
    cuke_db[cuke_db$Sample == "9166",
           c("decimalLatitude", "decimalLongitude")] <- data.frame(-22.33917, 40.3388) ## prob cont.
    cuke_db[cuke_db$Sample == "MOLAF_0108",
           c("decimalLatitude", "decimalLongitude")] <- data.frame(-9.536, 147.289)
    cuke_db[cuke_db$Sample == "8931",
           c("decimalLatitude", "decimalLongitude")] <- data.frame(-6.14623, 39.13786)
    cuke_db[cuke_db$Sample == "8932",
           c("decimalLatitude", "decimalLongitude")] <- data.frame(-6.14623, 39.13786)
    cuke_db[cuke_db$Sample == "9190",
           c("decimalLatitude", "decimalLongitude")] <- data.frame(-22.34657, 40.33203)
    cuke_db[cuke_db$Sample == "8937",
           c("decimalLatitude", "decimalLongitude")] <- data.frame(-6.38733, 39.28727)
    cuke_db[cuke_db$Sample == "Hickman_needed3",
           c("decimalLatitude", "decimalLongitude")] <- data.frame(-0.41, -91.48)
    cuke_db[cuke_db$Sample == "8923",
           c("decimalLatitude", "decimalLongitude")] <- data.frame(-6.146230, 39.13786)
    cuke_db[cuke_db$Sample == "8933",
           c("decimalLatitude", "decimalLongitude")] <- data.frame(-6.38733, 39.28727)
    cuke_db[cuke_db$Sample == "8930",
           c("decimalLatitude", "decimalLongitude")] <- data.frame(-6.146230, 39.13786)
    cuke_db[cuke_db$Sample == "8938",
           c("decimalLatitude", "decimalLongitude")] <- data.frame(-6.38733, 39.28727)
    cuke_db[cuke_db$Sample == "6355",
           c("decimalLatitude", "decimalLongitude")] <- data.frame(-21.1008, 55.2437)
    cuke_db[cuke_db$Sample == "6322",
           c("decimalLatitude", "decimalLongitude")] <- data.frame(-21.0, 55.2437)

    invisible(cuke_db)
}

load_cuke_seqs <- function(cuke_db) {
    cuke_db[, c("Sample", "Sequence")]
}

clean_labels <- function(lbls)  {
    gsub(":|,|\\(|\\)|;|\\[|\\]|\\'|\\s|\t", "", lbls)
}


load_cuke_clusters <- function(cuke_db, cuke_tree) {

    dataLbls <- character(nrow(cuke_db))
    for (i in 1:nrow(cuke_db)) {
        dataLbls[i] <- genLabel(cuke_db[i, ])
    }

    ## using the same pattern as in chopper::cleanSeqLabels
    cuke_db$Labels <- clean_labels(dataLbls)

    treeTips <- data.frame(Labels_withAmb = tipLabels(cuke_tree),
                           Labels = gsub("_\\d+amb$", "", tipLabels(cuke_tree)),
                           stringsAsFactors=FALSE)

    stopifnot(all(treeTips$treeLabels %in% cuke_db$Labels))
    cuke_db_lbls <- merge(cuke_db, treeTips, by="Labels")

    invisible(cuke_db_lbls)
}
