load_echino_db <- function(cuke_file) {
    echino_db <- read.csv(file = cuke_file, stringsAsFactors = FALSE)
    if (length(nzchar(echino_db$guid)) !=
        length(echino_db$guid))
        stop("some records are missing a guid")
    if (any(is.na(echino_db$guid)))
        stop("there are NAs in the guid column")
    if (any(duplicated(echino_db$guid)))
        stop("duplicated guid!")
    echino_db
}

load_cuke_db <- function(echino_db) {

    cuke_db <- subset(echino_db, class == "Holothuroidea")

    ## remove samples that don't pass the sequence checks
    cuke_db <- subset(cuke_db, ! pass.seq %in%  c("GenBank",
                                                  "fix",
                                                  "no_seq_yet",
                                                  "no",
                                                  "duplicate",
                                                  "16S_only_to_check"))
    cuke_db <- subset(cuke_db, Notes != "MH sequence")
    cuke_db <- cuke_db[nzchar(cuke_db$Sample), ]


    ## check for duplicated samples
    dup <- cuke_db[duplicated(cuke_db$Sample), "Sample"]
    if (length(dup) > 0) {
        stop("There are duplicated samples: ",
             paste(dup, collapse = ", "))
    }

    ## check pass:voucher column: no blank, no NAs, check values
    if (length(nzchar(cuke_db$Pass.voucher)) != nrow(cuke_db)) {
        stop("There are blank values in Pass.voucher column")
    }
    if (any(is.na(cuke_db$Pass.voucher))) {
        stop("There are missing vlaues in Pass.voucher column")
    }
    allowed_pass_voucher_values <- c("yes", "no", "check/no", "check", "need",
                                     "not ours", "catalog tissue?")
    if (any(!cuke_db$Pass.voucher %in% allowed_pass_voucher_values)) {
        stop("Only values allowed in pass:voucher column: ",
             paste(allowed_pass_voucher_values, collapse = ", "),
             ". Found: ",
             paste(cuke_db$Pass.voucher[!cuke_db$Pass.voucher %in%
                                        allowed_pass_voucher_values],
                   collapse = ", "))
    }

    ## no empty values in repository
    if (length(nzchar(cuke_db$Repository)) !=
        nrow(cuke_db) |
        any(is.na(cuke_db$Repository)))
        stop("Blank values in Repository column")

    ## Taxonomic check
    test_genera <- as.matrix(xtabs(~ genusorhigher + family, data = cuke_db,
                                   subset = family != "Uncertain"))
    res_genera <- apply(test_genera, 1, function(x) sum(x != 0))
    if(any(res_genera != 1)) {
        stop("These genera have more than 1 family assigned: ",
             paste(names(res_genera[res_genera != 1]), sep = ", "))
    }

    test_family <- as.matrix(xtabs(~ family + order, data = cuke_db,
                                   subset = family != "Uncertain"))
    res_family <- apply(test_family, 1, function(x) sum(x != 0))
    if (any(res_family != 1)) {
        stop("These families have more than 1 order assigned: ",
             paste(names(res_family[res_family != 1]),
                   sep = ", "))
    }

    ## only non-ambiguous bp
    #l_seq <- sapply(cuke_db$Sequence, function(x)
    #    length(gregexpr("[actgACTG]", x)[[1]]))

    ## all bp
    l_amb <- sapply(cuke_db$Sequence, function(x)
        length(gregexpr("[^-]", x)[[1]]))

    ## this also takes care of empty sequences (only -)
    cuke_db <- cuke_db[l_amb >= 500, ]


    ## These 3 sequences are not represented by other representative

    ##  it might be worth trying to figure out if we can clean up the
    ##  sequences to deal with the issues
    ##  - FRM-194
    ##  - NMV F112128
    ##  - NIWA 38032

    invisible(cuke_db)
}

load_cuke_seqs <- function(cuke_db) {
    cuke_db[, c("guid", "Sequence")]
}


load_cuke_clusters <- function(cuke_db, cuke_tree) {

    dataLbls <- character(nrow(cuke_db))
    for (i in 1:nrow(cuke_db)) {
        dataLbls[i] <- genLabel(cuke_db[i, ])
    }

    cuke_db$Labels <- clean_labels(dataLbls)

    treeTips <- data.frame(Labels_withAmb = tipLabels(cuke_tree),
                           Labels = gsub("_\\d+amb$", "", tipLabels(cuke_tree)),
                           stringsAsFactors=FALSE)

    stopifnot(all(treeTips$treeLabels %in% cuke_db$Labels))
    cuke_db_lbls <- merge(cuke_db, treeTips, by="Labels")

    invisible(cuke_db_lbls)
}
