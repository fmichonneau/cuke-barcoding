## fetch the labels found in the trees for a given taxon
fetch_labels_from_taxa <- function(taxonomy, cuke_db, taxa="all") {

    taxa <- match.arg(as.character(taxa), taxonomy$taxa)

    if (taxa == "all") {
        lbls <- cuke_db[, "Sample"]
    } else if (taxonomy[taxonomy$taxa == taxa, "rank"] == "Order") {
        lbls <- cuke_db[cuke_db$order == taxa, "Sample"]
    } else if (taxonomy[taxonomy$taxa == taxa, "rank"] == "Family") {
        lbls <- cuke_db[cuke_db$family == taxa, "Sample"]
    } else {
        stop("something is wrong with ", taxa)
    }
    clean_labels(lbls)
}

## labels from sample numbers
make_labels_from_sample <- function(cuke_db, samples,
                                    fields = c("family", "genusorhigher",
                                               "species", "Repository",
                                               "Catalog_number", "Sample")) {
    dt <- cuke_db[match(samples, cuke_db$Sample), fields]
    lbl <- apply(dt, 1, paste, collapse = "_")
    gsub("_{2,}", "", lbl)
}
