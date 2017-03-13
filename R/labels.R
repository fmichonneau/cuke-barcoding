## remove characters that are not allowed in RAxML
## using the same pattern as in chopper::cleanSeqLabels
clean_labels <- function(lbls)  {
    gsub(":|,|\\(|\\)|;|\\[|\\]|\\'|\\s|\t", "", lbls)
}


##' Given a taxon, return all the GUIDs it includes
##'
##' @param taxonony the taxonomic database
##' @param cuke_db the full database of samples
##' @param taxa the taxon name for which the GUIDs are wanted
##' @return a vector of GUIDs
##' @export
fetch_guids_from_taxa <- function(taxonomy, cuke_db, taxa="all") {

    taxa <- match.arg(as.character(taxa), taxonomy$taxa)

    if (taxa == "all") {
        lbls <- cuke_db[, "guid"]
    } else if (taxonomy[taxonomy$taxa == taxa, "rank"] == "Order") {
        lbls <- cuke_db[cuke_db$order == taxa, "guid"]
    } else if (taxonomy[taxonomy$taxa == taxa, "rank"] == "Family") {
        lbls <- cuke_db[cuke_db$family == taxa, "guid"]
    } else {
        stop("something is wrong with ", taxa)
    }
    lbls
}

##' Convert GUIDs from alignments or tip labels in trees to something
##' intelligible
##'
##' @title Make pretty and informative labels for trees and alignments
##' @param cuke_db the database
##' @param guids the GUIDs that need to be replaced by something
##'     informative
##' @param fields a vector listing the fields in \code{cuke_db} to use
##'     to build the labels
##' @return a vector of labels
##' @author Francois Michonneau
##' @export
make_labels_from_guids <- function(cuke_db, guids,
                                    fields = c(
                                        "family",
                                        "genusorhigher",
                                        "species",
                                        "modifier",
                                        "Repository",
                                        "Catalog_number",
                                        "Sample"
                                    )) {
    mtch <- match(guids, cuke_db$guid)
    if (any(is.na(mtch))) {
        stop("guids not found in cuke_db: ", paste(guids[is.na(mtch)], collapse = ", "))
    }
    chk_col <- match(fields, names(cuke_db))
    if (any(is.na(chk_col))) {
        stop("fields not found in cuke_db: ", paste(fields[is.na(chk_col)], collapse = ", "))
    }
    dt <- cuke_db[mtch, fields, drop = FALSE]
    lbl <- apply(dt, 1, paste, collapse = "_")
    gsub("_{2,}", "_", lbl)
}
