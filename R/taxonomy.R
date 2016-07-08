## returns a "tidy" data frame of the classification for the families
## and orders included in the database. The column names are taxa,
## rank, higher.

load_taxonomy <- function(cuke_db) {
    uniq_order <- unique(cuke_db$order)
    uniq_family <- unique(cuke_db$family)
    uniq_family <- uniq_family[! uniq_family %in%
                               c("Dactylochirotida", "?", "Uncertain")]
    uniq_taxa <- c(uniq_order, uniq_family)
    stopifnot(! any(duplicated(uniq_taxa)))
    taxonomy_df <- data.frame(rank = c(rep("Order", length(uniq_order) + 1),
                                       rep("Family", length(uniq_family))),
                              taxa = c("all", uniq_order, uniq_family),
                              stringsAsFactors = FALSE)
    stopifnot(! any(duplicated(taxonomy_df$taxa)))
    test_family <- as.matrix(xtabs(~ family + order, data = cuke_db,
                                   subset = !family %in%
                                       c("Dactylochirotida", "?", "Uncertain")))
    which_order <- apply(test_family, 1, function(x) which(x != 0))
    tmp_family <- data.frame(taxa = names(which_order),
                             higher = dimnames(test_family)[[2]][which_order],
                             stringsAsFactors = FALSE)
    taxonomy_df <- merge(taxonomy_df, tmp_family, all.x = TRUE)
    taxonomy_df
}
