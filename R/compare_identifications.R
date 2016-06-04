compare_idig_cuke_db <- function(idig = "data/idigbio/20160603/occurrence.csv",
                                 cuke_db) {

    idig <- read.csv(file = idig, stringsAsFactors = FALSE)
    cuke_db_uf <- cuke_db[cuke_db$Repository == "UF", ]
    cuke_db_uf$idig_id <- paste(cuke_db_uf$Catalog_number, "echinodermata", sep = "-")
    res <- dplyr::full_join(cuke_db_uf, idig, by = c("idig_id" = "dwc.catalogNumber"))
    res <- res[, c("Repository", "Catalog_number", "x_genus", "genusorhigher", "species", "modifier",
                   "dwc.genus", "dwc.specificEpithet")]
    res <- res[nzchar(res$Catalog_number), ]
    res <- res[!is.na(res$Catalog_number), ]
    res <- res[tolower(paste(res$genusorhigher, res$species)) !=
               tolower(paste(res$dwc.genus, res$dwc.specificEpithet)), ]
    write.csv(res, file = "tmp/compare_idigibio_idenfications.csv", row.names = FALSE)
}
