build_idigbio_ids <- function(cuke_db) {
    uf_spcm <- cuke_db[cuke_db$Repository == "UF" &
                       cuke_db$Pass.voucher == "yes",
                       c("Repository", "Catalog_number")]

    all_uf_specimens <- cuke_db[cuke_db$Repository == "UF" &
                                !cuke_db$Catalog_number %in%
                                 c("need", "no voucher", "not cataloged"), ]
    n_uf_specimens <- nrow(all_uf_specimens)
    not_counted <- n_uf_specimens - nrow(uf_spcm)

    if (not_counted > 0) {
        warning("UF specimens that don't pass voucher QC: ",
                paste(setdiff(all_uf_specimens$Catalog_number,
                                   uf_spcm$Catalog_number),
                      collapse = ", "))
    }

    idig_cat_number <- paste(uf_spcm$Catalog_number, "echinodermata", sep = "-")
    idig_cat_number
}



get_idigbio <- function(idig_cat_numbers) {
    idig_res <- idig_search_records(list(catalognumber = idig_cat_numbers,
                                         institutioncode = "flmnh"))

    idig_res
}

check_idigbio <- function(idig_cat_numbers, idigbio_res) {
    res <- setdiff(idig_cat_numbers, idigbio_res$catalognumber)
    warning("These specimens are not in iDigBio: ",
            paste(res, collapse = ", "))
}
