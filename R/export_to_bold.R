lookup_repository_table <- function(code) {
    tbl <- c(
        ## UF
        "UF" = "Florida Museum of Natural History",
        #"UF-Karim" = "Florida Museum of Natural History",
        #"Hickman" = "Florida Museum of Natural History",
        ## US museums
        "USNM" = "Smithsonian Institution, National Museum of Natural History",
        "LACM" = "Natural History Museum of Los Angeles County",
        "CAS" = "California Academy of Sciences",
        "BPBM" = "Bishop Museum",
        "SIO" = "Scripps Institution of Oceanography",
        ## Australia + NZ museums
        "NMV" = "Museum Victoria",
        "QM" = "Queensland Museum",
        "SAMA" = "South Australian Museum",
        "WAM" = "Western Australian Museum",
        "AM" = "Australian Museum",
        "MAGNT" = "Museum and Art Gallery of the Northern Territory",
        "NIWA" = "National Institute of Water and Atmospheric Research",
        ## European museums
        "GZG" = "Georg-August-University Göttingen Geoscience Centre",
        #"Samyn" = "Royal Museum of Natural History",
        "RBINS" = "Royal Belgian Institute of Natural Sciences",
        "MRAC" = "Royal Museum for Central Africa",
        "MNHN" = "Museum National d'Histoire Naturelle",
        "NHMUK" = "National History Museum, London",
        ## Other museums
        "UNAM" = "Universidad Nacional Autonoma de Mexico",
        "RUMF" = "University of the Ryukyus Museum"
        #"EMU" =  "EMU",
        ## Not museums
        #"Reunion" = "University of La Reunion",
        #"Masami Obuchi Collection" = "Masami Obuchi Collection",
        #"UP tissues" = "University of the Philippines",
        ## No voucher
        #"no voucher" = "no voucher",
        ## problems
        #"Werner tissues" = "Werner tissues",
        #"NCDRS_2005_01" = "NCDRS_2005_01",
        #"PHR_18" = "PHR_18",
        #"RSAKZN_2003_61" = "RSAKZN_2003_61",
        #"AAD BRC" = "AAD BRC",
        #"JBEL" = "JBEL"
    )
    res <- tbl[code]
    if (any(is.na(res)))
        stop("Repositories not found: ", paste(code[is.na(res)], collapse = ", "))
    res
}

#lookup_collection_id_table <- function()


bold_voucher_table <- function(cuke_db) {

    stopifnot(!any(duplicated(cuke_db$Sample)))
    cuke_db <- cuke_db[cuke_db$Pass.voucher == "yes", ]
    n_row <- nrow(cuke_db)

     res <- list(
        `Sample ID` = character(n_row),
        `Field ID` = character(n_row),
        `Museum ID` = character(n_row),
        `Collection Code` = character(n_row),
        `Institution Storing` = character(n_row)
     )

    museum_id <- paste(cuke_db$Repository,
                       cuke_db$collection_code,
                       ifelse(cuke_db$Repository == "UF",
                              paste(cuke_db$Catalog_number,
                                    "Echinodermata", sep = "-"),
                                 cuke_db$Catalog_number),
                              sep = ":")

    sample_id <- paste(museum_id, substr(cuke_db$guid, 1, 8),
                       sep = "_")

    res[["Sample ID"]] <- sample_id
    res[["Museum ID"]] <- museum_id

    res[["Collection Code"]] <- cuke_db$collection_code
    res[["Institution Storing"]] <- lookup_repository_table(cuke_db$Repository)

    as.data.frame(res, stringsAsFactors = FALSE)

}
