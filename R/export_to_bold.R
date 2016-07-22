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


bold_sample_id <- function(cuke_db) {
    museum_id <- bold_museum_id(cuke_db)
    paste(museum_id, substr(cuke_db$guid, 1, 8),
          sep = "_")
}

bold_museum_id <- function(cuke_db) {
    paste(cuke_db$Repository,
          cuke_db$collection_code,
          ifelse(cuke_db$Repository == "UF",
                 paste(cuke_db$Catalog_number,
                       "Echinodermata", sep = "-"),
                 cuke_db$Catalog_number),
          sep = ":")
}

bold_voucher_table <- function(cuke_db) {

    n_row <- nrow(cuke_db)

     res <- list(
        `Sample ID` = character(n_row),
        `Field ID` = character(n_row),
        `Museum ID` = character(n_row),
        `Collection Code` = character(n_row),
        `Institution Storing` = character(n_row)
     )

    museum_id <- bold_museum_id(cuke_db)
    sample_id <- bold_sample_id(cuke_db)

    res[["Sample ID"]] <- sample_id
    res[["Museum ID"]] <- museum_id

    res[["Collection Code"]] <- cuke_db$collection_code
    res[["Institution Storing"]] <- lookup_repository_table(cuke_db$Repository)

    as.data.frame(res, stringsAsFactors = FALSE, check.names = FALSE)

}


bold_taxonomy_table <- function(cuke_db) {
     n_row <- nrow(cuke_db)

     res <- list(
         `Sample ID` = bold_sample_id(cuke_db),
         `Phylum ID` = rep("Echinodermata", n_row),
         `Class` = cuke_db$class,
         `Order` = cuke_db$order,
         `Family` = cuke_db$family,
         `Subfamily` = character(n_row),
         `Genus` = cuke_db$genusorhigher,
         `Species` = cuke_db$species
     )
     as.data.frame(res, stringsAsFactors = FALSE, check.names = FALSE)

}

bold_collection_table <- function(cuke_db) {

    n_row <- nrow(cuke_db)

    res <- list(
        `Sample ID` = bold_sample_id(cuke_db),
	`Collectors` = character(n_row),
	`Collection Date` = character(n_row),
	`Country/Ocean` = cuke_db[["Country/Ocean"]],
        `State/Province` = character(n_row),
        `Region` = character(n_row),
        `Sector` = character(n_row),
        `Exact Site` = character(n_row),
        `Latitude` = cuke_db[["decimalLatitude"]],
        `Longitude` =  cuke_db[["decimalLongitude"]],
        `Elevation` = character(n_row)
    )

    if (any(is.na(res[["Country/Ocean"]])))
        stop("Missing values in ", sQuote("Country/Ocean"), " for ",
             paste(res[["Sample ID"]][is.na(res[["Country/Ocean"]])],
                   collapse = ", \n"))

    as.data.frame(res, stringsAsFactors = FALSE, check.names = FALSE)
}


write_bold_csv <- function(res, which_table, path = ".") {
    fnm <- file.path(path, "bold_uploads",
                     paste0(format(Sys.Date(), "%Y%m%d"),
                            "-", which_table, ".csv"))
    write_csv_(res, fnm)
}

write_csv_ <- function(res, fnm) {
    res <- lapply(res, function(x) {
        x[which(is.na(x))] <- ""
        x
    })
    write.csv(data.frame(res, check.names = FALSE, stringsAsFactors = FALSE),
              row.names = FALSE,
              file = fnm
              )
    fnm
}

generate_bold_sequences <- function(cuke_db) {
    cuke_seqs <- data.frame(
        seq_ids = bold_sample_id(cuke_db),
        seqs = gsub("-", "", cuke_db$Sequence)
    )
    fnm <- file.path("bold_uploads",
                     paste0(format(Sys.Date(), "%Y%m%d"),
                            "-", "bold-sequences.fas"))
    generate_unaligned_cuke_fasta(cuke_seqs,
                                  out = fnm)
    fnm
}


generate_bold_data <- function(cuke_db) {
    if (!exists("Country/Ocean", cuke_db)) {
        stop("The cuke_db object is the target add_geodata")
    }
    stopifnot(!any(duplicated(cuke_db$Sample)))
    cuke_db <- cuke_db[cuke_db$Pass.voucher == "yes" &
                       (!is.na(cuke_db$decimalLatitude) |
                        !is.na(cuke_db$decimalLongitude)), ]


    bold_data <- lapply(list(
        bold_voucher_table,
        bold_taxonomy_table,
        bold_collection_table
    ), function(f) f(cuke_db))

    res <- mapply(function(tbl, nm) {
        write_bold_csv(tbl, nm)
    }, bold_data, c("voucher", "taxonomy", "collection"))

    seq <- generate_bold_sequences(cuke_db)
    c(unlist(res), seq)
}



subset_bold_data <- function(date, sample_to_keep) {
    csv_files <- paste0(date, "-",
                       c("collection",
                         "taxonomy",
                         "voucher"), ".csv")
    csv_files <- file.path("bold_uploads", csv_files)

    if (!all(file.exists(csv_files))) {
        stop("Some files are missing: ",
             paste(csv_files[! file.exists(csv_files)],
                   collapse = ", "))
    }

    res <- lapply(csv_files, function(x) {
        tmp_csv <- read.csv(file = x, stringsAsFactors = FALSE,
                            check.names = FALSE)
        tmp_csv[match(sample_to_keep, tmp_csv$"Sample ID"), ]
    })

    subset_csv_files <- gsub("-", "-subset-", csv_files)

    fnm <- mapply(function(x, y) {
        write_csv_(x, y)
    }, res, subset_csv_files)

    alg_in <- file.path("bold_uploads", paste0(date, "-bold-sequences.fas"))
    alg <- ape::read.dna(file = alg_in, format = "fasta", as.character = TRUE)
    alg <- alg[sample_to_keep]
    alg_sub <- gsub("([0-9]{8})-", "\\1-subset-", alg_in)
    alg <- cbind(names(alg),  sapply(alg, function(x) paste(x, collapse = "")))
    generate_unaligned_cuke_fasta(alg, out = alg_sub)
    c(fnm, alg_sub)
}
