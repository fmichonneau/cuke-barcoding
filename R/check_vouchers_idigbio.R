check_coll_code <- function(coll_code) {
     if (length(coll_code) > 1) {
        stop("For ", sQuote(inst), " more than one collection code: ",
             paste(coll_code, collapse = ", "))
    }
}

## build_idigbio_ids_flmnh and build_idigbio_ids_generic
## create a data frame from cuke_db that contains 3 columns:
## - the catalog numbers (in idigbio's format with phylum name after catalog number)
## - the latitude
## - the longitude of the specimens
## We give FLMNH specimen a special treatment that warrant their own function.
## - we check for specimens that probably have a catalog number but don't have
##   "yes" in the "pass:voucher" column
## - the catalog numbers in iDigBio are different but it's not worth changing in
##   cuke_db, so do it here. They need to have "-echinodermata" after the number
build_idigbio_ids_flmnh <- function(cuke_db) {

    uf_spcm <- cuke_db[cuke_db$Repository == "UF" &
                       cuke_db$Pass.voucher == "yes",
                       c("Repository", "Catalog_number",
                         "collection_code", "decimalLatitude",
                         "decimalLongitude", "genusorhigher",
                         "species")]

    ## check the only collection_code is "invertebrate zoology"
    coll_code <- tolower(unique(uf_spcm$collection_code))
    check_coll_code(coll_code)

    ## for FLMNH specimens show specimens that have a number but are
    ## not marked as "yes" in Pass.voucher column
    all_uf_specimens <- cuke_db[cuke_db$Repository == "UF" &
                                !cuke_db$Catalog_number %in%
                                 c("need", "no voucher", "not cataloged"), ]
    n_uf_specimens <- nrow(all_uf_specimens)
    not_counted <- n_uf_specimens - nrow(uf_spcm)

    if (not_counted > 0) {
        warning("UF specimens with Catalog Number that don't pass voucher QC: \n",
                paste(" -", setdiff(all_uf_specimens$Catalog_number,
                              uf_spcm$Catalog_number),
                      collapse = "\n"),
                call. = FALSE)
    }

    ## add phylum info to catalog number
    idig_cat_number <- paste(uf_spcm$Catalog_number, "echinodermata",
                            sep = "-")

    list(institutioncode = "uf",
         collectioncode = "invertebrate zoology",
         records =  data.frame(
             catalognumber = idig_cat_number,
             scientificname = gsub("\\s+$", "", tolower(paste(uf_spcm$genusorhigher, uf_spcm$species, sep = " "))),
             decimalLatitude = uf_spcm$decimalLatitude,
             decimalLongitude = uf_spcm$decimalLongitude,
             stringsAsFactors = FALSE
         ))
}


build_idigbio_ids_generic <- function(inst, cuke_db) {

    stopifnot(length(inst) == 1)

    spcm <- cuke_db[cuke_db$Repository == inst &
                    cuke_db$Pass.voucher == "yes",
                    c("Repository", "collection_code",
                      "Catalog_number", "decimalLatitude",
                      "decimalLongitude", "genusorhigher",
                      "species")]

    if (nrow(spcm) < 1)
        stop("Invalid institution.", sQuote(inst),
             " No results found.")

    coll_code <- unique(tolower(spcm$collection_code))
    check_coll_code(coll_code)

    list(institutioncode = tolower(inst),
         collectioncode = coll_code,
         records = data.frame(
             catalognumber = spcm$Catalog_number,
             scientificname = gsub("\\s+$", "", tolower(paste(spcm$genusorhigher, spcm$species, sep = " "))),
             decimalLatitude = spcm$decimalLatitude,
             decimalLongitude = spcm$decimalLongitude,
             stringsAsFactors = FALSE
         ))
}


## This is the main function called by remake that builds the data frames
## in the format expected by get_idigbio (below). Most collections should
## be ok with build_idigbio_ids_generic,  but FLMNH specimens require
## special treatment (see above).
build_idigbio_ids <- function(cuke_db) {
    inst <- c("CAS", "MNHN")
    res <- lapply(inst, function(x) build_idigbio_ids_generic(x, cuke_db))
    flmnh <- build_idigbio_ids_flmnh(cuke_db)
    c(list(flmnh), res)
}



## Retrieve specimens from FLMNH, CAS and MNHN in iDigBio using
## catalog numbers found in cuke_db
get_idigbio_info <- function(inst_spcm) {
    lapply(inst_spcm, function(x) {
        res <- idig_search_records(list(catalognumber = x[["records"]]$catalognumber,
                                        `data.dwc:phylum` = "echinodermata",
                                        institutioncode = x[["institutioncode"]],
                                        collectioncode = x[["collectioncode"]]))
        if (nrow(res) < 1) {
            warning("No iDigBio results for ", x[["institutioncode"]], call. = FALSE)
        }
        res
    })
}

## Find catalogued UF's specimens that are in cuke_db but not in
## iDigBio
compare_idigbio_specimens <- function(sub_cuke_db, idigbio_res) {
    res <- mapply(function(cdb, idg) {
        tt <- setdiff(cdb[["records"]]$catalognumber,
                      idg$catalognumber)
        if (length(tt) > 1) {
            msg <- paste("For", cdb[["institutioncode"]],
                         ": \n  -",
                         paste(tt, collapse = ", "))
            return(msg)
        }
    }, sub_cuke_db, idigbio_res)
    warning("These specimens are not in iDigBio: \n",
            paste(res, collapse = "\n"), call. = FALSE)

    ## compare identifications
    res <- mapply(function(cdb, idg) {
        if (nrow(idg) < 1) return(NULL)
        dplyr::left_join(cdb[["records"]], idg, by = "catalognumber") %>%
            rename_("scientificname_cukedb" = "scientificname.x",
                    "scientificname_idig" = "scientificname.y") %>%
            mutate_(.dots = setNames(list(~if_else(is.na(scientificname_idig), "", scientificname_idig)), "scientificname_idig")) %>%
            filter_("scientificname_cukedb != scientificname_idig")

    }, sub_cuke_db, idigbio_res)
    names(res) <- vapply(sub_cuke_db, function(x) x$institutioncode, character(1))
    res <- dplyr::bind_rows(res, .id = "institutioncode")

    if (nrow(res) > 0) {
        f <- "tmp/identification_mismatches.csv"
        write.csv(res[, c("institutioncode", "catalognumber",
                          "scientificname_cukedb", "scientificname_idig")],
                  file = f, row.names = FALSE)
        message(nrow(res), " do not match. Mismatches written to: ", f)
    }
    res
}


## Compare the GPS coordinates from cuke_db and iDigBio
compare_idigbio_coordinates <- function(idigbio_ids, idigbio_res) {
    res <- mapply(function(ids, res) {
        if (nrow(res) < 1) return(NULL)
        tmp_res <- int_check_idigbio_coordinates(ids[["records"]], res)
        tmp_res$institioncode <- rep(ids[["institutioncode"]], nrow(tmp_res))
        tmp_res$collectioncode <- rep(ids[["collectioncode"]], nrow(tmp_res))
        tmp_res
    }, idigbio_ids, idigbio_res, SIMPLIFY = FALSE)
    res <- do.call("rbind", res)

    summary_coordinate_comparison(res)
    res
}

is_in_idigbio <- function(idig_lat, idig_lon) {
    !(is.na(idig_lat) & is.na(idig_lon))
}

is_in_uf <- function(uf_lat, uf_lon) {
    !(is.na(uf_lat) & is.na(uf_lon))
}

check_coord_status <- function(idig_lat, idig_lon,
                               uf_lat, uf_lon, ...) {
    if (is_in_idigbio(idig_lat, idig_lon) &&
        !is_in_uf(uf_lat, uf_lon))
        res <- "missing in UF"
    else if (!is_in_idigbio(idig_lat, idig_lon) &&
             is_in_uf(uf_lat, uf_lon)) {
        res <- "missing in iDigBio"
    } else if (!is_in_idigbio(idig_lat, idig_lon) &&
               !is_in_uf(uf_lat, uf_lon)){
        res <- "missing in both"
    } else if (identical(uf_lat, idig_lat) &&
               identical(uf_lon, idig_lon))
        res <- "identical coordinates"
    else {
        res <- "different coordinates"
    }
    res
}

calc_coord_distance <- function(status, idig_lat, idig_lon,
                                 uf_lat, uf_lon, ...) {
    if (identical(status, "different coordinates")) {
        gcd.hf(deg2rad(as.numeric(uf_lon)),
               deg2rad(as.numeric(uf_lat)),
               deg2rad(as.numeric(idig_lon)),
               deg2rad(as.numeric(idig_lat)),
               warn = FALSE)
    } else NA_real_
}


int_check_idigbio_coordinates <- function(idigbio_ids, idigbio_res) {
    res <- left_join(idigbio_ids,
                     idigbio_res[, c("catalognumber", "geopoint.lat", "geopoint.lon")]
                     ) %>%
        rename(idig_lat = geopoint.lat,
               idig_lon = geopoint.lon,
               uf_lat = decimalLatitude,
               uf_lon = decimalLongitude)

    res$status <- res %>%
        purrr::pmap_chr(check_coord_status)

    res$distance <- res %>%
        purrr::pmap_dbl(calc_coord_distance)

    res
}


summary_coordinate_comparison <- function(comp) {
    n_missing <- nrow(comp[comp$status == "missing in cuke_db", ])
    test_distance <- !is.na(comp$distance) & comp$distance > 1
    n_more_than_1km <- nrow(comp[test_distance, ])

    if (n_missing > 1)
        warning("Coordinates are missing in cuke_db for ", n_missing,
                " records.", call. = FALSE)
    if (n_more_than_1km > 1)
        warning("Coordinates differ by more than 1 km for ",
                n_more_than_1km, " records:\n",
                paste0("  - ", comp$catalognumber[test_distance],
                       " (", comp$distance[test_distance], ")", collapse = "\n"),
                call. = FALSE)
}


## Compare within ESU identifications useful to make sure that
## discrepencies are real (e.g., failure for DNA barcoding to
## recognize different morphological species)
get_id_within_esu <- function(cuke_db, clusters) {
    grps <- tdata(clusters, "tip")
    grps <- split(rownames(grps), grps$Groups)

    res <- lapply(grps, function(x) {
        .r <- make_labels_from_guids(cuke_db, x, fields = c("family", "genusorhigher", "species", "modifier"))
        tibble(
            guid = x,
            ids = .r
        )
    })
    res
}

compare_id_within_esu <- function(cuke_db, clusters) {
    ids <- get_id_within_esu(cuke_db, clusters)
    with_conflicts <- vapply(ids, function(x)
        length(unique(x$ids)) > 1, logical(1))
    res <- ids[with_conflicts]
    res <- lapply(res, function(x) {
        to_add <- tibble(repository = make_labels_from_guids(cuke_db = cuke_db, guids = x$guid, fields = "Repository"),
                         catalog_number = make_labels_from_guids(cuke_db, x$guid, "Catalog_number"),
                         sample = make_labels_from_guids(cuke_db, x$guid, "Sample"))
        dplyr::bind_cols(x, to_add)
    })
    dplyr::bind_rows(res, .id = "group")
}
