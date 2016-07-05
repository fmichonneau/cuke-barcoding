check_coll_code <- function(coll_code) {
     if (length(coll_code) > 1) {
        stop("For ", sQuote(inst), " more than one collection code: ",
             paste(coll_code, collapse = ", "))
    }
}

## create a data frame from cuke_db that contains 3 columns:
## - the catalog numbers (in idigbio's format with phylum name after catalog number)
## - the latitude
## - the longitude of the specimens
build_idigbio_ids_flmnh <- function(cuke_db) {

    uf_spcm <- cuke_db[cuke_db$Repository == "UF" &
                       cuke_db$Pass.voucher == "yes",
                       c("Repository", "Catalog_number",
                         "collection_code",
                         "decimalLatitude", "decimalLongitude")]

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
        warning("UF specimens that don't pass voucher QC: ",
                paste(setdiff(all_uf_specimens$Catalog_number,
                              uf_spcm$Catalog_number),
                      collapse = ", "))
    }

    ## add phylum info to catalog number
    idig_cat_number <- paste(uf_spcm$Catalog_number, "echinodermata",
                             sep = "-")

    list(institutioncode = "flmnh",
         collectioncode = "invertebrate zoology",
         records =  data.frame(
             catalognumber = idig_cat_number,
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
                      "decimalLongitude")]

    if (nrow(spcm) < 1)
        stop("Invalid institution.", sQuote(inst),
             " No results found.")

    coll_code <- unique(tolower(spcm$collection_code))
    check_coll_code(coll_code)


    list(institutioncode = tolower(inst),
         collectioncode = coll_code,
         records = data.frame(
             catalognumber = spcm$Catalog_number,
             decimalLatitude = spcm$decimalLatitude,
             decimalLongitude = spcm$decimalLongitude,
             stringsAsFactors = FALSE
         ))
}



build_idigbio_ids <- function(cuke_db) {
    inst <- c("CAS")
    res <- lapply(inst, function(x) build_idigbio_ids_generic(x, cuke_db))
    c(list(build_idigbio_ids_flmnh(cuke_db)), res)
}



## Retrieve specimens from FLMNH in iDigBio using catalog numbers
## found in cuke_db
get_idigbio <- function(inst_spcm) {
    lapply(inst_spcm, function(x) {
        idig_search_records(list(catalognumber = x[["records"]]$catalognumber,
                                 institutioncode = x[["institutioncode"]],
                                 collectioncode = x[["collectioncode"]]))
    })
}

## Find catalogued UF's specimens that are in cuke_db but not in
## iDigBio
check_idigbio <- function(sub_cuke_db, idigbio_res) {
    res <- mapply(function(cdb, idg) {
        tt <- setdiff(cdb[["records"]]$catalognumber,
                      idg$catalognumber)
        if (length(tt) > 1) {
            msg <- paste("For", cdb[["institutioncode"]],
                         ": \n - ",
                         paste(tt, collapse = ", "))
            return(msg)
        }
    }, sub_cuke_db, idigbio_res)
    warning("These specimens are not in iDigBio: \n",
            paste(res, collapse = "\n"), call. = FALSE)
}

## Compare the GPS coordinates from cuke_db and iDigBio
check_idigbio_coordinates <- function(idigbio_ids, idigbio_res) {
    res <- mapply(function(ids, res) {
        tmp_res <- int_check_idigbio_coordinates(ids[["records"]], res)
        tmp_res$institioncode <- rep(ids[["institutioncode"]], nrow(tmp_res))
        tmp_res$collectioncode <- rep(ids[["collectioncode"]], nrow(tmp_res))
        tmp_res
    }, idigbio_ids, idigbio_res, SIMPLIFY = FALSE)
    do.call("rbind", res)
}

int_check_idigbio_coordinates <- function(idigbio_ids, idigbio_res) {
    res <- left_join(idigbio_ids,
                     idigbio_res[, c("catalognumber", "geopoint.lat", "geopoint.lon")])

    res$status <- apply(res, 1, function(x) {
        if ((!is.na(x[2]) & !is.na(x[3])) &
            (is.na(x[4]) & is.na(x[5])))
            return("missing in iDigBio")
        if ((is.na(x[2]) & is.na(x[3])) &
            (!is.na(x[4]) & !is.na(x[5])))
            return("missing in UF")
        if ((is.na(x[2]) & is.na(x[3])) &
            (is.na(x[4]) & is.na(x[5])))
            return("missing in both")
        if (identical(x[2], x[4]) &
            identical(x[3],  x[5]))
            return("identical coordinates")
        return("different coordinates")
    })

    res$distance <- apply(res, 1, function(x) {
        if (x[6] == "different coordinates") {
            gcd.hf(deg2rad(as.numeric(x[2])),
                   deg2rad(as.numeric(x[3])),
                   deg2rad(as.numeric(x[4])),
                   deg2rad(as.numeric(x[5])),
                   warn = FALSE)
        } else {
            NA
        }
    })

    res
}
