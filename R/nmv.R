extract_nmv_records <- function(cuke_db) {
    cuke_db[cuke_db$Repository == "NMV" &
            cuke_db$Pass.voucher == "yes", ]
}


## this function is not triggered by remake, it should be done
## manually when there are new NMV specimens that are added.
get_nmv_records <- function(cuke_db) {
    nmv_spcm <- extract_nmv_records(cuke_db)

    res <- lapply(nmv_spcm$Catalog_number, function(x) {
        httr::GET("http://collections.museumvictoria.com.au",
                  path = "api/search",
                  query = list("query" = x,
                               "specimenscientificgroup" = "invertebrate zoology"))
    })
    saveRDS(res, file = "data/nmv_data.rds")
    res
}

## Extract from the GPS coordinates and the catalog numbers from the
## GET requests
get_nmv_info <- function(nmv_check_path) {
    nmv_check <- readRDS(file = nmv_check_path)
    res <- lapply(nmv_check, function(x) {
        con <- content(x)
        if (length(con) > 0)
            extract_nmv_info(con)
        else {
            cat_num <- gsub(".+query=([a-z]{1}[0-9]+)\\&.+", "\\1",
                            x$all_headers[[1]]$headers$location)
            c(cat_num, "", "", "not in database")
        }
    })
    res <- as.data.frame(do.call("rbind", res),
                         stringsAsFactors = FALSE)
    names(res) <- c("catalognumber", "latitude", "longitude", "status")
    res$catalognumber <- gsub("\\s", "", res$catalognumber)
    res
}


extract_nmv_info <- function(nmv_content) {
    catalog_number <- nmv_content[[1]]$registrationNumber
    latitude <- nmv_content[[1]]$collectionSite$latitudes
    longitude <- nmv_content[[1]]$collectionSite$longitudes
    if (length(latitude) > 2) {
        stop("More than two latitudes?! for: ", catalog_number)
    }
    else
        ## WARNING -- When there are 2 coordinates for a station
        ## (assuming trawl), we take the average of the 2 as the point
        ## of reference
        latitude <- mean(as.numeric(unlist(latitude)))

    if (length(longitude) > 2)
        stop("More than one longitude?! for: ", catalog_number)
    else
        longitude <- mean(as.numeric(unlist(longitude)))

    res <- c(catalog_number, latitude, longitude, "in database")

    res
}

compare_nmv_specimens <- function(nmv_res) {
    nmv_missing <- nmv_res[nmv_res$status == "not in database", "catalognumber"]
    if (length(nmv_missing) > 0)
        warning("These specimens are not in the NMV database: ",
                paste(nmv_missing, collapse = ", "), call. = FALSE)
}


compare_nmv_coordinates <- function(nmv_res, cuke_db) {
    cdb_nmv <- extract_nmv_records(cuke_db)[, c("Catalog_number",
                                                "decimalLatitude",
                                                "decimalLongitude")]

    ## TODO: issue #1 -- check all the NMV records are included in
    ## nmv_res cache.

    comp_nmv <- left_join(cdb_nmv, nmv_res[, c("catalognumber", "latitude", "longitude")],
                          by = c("Catalog_number" = "catalognumber"))

    comp_nmv$status <- apply(comp_nmv, 1, function(x) {
        if ((!is.na(x[2]) & !is.na(x[3])) &
            (is.na(x[4]) & is.na(x[5])))
            return("missing in NMV database")
        if ((is.na(x[2]) & is.na(x[3])) &
            (!is.na(x[4]) & !is.na(x[5])))
            return("missing in cuke_db")
        if ((is.na(x[2]) & is.na(x[3])) &
            (is.na(x[4]) & is.na(x[5])))
            return("missing in both")
        if (identical(x[2], x[4]) &
            identical(x[3],  x[5]))
            return("identical coordinates")
        return("different coordinates")
    })

    comp_nmv$distance <- apply(comp_nmv, 1, function(x) {
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

    summary_coordinate_comparison(comp_nmv)

    comp_nmv

}
