
### Google Maps-----------------------------------------------------------------

gmap_get_countries <- function(x) {
    if (x$status == "ZERO_RESULTS") {
        c(country = NA,
          first_div = NA)
    } else if (x$status == "UNKNOWN_ERROR") {
        c(country = "unknown_error",
          first_div = "unknown_error")
    }        else {
        x <- x[["results"]][[1]]$address_components
        country_i <- sapply(x, function(i)
            any(grepl("country", i)))
        admin_i <- sapply(x, function(i)
            any(grepl("administrative_area_level_1", i)))
        country <- ifelse(any(country_i), x[country_i][[1]]$long_name, NA)
        admin <- ifelse(any(admin_i), x[admin_i][[1]]$long_name, NA)
        c(country = country,
          first_div = admin)
    }
}


gmap_store <- function(path = "data/googlemap_store",
                       hook = gmap_hook) {
    storr::storr_external(storr::driver_rds(path = path,
                                            mangle_key = TRUE),
                          hook)
}

gmap_hook <- function(key, namespace) {
    res <- httr::GET(url = "https://maps.googleapis.com/",
                     path = "maps/api/geocode/json",
                     query = list(latlng = key,
                                  key = Sys.getenv("GOOGLE_API_KEY")))
    content(res)
}

fill_gmap_store <- function(cuke_db, store = gmap_store()) {
    u_gps <- unique(paste(cuke_db$decimalLatitude,
                          cuke_db$decimalLongitude, sep = ","))
    lapply(u_gps, function(x) store$get(x))
}

gmap_country <- function(key, store = gmap_store()) {
    res <- store$get(key)
    gmap_get_countries(res)
}


### openstreet map -------------------------------------------------------------

openstreetmap_store <- function(path = "data/openstreetmap_store",
                                hook = openstreetmap_hook) {
    storr::storr_external(storr::driver_rds(path = path,
                                            mangle_key = TRUE),
                          hook)
}

openstreetmap_hook <- function(key, namespace) {
    lat_lon <- unlist(strsplit(key, ","))
    res <- httr::GET(url = "http://nominatim.openstreetmap.org",
                     path = "reverse",
                     query = list(format = "json", lat = lat_lon[1],
                                  lon = lat_lon[2], `accept-language` = "en",
                                  email = "francois.michonneau@gmail.com"))
    Sys.sleep(.7)
    if (length(content(res)) < 1) stop("something is wrong with: ", key)
    content(res)
}


fill_openstreetmap_store <- function(cuke_db, store = openstreetmap_store()) {
    u_gps <- unique(paste(cuke_db$decimalLatitude,
                          cuke_db$decimalLongitude, sep = ","))
    lapply(u_gps, function(x) store$get(x))
}

openstreetmap_country <- function(key, store = openstreetmap_store()) {
    res <- store$get(key)
    c(country = res$address$country,
      first_div = res$address$state)
}


### Marineregions --------------------------------------------------------------

## We use 2 data stores for marine regions:
## - one to store the names of the location based on the GPS coordinates
## - one to store the hierarchical geographical relations for the place names
##   based on the MRGID (marine region guids)

mregion_store <- function(path = "data/mregion_store", hook = mregions_hook) {
    storr::storr_external(storr::driver_rds(path = path,
                                            mangle_key = TRUE),
                          hook)
}

## If there are no match in a given radius, we increase the radius hoping we'll
## find a match
mregions_hook <- function(key, namespace) {
    if (key == "NA,NA") return(NA)
    lat_lon <- unlist(strsplit(key, ","))
    res <- mregions::rev_geo_code(lat = lat_lon[1], lon = lat_lon[2],
                                  lat_radius = .1, lon_radius = .1)
    if (length(res) == 0) {
        res <- mregions::rev_geo_code(lat = lat_lon[1], lon = lat_lon[2],
                                      lat_radius = .5, lon_radius = .5)
    }
    if (length(res) == 0) {
        res <- mregions::rev_geo_code(lat = lat_lon[1], lon = lat_lon[2],
                                      lat_radius = 1, lon_radius = 1)
    }
    if (length(res) == 0) {
        res <- mregions::rev_geo_code(lat = lat_lon[1], lon = lat_lon[2],
                                      lat_radius = 2, lon_radius = 2)
    }
     if (length(res) == 0) {
        res <- mregions::rev_geo_code(lat = lat_lon[1], lon = lat_lon[2],
                                      lat_radius = 3, lon_radius = 3)
    }
    res
}

## convenience function to fill the store, it's not needed in theory as
##  the calls to $get() will do the work, but it can be useful to cache
##  all the data at once instead.
fill_mregion_store <- function(cuke_db, store = mregion_store()) {
    u_gps <- unique(paste(cuke_db$decimalLatitude,
                          cuke_db$decimalLongitude, sep = ","))
    lapply(u_gps, function(x)
           store$get(x))
}

mregion_names_store <- function(path = "data/mregion_names_store",
                                hook = mregion_get_name) {
    storr::storr_external(storr::driver_rds(path = path,
                                            mangle_key = TRUE),
                          hook)
}

mregion_get_name <- function(key, namespace) {
    mregions::place_relations(key, direction = "upper", type = "partof")
}

## We go up in the hierarchy until we try to look the parent of "World"
get_higher_names_complete <- function(latlon_key, store = mregion_names_store()) {
    place <- mregion_store()$get(latlon_key)
    tmp_res <- list()
    if (length(place) > 0 && !is.na(place)) {
        curr_mrgid <- place$MRGID[1]
        while (curr_mrgid != 5218 && !is.null(curr_mrgid)) {
            curr_res <- store$get(as.character(curr_mrgid))
            tmp_res <- c(list(curr_res), tmp_res)
            curr_mrgid <- curr_res$MRGID[1]
        }
    }
    do.call("bind_rows", tmp_res)
}


### Summarize ------------------------------------------------------------------

get_external_geodata <- function(cuke_db) {
    u_gps <- unique(paste(cuke_db$decimalLatitude,
                          cuke_db$decimalLongitude, sep = ","))

    res <- lapply(u_gps, function(x) {
        gmap_res <- gmap_country(x)
        osm_res <- openstreetmap_country(x)
        data.frame(coords = x,
                   gmap_country = ifelse(is.character(gmap_res[1]), gmap_res[1], NA),
                   gmap_first_div = ifelse(is.character(gmap_res[2]), gmap_res[2], NA),
                   osm_country = ifelse(is.character(osm_res[1]), osm_res[1], NA),
                   osm_first_div = ifelse(is.character(osm_res[2]), osm_res[2], NA),
                   stringsAsFactors = FALSE
                   )
    })
    res <- do.call("bind_rows", res)
    to_add <- lapply(seq_len(nrow(res)), function(i) {
        x <- res[i, ]
        tmp_res <- data.frame(coords = x$coords,
                              mreg_sea_area = NA,
                              mreg_ocean = NA,
                              stringsAsFactors = FALSE)
        mreg <- get_higher_names_complete(x$coords)
        if (length(mreg) > 0) {
            ## First we attempt to look for ocean names
            mreg_ocean <- mreg[mreg$placeType == "Ocean", ]$preferredGazetteerName
            mreg_sea_area <- mreg[grepl("sea area", mreg$placeType, ignore.case = TRUE), ]$preferredGazetteerName

            ## but sometimes the point doesn't return an ocean, so we need to get the
            ## country information
            if (length(mreg_ocean) == 0)
                mreg_ocean <- mreg[mreg$placeType == "Nation", ]$preferredGazetteerName

            if (length(mreg_sea_area) == 0)
                mreg_sea_area <- mreg[mreg$placeType == "State", ]$preferredGazetteerName

            tmp_res <- data.frame(coords = x$coords,
                                  mreg_ocean = ifelse(length(mreg_ocean) >  0,
                                                      mreg_ocean,
                                                      NA),
                                  mreg_sea_area = ifelse(length(mreg_sea_area) > 0,
                                                         mreg_sea_area,
                                                         NA),
                                  stringsAsFactors = FALSE
                                  )
        }
        tmp_res
    })
    to_add <- do.call("bind_rows", to_add)
    left_join(res, to_add, by = "coords")
}

### Compare with cuke_db data --------------------------------------------------

get_geodata_country <- function(coords, geodata) {
    tmp_coords <- geodata[geodata$coords == coords, ]
    if (nrow(tmp_coords) == 0) stop("invalid coordinates")


    if (!is.na(tmp_coords$osm_country)) {
         if (tmp_coords$osm_country == "France" &&
                   is.na(tmp_coords$osm_first_div)) {
            res <- tmp_coords$gmap_country
        } else if (tmp_coords$osm_country == "France" &&
            (tmp_coords$gmap_country != "France" || is.na(tmp_coords$gmap_country))) {
            res <- tmp_coords$osm_first_div
        } else if (grepl("Polynésie", tmp_coords$osm_country)) {
            res <- tmp_coords$osm_first_div
        } else if (tmp_coords$osm_country == "United States of America" &&
                   tmp_coords$osm_first_div %in% c("Guam", "Hawaii",
                                                   "United States Virgin Islands",
                                                   "Northern Mariana Islands",
                                                   "United States Minor Outlying Islands")) {
            res <- tmp_coords$osm_first_div

        } else if (tmp_coords$osm_country == "Norway" &&
                   tmp_coords$osm_first_div == "Bouvet Island") {
            res <- tmp_coords$osm_first_div
        } else if (tmp_coords$osm_country == "RSA") {
            res <- "South Africa"
        } else {
            res <- gsub("\\s\\(eaux territoriales\\)", "", tmp_coords$osm_country)
        }
    } else if (!is.na(tmp_coords$gmap_country)) {
        res <- tmp_coords$gmap_country
    } else if (!is.na(tmp_coords$mreg_ocean)) {
        res <- tmp_coords$mreg_ocean
    } else {
        res <- tmp_coords$mreg_sea_area
    }

    res <- gsub("Nouvelle-Calédonie", "New Caledonia", res)
    res <- gsub("Réunion", "Reunion", res)
    res
}

add_geodata_to_cuke_db <- function(cuke_db, geodata) {
    tmp_cuke_db <- cuke_db[, c("guid", "Loc", "decimalLatitude", "decimalLongitude")]

    public_countries <- lapply(seq_len(nrow(tmp_cuke_db)), function(i) {
        x <- tmp_cuke_db[i, ]
        get_geodata_country(paste(x$decimalLatitude, x$decimalLongitude, sep = ","),
                           geodata)
    })
    public_countries <- unlist(public_countries)
    res <- data.frame(tmp_cuke_db,
                      "Country/Ocean" = public_countries,
                      stringsAsFactors=FALSE,
                      check.names = FALSE)
    left_join(cuke_db, select(res, guid, `Country/Ocean`),  by = "guid")
}
