
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

mregion_store <- function(path = "data/mregion_store", hook = mregions_hook) {
    storr::storr_external(storr::driver_rds(path = path,
                                            mangle_key = TRUE),
                          hook)
}

mregions_hook <- function(key, namespace) {
    if (key == "NA,NA") return(NA)
    lat_lon <- unlist(strsplit(key, ","))
    res <- mregions::rev_geo_code(lat = lat_lon[1], lon = lat_lon[2],
                                  lat_radius = .1, lon_radius = .1)
    res
}

fill_mregion_store <- function(cuke_db, store) {
    u_gps <- unique(paste(cuke_db$decimalLatitude,
                          cuke_db$decimalLongitude, sep = ","))
    lapply(u_gps, function(x)
           store$get(x))
}

mregion_store_names <- function(path = "data/mregion_store_names",
                                hook = mregion_get_name) {
    storr::storr_external(storr::driver_rds(path = path,
                                            mangle_key = TRUE),
                          hook)
}

mregion_get_name <- function(key, namespace) {
    mregions::place_relations(key, direction = "upper", type = "partof")
}

get_higher_names_complete <- function(latlon_key, store = mregion_store_names()) {
    place <- mregion_store()$get(latlon_key)
    curr_mrgid <- place$MRGID[1]
    tmp_res <- list()
    while (curr_mrgid != 5218) {
        curr_res <- store$get(as.character(curr_mrgid))
        tmp_res <- c(list(curr_res), tmp_res)
        curr_mrgid <- curr_res$MRGID[1]
    }
    tmp_res
}


### Summarize ------------------------------------------------------------------

compare_gmap_openstreet <- function(cuke_db) {
    u_gps <- unique(paste(cuke_db$decimalLatitude,
                          cuke_db$decimalLongitude, sep = ","))

    gmap <- lapply(u_gps, function(x) gmap_country(x))
    gmap <- do.call("rbind", gmap)
    gmap <- cbind(rep("gmap", nrow(gmap)), gmap)
    osm <- lapply(u_gps, function(x) openstreetmap_country(x))
    osm <- do.call("rbind", osm)
    osm <- cbind(rep("osm", nrow(osm)), osm)
    rbind(gmap, osm)
}
