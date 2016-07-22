check_gps_coordinates <- function(cuke_db, column = "Loc") {

    globalMap <- map_data("world2")
    cuke_db <- cuke_db[!is.na(cuke_db$decimalLatitude), ]
    cuke_db <- cuke_db[!is.na(cuke_db$decimalLongitude), ]
    cuke_db <- cuke_db[nzchar(cuke_db$decimalLatitude) &
                       nzchar(cuke_db$decimalLongitude), ]

    uniq_loc <- sort(unique(cuke_db[[column]]))

    cuke_maps <- mclapply(uniq_loc, function(uloc) {
        tmp_coord <- cuke_db[cuke_db[[column]] == uloc, ]
        if (nrow(tmp_coord) == 0) {
            warning("no data for: ", uloc)
            return(NULL)
        }
        x_lims <- c(min(tmp_coord$decimalLongitude, na.rm = TRUE) - 5,
                    max(tmp_coord$decimalLongitude, na.rm = TRUE) + 5)
        y_lims <- c(min(tmp_coord$decimalLatitude, na.rm = TRUE) - 5,
                    max(tmp_coord$decimalLatitude, na.rm = TRUE) + 5)

        p <- ggplot(tmp_coord) + annotation_map(globalMap, fill="gray40", colour="gray40") +
            geom_point(aes(x=decimalLongitude, y=decimalLatitude, colour=species), size = 3, alpha = .6) +
            coord_map(projection = "mercator", orientation=c(90, 160, 0), xlim=x_lims,
                      ylim=y_lims) +
            theme(panel.background = element_rect(fill="aliceblue"),
                  legend.position = "none") +
            ggtitle(uloc) +
            xlab("Longitude") + ylab("Latitude")
    }, mc.cores = getOption("mc.cores"))

    pdf(file = "tmp/check_loc_coordinates.pdf")
    for (i in seq_along(uniq_loc)) {
        print(cuke_maps[[i]])
    }
    dev.off()
}


global_map <- function(cuke_db) {

    globalMap <- borders("world") #map_data("world2")
    cuke_db <- cuke_db[!is.na(cuke_db$decimalLatitude), ]
    cuke_db <- cuke_db[!is.na(cuke_db$decimalLongitude), ]
    cuke_db <- cuke_db[nzchar(cuke_db$decimalLatitude) &
                       nzchar(cuke_db$decimalLongitude), ]

    world <- cshapes::cshp(date = as.Date("2008-1-1"))
    world.points <- fortify(world, region = "COWCODE")

    ggplot(cuke_db) + geom_polygon(data = world.points, aes(long, lat,  group = group),
                                   fill = "gray60") +
        geom_point(aes(x=decimalLongitude, y=decimalLatitude), colour = "red",
                   position = position_jitter(2, 2),
                   size = 2, alpha = .6) +
        coord_map(projection = "mercator", orientation=c(90, 0, 0)) +
            theme(panel.background = element_rect(fill="#eaedff"),
                  legend.position = "none") +
        xlab("") + ylab("")

}
