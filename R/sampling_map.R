make_sampling_map <- function(cuke_db) {

    summ_gps <- data.frame(uniq=paste(cuke_db$decimalLatitude, cuke_db$decimalLongitude, sep="/"))
    tabu_gps <- table(summ_gps)
    summ_gps <- data.frame(latitude=sapply(names(tabu_gps), function(x) strsplit(x, "/")[[1]][1]),
                          longitude=sapply(names(tabu_gps), function(x) strsplit(x, "/")[[1]][2]),
                          n_ind=as.numeric(tabu_gps), row.names=1:length(tabu_gps), stringsAsFactors=FALSE)
    summ_gps <- summ_gps[-c(1, nrow(summ_gps)), ]
    summ_gps$latitude <- as.numeric(summ_gps$latitude)
    summ_gps$longitude <- as.numeric(summ_gps$longitude)

    center <- 0
    summ_gps$long.recenter <- ifelse(summ_gps$longitude < center - 180,
                                     summ_gps$longitude + 360, summ_gps$longitude)

    global_map <- ggplot2::map_data("world")

    pacificmap <- ggplot(summ_gps) +
        annotation_map(global_map, fill="gray40", colour="gray40", size = 0) +
        geom_point(aes(x = long.recenter, y = latitude, size= n_ind), colour="red", data=summ_gps) +
        coord_map(projection = "mercator", orientation=c(90, 160, 0)) +
        theme(panel.background = element_rect(fill="aliceblue"),
              legend.position="top",
              axis.title.x = element_blank(),
              axis.title.y = element_blank(),
              plot.margin=unit(rep(0, 4), "mm")) +
        scale_size_continuous(name="Number of individuals", range=c(1, 4)) +
        ylim(c(-45,45))

    southmap <- ggplot(summ_gps) +
        annotation_map(global_map, fill="gray40", colour="gray40", size = 0) +
        geom_point(aes(x = long.recenter, y = latitude, size= n_ind), colour="red", data=summ_gps) +
        coord_map(projection = "ortho", orientation=c(-90, 0, 0)) +
        theme(panel.background = element_rect(fill="aliceblue"),
              legend.position = "none",
              axis.title.x = element_blank(),
              axis.title.y = element_blank(),
              plot.margin=unit(rep(0, 4), "mm")) +
        ylim(c(-90, -45))

    northmap <- ggplot(summ_gps) +
        annotation_map(global_map, fill="gray40", colour="gray40", size = 0) +
        geom_point(aes(x = long.recenter, y = latitude, size= n_ind), colour="red", data=summ_gps) +
        coord_map(projection = "ortho", orientation=c(90, 0, 0)) +
        theme(panel.background = element_rect(fill="aliceblue"),
              legend.position="none",
              axis.title.x = element_blank(),
              axis.title.y = element_blank(),
              plot.margin=unit(rep(0, 4), "mm")) +
        ylim(c(45, 90))

    multiplot(pacificmap, northmap, southmap, layout=matrix(c(1,1,2,3), ncol=2, byrow=T))
}
