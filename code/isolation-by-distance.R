### ----- isolation-by-distance-data -----
source("R/test-allopatry-functions.R")
sppGrps <- load_species_clusterGrps(distance="raw", taxa="all", threshold=0.02)
cukeDB <- load_cukeDB()
cukeDistRaw <- load_cukeDist_raw()

nInd <- sapply(sppGrps, length)
geoDist <- lapply(sppGrps, geoDistESU, cukeDB)
genDist <- lapply(sppGrps, intraESUDist, cukeDistRaw)

families <- sapply(sppGrps, function(x) {
    fams <- sapply(strsplit(x, "_"), function(y) y[1])
    ## TODO -- if (length(fams) > 1) warning("Group", x, "more than one fam")
    fams[1]
})

species <- sapply(sppGrps, function(x) {
    spp <- sapply(strsplit(x, "_"), function(y) paste(y[2:3], collapse="_"))
    names(which.max(table(spp)))
})

maxGeoDist <- sapply(geoDist, function(x) x$max)
maxGenDist <- sapply(genDist, function(x) x$max)

distBySpecies <- data.frame(species = species,
                            maxGenDist = sapply(genDist, function(x) x$max),
                            maxGeoDist = sapply(geoDist, function(x) x$max),
                            nInd = nInd,
                            family = families)

taxonomyDf <- load_taxonomyDf()
distBySpecies <- merge(distBySpecies, taxonomyDf, by.x="family", by.y="taxa", all.x=TRUE)
distBySpecies <- subset(distBySpecies, nInd >= 3)
distBySpecies <- distBySpecies[complete.cases(distBySpecies), ]
names(distBySpecies)[match("higher", names(distBySpecies))] <- "Order"

orderToInclude <- c("Apodida", "Aspidochirotida", "Dendrochirotida")

medRgSize <- with(subset(distBySpecies, Order %in% orderToInclude),
                  tapply(maxGeoDist, Order, median, na.rm=TRUE))
medRgSizeAsp <- medRgSize["Aspidochirotida"]
medRgSizeApo <- medRgSize["Apodida"]
medRgSizeDen <- medRgSize["Dendrochirotida"]

manGrps <- load_species_manualGrps()
maxGeoDistHol <- sapply(manGrps, function(x) geoDistESU(x, cukeDB)$max)
genDistHol <- lapply(manGrps, function(x) intraESUDist(x, cukeDistRaw))
maxGenDistHol <- sapply(genDistHol, function(x) x$max)
distBySpeciesHol <- data.frame(maxGenDist=maxGenDistHol, maxGeoDist=maxGeoDistHol)
distBySpeciesHol <- distBySpeciesHol[is.finite(distBySpeciesHol$maxGenDist) &
                                     is.finite(distBySpeciesHol$maxGeoDist), ]
ibdHol <- summary(lm(maxGenDist ~ maxGeoDist, data=distBySpeciesHol))

### ---- isolation-by-distance-stats ----
leveneIbd <- car::leveneTest(lm(maxGenDist ~ Order, data=distBySpecies,
                                subset=Order %in% orderToInclude))
leveneIbdFstat <- paste("$F(", paste0(leveneIbd$Df, collapse=", "), ")=",
                        formatC(leveneIbd$"F value"[1]), "$", sep="")
leveneIbdP <- paste("$P >", formatC(leveneIbd$"Pr(>F)"[1], digits=2),
                    "$", sep="")

ibdAncova <- lm(maxGenDist ~ maxGeoDist*Order, data=distBySpecies,
                subset=Order %in% orderToInclude)

finalIbdAncova <- update(ibdAncova, . ~ . - Order - maxGeoDist:Order)

summLmIbdFinal <- summary.lm(finalIbdAncova)
summAovIbdFinal <- summary.aov(finalIbdAncova)

summLmIbd <- summary.lm(ibdAncova)
summAovIbd <- summary.aov(ibdAncova)

ancovaFstat <- paste("$F(",
                     paste(summLmIbdFinal$fstatistic[2],
                           summLmIbdFinal$fstatistic[3], sep=", "),
                     ")=",
                     formatC(summLmIbdFinal$fstatistic[1], digits=3), "$",
                     sep="")

ancovaP <- paste("$P <",
                 ifelse(pf(summLmIbdFinal$fstatistic[1],
                           summLmIbdFinal$fstatistic[2],
                           summLmIbdFinal$fstatistic[3], lower=F) < 0.001,
                        0.001,
                        pf(summLmIbdFinal$fstatistic[1],
                           summLmIbdFinal$fstatistic[2],
                           summLmIbdFinal$fstatistic[3], lower=F)),
                 "$")

interactionFstat <- paste("$F(",
                          paste(summAovIbd[[1]]$Df[3], summAovIbd[[1]]$Df[4], sep=", "),
                          ")=",
                          formatC(summAovIbd[[1]]$"F value"[3], digits=3), "$", sep="")

interactionP <- paste("$P =",
                      formatC(summAovIbd[[1]]$"Pr(>F)"[3], digits=2),
                      "$")

interceptFstat <- paste("$F(",
                        paste(summAovIbd[[1]]$Df[2], summAovIbd[[1]]$Df[4], sep=", "),
                        ")=",
                        formatC(summAovIbd[[1]]$"F value"[2], digits=3), "$", sep="")

interceptP <- paste("$P =", formatC(summAovIbd[[1]]$"Pr(>F)"[2], digits=2), "$")

ibdDendro <- lm(maxGenDist ~ maxGeoDist,
                data=distBySpecies,
                subset=Order == "Dendrochirotida")

### ---- isolation-by-distance-table ----
print(xtable(summary(finalIbdAncova), display=rep("g", 5), caption=c(paste("Coefficients of the regression",
             "between maximum genetic distances and",
              "maximum geographic distances for all ESUs identified with the clustering method",
              "with a threshold of 4\\%, represented by 3 or more individuals.",
              "See Fig.~\\ref{fig:isolation-by-distance-plot}."),
              "Coefficients of the regression between maximum genetic and maximum geographic distances."),
             label="tab:ibd-stats"),
             caption.placement="top")

### ---- isolation-by-distance-plot ----
ggplot(subset(distBySpecies, Order %in% orderToInclude),
       aes(x=maxGeoDist, y=maxGenDist, colour=Order)) +
    geom_point() + geom_abline(intercept=finalIbdAncova$coefficients[1],
                               slope=finalIbdAncova$coefficients[2], aes(colour=Order),
                               linetype=2) +
    facet_wrap(~ Order) + ylab("Maximum genetic distance (K2P)") +
    xlab("Maximum geographic distance (km)") +
    theme(legend.position="none")
