source("R/load.R")
source("R/CalcGeoDist.R")
source("R/pairwise-groups-functions.R")

spatialFromSpecies <- function(listSpecies, cukeDB) {

    allHll <- allSpatial <- vector("list", length(listSpecies))

    species <- sapply(listSpecies, function(x) {
        sapply(strsplit(x, "_"), function(y) paste(y[2:3], collapse="_"))
    })
    nmAllHll <- sapply(species, function(x) names(which.max(table(x))))

    names(allSpatial) <- names(allHll) <- paste(names(listSpecies), nmAllHll, sep="-")

    listSpecies <- lapply(listSpecies, function(x) gsub("\\\"", "", x))

    ## Check that all species match labels in the database
    tipLbl <- cukeDB$Labels_withAmb
    chkLbl <- all(sapply(listSpecies, function(x) all(x %in% tipLbl)))
    stopifnot(chkLbl)

    ## Convex Hull for each group
    for (i in 1:length(listSpecies)) {
        tmpSpp <- listSpecies[[i]]
        tmpCoords <- cukeDB[match(tmpSpp, cukeDB$Labels_withAmb),
                            c("decimalLatitude", "decimalLongitude")]
        center <- 180
        tmpCoords$long.recenter <- ifelse(tmpCoords$decimalLongitude < center - 180,
                                        tmpCoords$decimalLongitude + 360,
                                        tmpCoords$decimalLongitude)
        tmpCoords <- tmpCoords[complete.cases(tmpCoords), ]
        hll <- try(chull(tmpCoords$long.recenter, tmpCoords$decimalLatitude),
                   silent=TRUE)
        if (inherits(hll, "try-error") ||
            length(unique(paste(tmpCoords$decimalLatitude,
                                tmpCoords$decimalLongitude))) < 3) {
            allHll[[i]] <- tmpCoords
            attr(allHll[[i]], "type-coords") <- "points"
        }
        else {
            hll <- tmpCoords[hll, ]
            allHll[[i]] <- hll
            attr(allHll[[i]], "type-coords") <- "polygon"
        }
    }

    ## Spatial format for each group
    crs <- CRS("+proj=longlat +datum=WGS84 +units=m")
    for (i in 1:length(listSpecies)) {
        if (attr(allHll[[i]], "type-coords") == "polygon") {
            tmpHll <- allHll[[i]]
            tmpPoly <- Polygon(cbind(tmpHll$long.recenter, tmpHll$decimalLatitude))
            tmpPolys <- Polygons(list(tmpPoly), names(allSpatial)[i])
            tmpSpPoly <- SpatialPolygons(list(tmpPolys), proj4string=crs)
            allSpatial[[i]] <- tmpSpPoly
        } else if (attr(allHll[[i]], "type-coords") == "points") {
            tmpCoords <- allHll[[i]][complete.cases(allHll[[i]]), ]
            if (nrow(tmpCoords) < 1) {
                allSpatial[[i]] <- NA
            }
            else {
                allSpatial[[i]] <- SpatialPoints(cbind(tmpCoords$long.recenter,
                                                       tmpCoords$decimalLatitude),
                                                 proj4string=crs)
            }
        }
    }
    list(allHll, allSpatial)
}

rangeType <- function(i, j, poly, percentOverlap=10) {
    r1 <- poly[[i]]
    r2 <- poly[[j]]

    if (inherits(r1, "SpatialPolygons") && inherits(r2, "SpatialPolygons")) {
        r1Area <- suppressWarnings(gArea(r1))
        r2Area <- suppressWarnings(gArea(r2))
        if (gIntersects(r1, r2)) {
            rIArea <- suppressWarnings(gArea(gIntersection(r1, r2)))
            if (rIArea < (min(r1Area, r2Area)/percentOverlap)) {
                res <- "parapatric"
            }
            else {
                res <- "sympatric"
            }
        }
        else {
            res <- "allopatric"
        }
    } else if ((inherits(r1, "SpatialPoints") && inherits(r2, "SpatialPolygons")) ||
               (inherits(r1, "SpatialPolygons") && inherits(r2, "SpatialPoints"))) {
        if (inherits(r1, "SpatialPoints") && nrow(r1@coords) < 2)
            res <- NA
        else if (inherits(r2, "SpatialPoints") && nrow(r2@coords) < 2)
            res <- NA
        else {
            res <- ifelse(gIntersects(r1, r2), "sympatric", "allopatric")
        }
    } else {
        res <- NA
    }

    list(rangeType=res, species=paste(names(poly)[i], names(poly)[j], sep="/"))
}

esuPairs <- function(tr) {
    grps <- tdata(tr, "tip")[, "Groups", drop=FALSE]
    sppGrps <- split(rownames(grps), grps$Groups)
    BS <- tdata(tr, "internal")[, "bs", drop=FALSE]

    mrca <- lapply(sppGrps, function(x) {
        ancestor(tr, MRCA(tr, x))
    })

    res <- lapply(mrca, function(tmpMrca) {
        tmpDesc <- descendants(tr, tmpMrca)
        tmpGrps <- unique(grps[tmpDesc, ])
        if (length(tmpGrps) > 2) {
            NA
        }
        else {
            tmpGrps
        }
    })
    toKeep1 <- !sapply(res, is.null)
    res <- res[toKeep1]
    toKeep2 <- !duplicated(res)
    res <- res[toKeep2]

    mrca <- mrca[toKeep1]
    mrca <- mrca[toKeep2]
    mrca <- unname(unlist(mrca))

    attr(res, "bootstrap") <- BS[as.character(mrca), ]
    res
}

interESUDist <- function(ind1, ind2, distMat, distance="raw") {
    ind1 <- gsub("\\\"", "", ind1)
    ind2 <- gsub("\\\"", "", ind2)

    if (length(ind1) && length(ind2)) {
        iMat <- match(c(ind1, ind2), dimnames(distMat)[[1]])
        if (any(is.na(iMat)))
            stop("problem")

        tmpDistMat <- distMat[iMat, iMat]

        tmpDistInter <- tmpDistMat[ind1, ind2]
        list(mean=mean(tmpDistInter),
             max=max(tmpDistInter),
             min=min(tmpDistInter))
    }
    else {
        list(mean=NA, max=NA, min=NA)
    }
}

intraESUDist <- function(listSpecies, alg, distance="raw") {
    listSpecies <- gsub("\\\"", "", listSpecies)
    sppLbl <- match(listSpecies, dimnames(alg)[[1]])
    if (any(is.na(sppLbl))) stop("problem with the labels")
    if (length(sppLbl) > 1) {
        tmpAlg <- alg[sppLbl, ]
        tmpDist <- ape::dist.dna(tmpAlg, model=distance)
        list(mean=mean(tmpDist), max=max(tmpDist), min=min(tmpDist))
    } else {
        list(mean=NA, max=NA, min=NA)
    }
}


geoDistESU <- function(spp, cukeDB) {
    tmpCoords <- cukeDB[match(spp, cukeDB$Labels_withAmb),
                        c("decimalLatitude", "decimalLongitude")]
    tmpCoords <- tmpCoords[complete.cases(tmpCoords), ]
    if (nrow(tmpCoords) > 1) {
        matGeoDistTmp <- CalcGeoDists(cbind(deg2rad(tmpCoords$decimalLongitude),
                                            deg2rad(tmpCoords$decimalLatitude)))
        list(mean=mean(matGeoDistTmp, na.rm=TRUE),
             max=max(matGeoDistTmp, na.rm=TRUE),
             min=min(matGeoDistTmp, na.rm=TRUE))
    } else {
        list(mean=NA, max=NA, min=NA)
    }
}

testRangeType <- function(tr, distMat, cukeDB, percentOverlap=10) {

    grps <- tdata(tr, "tip")[, "Groups", drop=FALSE]
    sppGrps <- split(rownames(grps), grps$Groups)

    polygons <- spatialFromSpecies(sppGrps, cukeDB)[[2]]

    esuPrs <- esuPairs(tr)

    rgType <- lapply(esuPrs, function(x) {
        if (! all(is.na(x)) ) {
            rangeType(x[1], x[2], polygons)
        } else {
            list(rangeType=NA, species=NA)
        }})

    interDist <- lapply(esuPrs, function(x) {
        ind1 <- sppGrps[[x[1]]]
        ind2 <- sppGrps[[x[2]]]
        interESUDist(ind1, ind2, distMat)
    })

    data.frame(species = unlist(do.call("rbind", lapply(rgType, function(x) x$species))),
               rangeType = unlist(do.call("rbind", lapply(rgType, function(x) x$rangeType))),
               meanInterDist = unlist(do.call("rbind", lapply(interDist, function(x) x$mean))),
               maxInterDist = unlist(do.call("rbind", lapply(interDist, function(x) x$max))),
               minInterDist = unlist(do.call("rbind", lapply(interDist, function(x) x$min))),
               bootstrap = attr(esuPrs, "bootstrap")
               )
}

build_species_overlap <- function() {
    ## constants for this analysis
    thres <- 0.04
    dst <- "raw"
    tax <- "all"

    taxonomyDf <- load_taxonomyDf()
    ##tax <- taxonomyDf$taxa ## need to double check but I don't think the code below works if used on other taxa

    cukeDB <- load_cukeDB()
    cukeAlg <- load_cukeAlg()

    method <- c("cluster", "pairwise")

    res <- vector("list", length(method) * length(tax))
    k <- 1

    tmpNm <- expand.grid(tax, method)
    names(res) <- apply(tmpNm, 1, function(x) paste(x[2], x[1], gsub("\\.", "", thres), sep="-"))

    for (i in 1:length(method)) {
        for (j in 1:length(tax)) {
            meth <- method[i]
            if (method[i] == "cluster") {
                tree <- load_tree_clusterGrps(distance=dst, taxa=tax[j], threshold=thres/2)
            }
            else {
                tree <- load_tree_pairwiseGrps(distance=dst, taxa=tax[j], threshold=thres)
            }
            tmpSpatial <- spatialFromSpecies(tree, cukeDB)
            tmpGeoCtxt <- testRangeType(tree, tmpSpatial[[2]], cukeDistRaw)
            res[[k]] <- tmpGeoCtxt
            k <- k+1
        }
    }
    res
}

load_species_overlap <- function(overwrite=FALSE) {
    fnm <- "data/species-overlap.rds"
    if (file.exists(fnm) && !overwrite) {
        spOver <- readRDS(file=fnm)
    }
    else {
        spOver <- build_species_overlap()
        saveRDS(spOver, file=fnm)
    }
    invisible(spOver)
}

load_species_overlap_comparison <- function(taxa="all", threshold=0.04) {
     taxonomyDf <- load_taxonomyDf()
     taxa <- match.arg(as.character(taxa), taxonomyDf$taxa)
     spOver <- load_species_overlap()
     prwseStr <- paste("pairwise", taxa, gsub("\\.", "", threshold), sep="-")
     clstrStr <- paste("cluster", taxa, gsub("\\.", "", threshold), sep="-")
     mPrwseStr <- match(prwseStr, names(spOver))
     mClstrStr <- match(clstrStr, names(spOver))
     stopifnot(length(mPrwseStr) == 1 && length(mClstrStr) == 1)
     pairwiseData <- spOver[[mPrwseStr]]
     clusterData <- spOver[[mClstrStr]]
     pairwiseData$method <- "pairwise"
     clusterData$method <- "cluster"
     rbind(pairwiseData, clusterData)
}
