spatialFromSpecies <- function(tree, cukeDB) {
    
    grps <- tdata(tree, "tip")[, "Groups", drop=FALSE]
    uniqGrps <- unique(grps$Groups)

    species <- vector("list", length(uniqGrps))
    for (i in 1:length(uniqGrps)) {
        tmpSpp <- rownames(subset(grps, Groups == uniqGrps[i]))
        species[[i]] <- sapply(strsplit(tmpSpp, "_"), function(x) paste(x[2:3], collapse="_"))
    }
    
    allHll <- allSpatial <- vector("list", length(uniqGrps))
    
    nmAllHll <- sapply(species, function(x) x[1])
    names(allHll) <- paste(uniqGrps, nmAllHll, sep="-")
    names(allSpatial) <- names(allHll)

    
    ## Convex Hull for each group
    for (i in 1:length(uniqGrps)) {
        tmpSpp <- rownames(subset(grps, Groups == uniqGrps[i]))
        tmpCoords <- cukeDB[match(tmpSpp, cukeDB$Labels_withAmb),
                            c("decimalLatitude", "decimalLongitude")]
        
        center <- 180
        tmpCoords$long.recenter <- ifelse(tmpCoords$decimalLongitude < center - 180,
                                        tmpCoords$decimalLongitude + 360,
                                        tmpCoords$decimalLongitude)
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
    for (i in 1:length(uniqGrps)) {
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

rangeTypePolygon <- function(i, j, poly, threshold=10) {
    r1 <- poly[[i]]
    r2 <- poly[[j]]
    r1Area <- suppressWarnings(gArea(r1))
    r2Area <- suppressWarnings(gArea(r2))
    if (gIntersects(r1, r2)) {
        rIArea <- suppressWarnings(gArea(gIntersection(r1, r2)))
        if (rIArea < (min(r1Area, r2Area)/threshold)) {
            "parapatric"
        }
        else {
            "sympatric"
        }
    }
    else {
        "allopatric"
    }
}

testRangeType <- function(tr, polygons, alg, threshold) {
    grps <- tdata(tr, "tip")[, "Groups", drop=FALSE]
    rownames(grps) <- gsub("\"$", "", rownames(grps))
    tipLabels(tr) <- gsub("\"$", "", tipLabels(tr))
    
    uniqGrps <- unique(grps$Groups)
    res <- vector("list", length(uniqGrps))
    for (i in 1:length(uniqGrps)) {
        tmpMrca <- ancestor(tr, MRCA(tr, rownames(grps)[grps$Groups == uniqGrps[i]]))
        tmpDesc <- descendants(tr, tmpMrca)
        tmpGrps <- unique(grps[tmpDesc, ])
        if (length(tmpGrps) > 2) {
            next
        }
        else {
            res[[i]] <- tmpGrps   
        }
    }
    res <- res[!sapply(res, is.null)]
    res <- res[!duplicated(res)]
    
    rgType <- lapply(res, function(x) {
        if (inherits(polygons[[x[1]]], "SpatialPolygons") &&
            inherits(polygons[[x[2]]], "SpatialPolygons")) {
            list(rangeTypePolygon(x[1], x[2], polygons, threshold=threshold),
                 species=paste(names(polygons)[x[1]],
                     names(polygons)[x[2]], sep="/"))
        } else if ((inherits(polygons[[x[1]]], "SpatialPoints") &&
                    inherits(polygons[[x[2]]], "SpatialPolygons"))) {
            isSympatric <- ifelse(gContains(polygons[[x[2]]], polygons[[x[1]]]),
                                  "sympatric", "allopatric")            
            list(isSympatric, species=paste(names(polygons)[x[1]],
                                  names(polygons)[x[2]], sep="/"))            
        } else if ((inherits(polygons[[x[1]]], "SpatialPolygons") &&
                    inherits(polygons[[x[2]]], "SpatialPoints"))) {
            isSympatric <- ifelse(gContains(polygons[[x[1]]], polygons[[x[2]]]),
                                  "sympatric", "allopatric")
            list(isSympatric, species=paste(names(polygons)[x[1]],
                                  names(polygons)[x[2]], sep="/"))       
        } else {
            list(NA, NA)
        }
    })
    
    interDist <- lapply(res, function(x) {
        ind1 <- rownames(grps)[grps$Groups == x[1]]
        ind2 <- rownames(grps)[grps$Groups == x[2]]
        if (length(ind1) && length(ind2)) {
            if (any(is.na(match(c(ind1, ind2), dimnames(alg)[[1]]))))
                stop("problem")
            tmpAlg <- alg[c(ind1, ind2), ]
            tmpDist <- dist.dna(tmpAlg, model="raw", as.matrix=TRUE)
            list(mean=mean(tmpDist[ind1, ind2]), max=max(tmpDist[ind1, ind2]),
                 min=min(tmpDist[ind1, ind2]))
        }
        else {
            NA
        }
    })
    data.frame(species = unlist(do.call("rbind", lapply(rgType, function(x) x[2]))),
               rangeType = unlist(do.call("rbind", lapply(rgType, function(x) x[1]))),
               meanInterDist = unlist(do.call("rbind", lapply(interDist, function(x) x$mean))),
               maxInterDist = unlist(do.call("rbind", lapply(interDist, function(x) x$max))),
               minInterDist = unlist(do.call("rbind", lapply(interDist, function(x) x$min))))
}
