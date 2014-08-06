
### ---- find-group-full-tree ----
## library(ape)
## library(phylobase)
## library(doMC)
## registerDoMC()
## source("R/findGroups.R")
## phylobase.options(allow.duplicated.labels="ok")

## ## based on tree calculated from raw distances
## treeH <- readRDS(file="data/cukeTree-raw.rds")

## treeH$edge.length[treeH$edge.length < 0] <- 1e-6
## treeHr <- ape::root(treeH, 1, resolve.root=TRUE)
## treeH4 <- as(treeHr, "phylo4")
## bs <- nodeLabels(treeH4)
## is.na(bs[as.character(rootNode(treeH4))]) <- TRUE
## bs <- data.frame(bs, stringsAsFactors=FALSE)
## bs$bs <- as.numeric(bs$bs)
## treeH4 <- phylo4d(treeH4, node.data=bs)
## treeH4 <- removeNodeLabels(treeH4)

## treeHgr <- findGroups(treeH4, experimental=FALSE, parallel=TRUE)

## saveRDS(treeHgr, file="data/cukeTree-raw-withGroups.rds")

## ## based on tree calculated with K2P distances
## treeHk2p <- readRDS(file="data/cukeTree-k2p.rds")

## treeHk2p$edge.length[treeHk2p$edge.length < 0] <- 1e-6
## treeHk2pr <- ape::root(treeHk2p, 1, resolve.root=TRUE)
## treeHk2p4 <- as(treeHk2pr, "phylo4")
## bs <- nodeLabels(treeHk2p4)
## is.na(bs[as.character(rootNode(treeHk2p4))]) <- TRUE
## bs <- data.frame(bs, stringsAsFactors=FALSE)
## bs$bs <- as.numeric(bs$bs)
## treeHk2p4 <- phylo4d(treeHk2p4, node.data=bs)
## treeHk2p4 <- removeNodeLabels(treeHk2p4)

## treHk2pgr <- findGroups(treeHk2p4, experimental=FALSE, parallel=TRUE)

## saveRDS(treHgr, file="data/cukeTree-k2p-withGroups.rds")

### 
source("make/build_cukeTree-withGroups.R")
source("make/build_cukeTree-byTaxa-withGroups.R")

### ---- summary-bPTP-results ----
resPTP <- readLines("data/bPTP/bPTP_fullTree.PTPPartitonSummary.txt")

getSpp <- lapply(resPTP, function(x) {
    pp <- strsplit(x, ": ")[[1]][2]
    x <- gsub(":\\s?\\d\\.\\d+$", "", x)
    x <- gsub("^\\s?", "", x)
    list(species = unlist(strsplit(x, ", ")),
         posterior = as.numeric(pp))
})

thres <- 9:0/10
res <- vector("list", length(thres))
for (i in 1:length(thres)) {
    subPP <- sapply(getSpp, function(x) x$posterior > thres[i])
    sppPP <- sapply(getSpp[subPP], function(x) x$species)
    res[[i]]$nClusters <- length(sppPP)
    res[[i]]$nSpecies <- length(unlist(sppPP))
    res[[i]]$anyDup <- any(duplicated(unlist(sppPP)))
    res[[i]]$species <- unlist(sppPP)
}
    
### ---- add-labels-to-database ----
treeH <- readRDS("data/cukeTree-raw-withGroups.rds")

source("R/genFasta.R")
echinoDB <- readRDS("data/raw/cukeBarcodes.csv.rds")
cukeDB <- subset(echinoDB, class_ == "Holothuroidea")

dataLbls <- character(nrow(cukeDB))
for (i in 1:nrow(cukeDB)) {
    dataLbls[i] <- genLabel(cukeDB[i, ])
}

## using the same pattern as in seqManagement::cleanSeqLabels
cukeDB$Labels <- gsub(":|,|\\(|\\)|;|\\[|\\]|\\'|\\s|\t", "", dataLbls)

treeTips <- data.frame(Labels_withAmb = tipLabels(treeH),
                       Labels = gsub("_\\d+amb$", "", tipLabels(treeH)),
                       stringsAsFactors=FALSE)

stopifnot(all(treeTips$treeLabels %in% cukeDB$Labels))

cukeDB_lbls <- merge(cukeDB, treeTips, by="Labels")
saveRDS(cukeDB_lbls, "data/cukeDB_withLabels.rds")

### ---- test-thresholds----
library(ggplot2)
taxonomyDf <- readRDS(file="data/taxonomyDf.rds")
treeGrpsFiles <- list.files(pattern="cukeTree-.+-\\d+\\.rds$",
                            path="data", full.names=TRUE)

treeGrps <- lapply(treeGrpsFiles, readRDS)
nGrpsVec <- sapply(treeGrps, function(tr) max(tdata(tr, "tip")$Groups))
pSnglVec <- sapply(treeGrps, function(tr) {
    tabGrp <- table(tdata(tr, "tip")$Groups) == 1
    sum(tabGrp)/length(tabGrp)
})

factStr <- lapply(treeGrpsFiles, function(x) {
    tmpStr <- gsub("^data/", "", x)
    tmpStr <- gsub("\\.rds$", "", tmpStr)
    ## to fix the ones done on all species
    tmpStr <- gsub("(cukeTree-(raw|k2p)-)(\\d+)", "\\1all-\\3", tmpStr)
    tmpStr <- strsplit(tmpStr, "-")
    data.frame(distance=tmpStr[[1]][2], taxa=tmpStr[[1]][3],
               threshold=tmpStr[[1]][4], stringsAsFactors=FALSE)
})
factDf <- do.call("rbind", factStr)
factDf$threshold <- as.numeric(gsub("^0", "0.", factDf$threshold))
factDf$threshold <- factDf$threshold*2

nGrpsClustersDf <- cbind(factDf, nGrps=nGrpsVec, pSngl=nSnglVec)
nGrpsClustersDf <- merge(nGrpsClustersDf, taxonomyDf, all.x=TRUE)
nGrpsClustersDf$taxa <- factor(nGrpsClustersDf$taxa)
nGrpsClustersDf$distance <- factor(nGrpsClustersDf$distance)
nGrpsClustersDf$method <- "Clusters"


### ---- pairwise-groups ----
library(ape)
library(igraph)

pairwiseGrps <- function(d, threshold) {
    iGrp <- lapply(1:ncol(d), function(j) dimnames(d)[[1]][d[, j] < threshold])
    sngl <- sapply(iGrp, function(x) length(x) == 1)
    edg <- do.call("rbind", lapply(iGrp, function(x) {
        if (length(x) > 1) cbind(head(x, -1), tail(x, -1)) else NULL
    }))
    g <- graph.data.frame(edg)
    c(split(V(g)$name, clusters(g)$membership), iGrp[sngl])
}

getPairwiseGrp <- function(alg=cukeAlg, model=c("raw", "K80"),
                           db=cukeDB,
                           taxonomySubset, threshold=thresVec,
                           taxonomyData=taxonomyDf) {
    model <- match.arg(model)
    if (missing(taxonomySubset) || taxonomySubset == "all")
        distAlg <- ape::dist.dna(alg, model=model, as.matrix=TRUE)
    else {
        taxLvl <- taxonomyData[taxonomyData$taxa == taxonomySubset, "rank"]
        if (length(taxLvl) > 1)
            stop("Something is wrong with this taxonomic name.")
        if (taxLvl == "Order") {
            distLbl <- subset(db, order == taxonomySubset)$Labels_withAmb

        } else if (taxLvl == "Family") {
            distLbl <- subset(db, family == taxonomySubset)$Labels_withAmb

        }
        else stop("Houston, we have a problem.")
        distLbl <- gsub("\\\"", "", distLbl)
        indexLbl <- match(distLbl, dimnames(alg)[[1]])
        stopifnot(! any(is.na(indexLbl)))
        stopifnot(length(indexLbl) > 3)
        distAlg <- dist.dna(alg[indexLbl, ], model=model, as.matrix=TRUE)
    }
    lapply(threshold, function(thres) pairwiseGrps(distAlg, thres))
}

getPropSngl <- function(pairwiseGrp) {
    sapply(pairwiseGrpRawAll, function(x) {
        tt <- sapply(x, length)
        sum(tt == 1)/length(tt)
    })
}

thresVec <- c(seq(1, 5, by=.5), 6:8)/100
cukeDB <- readRDS(file="data/cukeDB_withLabels.rds")

cukeAlg <- ape::read.dna(file="data/cukeBarcodes-flagAmb.phy", format="seq")
dimnames(cukeAlg)[[1]] <- gsub("\\\"", "", dimnames(cukeAlg)[[1]])

pairwiseGrpMdl <- c("raw", "K80")
pairwiseGrpTax <- c("all", "Aspidochirotida", "Holothuriidae", "Apodida", "Dendrochirotida")
pairwiseGrpFactors <- expand.grid(threshold=thresVec, taxa=pairwiseGrpTax, distance=pairwiseGrpMdl)
pairwiseGrpFactors <- merge(pairwiseGrpFactors, taxonomyDf)

getPairwiseGrpRes <- function() {
    pairwiseGrpRes <- vector("list", length(pairwiseGrpMdl) * length(pairwiseGrpTax))
    nmGrpRes <- character(length(pairwiseGrpMdl) * length(pairwiseGrpTax))
    i <- 1
    for (eachMdl in 1:length(pairwiseGrpMdl)) {
        for (eachTax in 1:length(pairwiseGrpTax)) {
            pairwiseGrpRes[[i]] <- getPairwiseGrp(model=pairwiseGrpMdl[eachMdl],
                                                  taxonomySubset=pairwiseGrpTax[eachTax])
            nmGrpRes[i] <- paste(pairwiseGrpMdl[eachMdl], pairwiseGrpTax[eachTax], sep="-")
            i <- i + 1
        }
    }
    names(pairwiseGrpRes) <- nmGrpRes
    pairwiseGrpRes
}

pairwiseGrpRes <- getPairwiseGrpRes()
saveRDS(pairwiseGrpRes, file="data/pairwiseGrpRes.rds")

nGrpsPairwise <- lapply(pairwiseGrpRes, function(x) sapply(x, length))
pSnglPairwise <- lapply(pairwiseGrpRes, function(x) sapply(x, getPropSngl))

nGrpsPairwiseDf <- data.frame(pairwiseGrpFactors, nGrps=unlist(nGrpsPairwise),
                              pSngl=unlist(pSnglPairwise))
nGrpsPairwiseDf$method <- "Pairwise"

nGrpsDf <- rbind(nGrpsClustersDf, nGrpsPairwiseDf)


ggplot(nGrpsDf,
              aes(x=threshold, y=nGrps, colour=interaction(distance, method), shape=distance)) + geom_line() +
    geom_point() + facet_wrap( ~ taxa)

ggplot(subset(nGrpsDf, taxa %in% c("all", "Holothuriidae", "Aspidochirotida")),
              aes(x=threshold, y=pSngl, colour=interaction(distance, method), shape=distance)) + geom_line() +
    geom_point() + facet_wrap( ~ taxa)


### ----- compare-genetic-geo-distances ----
source("R/CalcGeoDist.R")
treeH <- readRDS("data/cukeTree-raw-withGroups.rds")
cukeDB <- readRDS("data/cukeDB_withLabels.rds")
grps <- tdata(treeH, "tip")[, "Groups", drop=FALSE]
cukeAlg <- ape::read.dna(file="data/cukeBarcodes-flagAmb.phy", format="sequential")

uniqGrps <- unique(grps$Groups)

matGeoDist <- vector("list", length(uniqGrps))
latCat <- character(length(uniqGrps))
nInd <- numeric(length(uniqGrps))
matGenDist <- vector("list", length(uniqGrps))
families <- character(length(uniqGrps))
species <- vector("list", length(uniqGrps))
listCoords <- vector("list", length(uniqGrps))
for (i in 1:length(uniqGrps)) {
    tmpSpp <- rownames(subset(grps, Groups == uniqGrps[i]))
    tmpCoords <- cukeDB[match(tmpSpp, cukeDB$Labels_withAmb),
                        c("decimalLatitude", "decimalLongitude")]
    listCoords[[i]] <- tmpCoords
    matGeoDistTmp <- CalcGeoDists(cbind(deg2rad(tmpCoords$decimalLongitude),
                                        deg2rad(tmpCoords$decimalLatitude)))
    latCat[i] <- ifelse(mean(abs(tmpCoords$decimalLatitude)) < 22, "tropical", "polar")
    nInd[i] <- nrow(tmpCoords)
    matGeoDist[[i]] <- matGeoDistTmp
    tmpAlg <- cukeAlg[match(tmpSpp, dimnames(cukeAlg)[[1]]), ]
    matGenDist[[i]] <- dist.dna(tmpAlg, model="raw")
    tmpFamilies <- sapply(strsplit(tmpSpp, "_"), function(x) x[1])
    species[[i]] <- sapply(strsplit(tmpSpp, "_"), function(x) paste(x[2:3], collapse="_"))
    #stopifnot(length(unique(tmpFamilies)) < 2)
    if (length(unique(tmpFamilies)) >= 2 ) {
        warning("Group ", i, " had more than 1 family.")
        tmpFamilies <- NA
    }
    families[i] <- unique(tmpFamilies)
}

meanGeoDist <- lapply(matGeoDist, mean)
maxGeoDist <- lapply(matGeoDist, max)
fullGeoDist <- unlist(lapply(matGeoDist, as.numeric))
meanGenDist <- lapply(matGenDist, mean)
fullGenDist <- unlist(lapply(matGenDist, as.numeric))

dat <- data.frame(genDist = unlist(meanGenDist), geoDist = unlist(meanGeoDist),
                  maxGeoDist = unlist(maxGeoDist), family = families,
                 latCat = latCat, nInd = nInd)

plot(genDist ~ maxGeoDist, data=dat, subset = nInd >= 5 & geoDist > 100)
summary(lm(genDist ~ geoDist, data=dat, subset=nInd >= 5 & geoDist > 100))

ggplot(data=dat, aes(x=geoDist, y=genDist, colour=latCat, size=nInd)) + geom_point()

nGeoData <- tapply(dat$geoDist, dat$family, function(x) sum(!is.na(x)))
keepFamilies <- names(nGeoData)[nGeoData >= 9]

ggplot(subset(dat, nInd >= 3 & family %in% keepFamilies & latCat == "tropical")) +
    geom_point(aes(x=family, y=geoDist, colour=family), position=position_jitter(width=.05))



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

testRangeType <- function(tr, poly, alg, threshold) {
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
    rType <- lapply(res, function(x) {
        #iPoly1 <- grep(paste0("^", x[1], "-"), names(poly))
        #iPoly2 <- grep(paste0("^", x[2], "-"), names(poly))
        if (inherits(allSpatial[[x[1]]], "SpatialPolygons") &&
            inherits(allSpatial[[x[2]]], "SpatialPolygons")) {
            list(rangeTypePolygon(x[1], x[2], poly, threshold=threshold),
                 species=paste(names(poly)[x[1]],
                     names(poly)[x[2]], sep="/"))
        } else if ((inherits(allSpatial[[x[1]]], "SpatialPoints") &&
                    inherits(allSpatial[[x[2]]], "SpatialPolygons"))) {
            isSympatric <- ifelse(gContains(allSpatial[[x[2]]], allSpatial[[x[1]]]),
                                  "sympatric", "allopatric")            
            list(isSympatric, species=paste(names(poly)[x[1]],
                                  names(poly)[x[2]], sep="/"))            
        } else if ((inherits(allSpatial[[x[1]]], "SpatialPolygons") &&
                    inherits(allSpatial[[x[2]]], "SpatialPoints"))) {
            isSympatric <- ifelse(gContains(allSpatial[[x[1]]], allSpatial[[x[2]]]),
                                  "sympatric", "allopatric")
            list(isSympatric, species=paste(names(poly)[x[1]],
                                  names(poly)[x[2]], sep="/"))       
        } else {
            list(NA, NA)
        }
    })
    interDist <- lapply(res, function(x) {
        ind1 <- rownames(grps)[grps$Groups == x[1]]
        ind2 <- rownames(grps)[grps$Groups == x[2]]
        if (length(ind1) && length(ind2)) {
            if (any(is.na(match(c(ind1, ind2), dimnames(alg)[[1]]))))
                browser()
            tmpAlg <- alg[c(ind1, ind2), ]
            tmpDist <- dist.dna(tmpAlg, model="raw", as.matrix=TRUE)
            list(mean=mean(tmpDist[ind1, ind2]), max=max(tmpDist[ind1, ind2]),
                 min=min(tmpDist[ind1, ind2]))
        }
        else {
            NA
        }
    })
    data.frame(species = unlist(do.call("rbind", lapply(rType, function(x) x[2]))),
               rangeType = unlist(do.call("rbind", lapply(rType, function(x) x[1]))),
               meanInterDist = unlist(do.call("rbind", lapply(interDist, function(x) x$mean))),
               maxInterDist = unlist(do.call("rbind", lapply(interDist, function(x) x$max))),
               minInterDist = unlist(do.call("rbind", lapply(interDist, function(x) x$min))))
}

spatialTreeH <- spatialFromSpecies(treeH, cukeDB)

actTree <- subset(treeH, tips.include=grep("Actinopyga", tipLabels(treeH)))
bohTree <- subset(treeH, tips.include=grep("Bohadschia", tipLabels(treeH)))

testRangeType(actTree, spatialTreeH[[2]], cukeAlg, 10)
testRangeType(bohTree, allSpatial, cukeAlg, 10)

geoSpe <- testRangeType(treeH, spatialTreeH[[2]], cukeAlg, 10)

ggplot(geoSpe) +
    geom_boxplot(data=(geoSpe[!is.na(geoSpe$rangeType),]),
                 aes(x=rangeType, y=minInterDist))
    
       geom_bar(data=(geoSpe[!is.na(geoSpe$rangeType),]),
                aes(x=meanInterDist, fill=rangeType)) #+ geom_bar(position="fill")


library(maps)
library(ggplot2)

globalMap <- map_data("world2")

pdf(file="distributionMaps.pdf", paper="letter")
for (i in 1:length(uniqGrps)) {
    tmpCoords <- listCoords[[i]]
    if (nrow(tmpData) < 4) next
    if (length(unique(paste(tmpData$decimalLatitude, tmpData$decimalLongitude))) < 4) next
    center <- 150
    tmpData$long.recenter <- ifelse(tmpData$decimalLongitude < center - 180,
                                    tmpData$decimalLongitude + 360,
                                    tmpData$decimalLongitude)
    hll <- try(chull(tmpData$long.recenter, tmpData$decimalLatitude), silent=FALSE)
    if (inherits(hll, "try-error")) next
    hll <- tmpData[hll, ]
    xx <- ggplot(tmpData) + annotation_map(globalMap, fill="gray40", colour="gray40") +
        geom_point(aes(x = long.recenter, y = decimalLatitude), data=tmpData) +
        geom_polygon(data=hll, aes(x=long.recenter, y=decimalLatitude), alpha=.2) + 
        coord_map(projection = "mercator", orientation=c(90, 160, 0)) +
        theme(panel.background = element_rect(fill="aliceblue")) +
         ggtitle(paste(i, species[[i]][1]))
        ##ylim(c(-45,45))
    print(xx)
}
dev.off()



tmpHll <- allHll[sapply(allHll, function(x) attr(x, "type-coords") == "polygon")]
tmpSpp <- rep(names(tmpHll), sapply(tmpHll, function(x) nrow(x[complete.cases(x), ])))
allHllDf <- do.call("rbind", tmpHll)
allHllDf <- cbind(allHllDf, species=tmpSpp)

ggplot(allHllDf) + annotation_map(globalMap, fill="gray40", colour="gray40") +
    ##geom_point(aes(x = long.recenter, y = decimalLatitude), data=tmpData) +
   ## geom_polygon(data=(allHllDf[grep("^21-|^9-", allHllDf$species), ]),
    ##             aes(x=long.recenter, y=decimalLatitude, colour=species), fill="blue",
      ##           linetype=1, size=.01, alpha=.05) +
    geom_point(data=(allHllDf[grep("^21-|^9-", allHllDf$species), ]),
                aes(x=long.recenter, y=decimalLatitude, colour=species)) +
    coord_map(projection = "mercator", orientation=c(90, 160, 0)) +
    theme(panel.background = element_rect(fill="aliceblue")) +   
    ylim(c(-30,30)) + xlim(c(0, 300))


### ---------  Make trees for each family
## Get the families

getFam <- sapply(dimnames(seqHol)[[1]], function(x) { unlist(strsplit(x, "_"))[1] })
uniqFam <- unique(getFam)
missing <- dimnames(seqHol)[[1]][which(getFam == "")]
stopifnot(length(missing) == 0)

treeForEachFamily <- function(uniqFam, alg, drawTrees=TRUE) {
    res <- vector("list", length(uniqFam))
    for (i in 1:length(uniqFam)) {
        fam <- uniqFam[i]
        message(fam)
        if (nchar(fam) == 0 || fam == "?") next
        selSeqI <- grep(paste("^", fam, "_", sep=""), dimnames(seqHol)[[1]])
        selSeq <- alg[selSeqI, ]
        if (nrow(selSeq) < 3) {
            message("not enough sequences for ", fam, " to make a tree.")
            next
        }
        dimnames(selSeq)[[1]] <- gsub(paste("^", fam, "_", sep=""), "", dimnames(selSeq)[[1]])
        treTmp <- nj(dist.dna(selSeq))
        res[[i]] <- treTmp
        if (drawTrees) {
            h <- (dim(selSeq)[1]/10) + 5
            pdf(file=paste(fam, ".pdf", sep=""), height=h)
            plot(ladderize(treTmp), no.margin=TRUE, cex=.7)
            dev.off()
        }
        message("Done.")
    }
    res
}

treeForEachFamily(uniqFam, seqHol)

### make findGroup more efficient

synTr <- treeForEachFamily("Synaptidae", seqHol, drawTrees=FALSE)[[1]]
synTr$edge.length[synTr$edge.length < 0] <- 1e-7
synTr$tip.label <- make.unique(synTr$tip.label)
synTr <- root(synTr, outgroup=grep("Leptosynapta|Patinapta", synTr$tip.label), resolve.root=TRUE)
synTr4 <- as(synTr, "phylo4")

system.time(synGr <- findGroups(synTr4, experimental=FALSE, parallel=TRUE))

library(microbenchmark)
microbenchmark(
    sapply(nodeId(synTr4, "internal"), function(x)
           distFromTip(synTr4, x, parallel=FALSE)),
    sapply((nTips(synTr4)+1):(nEdges(synTr4)), function(x)
           distFromTip(synTr4, x, parallel=FALSE)),
    times=50L
    )

microbenchmark(
    nodeId(synTr4, "all"),
    nodeId(synTr4, "all2"))

microbenchmark(
    findGroups(synTr4, experimental=FALSE, parallel=FALSE),
    findGroups(synTr4, experimental=TRUE, parallel=FALSE))

library(lineprof)
lp <- lineprof(findGroups(synTr4, experimental=FALSE, parallel=FALSE))

### Summary coords
library(maps)
library(ggplot2)

summGPS <- data.frame(uniq=paste(holDB$decimalLatitude, holDB$decimalLongitude, sep="/"))
tabuGPS <- table(summGPS)
summGPS <- data.frame(latitude=sapply(names(tabuGPS), function(x) strsplit(x, "/")[[1]][1]),
                      longitude=sapply(names(tabuGPS), function(x) strsplit(x, "/")[[1]][2]),
                      nInd=as.numeric(tabuGPS), row.names=1:length(tabuGPS), stringsAsFactors=FALSE)
summGPS <- summGPS[-c(1, nrow(summGPS)), ]
summGPS$latitude <- as.numeric(summGPS$latitude)
summGPS$longitude <- as.numeric(summGPS$longitude)

center <- 0
summGPS$long.recenter <- ifelse(summGPS$longitude < center - 180, summGPS$longitude + 360, summGPS$longitude)
globalMap <- map_data("world")

pacificmap <- ggplot(summGPS) + annotation_map(globalMap, fill="gray70", colour="gray70") +
    geom_point(aes(x = long.recenter, y = latitude, size= nInd), colour="red", data=summGPS) +
    coord_map(projection = "mercator", orientation=c(90, 160, 0)) +
    theme(panel.background = element_rect(fill="aliceblue")) +
    ylim(c(-45,45))
pacificmap

southmap <- ggplot(summGPS) + annotation_map(globalMap, fill="gray70", colour="gray70") +
    geom_point(aes(x = long.recenter, y = latitude, size= nInd), colour="red", data=summGPS) +
    coord_map(projection = "ortho", orientation=c(-90, 0, 0)) +
    theme(panel.background = element_rect(fill="aliceblue"))
southmap

northmap <- ggplot(summGPS) + annotation_map(globalMap, fill="gray70", colour="gray70") +
    geom_point(aes(x = long.recenter, y = latitude, size= nInd), colour="red", data=summGPS) +
    coord_map(projection = "ortho", orientation=c(90, 0, 0)) +
    theme(panel.background = element_rect(fill="aliceblue"))
northmap

pdf(file="cukebarcodingmaps.pdf", paper="USr", width=0, height=0)
print(pacificmap)
print(southmap)
print(northmap)
dev.off()



################# below is old code

treH$edge.length[treH$edge.length < 0] <- 0
bootH <- boot.phylo(treH, seqH, function(xx) nj(dist.dna(xx)), B=100)
treH$node.label <- bootH
treHr <- root(treH, grep("^Pseudostichopus", treH$tip.label), resolve.root=TRUE)

treH4 <- phylo4d(treHr, check.node.label="asdata")

grpH1 <- findGroups(treH4, threshold=.005)
grpH2 <- findGroups(treH4, threshold=.010)
grpH3 <- findGroups(treH4, threshold=.015)
grpH4 <- findGroups(treH4, threshold=.020)
grpH5 <- findGroups(treH4, threshold=.025)
grpH6 <- findGroups(treH4, threshold=.030)









################## maybe for a later time, if doing analyses by family/orders.


## Dendro
dbDen <- subset(dbH, order_ == "Dendrochirotida")
lblDen <- sapply(1:nrow(dbDen), function(i) genLabel(dbDen[i, ]))
seqDen <- allEch[dimnames(allEch)[[1]] %in% lblDen, ]
dimnames(seqDen)[[1]] <- make.unique(dimnames(seqDen)[[1]])
treDen <- nj(dist.dna(seqDen))
treDen$edge.length[treDen$edge.length < 0] <- 0
bootDen <- boot.phylo(treDen, seqDen, function(xx) nj(dist.dna(xx)), B=100)
treDen$node.label <- bootDen
treDenr <- root(treDen, grep("^Lissothuria", treDen$tip.label), resolve.root=TRUE)
treDen4 <- phylo4d(treDenr, check.node.label="asdata")

grpDen1 <- findGroups(treDen4, threshold=.005)
grpDen2 <- findGroups(treDen4, threshold=.010)
grpDen3 <- findGroups(treDen4, threshold=.015)
grpDen4 <- findGroups(treDen4, threshold=.020)
grpDen5 <- findGroups(treDen4, threshold=.025)
grpDen6 <- findGroups(treDen4, threshold=.030)

## Apodid
dbApo <- subset(dbH, order_ == "Apodida")
lblApo <- sapply(1:nrow(dbApo), function(i) genLabel(dbApo[i, ]))
seqApo <- allEch[dimnames(allEch)[[1]] %in% lblApo, ]
dimnames(seqApo)[[1]] <- make.unique(dimnames(seqApo)[[1]])
treApo <- nj(dist.dna(seqApo))
treApo$edge.length[treApo$edge.length < 0] <- 0
bootApo <- boot.phylo(treApo, seqApo, function(xx) nj(dist.dna(xx)), B=100)
treApo$node.label <- bootApo
treApor <- root(treApo, grep("^Polycheira", treApo$tip.label), resolve.root=TRUE)
treApo4 <- phylo4d(treApor, check.node.label="asdata")

grpApo1 <- findGroups(treApo4, threshold=.005)
grpApo2 <- findGroups(treApo4, threshold=.010)
grpApo3 <- findGroups(treApo4, threshold=.015)
grpApo4 <- findGroups(treApo4, threshold=.020)
grpApo5 <- findGroups(treApo4, threshold=.025)
grpApo6 <- findGroups(treApo4, threshold=.030)

## Aspido
dbAsp <- subset(dbH, order_ == "Aspidochirotida")
lblAsp <- sapply(1:nrow(dbAsp), function(i) genLabel(dbAsp[i, ]))
seqAsp <- allEch[dimnames(allEch)[[1]] %in% lblAsp, ]
dimnames(seqAsp)[[1]] <- make.unique(dimnames(seqAsp)[[1]])
treAsp <- nj(dist.dna(seqAsp))
treAsp$edge.length[treAsp$edge.length < 0] <- 0
bootAsp <- boot.phylo(treAsp, seqAsp, function(xx) nj(dist.dna(xx)), B=100)
treAsp$node.label <- bootAsp
treAspr <- root(treAsp, grep("^Pseudostichopus", treAsp$tip.label), resolve.root=TRUE)
treAsp4 <- phylo4d(treAspr, check.node.label="asdata")

grpAsp1 <- findGroups(treAsp4, threshold=.005)
grpAsp2 <- findGroups(treAsp4, threshold=.010)
grpAsp3 <- findGroups(treAsp4, threshold=.015)
grpAsp4 <- findGroups(treAsp4, threshold=.020)
grpAsp5 <- findGroups(treAsp4, threshold=.025)
grpAsp6 <- findGroups(treAsp4, threshold=.030)

## Holoth
dbHol <- subset(dbH, family == "Holothuriidae")
lblHol <- sapply(1:nrow(dbHol), function(i) genLabel(dbHol[i, ]))
seqHol <- allEch[dimnames(allEch)[[1]] %in% lblHol, ]
dimnames(seqHol)[[1]] <- make.unique(dimnames(seqHol)[[1]])
treHol <- nj(dist.dna(seqHol))
treHol$edge.length[treHol$edge.length < 0] <- 0
bootHol <- boot.phylo(treHol, seqHol, function(xx) nj(dist.dna(xx)), B=100)
treHol$node.label <- bootHol
treHolr <- root(treHol, grep("platei", treHol$tip.label), resolve.root=TRUE)
treHol4 <- phylo4d(treHolr, check.node.label="asdata")

grpHol1 <- findGroups(treHol4, threshold=.005)
grpHol2 <- findGroups(treHol4, threshold=.010)
grpHol3 <- findGroups(treHol4, threshold=.015)
grpHol4 <- findGroups(treHol4, threshold=.020)
grpHol5 <- findGroups(treHol4, threshold=.025)
grpHol6 <- findGroups(treHol4, threshold=.030)


nrow(seqH)
c(getUniqueSeq(seqH))
max(tipData(grpH3))
table(table(tipData(grpH3)$Group))

## holnm <- sapply(1:nrow(dbH), function(i) genSp(dbH[i, ]))
## holnm <- cbind(lblH, holnm)
## holnm <- holnm[lblH %in% dimnames(seqH)[[1]], ]
## write.csv(unique(holnm[,2]), file="~/Documents/echinoBarcode/uniqueHolNm.csv")



distH1 <- getIntraInterDist(grpH1, dist.dna(seqH, as.matrix=T), check.boot=80)
distH2 <- getIntraInterDist(grpH2, dist.dna(seqH, as.matrix=T), check.boot=80)
distH3 <- getIntraInterDist(grpH3, dist.dna(seqH, as.matrix=T), check.boot=80)
distH4 <- getIntraInterDist(grpH4, dist.dna(seqH, as.matrix=T), check.boot=80)
distH5 <- getIntraInterDist(grpH5, dist.dna(seqH, as.matrix=T), check.boot=80)
distH6 <- getIntraInterDist(grpH6, dist.dna(seqH, as.matrix=T), check.boot=80)

distDen1 <- getIntraInterDist(grpDen1, dist.dna(seqDen, as.matrix=T), check.boot=80)
distDen2 <- getIntraInterDist(grpDen2, dist.dna(seqDen, as.matrix=T), check.boot=80)
distDen3 <- getIntraInterDist(grpDen3, dist.dna(seqDen, as.matrix=T), check.boot=80)
distDen4 <- getIntraInterDist(grpDen4, dist.dna(seqDen, as.matrix=T), check.boot=80)
distDen5 <- getIntraInterDist(grpDen5, dist.dna(seqDen, as.matrix=T), check.boot=80)
distDen6 <- getIntraInterDist(grpDen6, dist.dna(seqDen, as.matrix=T), check.boot=80)

distApo1 <- getIntraInterDist(grpApo1, dist.dna(seqApo, as.matrix=T), check.boot=80)
distApo2 <- getIntraInterDist(grpApo2, dist.dna(seqApo, as.matrix=T), check.boot=80)
distApo3 <- getIntraInterDist(grpApo3, dist.dna(seqApo, as.matrix=T), check.boot=80)
distApo4 <- getIntraInterDist(grpApo4, dist.dna(seqApo, as.matrix=T), check.boot=80)
distApo5 <- getIntraInterDist(grpApo5, dist.dna(seqApo, as.matrix=T), check.boot=80)
distApo6 <- getIntraInterDist(grpApo6, dist.dna(seqApo, as.matrix=T), check.boot=80)

distAsp1 <- getIntraInterDist(grpAsp1, dist.dna(seqAsp, as.matrix=T), check.boot=80)
distAsp2 <- getIntraInterDist(grpAsp2, dist.dna(seqAsp, as.matrix=T), check.boot=80)
distAsp3 <- getIntraInterDist(grpAsp3, dist.dna(seqAsp, as.matrix=T), check.boot=80)
distAsp4 <- getIntraInterDist(grpAsp4, dist.dna(seqAsp, as.matrix=T), check.boot=80)
distAsp5 <- getIntraInterDist(grpAsp5, dist.dna(seqAsp, as.matrix=T), check.boot=80)
distAsp6 <- getIntraInterDist(grpAsp6, dist.dna(seqAsp, as.matrix=T), check.boot=80)

distHol1 <- getIntraInterDist(grpHol1, dist.dna(seqHol, as.matrix=T), check.boot=80)
distHol2 <- getIntraInterDist(grpHol2, dist.dna(seqHol, as.matrix=T), check.boot=80)
distHol3 <- getIntraInterDist(grpHol3, dist.dna(seqHol, as.matrix=T), check.boot=80)
distHol4 <- getIntraInterDist(grpHol4, dist.dna(seqHol, as.matrix=T), check.boot=80)
distHol5 <- getIntraInterDist(grpHol5, dist.dna(seqHol, as.matrix=T), check.boot=80)
distHol6 <- getIntraInterDist(grpHol6, dist.dna(seqHol, as.matrix=T), check.boot=80)


save.image(file="20120630.echinoBarcode.RData")

### number of ESUs vs threshold
nbEsu <- data.frame(Threshold=rep(1:6, 7),
                    Class=c(rep("Asteroidea", 6),
                      rep("Echinoidea", 6),
                      rep("Holothuroidea", 6),
                      rep("Dendrochirotida", 6),
                      rep("Apodida", 6),
                      rep("Aspidochirotida", 6),
                      rep("Holothuriidae", 6)),                    
                    NbESU=c(
                      max(tipData(grpA1)$Groups),
                      max(tipData(grpA2)$Groups),
                      max(tipData(grpA3)$Groups),
                      max(tipData(grpA4)$Groups),
                      max(tipData(grpA5)$Groups),
                      max(tipData(grpA6)$Groups),
                      max(tipData(grpE1)$Groups),
                      max(tipData(grpE2)$Groups),
                      max(tipData(grpE3)$Groups),
                      max(tipData(grpE4)$Groups),
                      max(tipData(grpE5)$Groups),
                      max(tipData(grpE6)$Groups),
                      max(tipData(grpH1)$Groups),
                      max(tipData(grpH2)$Groups),
                      max(tipData(grpH3)$Groups),
                      max(tipData(grpH4)$Groups),
                      max(tipData(grpH5)$Groups),
                      max(tipData(grpH6)$Groups),
                      max(tipData(grpDen1)$Groups),
                      max(tipData(grpDen2)$Groups),
                      max(tipData(grpDen3)$Groups),
                      max(tipData(grpDen4)$Groups),
                      max(tipData(grpDen5)$Groups),
                      max(tipData(grpDen6)$Groups),
                      max(tipData(grpApo1)$Groups),
                      max(tipData(grpApo2)$Groups),
                      max(tipData(grpApo3)$Groups),
                      max(tipData(grpApo4)$Groups),
                      max(tipData(grpApo5)$Groups),
                      max(tipData(grpApo6)$Groups),
                      max(tipData(grpAsp1)$Groups),
                      max(tipData(grpAsp2)$Groups),
                      max(tipData(grpAsp3)$Groups),
                      max(tipData(grpAsp4)$Groups),
                      max(tipData(grpAsp5)$Groups),
                      max(tipData(grpAsp6)$Groups),
                      max(tipData(grpHol1)$Groups),
                      max(tipData(grpHol2)$Groups),
                      max(tipData(grpHol3)$Groups),
                      max(tipData(grpHol4)$Groups),
                      max(tipData(grpHol5)$Groups),
                      max(tipData(grpHol6)$Groups)
                      ))

pdf(file="nbESU.pdf", width=8, height=3)
ggplot(subset(nbEsu, Class %in% c("Echinoidea", "Asteroidea", "Holothuroidea")),
              aes(x=Threshold, y=NbESU, group=Class, colour=Class)) + geom_line()
dev.off()

pdf(file="nbESU_hol.pdf", width=8, height=3)
ggplot(subset(nbEsu, Class %in% c("Aspidochirotida", "Dendrochirotida", "Apodida", "Holothuriidae")),
              aes(x=Threshold, y=NbESU, group=Class, colour=Class)) + geom_line()
dev.off()


## ### visualize sequence differences
## varSA <- apply(seqA, 2, function(x) max(base.freq(x)))
## varSH <- apply(seqH, 2, function(x) max(base.freq(x)))
## varSE <- apply(seqE, 2, function(x) max(base.freq(x)))

## par(mfrow=c(3,1))
## plot(1-varSA, col="blue", type="h", pch=15, cex=.5)
## plot(1-varSH, col="red", type="h", pch=15, cex=.5)
## plot(1-varSE, col="green", type="h", pch=15, cex=.5)

## matVar <- array(, dim=c(3,max(c(length(varSA), length(varSE), length(varSH)))))
## matVar[1,] <- 1-varSA
## matVar[2,] <- 1-varSE
## matVar[3,] <- 1-varSH
## barplot(matVar[,1:100], col=c("blue", "red", "green"), beside=T)

allDist <- rbind(cbind(Order=rep("Asteroidea", nrow(distA1)), Threshold=rep(1, nrow(distA1)), distA1),
        cbind(Order=rep("Asteroidea", nrow(distA2)), Threshold=rep(2, nrow(distA2)), distA2),
        cbind(Order=rep("Asteroidea", nrow(distA3)), Threshold=rep(3, nrow(distA3)), distA3),
        cbind(Order=rep("Asteroidea", nrow(distA4)), Threshold=rep(4, nrow(distA4)), distA4),
        cbind(Order=rep("Asteroidea", nrow(distA5)), Threshold=rep(5, nrow(distA5)), distA5),
        cbind(Order=rep("Asteroidea", nrow(distA6)), Threshold=rep(6, nrow(distA6)), distA6),        
        cbind(Order=rep("Echinoidea", nrow(distE1)), Threshold=rep(1, nrow(distE1)), distE1),
        cbind(Order=rep("Echinoidea", nrow(distE2)), Threshold=rep(2, nrow(distE2)), distE2),
        cbind(Order=rep("Echinoidea", nrow(distE3)), Threshold=rep(3, nrow(distE3)), distE3),
        cbind(Order=rep("Echinoidea", nrow(distE4)), Threshold=rep(4, nrow(distE4)), distE4),
        cbind(Order=rep("Echinoidea", nrow(distE5)), Threshold=rep(5, nrow(distE5)), distE5),
        cbind(Order=rep("Echinoidea", nrow(distE6)), Threshold=rep(6, nrow(distE6)), distE6),
        cbind(Order=rep("Holothuroidea", nrow(distH1)), Threshold=rep(1, nrow(distH1)), distH1),
        cbind(Order=rep("Holothuroidea", nrow(distH2)), Threshold=rep(2, nrow(distH2)), distH2),
        cbind(Order=rep("Holothuroidea", nrow(distH3)), Threshold=rep(3, nrow(distH3)), distH3),
        cbind(Order=rep("Holothuroidea", nrow(distH4)), Threshold=rep(4, nrow(distH4)), distH4),
        cbind(Order=rep("Holothuroidea", nrow(distH5)), Threshold=rep(5, nrow(distH5)), distH5),
        cbind(Order=rep("Holothuroidea", nrow(distH6)), Threshold=rep(6, nrow(distH6)), distH6),
        cbind(Order=rep("Dendrochirotida", nrow(distDen1)), Threshold=rep(1, nrow(distDen1)), distDen1),
        cbind(Order=rep("Dendrochirotida", nrow(distDen2)), Threshold=rep(2, nrow(distDen2)), distDen2),
        cbind(Order=rep("Dendrochirotida", nrow(distDen3)), Threshold=rep(3, nrow(distDen3)), distDen3),
        cbind(Order=rep("Dendrochirotida", nrow(distDen4)), Threshold=rep(4, nrow(distDen4)), distDen4),
        cbind(Order=rep("Dendrochirotida", nrow(distDen5)), Threshold=rep(5, nrow(distDen5)), distDen5),
        cbind(Order=rep("Dendrochirotida", nrow(distDen6)), Threshold=rep(6, nrow(distDen6)), distDen6),
        cbind(Order=rep("Apodida", nrow(distApo1)), Threshold=rep(1, nrow(distApo1)), distApo1),
        cbind(Order=rep("Apodida", nrow(distApo2)), Threshold=rep(2, nrow(distApo2)), distApo2),
        cbind(Order=rep("Apodida", nrow(distApo3)), Threshold=rep(3, nrow(distApo3)), distApo3),
        cbind(Order=rep("Apodida", nrow(distApo4)), Threshold=rep(4, nrow(distApo4)), distApo4),
        cbind(Order=rep("Apodida", nrow(distApo5)), Threshold=rep(5, nrow(distApo5)), distApo5),
        cbind(Order=rep("Apodida", nrow(distApo6)), Threshold=rep(6, nrow(distApo6)), distApo6),
        cbind(Order=rep("Aspidochirotida", nrow(distAsp1)), Threshold=rep(1, nrow(distAsp1)), distAsp1),
        cbind(Order=rep("Aspidochirotida", nrow(distAsp2)), Threshold=rep(2, nrow(distAsp2)), distAsp2),
        cbind(Order=rep("Aspidochirotida", nrow(distAsp3)), Threshold=rep(3, nrow(distAsp3)), distAsp3),
        cbind(Order=rep("Aspidochirotida", nrow(distAsp4)), Threshold=rep(4, nrow(distAsp4)), distAsp4),
        cbind(Order=rep("Aspidochirotida", nrow(distAsp5)), Threshold=rep(5, nrow(distAsp5)), distAsp5),
        cbind(Order=rep("Aspidochirotida", nrow(distAsp6)), Threshold=rep(6, nrow(distAsp6)), distAsp6),
        cbind(Order=rep("Holothuriidae", nrow(distHol1)), Threshold=rep(1, nrow(distHol1)), distHol1),
        cbind(Order=rep("Holothuriidae", nrow(distHol2)), Threshold=rep(2, nrow(distHol2)), distHol2),
        cbind(Order=rep("Holothuriidae", nrow(distHol3)), Threshold=rep(3, nrow(distHol3)), distHol3),
        cbind(Order=rep("Holothuriidae", nrow(distHol4)), Threshold=rep(4, nrow(distHol4)), distHol4),
        cbind(Order=rep("Holothuriidae", nrow(distHol5)), Threshold=rep(5, nrow(distHol5)), distHol5),
        cbind(Order=rep("Holothuriidae", nrow(distHol6)), Threshold=rep(6, nrow(distHol6)), distHol6))


svg(file="compareIntraInter.svg", width=9, height=5)
ggplot(subset(allDist, Order %in% c("Asteroidea", "Echinoidea", "Holothuroidea") &
              Threshold > 1),
       aes(x=dist, fill=typeDist, alpha=.5)) + geom_density() +
  facet_grid(Order ~ Threshold) + scale_fill_manual(values=c("yellow", "red"))
dev.off()

svg(file="compareIntraInterHol.svg", width=9, height=5)
ggplot(subset(allDist, Order %in% c("Dendrochirotida", "Apodida", "Aspidochirotida", "Holothuriidae") &
              Threshold > 1),
       aes(x=dist, fill=typeDist, alpha=.5)) + geom_density() +
  facet_grid(Order ~ Threshold) + scale_fill_manual(values=c("yellow", "red"))
dev.off()

### Geographic barriers
treE3 <- grpE3
tipLabels(treE3) <- paste(tipData(grpE3)$Group, tipLabels(treE3), sep="_")
bs <- nodeData(treE3)[,1, drop=F]
bsToShow <- rownames(bs)[bs > 80]
pdf(file="echinoids-015.pdf", height=100)
plot(as(treE3, "phylo"))
nodelabels(rep("", length(bsToShow)), as.integer(bsToShow), frame="circ", cex=.5)
dev.off()

treA3 <- grpA3
tipLabels(treA3) <- paste(tipData(grpA3)$Group, tipLabels(treA3), sep="_")
bs <- nodeData(treA3)[,1, drop=F]
bsToShow <- rownames(bs)[bs > 80]
pdf(file="asteroids-015.pdf", height=100)
plot(as(treA3, "phylo"))
nodelabels(rep("", length(bsToShow)), as.integer(bsToShow), frame="circ", cex=.5)
dev.off()

treH3 <- grpH3
tipLabels(treH3) <- paste(tipData(grpH3)$Group, tipLabels(treH3), sep="_")
bs <- nodeData(treH3)[,1, drop=F]
bsToShow <- rownames(bs)[bs > 80]
pdf(file="holothuroids-015.pdf", height=100)
plot(as(treH3, "phylo"), cex=.25)
nodelabels(rep("", length(bsToShow)), as.integer(bsToShow), frame="circ", cex=.125)
dev.off()

geo <- read.csv(file="geographicBarriers.csv")
propGeo <- with(geo, tapply(Allopatric, Class, function(x) {x <- x[!is.na(x)]; sum(x)/length(x)}))

svg(file="rateGeo.svg")
barplot(propGeo, ylim=c(0,1))
dev.off()

#############

stdSeqLength <- function(alg) {

}




## ## if only ufid

## ## yy <- sapply(xx, function(x) {
## ##   xt <- unlist(strsplit(x, "_"))
## ##   xt <- xt[length(xt)]
## ##   if (length(grep("^[0-9]+$", xt)) == 1)
## ##     TRUE
## ##   else FALSE
## ## })


## ## zz <- sapply(xx[yy], function(x) {
## ##   t <- unlist(strsplit(x, "_"))
## ##   t[length(t)]
## ##   })
