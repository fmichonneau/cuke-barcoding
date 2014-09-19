### ---- load-packages ---
source("R/packages.R")
source("R/load.R")
library(xtable)
library(car)
library(wesanderson)
library(tikzDevice)
library(reshape2)

### ---- sampling-maps ----
cukeDB <- load_cukeDB()

summGPS <- data.frame(uniq=paste(cukeDB$decimalLatitude, cukeDB$decimalLongitude, sep="/"))
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

pacificmap <- ggplot(summGPS) + annotation_map(globalMap, fill="gray40", colour="gray40") +
    geom_point(aes(x = long.recenter, y = latitude, size= nInd), colour="red", data=summGPS) +
    coord_map(projection = "mercator", orientation=c(90, 160, 0)) +
    theme(panel.background = element_rect(fill="aliceblue"),
          legend.position="top",
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          plot.margin=unit(rep(0, 4), "mm")) +
    scale_size_continuous(name="Number of individuals", range=c(1, 4)) +
    ylim(c(-45,45))

southmap <- ggplot(summGPS) + annotation_map(globalMap, fill="gray40", colour="gray40") +
    geom_point(aes(x = long.recenter, y = latitude, size= nInd), colour="red", data=summGPS) +
    coord_map(projection = "ortho", orientation=c(-90, 0, 0)) +
    theme(panel.background = element_rect(fill="aliceblue"),
          legend.position = "none",
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          plot.margin=unit(rep(0, 4), "mm")) +
    ylim(c(-90, -45))

northmap <- ggplot(summGPS) + annotation_map(globalMap, fill="gray40", colour="gray40") +
    geom_point(aes(x = long.recenter, y = latitude, size= nInd), colour="red", data=summGPS) +
    coord_map(projection = "ortho", orientation=c(90, 0, 0)) +
    theme(panel.background = element_rect(fill="aliceblue"),
          legend.position="none",
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          plot.margin=unit(rep(0, 4), "mm")) +
    ylim(c(45, 90))

multiplot(pacificmap, northmap, southmap, layout=matrix(c(1,1,2,3), ncol=2, byrow=T))

### ---- morphospecies-stats ----
taxonomyDf <- load_taxonomyDf()
cukeDB <- load_cukeDB()
cukeAlg <- load_cukeAlg()

cukeDBTmp <- cukeDB[match(dimnames(cukeAlg)[[1]], cukeDB$Labels_withAmb), ]
genera <- unique(cukeDBTmp$genusorhigher)
genera <- genera[-grep("\\?", genera)]

allSpp <- cukeDBTmp[cukeDBTmp$genusorhigher %in% genera,
                    c("family", "genusorhigher", "species")]
allSpp <- allSpp[-which(is.na(allSpp$family)), ]
uniqSpp <- unique(paste(allSpp$family, allSpp$genusorhigher, allSpp$species, sep="_"))
uniqSpp <- uniqSpp[-grep("\\d", uniqSpp)]
uniqSpp <- uniqSpp[-grep("n_?sp|new\\s?species|unique\\s?species", uniqSpp)]
uniqSpp <- uniqSpp[-grep("(cf|aff)\\.?\\s", uniqSpp)]
uniqSpp <- uniqSpp[-grep("_.{1,3}$", uniqSpp)] # too generous, removes H. bo
uniqSpp <- uniqSpp[-grep("pink|gray", uniqSpp)]
uniqSpp <- uniqSpp[-grep("_$", uniqSpp)]

fams <- sapply(strsplit(uniqSpp, "_"), function(x) x[1])
uniqSpp <- data.frame(family=fams, species=uniqSpp)
uniqSpp <- merge(uniqSpp, taxonomyDf, by.x="family", by.y="taxa", all.x=TRUE)

nSppAll <- nrow(uniqSpp)
nSppAsp <- nrow(subset(uniqSpp, higher=="Aspidochirotida"))
nSppHol <- nrow(subset(uniqSpp, family=="Holothuriidae"))

nSppApo <- nrow(subset(uniqSpp, higher=="Apodida"))
nSppDen <- nrow(subset(uniqSpp, higher=="Dendrochirotida"))
nSppEla <- nrow(subset(uniqSpp, higher=="Elasipodida"))

undescSpp <- unique(paste(allSpp$family, allSpp$genusorhigher, allSpp$species, sep="_"))
undescSpp <- grep("(n)?_sp(_|\\.|\\s|nov)?\\d?", undescSpp, value=T)
undescSpp <- undescSpp[-grep("spic|spin|spect", undescSpp)]

propHol <- 100*length(grep("^Holothuriidae", dimnames(cukeAlg)[[1]]))/dim(cukeAlg)[1]


### ---- esu-stats ----
isMonophyletic <- function(lbl, tree) {
    if (length(lbl) < 2) {
        NA
    } else {
        desc <- descendants(tree, MRCA(tree, lbl), "tips")
        if (length(desc) == length(lbl))
            TRUE
        else FALSE
    }
}

source("R/test-allopatry-functions.R")
manESU <- load_manESU()
cukeDistRaw <- load_cukeDist_raw()
localGap <- load_localGap()

manGrps <- load_species_manualGrps()

nESUs <- length(unique(manESU$ESU_noGeo))

esus <- strsplit(unique(manESU$ESU_noGeo), "_")
hasCryptic <- sapply(esus, function(x) length(x) > 2 & length(grep("nsp", x)) < 1)
esuCryptic <- esus[hasCryptic]
esuCryptic <- sapply(esuCryptic, function(x) paste0(x[1:2], collapse="_"))
nCryptic <- length(unique(esuCryptic))

newSpp <- sum(sapply(esus, function(x) length(grep("nsp", x)) > 0))

percentSinglHol <- 100*sum(sapply(manGrps, function(x) length(x) == 1))/length(manGrps)

percentGap   <- 100*sum(is.na(localGap$species))/nrow(localGap)
nGap <- sum(is.na(localGap$species))

nGreater02 <- sum(localGap$minInter > 0.02)
percentGreater02 <- 100*nGreater02/nrow(localGap)

holTree <- load_cukeTree_raw_phylo4()
toKeep <- intersect(manESU$Labels, tipLabels(holTree))
holTree <- subset(holTree, tips.include=toKeep)

taxLbl <- split(manESU$Labels, manESU$ESU_taxonomy)

taxMono <- sapply(taxLbl, function(x) isMonophyletic(x, holTree))
percentTaxMono <- 100*(1 - sum(taxMono, na.rm=TRUE)/sum(!is.na(taxMono)))

### ---- cluster-groups-data ----
taxonomyDf <- load_taxonomyDf()
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

nGrpsClustersDf <- cbind(factDf, nGrps=nGrpsVec, pSngl=pSnglVec)
nGrpsClustersDf <- merge(nGrpsClustersDf, taxonomyDf, all.x=TRUE)
nGrpsClustersDf$taxa <- factor(nGrpsClustersDf$taxa)
nGrpsClustersDf$distance <- factor(nGrpsClustersDf$distance)
nGrpsClustersDf$method <- "Clusters"

### ---- pairwise-groups-data ----
source("R/pairwise-groups-functions.R")
pairwiseGrpRes <- load_pairwiseGrpRes()
thresVec <- load_thresholdPairwise()
pairwiseGrpFactors <- rep(names(pairwiseGrpRes), each=length(thresVec))
pairwiseGrpFactors <- strsplit(pairwiseGrpFactors, "-")
pairwiseGrpFactors <- do.call("rbind", pairwiseGrpFactors)
pairwiseGrpFactors <- data.frame(distance=pairwiseGrpFactors[, 1],
                                 taxa=pairwiseGrpFactors[, 2],
                                 threshold=rep(thresVec, length(pairwiseGrpRes)))

nGrpsPairwise <- lapply(pairwiseGrpRes, function(x) sapply(x, length))
pSnglPairwise <- lapply(pairwiseGrpRes, function(x) sapply(x, getPropSngl))

nGrpsPairwiseDf <- data.frame(pairwiseGrpFactors, nGrps=unlist(nGrpsPairwise),
                              pSngl=unlist(pSnglPairwise))
nGrpsPairwiseDf$method <- "Pairwise"

## merging earlier messes up the ordering
nGrpsPairwiseDf <- merge(nGrpsPairwiseDf, taxonomyDf)

### ---- sampled-species ----
worms <- read.csv(file="data/raw/cuke-worms.csv", stringsAsFactors=FALSE)
worms <- subset(worms, worms$"Status..WoRMS." == "accepted")
wormsTab <- data.frame(xtabs(~ Order.current + Family.current + Genus.current, data=worms))
wormsTab <- wormsTab[wormsTab$Freq != 0, ]
names(wormsTab) <- c("Order", "Family", "Genus", "nSpp")

wormsSumOrder <- with(wormsTab, tapply(nSpp, Order, sum))
wormsSumFam <- with(wormsTab, tapply(nSpp, Family, sum))

wormsSum <- data.frame(Family = names(wormsSumFam), "# accepted"=wormsSumFam, check.names=FALSE)
wormsSum <- merge(wormsSum, wormsTab, by="Family", all.x=TRUE)
wormsSum <- wormsSum[, -match(c("Genus", "nSpp"), names(wormsSum))]
wormsSum <- wormsSum[!duplicated(wormsSum), ]
wormsSum <- rbind(wormsSum, data.frame(Family="", Order=names(wormsSumOrder),
                                       "# accepted"=wormsSumOrder, check.names=FALSE))

sampTab <- data.frame(xtabs(~ higher + family, data=uniqSpp))
sampTabOrder <- data.frame(xtabs(~ higher, data=uniqSpp), family = "")
sampTab <- rbind(sampTab, sampTabOrder)
sampTab <- sampTab[sampTab$Freq != 0, ]
names(sampTab) <- c("Order", "Family", "# sampled")

uniqFam <- as.character(taxonomyDf$taxa[taxonomyDf$rank == "Family"])
nESUsPerFam <- sapply(uniqFam, function(x) {
    length(load_species_clusterGrps(distance="raw", threshold=0.0225, taxa=x))
})
nESUsPerFam <- data.frame(taxa=names(nESUsPerFam), mtLineages=nESUsPerFam,
                          stringsAsFactors=FALSE)
nESUsPerFam <- merge(taxonomyDf, nESUsPerFam)

uniqOrder <- as.character(taxonomyDf$taxa[taxonomyDf$rank == "Order"])
uniqOrder <- uniqOrder[-match("all", uniqOrder)]
nESUsPerOrder <- sapply(uniqOrder, function(x) {
    length(load_species_clusterGrps(distance="raw", threshold=0.0225, taxa=x))
})
nESUsPerOrder <- data.frame(taxa=names(nESUsPerOrder), mtLineages=nESUsPerOrder,
                            stringsAsFactors=FALSE)
nESUsPerOrder <- merge(taxonomyDf, nESUsPerOrder)

nESUsPerTaxa <- rbind(nESUsPerFam, nESUsPerOrder)
nESUsPerTaxa$taxa <- as.character(nESUsPerTaxa$taxa)
nESUsPerTaxa$higher <- as.character(nESUsPerTaxa$higher)
nESUsPerTaxa <- nESUsPerTaxa[, -match("rank", colnames(nESUsPerTaxa))]
toSwitch <- is.na(nESUsPerTaxa$higher)
nESUsPerTaxa$higher[toSwitch] <- nESUsPerTaxa$taxa[toSwitch]
nESUsPerTaxa$taxa[toSwitch] <- ""
names(nESUsPerTaxa) <- c("Family", "Order", "# mtLineages")

sampTab <- merge(sampTab, wormsSum)
sampTab <- merge(sampTab, nESUsPerTaxa)
sampTab$"# mtLineages"[is.na(sampTab$"# mtLineages")] <- nESUsPerOrder
sampTab <- sampTab[, c("Order", "Family", "# accepted", "# sampled", "# mtLineages")]

sampTab$Family <- factor(sampTab$Family,
                         levels=c(levels(sampTab$Family)[nlevels(sampTab$Family)],
                         levels(sampTab$Family)[-nlevels(sampTab$Family)]))
sampTab <- sampTab[order(sampTab$Family), ]
sampTab <- sampTab[order(sampTab$Order), ]
hlinePos <- cumsum(table(sampTab$Order))
sampTab$Order <- as.character(sampTab$Order)
sampTab$Order <- paste("\\textbf{", sampTab$Order, "}", sep="")
sampTab$Order[duplicated(sampTab$Order)] <- ""
sampTab$"# sampled"[nzchar(sampTab$Order)] <-
    paste("\\textbf{", sampTab$"# sampled"[nzchar(sampTab$Order)], "}", sep="")
sampTab$"# accepted"[nzchar(sampTab$Order)] <-
    paste("\\textbf{", sampTab$"# accepted"[nzchar(sampTab$Order)], "}", sep="")
sampTab$"# mtLineages"[nzchar(sampTab$Order)] <-
    paste("\\textbf{", sampTab$"# mtLineages"[nzchar(sampTab$Order)], "}", sep="")


sampTab <- rbind(sampTab,
                 cbind(Order = "\\textbf{Total}", Family = "",
                       "# sampled" = paste0("\\textbf{", nrow(uniqSpp), "}"),
                       "# accepted" = paste0("\\textbf{", sum(subset(wormsSum, Family != "")$acceptedSpp), "}"),
                       "# mtLineages" = paste0("\\textbf{", length(load_species_clusterGrps(distance="raw",
                           threshold=0.0225,
                           taxa="all")), "}")))

colnames(sampTab)[3:5] <- paste0("\\", colnames(sampTab)[3:5])

sampXtab <- xtable(sampTab,
                   caption=paste("Number of named morpho-species sampled (\\# sampled)",
                       "number of accepted species (\\# accepted) and",
                       "number of mtLineages (\\# mtLineages) estimated with the clustering",
                       "method and a 4.5\\% threshold for each",
                       "family and each order of sea cucumber. There were", nESUs, "ESUs",
                       "delineated for the Holothuriidae. Not all families were",
                       "sampled, thus totals in some orders are more than sum of family",
                       "diversities. This classification does not include modifications",
                       "proposed by Smirnov \\cite{Smirnov2012}. Estimations of the number",
                       "of mtLineages is based on different datasets for families and orders,",
                       "thus totals for some orders may differ because of taxonomic uncertainty",
                       "(samples identified at the order level but not at the family level), or",
                       "differences in lineages delineation when the entire order is considered"),
                   label="tab:sampled-species", display=c("s", "s", "s", "d", "d", "d"),
                   align=rep("l", 6))

print(sampXtab, include.rownames=FALSE, hline.after=c(-1, 0, hlinePos, nrow(sampXtab)),
      sanitize.text.function=function(x) {x},
      table.placement="ht!")




### ---- compare-manual-cluster-ESUs-data ----
accuracyGrps <- function(tree) {
    clstrGrps <- tdata(tree, "tip")[, "Groups", drop=FALSE]
    clstrGrps <- cbind(Labels=rownames(clstrGrps), clstrGrps=clstrGrps)
    manGrps <- manESU[, c("Labels", "ESU_noGeo")]
    manGrps$manGrps <- as.numeric(factor(manESU$ESU_noGeo))
    compGrps <- merge(clstrGrps, manGrps, by="Labels")

    uniqGrps <- unique(compGrps$manGrps)
    compRes <- vector("list", length(uniqGrps))
    for (i in 1:length(uniqGrps)) {
        tmpDt <- subset(compGrps, manGrps == uniqGrps[i])
        if (length(unique(tmpDt$Groups)) > 1)
            splits <- unique(tmpDt$Groups)
        else splits <- NA
        tmpDt2 <- subset(compGrps, Groups %in% unique(tmpDt$Groups))
        if( length(unique(tmpDt2$manGrps)) > 1)
            lumps <- unique(tmpDt2$manGrps)
        else lumps <- NA
        compRes[[i]] <- list(splits=splits, lumps=lumps)
    }
    nLumps <- sapply(compRes, function(x) x$lumps)
    nLumps <- unlist(nLumps)
    nLumps <- length(unique(nLumps[!is.na(nLumps)]))
    nSplits <- sapply(compRes, function(x) x$splits)
    nSplits <- unlist(nSplits)
    nSplits <- length(unique(nSplits[!is.na(nSplits)]))
    list(nLumps=nLumps, nSplits=nSplits)
}

pwiseTreeRaw <- lapply(load_thresholdPairwise(), function(x) load_tree_pairwiseGrps("raw", taxa="Holothuriidae", x))
pwiseTreeK80 <- lapply(load_thresholdPairwise(), function(x) load_tree_pairwiseGrps("K80", taxa="Holothuriidae", x))
clstrTreeRaw <- lapply(load_thresholdClusters(), function(x) load_tree_clusterGrps("raw", taxa="Holothuriidae", x))
clstrTreeK80 <- lapply(load_thresholdClusters(), function(x) load_tree_clusterGrps("K80", taxa="Holothuriidae", x))

pwiseGrpsRaw <- lapply(pwiseTreeRaw, accuracyGrps)
pwiseGrpsK80 <- lapply(pwiseTreeK80, accuracyGrps)
clstrGrpsRaw <- lapply(clstrTreeRaw, accuracyGrps)
clstrGrpsK80 <- lapply(clstrTreeK80, accuracyGrps)

pwiseDatRaw <- data.frame(method="pairwise", distance="raw", do.call("rbind", lapply(pwiseGrpsRaw, function(x) c(x[1], x[2]))))
pwiseDatK80 <- data.frame(method="pairwise", distance="K2P", do.call("rbind", lapply(pwiseGrpsK80, function(x) c(x[1], x[2]))))
clstrDatRaw <- data.frame(method="clustering", distance="raw", do.call("rbind", lapply(clstrGrpsRaw, function(x) c(x[1], x[2]))))
clstrDatK80 <- data.frame(method="clustering", distance="K2P", do.call("rbind", lapply(clstrGrpsK80, function(x) c(x[1], x[2]))))

compareManESUs <- rbind(pwiseDatRaw, pwiseDatK80, clstrDatRaw, clstrDatK80)

nGrps <- sapply(c(pwiseTreeRaw, pwiseTreeK80, clstrTreeRaw, clstrTreeK80),
                function(x) max(tdata(x, "tip")[, "Groups"]))

compareManESUs <- data.frame(threshold=rep(load_thresholdPairwise(), 4), compareManESUs, nGrps=nGrps)
compareManESUs$nSplits <- as.numeric(compareManESUs$nSplits)
compareManESUs$nLumps <- as.numeric(compareManESUs$nLumps)
allErrors <- compareManESUs$nSplits + compareManESUs$nLumps
compareManESUs <- cbind(compareManESUs, allErrors=allErrors)
minError <- compareManESUs[which.min(allErrors), ]

compareManESUs$pLumps <- compareManESUs$nLumps/compareManESUs$nGrps
compareManESUs$pSplits <- compareManESUs$nSplits/compareManESUs$nGrps

saveRDS(compareManESUs, file="tmp/compareManESUs.rds")

### ---- compare-manual-cluster-ESUs-plot ----
compareManESUs <- melt(compareManESUs, id.vars=c("threshold", "method", "distance"),
                       measure.vars=c("pLumps", "pSplits"))
levels(compareManESUs$distance)[levels(compareManESUs$distance) == "raw"] <- "Uncorrected"

ggplot(compareManESUs, aes(x=threshold, y=value, fill=variable)) + geom_bar(stat="identity") +
    facet_wrap(~ distance + method) + ylab("Proportion of misassigned mtLineages") +
    xlab("Genetic distance threshold") +
    scale_fill_discrete(labels=c("lumped", "oversplit")) +
    guides(fill=guide_legend(title=NULL)) +
    theme(legend.position="top")

### ---- local-gap-plot ----
localGap <- load_localGap()
esuNm <- read.csv(file="data/raw/ESU_names.csv", stringsAsFactors=FALSE, header=FALSE)

localGap$species <- gsub("\\..+$", "", localGap$species)
for (i in 1:nrow(esuNm)) { localGap$species <- gsub(esuNm[i, 1], esuNm[i, 2], localGap$species)}
localGap$species <- gsub("_", " ", localGap$species)

ggplot(localGap) + geom_point(aes(x=maxIntra, y=minInter, colour=species)) +
    geom_abline(slope=1, linetype=3, colour="gray40") + coord_fixed() +
    xlim(c(0, 0.18)) + ylim(c(0, 0.18)) +
    xlab("Maximum intra-ESU distance")  + ylab("Minimum distance to the closest neighboring ESU") +
    scale_colour_discrete(name="ESUs")

### ---- geography-diversification ----
source("R/test-allopatry-functions.R")
holTree <- load_tree_manualGrps()
esuRange <- testRangeType(holTree, load_cukeDist_raw(), load_cukeDB())

esuRange <- esuRange[complete.cases(esuRange), ]

esuRange$species <- as.character(esuRange$species)

esuRange$rangeType <- factor(esuRange$rangeType, levels=c("allopatric", "sympatric", "parapatric"))

esuRange <- esuRange[-match("54-Holothuria_excellens/192-Pearsonothuria_graeffei", esuRange$species), ]
esuRange[match("177-Holothuria_unicolor/184-Holothuria_zihuatanensis", esuRange$species), "rangeType"] <- "allopatric"

tabRangeType <- table(esuRange$rangeType)

percentAllo <- 100*tabRangeType["allopatric"]/sum(tabRangeType)
percentSymp <- 100*tabRangeType["sympatric"]/sum(tabRangeType)
percentPara <- 100*tabRangeType["parapatric"]/sum(tabRangeType)

nSymp <- tabRangeType["sympatric"]

medianDist <- median(subset(esuRange, rangeType == "allopatric")$meanInterDist)
nBelowMedian <- sum(subset(esuRange, rangeType=="sympatric")$meanInterDist < medianDist)

ggplot(esuRange) +
    geom_point(data=esuRange[!is.na(esuRange$rangeType), ],
               aes(x=rangeType, y=meanInterDist, colour=rangeType),
               position=position_jitter(width=.07, height=0)) +
    coord_flip() + xlab("") + ylab("Mean genetic distances between sister ESUs (uncorrected)") +
    theme(legend.position="none")

### ---- barriers ----
manESU <- load_manESU()
cukeDistRaw <- load_cukeDist_raw()
iopoDist <- manESU[grep("_(IO|PO)$", manESU$ESU_genetic), ]
iopoESU <- gsub("_[A-Z]{2}$", "", iopoDist$ESU_genetic)
iopoESU <- gsub("(\\d{1})[a-z]{1}", "\\1", iopoESU)
iopoDist$iopoESU <- iopoESU

iopoGrps <- split(iopoDist$ESU_genetic, iopoDist$iopoESU)
iopoGrps <- lapply(iopoGrps, function(x) unique(x))
iopoGrps <- iopoGrps[sapply(iopoGrps, function(x) length(x) == 2)]

iopoGenDist <- sapply(iopoGrps, function(x) {
    ind1 <- subset(iopoDist, iopoDist$ESU_genetic== x[1])$Labels
    ind2 <- subset(iopoDist, iopoDist$ESU_genetic== x[2])$Labels
    interESUDist(ind1, ind2, cukeDistRaw)$mean
})

iopoESUcat <- sapply(iopoGrps, function(x) {
    if (length(grep("\\d{1}[a-z]{1}", x)) == 2) "inter"
    else "intra"
})

iopoDat <- data.frame(barrier="IO/PO", genDist=iopoGenDist, esuCat=iopoESUcat)

##

rspoDist <- manESU[grep("_(RS|PO)$", manESU$ESU_genetic), ]
rspoESU <- gsub("_[A-Z]{2}$", "", rspoDist$ESU_genetic)
rspoESU <- gsub("(\\d{1})[a-z]{1}", "\\1", rspoESU)
rspoDist$rspoESU <- rspoESU

rspoGrps <- split(rspoDist$ESU_genetic, rspoDist$rspoESU)
rspoGrps <- lapply(rspoGrps, function(x) unique(x))
rspoGrps <- rspoGrps[sapply(rspoGrps, function(x) length(x) == 2)]

rspoGenDist <- sapply(rspoGrps, function(x) {
    ind1 <- subset(rspoDist, rspoDist$ESU_genetic== x[1])$Labels
    ind2 <- subset(rspoDist, rspoDist$ESU_genetic== x[2])$Labels
    interESUDist(ind1, ind2, cukeDistRaw)$mean
})

rspoESUcat <- sapply(rspoGrps, function(x) {
    if (length(grep("\\d{1}[a-z]{1}", x)) == 2) "inter"
    else "intra"
})

rspoDat <- data.frame(barrier="RS/PO", genDist=rspoGenDist, esuCat=rspoESUcat)

##

hipoDist <- manESU[grep("_(HI|PO|IW)$", manESU$ESU_genetic), ]
hipoESU <- gsub("_[A-Z]{2}$", "", hipoDist$ESU_genetic)
hipoESU <- gsub("(\\d{1})[a-z]{1}", "\\1", hipoESU)
hipoDist$hipoESU <- hipoESU

hipoGrps <- split(hipoDist$ESU_genetic, hipoDist$hipoESU)
hipoGrps <- lapply(hipoGrps, function(x) unique(x))
hipoGrps <- hipoGrps[sapply(hipoGrps, function(x) length(x) == 2)]

hipoGenDist <- sapply(hipoGrps, function(x) {
    ind1 <- subset(hipoDist, hipoDist$ESU_genetic== x[1])$Labels
    ind2 <- subset(hipoDist, hipoDist$ESU_genetic== x[2])$Labels
    interESUDist(ind1, ind2, cukeDistRaw)$mean
})

hipoESUcat <- sapply(hipoGrps, function(x) {
    if (length(grep("\\d{1}[a-z]{1}", x)) == 2) "inter"
    else "intra"
})

hipoDat <- data.frame(barrier="HI/PO", genDist=hipoGenDist, esuCat=hipoESUcat)

geoDat <- rbind(iopoDat, rspoDat, hipoDat)
geoDat <- subset(geoDat, esuCat == "inter")

percentIOPO <- 100*sum(geoDat$barrier == "IO/PO")/tabRangeType["allopatric"]
percentRSPO <- 100*sum(geoDat$barrier == "RS/PO")/tabRangeType["allopatric"]
percentHIPO <- 100*sum(geoDat$barrier == "HI/PO")/tabRangeType["allopatric"]

ggplot(geoDat) + geom_point(aes(y=barrier, x=genDist, colour=barrier)) +
    ylab("") + xlab("Genetic distances (uncorrected)") +
    theme(legend.position="none")






### ---- groups-comparison-plot ----
manGrps <- load_species_manualGrps()
nGrpsDf <- rbind(nGrpsClustersDf, nGrpsPairwiseDf)
levels(nGrpsDf$distance)[levels(nGrpsDf$distance) == "K80"] <- "K2P"
levels(nGrpsDf$distance)[levels(nGrpsDf$distance) == "k2p"] <- "K2P"

zisPal <- wes.palette(5, "Zissou")[c(1,2,3,5)]

nSpp <- data.frame(taxa=c("all", "Aspidochirotida", "Holothuriidae"),
                   nspp=c(nSppAll, nSppAsp, nSppHol))

nHol <- data.frame(taxa=c("all", "Aspidochirotida", "Holothuriidae"),
                   nspp=c(NA, NA, length(manGrps)))

pSinglHol <- data.frame(taxa=c("all", "Aspidochirotida", "Holothuriidae"),
                        psngl=c(NA, NA, sum(sapply(manGrps,
                            function(x) length(x) == 1))/length(manGrps)))

tmpDt <- subset(nGrpsDf, taxa %in% c("all", "Aspidochirotida", "Holothuriidae"))
nESU <- ggplot(tmpDt, aes(x=threshold, y=nGrps, colour=interaction(distance, method),
                          shape=interaction(distance, method))) +
    geom_line() + geom_point() + facet_wrap( ~ taxa) +
    geom_hline(data=nSpp, aes(yintercept=nspp), colour="gray40", linetype=2) +
    geom_hline(data=nHol, aes(yintercept=nspp), colour="darkgreen", linetype=3) +
    ylab("Number of mtLineages") +
    theme(legend.position="top", axis.title.x=element_blank()) +
    scale_colour_manual(name=element_blank(),
                        values=zisPal,
                        breaks=c("K2P.Clusters", "raw.Clusters",
                            "K2P.Pairwise", "raw.Pairwise"),
                        labels=c("K2P - Clustering",
                            "Uncorrected - Clustering",
                            "K2P - Pairwise",
                            "Uncorrected - Pairwise")) +
    scale_shape_manual(name=element_blank(),
                        values=seq(from=15, length.out=4),
                        breaks=c("K2P.Clusters", "raw.Clusters",
                            "K2P.Pairwise", "raw.Pairwise"),
                        labels=c("K2P - Clustering",
                            "Uncorrected - Clustering",
                            "K2P - Pairwise",
                            "Uncorrected - Pairwise"))

pSngl <- ggplot(subset(nGrpsDf, taxa %in% c("all", "Holothuriidae", "Aspidochirotida")),
              aes(x=threshold, y=pSngl, colour=interaction(distance, method),
                  shape=interaction(distance, method))) +
    geom_line() +  geom_point() + facet_wrap( ~ taxa) +
    geom_hline(data=pSinglHol, aes(yintercept=psngl), colour="darkgreen", linetype=3) +
    ylab("Proportion of singletons") + xlab("Distance threshold") +
    scale_colour_manual(values=zisPal) + scale_shape_manual(values=seq(from=15, length.out=4)) +
    theme(legend.position="none")

multiplot(nESU, pSngl, layout=matrix(c(rep(1, 4), rep(2, 3)), ncol=1))

### ---- sm-groups-comparison-plot ----
tmpDtSm <- subset(nGrpsDf, taxa %in% c("Dendrochirotida", "Apodida", "Elasipodida"))

nSppSm <- data.frame(taxa=c("Dendrochirotida", "Apodida", "Elasipodida"),
                   nspp=c(nSppDen, nSppApo, nSppEla))

nESUSm <- ggplot(tmpDtSm, aes(x=threshold, y=nGrps, colour=interaction(distance, method),
                          shape=interaction(distance, method))) +
    geom_line() + geom_point() + facet_wrap( ~ taxa) +
    geom_hline(data=nSppSm, aes(yintercept=nspp), colour="gray40", linetype=2) +
    ylab("Number of mtLineages") +
    theme(legend.position="top", axis.title.x=element_blank()) +
    scale_colour_manual(name=element_blank(),
                        values=zisPal,
                        breaks=c("K2P.Clusters", "raw.Clusters",
                            "K2P.Pairwise", "raw.Pairwise"),
                        labels=c("K2P - Clustering",
                            "Uncorrected - Clustering",
                            "K2P - Pairwise",
                            "Uncorrected - Pairwise")) +
    scale_shape_manual(name=element_blank(),
                        values=seq(from=15, length.out=4),
                        breaks=c("K2P.Clusters", "raw.Clusters",
                            "K2P.Pairwise", "raw.Pairwise"),
                        labels=c("K2P - Clustering",
                            "Uncorrected - Clustering",
                            "K2P - Pairwise",
                            "Uncorrected - Pairwise"))

pSnglSm <- ggplot(tmpDtSm,
              aes(x=threshold, y=pSngl, colour=interaction(distance, method),
                  shape=interaction(distance, method))) +
    geom_line() +  geom_point() + facet_wrap( ~ taxa) +
    ylab("Proportion of singletons") + xlab("Distance threshold") +
    scale_colour_manual(values=zisPal) + scale_shape_manual(values=seq(from=15, length.out=4)) +
    theme(legend.position="none")

multiplot(nESUSm, pSnglSm, layout=matrix(c(rep(1, 4), rep(2, 3)), ncol=1))



### ---- map-orders ----
globalMap <- map_data("world2")

dendroMap <- subset(load_cukeDB(), order=="Dendrochirotida")[, c("order", "genusorhigher", "species", "decimalLatitude", "decimalLongitude")]
center <- 150
dendroMap$long.recenter <- ifelse(dendroMap$decimalLongitude < center - 180,
                                  dendroMap$decimalLongitude + 360,
                                  dendroMap$decimalLongitude)
dendroMap <- dendroMap[complete.cases(dendroMap), ]

den <- ggplot(dendroMap) + annotation_map(globalMap, fill="gray40", colour="gray40") +
    geom_point(aes(x = long.recenter, y = decimalLatitude), data=dendroMap, color="red") +
    coord_map(projection = "mercator", orientation=c(90, 160, 0)) +
    theme(panel.background = element_rect(fill="aliceblue"),
          legend.position="none") + xlab("Longitude") + ylab("Latitude") +
    xlim(c(0, 300)) + ylim(c(-45, 30))

aspidoMap <- subset(load_cukeDB(), order=="Aspidochirotida")[, c("order", "genusorhigher", "species", "decimalLatitude", "decimalLongitude")]
center <- 150
aspidoMap$long.recenter <- ifelse(aspidoMap$decimalLongitude < center - 180,
                                  aspidoMap$decimalLongitude + 360,
                                  aspidoMap$decimalLongitude)
aspidoMap <- aspidoMap[complete.cases(aspidoMap), ]

asp <- ggplot(aspidoMap) + annotation_map(globalMap, fill="gray40", colour="gray40") +
    geom_point(aes(x = long.recenter, y = decimalLatitude), data=aspidoMap, color="red") +
    coord_map(projection = "mercator", orientation=c(90, 160, 0)) +
    theme(panel.background = element_rect(fill="aliceblue"),
          legend.position="none") + xlab("Longitude") + ylab("Latitude") +
    xlim(c(0, 300)) + ylim(c(-45, 30))

apodiMap <- subset(load_cukeDB(), order=="Apodida")[, c("order", "genusorhigher", "species", "decimalLatitude", "decimalLongitude")]
center <- 150
apodiMap$long.recenter <- ifelse(apodiMap$decimalLongitude < center - 180,
                                  apodiMap$decimalLongitude + 360,
                                  apodiMap$decimalLongitude)
apodiMap <- apodiMap[complete.cases(apodiMap), ]

apo <- ggplot(apodiMap) +  annotation_map(globalMap, fill="gray40", colour="gray40") +
    geom_point(aes(x = long.recenter, y = decimalLatitude), data=apodiMap, color="red") +
    coord_map(projection = "mercator", orientation=c(90, 160, 0)) +
    theme(panel.background = element_rect(fill="aliceblue"),
          legend.position="none") + xlab("Longitude") + ylab("Latitude") +
    xlim(c(0, 300)) + ylim(c(-45, 30))

allOrderMap <- rbind(apodiMap, aspidoMap, dendroMap)

ggplot(allOrderMap) +  annotation_map(globalMap, fill="gray40", colour="gray40") +
    geom_point(aes(x = long.recenter, y = decimalLatitude, colour=order), data=allOrderMap,
               position=position_jitter(height=1, width=1)) +
    coord_map(projection = "mercator", orientation=c(90, 160, 0)) +
    theme(panel.background = element_rect(fill="aliceblue"),
          legend.position="none") + xlab("Longitude") + ylab("Latitude") +
    xlim(c(0, 300)) + ylim(c(-45, 30))

### ---- mantel-test ----
### Doesn't really work given too much data, could instead do by species
###   but then issues with multiple comparisons...
## cukeAlg <- load_cukeAlg()
## cukeDB <- load_cukeDB()
## tmpCoords <- cukeDB[match(dimnames(cukeAlg)[[1]], cukeDB$Labels_withAmb), ]
## stopifnot(nrow(cukeAlg) == nrow(tmpCoords))
## genDistMat <- ape::dist.dna(cukeAlg, model="K80", as.matrix=TRUE)
## geoDist <- CalcGeoDists(cbind(deg2rad(tmpCoords$decimalLongitude),
##                               deg2rad(tmpCoords$decimalLatitude)))


### ---- esu-table ----
manESU <- load_manESU()
esuNm <- read.csv(file="data/raw/ESU_names.csv", stringsAsFactors=FALSE, header=FALSE)
manGrps <- load_species_manualGrps()

uniqESU <- unique(manESU$ESU_noGeo)
for (i in 1:nrow(esuNm)) {
    uniqESU <- gsub(esuNm[i, 1], paste0("\\\\textit{", esuNm[i, 2], "}"), uniqESU)
}
uniqESU <- gsub("_", " ", uniqESU)

esuVouch <- sapply(manGrps, function(x) {
    tmpDB <- cukeDB[match(x, cukeDB$Labels_withAmb), ]
    paste(tmpDB$Collection.Code, tmpDB$Catalog_number, sep="", collapse=", ")
})

esuVouch <- gsub("_", " ", esuVouch)

headerTable <- paste("\\hline", paste(c("ESUs", "Vouchers"), collapse=" & "), "\\\\ \\hline \\endfirsthead \n",
                     "\\caption{(continued) specimen information} \n")

print(xtable(cbind(ESUs=uniqESU, Vouchers=esuVouch),
             caption="List of manually delineated ESUs and associated vouchers",
             label="tab:esu-info", align=c("llp{4in}")),
      tabular.environment="longtable", floating=FALSE,
      hline.after = c(-1, length(esuVouch)),
      add.to.row = list(pos = list(-1, 0),
          command = c(headerTable, "\\hline \\endhead \n")),
      sanitize.text.function=function(x) {x}, caption.placement="top",
      table.placement="ht!" include.rownames=FALSE)

### ---- figure-test ----
aa <- cukeAlg[grep("cinerascens", dimnames(cukeAlg)[[1]]), ]
aa <- aa[order(tdata(ptr, "tip")[grep("cinerascens", tipLabels(ptr)), "Groups"]), ]
aa <- aa[-grep("Ningaloo|Moorea", dimnames(aa)[[1]]), ]
dimnames(aa)[[1]] <- seq_len(dim(aa)[1])
write.csv(file="tmp/cine.csv", dist.dna(aa, model="raw", as.matrix=TRUE))

tt <- load_tree_clusterGrps("raw", "Holothuriidae", 0.02)
tt <- subset(tt, tips.include=grep("cinerascens", tipLabels(tt)))

svg(file="tmp/cine.svg")
plot(ladderize(as(tt, "phylo")), show.tip.label=FALSE)
add.scale.bar()
dev.off()


########################       below is old code      ##########################

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
