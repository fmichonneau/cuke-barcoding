#######  Add GPS coordinates
## load("20131115.guidm.RData") ## this file contains all the UF database
## guidmHol <- subset(guidm, class == "Holothuroidea") ## keep only Cukes
## save(guidmHol, file="guidmHol.RData") 
## load("guidmHol.RData")
## guidmHol$UFnumber <- sapply(guidmHol$catalogNumber, function(x) gsub("-.+$", "", x))
## guidmHol$UFnumber <- paste("UF", guidmHol$UFnumber, sep="")
## allDB$UFnumber <- paste(allDB$Collection.Code, allDB$Catalog_number, sep="")
## holDBtmp <- merge(allDB, guidmHol[, c("UFnumber", "decimalLatitude", "decimalLongitude")],
##                   by="UFnumber", all.x=TRUE, all.y=FALSE)

formatCoords <- function(x) {
    x <- gsub("\\s+$", "", x)
    x <- gsub("º", "°", x)
    x <- gsub("°$", "", x)
    x <- gsub("(.+)[SW]$", "-\\1", x)
    x <- gsub("[NE]$", "", x)
    x <- gsub("\\s+$", "", x)
    x <- gsub("°$", "", x)
    ### klunky -- could be much improved
    if(!nzchar(x)) {        
        ""
    }
    ## format (X)X°( )(X)X
    else if (length(grep("[0-9]{1,2}°\\s?[0-9]{1,2}\\.[0-9]{1,2}'?", x)) > 0) {
        x <- gsub("'$", "", x)
        tCoord <- unlist(strsplit(x, "°"))
        tCoord <- as.numeric(tCoord)
        -(tCoord[1] + tCoord[2]/60)
    }
    else if (length(grep("[0-9]{1,2}°\\s?[0-9]{1,2}'\\s?[0-9]{,2}\\\"?", x)) > 0) {
        x <- gsub(" ", "", x)
        tCoord <- unlist(strsplit(x, "[°'\"]"))
        tCoord <- as.numeric(tCoord)
        if (tCoord[1] < 0) {
            tCoord[1] - tCoord[2]/60 - tCoord[3]/3600
        }
        else {
            tCoord[1] + tCoord[2]/60 + tCoord[3]/3600
        }
    }
    else if (length(grep("[0-9]{2}\\s{1}[0-9]{2}\\s{1}[0-9]{2}", x)) > 0) {
        tCoord <- unlist(strsplit(x, " "))
        tCoord <- as.numeric(tCoord)
        if (tCoord[1] < 0)  tCoord[1] - tCoord[2]/60 - tCoord[3]/3600
        else tCoord[1] + tCoord[2]/60 + tCoord[3]/3600            
    }
    else if (length(grep("[0-9]{1,2}\\.?[0-9]?", x, perl=TRUE)) > 0) {
        x
    }
    else {
        message("pb")
        "pb"
    }
}

molaf <- read.csv(file="MOL_data.csv", stringsAsFactors=FALSE)
molaf$lat <- sapply(molaf$lat, formatCoords)
molaf$long <- sapply(molaf$long, formatCoords)
stopifnot(length(grep("°", molaf$lat)) == 0 &&
          length(grep("°", molaf$long)) == 0 &&
          !any(as.numeric(molaf$lat) > 90, na.rm=T) &&
          !any(as.numeric(molaf$lat) < -90, na.rm=T) &&
          !any(as.numeric(molaf$long) > 180, na.rm=T) &&
          !any(as.numeric(molaf$long) < -180, na.rm=T))
## write.csv(molaf, file="data/molaf_withDecCoord.csv")

alltissues <- allDB$Sample[! allDB$Sample %in% molaf$Tissue]
moltissues <- grep("^MOL|^NDMQ|^NIWAMOL", alltissues, value=T)
moltissues

indexMolTissues <- match(molaf$Tissue, allDB$Sample)
allDB$decimalLatitude[indexMolTissues[!is.na(indexMolTissues)]] <- molaf$lat[which(!is.na(indexMolTissues))]
allDB$decimalLongitude[indexMolTissues[!is.na(indexMolTissues)]] <- molaf$long[which(!is.na(indexMolTissues))]
