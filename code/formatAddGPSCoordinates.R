
##### ------- Need to load allDB in memory before using the code below -------- 

#######  Add GPS coordinates
## load("20131115.guidm.RData") ## this file contains all the UF database
## guidmHol <- subset(guidm, class == "Holothuroidea") ## keep only Cukes
## save(guidmHol, file="guidmHol.RData") 
## load("guidmHol.RData")
## guidmHol$UFnumber <- sapply(guidmHol$catalogNumber, function(x) gsub("-.+$", "", x))
## guidmHol$UFnumber <- paste("UF", guidmHol$UFnumber, sep="")
## allDB$UFnumber <- paste(allDB$Collection.Code, allDB$Catalog_number, sep="")
## allDBtmp <- merge(allDB, guidmHol[, c("UFnumber", "decimalLatitude", "decimalLongitude")],
##                   by="UFnumber", all.x=TRUE, all.y=FALSE)

formatCoords <- function(x) {
    x <- gsub("’", "'", x)
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
    ## format (X)X°( )(X)X.(X)X(')(’)
    else if (length(grep("[0-9]{1,2}°\\s?[0-9]{1,2}(\\.[0-9]{1,2}'?)?", x)) > 0) {
        x <- gsub("'$", "", x)
        tCoord <- unlist(strsplit(x, "°"))
        tCoord <- as.numeric(tCoord)
        if (tCoord[1] < 0) {
            tCoord[1] - tCoord[2]/60
        }
        else {
            tCoord[1] + tCoord[2]/60
        }            
    }
    ## format (X)X°( )(X)X'( )XX(")
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
    ## format XX XX XX
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

######

xx <- read.csv(file="/tmp/redSeaCoords.csv", stringsAsFactors=FALSE)

xx$decLat <- sapply(xx$lat, formatCoords)
xx$decLong <- sapply(xx$long, formatCoords)

write.csv(xx, file="/tmp/redSeaCoords-fixed.csv")

######
##### Get GPS coordinates that are missing for UF specimens
#####

allUF <- read.csv(file="data/20131209.AllCukes.csv", stringsAsFactors=FALSE)

missing <- subset(allUF, is.na(latitude) | is.na(longitude))
noDup <- missing[!duplicated(missing$locality), ]

write.csv(noDup, file="data/noDuplicate.csv")

uniqSt <- read.csv(file="data/noDuplicate.csv", stringsAsFactors=FALSE)


###### Got this from -- http://www.r-bloggers.com/using-google-maps-api-and-r/
#### This script uses RCurl and RJSONIO to download data from Google's API:
#### Latitude, longitude, location type (see explanation at the end), formatted address
#### Notice there is a limit of 2,500 calls per day
 
library(RCurl)
library(RJSONIO)
library(plyr)

higherToUse <- ifelse(nzchar(uniqSt$secondary_subdivision), uniqSt$secondary_subdivision,
                      ifelse(nzchar(uniqSt$primary_subdivision), uniqSt$primary_subdivision,
                             uniqSt$country_archipelago))

higherToUse <- uniqSt$country_archipelago

uniqSt$address <- paste(uniqSt$localityToUse, higherToUse, sep=",")
for (i in 1:nrow(uniqSt)) {
    if (length(grep("°", uniqSt$address[i]) > 0)) {
        uniqSt$address[i] <- uniqSt$localityToUse[i]
        next
    }
    ## if (length(grep("island|atoll", uniqSt$secondary_subdivision[i], ignore.case=TRUE)) == 0) {
    ##     uniqSt$address[i] <- paste(uniqSt$address[i], uniqSt$primary_subdivision[i], sep=",")
    ## }
    ## else next
    ## if (uniqSt$primary_subdivision[i] == uniqSt$country_archipelago[i]) {
    ##     next
    ## }
    ## else {
    ##     uniqSt$address[i] <- paste(uniqSt$address[i], uniqSt$country_archipelago[i], sep=",")
    ## }
}
uniqSt$address <- gsub("\\s+,", ", ", uniqSt$address)
uniqSt$address <- gsub("\\,{2,}", ",", uniqSt$address)
uniqSt$address <- gsub("\\s", "+", uniqSt$address)
uniqSt$address <- gsub("\\+{2,}", "+", uniqSt$address)
uniqSt$address <- gsub(",\\+", ",", uniqSt$address)
uniqSt$address <- gsub("^,", "", uniqSt$address)

url <- function(address, return.call = "json", sensor = "false") {
    root <- "http://maps.google.com/maps/api/geocode/"
    u <- paste(root, return.call, "?address=", address, "&sensor=", sensor, sep = "")
    return(URLencode(u))
}
 
geoCode <- function(address, verbose=FALSE) {
    if(verbose) cat(address,"\n")
    u <- url(address)
    doc <- getURL(u)
    x <- fromJSON(doc, simplify = FALSE)
    if(x$status=="OK") {
        lat <- x$results[[1]]$geometry$location$lat
        lng <- x$results[[1]]$geometry$location$lng
        location_type <- x$results[[1]]$geometry$location_type
        formatted_address <- x$results[[1]]$formatted_address
        return(c(lat, lng, location_type, formatted_address))
    } else {
        browser()
        return(c(NA,NA,NA, NA))
    }
}


uniqStGPS <- vector("list", length(uniqSt$address))
for (i in 1:length(uniqSt$address)) {
    uniqStGPS[[i]] <- geoCode(uniqSt$address[i])
    system("sleep 0.5")
}

uniqSt$gpslat <- sapply(uniqStGPS, function(x) x[1])
uniqSt$gpslong <- sapply(uniqStGPS, function(x) x[2])
uniqSt$formattedAddress <- sapply(uniqStGPS, function(x) x[4])

write.csv(uniqSt, file="/tmp/testgeo.csv")

###### after cleaning up this output and adding the gps coordinates that were missing
### reload and copy the gps coordinates back to the original dataset (allUF)

uniqSt <- read.csv(file="data/noDuplicate_withGPS.csv", stringsAsFactors=FALSE)

indexToMatch <- match(allUF$locality, uniqSt$locality)
itm <- match(uniqSt$locality, allUF$locality)

uniqStNew <- uniqSt

for (i in 1:nrow(uniqSt)) {
    message(i)
    tmpLat <- allUF[allUF$locality %in% uniqSt$locality[i], ]$latitude
    if (uniqSt$locality[i] == "")
        next
    if (any(is.na(tmpLat))) {
        if (all(is.na(tmpLat))) {
            allUF[allUF$locality %in% uniqSt$locality[i], ]$latitude <- uniqSt$gpslat[i]
            allUF[allUF$locality %in% uniqSt$locality[i], ]$longitude <- uniqSt$gpslong[i]
        }
        else {
            allUF[allUF$locality %in% uniqSt$locality[i], ]$latitude[is.na(allUF[allUF$locality %in% uniqSt$locality[i], ]$latitude)] <- uniqSt$gpslat[i]
            allUF[allUF$locality %in% uniqSt$locality[i], ]$longitude[is.na(allUF[allUF$locality %in% uniqSt$locality[i], ]$longitude)] <- uniqSt$gpslong[i]
        }
    }
    else next
}


write.csv(allUF, file="/tmp/alluf_withgps.csv")

## after that I had to fix a few things by hands. But this file now
## has coordinates for all specimens except 1.

allUF <- read.csv(file="data/alluf_withgps.csv", stringsAsFactors=FALSE)
indexToMatch <- with(allDB, which(Collection.Code == "UF" & is.na(decimalLatitude) &
                                  nzchar(Catalog_number)))


for (i in 1:length(indexToMatch)) {
    xx <- subset(allUF, uf_id == allDB$Catalog_number[indexToMatch[i]])
    if (nrow(xx) == 0) {
        message("no match for UFID", allDB$Catalog_number[indexToMatch[i]])
    }
    else if (nrow(xx) > 1) {
        browser()
    }
    else {
        testNA <- allDB$decimalLatitude[indexToMatch[i]]
        if(is.na(testNA)) {
            allDB$decimalLatitude[indexToMatch[i]] <- xx$latitude
            allDB$decimalLongitude[indexToMatch[i]] <- xx$longitude
        }
        else {
            message("already filled")
            brower()
        }
    }
}

write.csv(allDB, file="/tmp/allDB-withcoords.csv")

#### Filling out the blanks
### for this last pass, I'm going to look for unique locations in the database
### create a little spreadsheet with coordinates for these locations
### fill them out in the marbol spreadsheet

### first load the latest version of the spreadsheet

allDB <- read.csv(file="data/MARBoL_Echinos_VIII_2013.csv", stringsAsFactors=FALSE) # nrow = 7017

missingcoords <- subset(allDB, class_ == "Holothuroidea" & is.na(decimalLatitude))

tmpUniq <- missingcoords[, c("Region", "Loc", "further.loc")]
dupLoc <- duplicated(paste(tmpUniq$Region, tmpUniq$Loc, tmpUniq$further.loc))

write.csv(tmpUniq[!dupLoc, ], file="/tmp/uniqLoc.csv")
