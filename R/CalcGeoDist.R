
####
### By Scott Chamberlain http://r.789695.n4.nabble.com/Geographic-distance-between-lat-long-points-in-R-td3442338.html
## Convert degrees to radians
deg2rad <- function(deg) return(deg*pi/180)

## Calculates the geodesic distance between two points specified by
## radian latitude/longitude using the Haversine formula
gcd.hf <- function(long1, lat1, long2, lat2, warn = TRUE) {
    if (warn)
        warning("Did you convert the latitudes in radian first?!")
    R <- 6371 # Earth mean radius [km]
    delta.long <- (long2 - long1)
    delta.lat <- (lat2 - lat1)
    a <- sin(delta.lat/2)^2 + cos(lat1) * cos(lat2) * sin(delta.long/2)^2
    c <- 2 * asin(min(1,sqrt(a)))
    d = R * c
    return(d) # Distance in km
}

## Fxn to calculate matrix of distances between each two sites
CalcGeoDists <- function(latlongs) {
    name <- list(rownames(latlongs), rownames(latlongs))
    n <- nrow(latlongs)
    z <- matrix(0, n, n, dimnames = name)
    for (i in 1:n) {
        for (j in 1:n) z[i, j] <- gcd.hf(long1 = latlongs[i, 1],
                                         lat1 = latlongs[i, 2], long2 = latlongs[j, 1], lat2 = latlongs[j,2])
    }
    z <- as.dist(z)
    return(z)
}
