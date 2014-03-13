
setwd("~/Documents/CukeBarcoding/")
source("code/data_preparation.R")

ophDB <- subset(allDB, class_ == "Ophiuroidea")
ophDB <- subset(ophDB, pass.seq == "yes")
lSeq <- sapply(ophDB$Sequence, function(x) length(gregexpr("[actgACTG]", x)[[1]]))
lAmb <- sapply(ophDB$Sequence, function(x) length(gregexpr("[^-]", x)[[1]]))
##sum(table(lAmb)[as.numeric(names(table(lAmb))) > 500])
ophDB <- ophDB[lAmb > 500, ]

### Taxonomic check
testGenera <- as.matrix(xtabs(~ genusorhigher + family, data=ophDB, subset=family != "Uncertain"))
resGenera <- apply(testGenera, 1, function(x) sum(x != 0))
stopifnot(all(resGenera == 1))
testFamily <- as.matrix(xtabs(~ family + order, data=ophDB, subset=family != "Uncertain"))
resFamily <- apply(testFamily, 1, function(x) sum(x != 0))
stopifnot(all(resFamily == 1))
