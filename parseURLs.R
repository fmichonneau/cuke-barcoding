#!/usr/bin/Rscript
bib <- readLines("cuke-barcoding.bib")
output <- "cuke-barcoding_nourl.bib"

allUrl <- grep("^url", bib)
allTitle <- grep("^title", bib)

nReplaced <- 0

for (i in 1:length(allUrl)) {
    tmpUrl <- bib[allUrl[i]]
    if (length(grep("github|r-project", tmpUrl))) {
        iTitle <- allTitle[which.min(abs(allTitle - allUrl[i]))]
        tmpTitle <- bib[iTitle]
        tmpUrl <- gsub("\\s+=\\s+", "", tmpUrl)
        tmpUrl <- gsub(",$", "", tmpUrl)
        tmpUrl <- gsub("^url", "url", tmpUrl)
        newTitle <- gsub("(}},)$", paste("  \\\\", tmpUrl, "\\1", sep=""), tmpTitle)
        bib[iTitle] <- newTitle
        nReplaced <- nReplaced + 1L
    }
}

bib <- bib[ - allUrl]

cat(bib, sep="\n", file=output)
message(nReplaced, " URL have been added to the title.")

