load_cukeAlg <- function(overwrite=FALSE) {
    fnm <- "data/cukeAlg-flagAmb.rds"
    if (file.exists(fnm) && !overwrite) {
        cukeAlg <- readRDS(file=fnm)
    }
    else {
        if (overwrite)  message("Re-creating", fnm)
        algFile <- "data/cukeAlg-cleaned.fas"

        ## identify sequences with ambiguities and rename them
        ambSeq <- checkAmbiguity(file=algFile)
        oldNm <- names(ambSeq)
        newNm <- paste(oldNm, "_", sapply(ambSeq, length), "amb", sep="")

        cukeAlg <- read.dna(file=algFile, format="fasta")
        dimnames(cukeAlg)[[1]][match(oldNm, dimnames(cukeAlg)[[1]])] <- newNm
        cukeAlg <- cleanSeqLabels(cukeAlg, software="RAxML")
        dimnames(cukeAlg)[[1]] <- gsub("\\\"", "", dimnames(cukeAlg)[[1]])
        saveRDS(cukeAlg, file="data/cukeAlg-flagAmb.rds")
    }
    invisible(cukeAlg)
}
