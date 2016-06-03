generate_unaligned_cuke_fasta <- function(cukeDB_noLabels) {
    genFasta(cukeDB_noLabels,
             out = file.path("data", "seq", "cukeAlg-unaligned.fas"))
}

generate_aligned_cuke_fasta <- function(unaligned) {

    aligned <- file.path("data", "seq", "cukeAlg-aligned.fas")
    mafftCmd <- paste("mafft --auto --op 10 --thread -1",
                      unaligned, ">", aligned)

    system(mafftCmd)
}

generate_cleaned_cuke_fasta <- function(aligned) {
    cleanFile <- file.path("data", "seq", "cukeAlg-cleaned.fas")

    ## Identify sequences with internal gaps
    seqHolC <- ape::read.dna(file=aligned, format="fasta", as.character=TRUE)
    seqHolC <- apply(seqHolC, 1, function(x) paste(x, sep="", collapse=""))
    intGap <- sapply(seqHolC, function(x)
        gregexpr("[actgnrmsykw]-+[actgnrmsykw]", x)[[1]][1] != -1)
    seqWithIntGap <- names(intGap[intGap])

    ## Identify sequences with stop codons
    seqHol <- read.dna(file=aligned, format="fasta")
    tranE <- foreach (i = 1:nrow(seqHol)) %dopar% {
        translate(as.character(seqHol[i, ]), frame=1, numcode=9)
    }
    seqWithStop <- dimnames(seqHol)[[1]][grep("\\*", tranE)]

    ## Remove sequences with internal gaps and stop codons
    toRm <- union(seqWithStop, seqWithIntGap)

    message("Removing: ", toRm)
    ## dimnames(seqHol)[[1]][match(toRm, dimnames(seqHol)[[1]])] <-
    ## paste("stop-intgap", toRm, sep="_")
    toRmInd <- match(toRm, dimnames(seqHol)[[1]])
    seqHol <- seqHol[-toRmInd, ]

    ## Write working copy of fasta file
    write.dna(seqHol, file=cleanFile, format="fasta", colsep="")

}

load_cukeAlg <- function(algFile) {

    ## identify sequences with ambiguities and rename them
    ambSeq <- checkAmbiguity(file=algFile, quiet=TRUE)
    oldNm <- names(ambSeq)
    newNm <- paste(oldNm, "_", sapply(ambSeq, length), "amb", sep="")

    cukeAlg <- read.dna(file=algFile, format="fasta")

    ## dimnames(cukeAlg)[[1]][match(oldNm, dimnames(cukeAlg)[[1]])] <- newNm

    cukeAlg <- cukeAlg[-match(oldNm, dimnames(cukeAlg)[[1]]), ]
    cukeAlg <- cleanSeqLabels(cukeAlg, software="RAxML")
    dimnames(cukeAlg)[[1]] <- gsub("\\\"", "", dimnames(cukeAlg)[[1]])
    invisible(cukeAlg)
}
