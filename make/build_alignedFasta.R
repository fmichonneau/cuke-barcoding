#!/usr/bin/env Rscript

source("R/genFasta.R")
source("R/load.R")
library(ape)
library(doMC)
registerDoMC()
library(seqinr)

tmpDir <- tempdir()
unaligned <- "cukeAlg-unaligned.fas"
aligned <- "cukeAlg-aligned.fas"
cleaned <- "cukeAlg-cleaned.fas"

algFile <- file.path("data", aligned)
cleanFile <- file.path("data", cleaned)

## Generate FASTA file
if (genFasta(load_cukeDB_noLabels(), out=file.path(tmpDir, unaligned))) {
    message("Unaligned fasta file generated from CSV")
}

mafftCmd <- paste("mafft --auto --op 10 --thread -1",
                  file.path(tmpDir, unaligned), ">", file.path(tmpDir, aligned),
                  "2>",
                  file.path("tmp", paste0(format(Sys.time(), "%Y%m%d-%H%M%S"),
                                          "-mafft.out")))
message(mafftCmd)
message("mafft output is written to ", file.path("tmp", "mafft.out"))
system(mafftCmd)

if (file.copy(file.path(tmpDir, aligned),
              algFile, overwrite=TRUE)) {
    message("Alignment generated")
}

## Identify sequences with internal gaps
seqHolC <- ape::read.dna(file=algFile, format="fasta", as.character=TRUE)
seqHolC <- apply(seqHolC, 1, function(x) paste(x, sep="", collapse=""))
intGap <- sapply(seqHolC, function(x)
                 gregexpr("[actgnrmsykw]-+[actgnrmsykw]", x)[[1]][1] != -1)
seqWithIntGap <- names(intGap[intGap])

## Identify sequences with stop codons
seqHol <- read.dna(file=algFile, format="fasta")
tranE <- foreach (i = 1:nrow(seqHol)) %dopar% {
    translate(as.character(seqHol[i, ]), frame=1, numcode=9)
}
seqWithStop <- dimnames(seqHol)[[1]][grep("\\*", tranE)]

## Remove sequences with internal gaps and stop codons
toRm <- union(seqWithStop, seqWithIntGap)
## dimnames(seqHol)[[1]][match(toRm, dimnames(seqHol)[[1]])] <- paste("stop-intgap", toRm, sep="_")
toRmInd <- match(toRm, dimnames(seqHol)[[1]])
seqHol <- seqHol[-toRmInd, ]

## Write working copy of fasta file
write.dna(seqHol, file=cleanFile, format="fasta", colsep="")

### These 3 sequences are not represented by other representative
###  it might be worth trying to figure out if we can clean up the
###  sequences to deal with the issues
###  - FRM-194
###  - NMV F112128
###  - NIWA 38032
