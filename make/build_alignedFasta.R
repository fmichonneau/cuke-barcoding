#!/usr/bin/env Rscript

allDB <- readRDS("data/raw/cukeBarcodes.csv.rds")
source("R/genFasta.R")

## Select sequences
holDB <- subset(allDB, class_ == "Holothuroidea")  # nrow = 4385
holDB <- subset(holDB, pass.seq != "GenBank")      # nrow = 4360
holDB <- subset(holDB, pass.seq != "fix")          # nrow = 4358
holDB <- subset(holDB, pass.seq != "no_seq_yet")   # nrow = 3466
holDB <- subset(holDB, pass.seq != "no")           # nrow = 3402
holDB <- subset(holDB, Notes != "MH sequence")     # nrow = 3379
holDB <- subset(holDB, pass.seq != "duplicate")    # 
lSeq <- sapply(holDB$Sequence, function(x) length(gregexpr("[actgACTG]", x)[[1]])) # only non-ambiguous bp
lAmb <- sapply(holDB$Sequence, function(x) length(gregexpr("[^-]", x)[[1]]))       # all bp
## sum(table(lSeq)[as.numeric(names(table(lSeq))) > 500 ])
holDB <- holDB[lAmb > 500, ] # nrow = 2894 -- this also takes care of empty sequences (only -)

## Taxonomic check
testGenera <- as.matrix(xtabs(~ genusorhigher + family, data=holDB, subset=family != "Uncertain"))
resGenera <- apply(testGenera, 1, function(x) sum(x != 0))
stopifnot(all(resGenera == 1))
testFamily <- as.matrix(xtabs(~ family + order, data=holDB, subset=family != "Uncertain"))
resFamily <- apply(testFamily, 1, function(x) sum(x != 0))
stopifnot(all(resFamily == 1))

## check for duplicated samples
dup <- holDB[duplicated(holDB$Sample), "Sample"]
stopifnot(length(dup) == 0)

## Generate FASTA file
tmpDir <- tempdir()
unaligned <- "cukeBarcodes-unaligned.fas"
aligned <- "cukeBarcodes-aligned.fas"
if (genFasta(holDB, out=file.path(tmpDir, unaligned))) {
    message("Unaligned fasta file generated from CSV")
}
mafftCmd <- paste("mafft --auto --op 10 --thread -1",
                  file.path(tmpDir, unaligned), ">", file.path(tmpDir, aligned),
                  "2>", file.path("tmp", paste0(format(Sys.time(), "%Y%m%d-%H%M%S"), "-mafft.out")))
message(mafftCmd)
message("mafft output is written to ", file.path("tmp", "mafft.out"))
system(mafftCmd)

if (file.copy(file.path(tmpDir, aligned), file.path("data", aligned), overwrite=TRUE)) {
    message("Alignment generated")
}
