

## Select sequences
holDB <- subset(allDB, class_ == "Holothuroidea")  # nrow = 4385
holDB <- subset(holDB, pass.seq != "GenBank")      # nrow = 4360
holDB <- subset(holDB, pass.seq != "fix")          # nrow = 4358
holDB <- subset(holDB, pass.seq != "no_seq_yet")   # nrow =h 3466
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
filename <- paste(format(Sys.time(), "%Y%m%d-%H%M%S"), "cukes.fas", sep="-")
aligned <- "cukeBarcodes-aligned.fas"
genFasta(holDB, out=file.path(tmpDir, filename))
system(paste("mafft --auto --op 10 --thread -1", file.path(tmpDir, filename), ">", file.path(tmpDir, aligned)))
file.copy(file.path(tmpDir, aligned), file.path("data", filename))

## Identify sequences with internal gaps
seqHolC <- read.dna(file="data/latestAlg.fas", format="fasta", as.character=TRUE)
seqHolC <- apply(seqHolC, 1, function(x) paste(x, sep="", collapse=""))
intGap <- sapply(seqHolC, function(x) gregexpr("[actgn]-+[actgn]", x)[[1]][1] != -1) # should I consider ambiguities here?
seqWithIntGap <- names(intGap[intGap])

## Identify sequences with stop codons
seqHol <- read.dna(file="data/latestAlg.fas", format="fasta")
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
write.dna(seqHol, file="data/workingAlg.fas", format="fasta", colsep="")


### These 3 sequences are not represented by other representative
###  it might be worth trying to figure out if we can clean up the
###  sequences to deal with the issues
###  - FRM-194
###  - NMV F112128
###  - NIWA 38032
