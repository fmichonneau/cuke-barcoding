### Compare with sequences submitted to genbank
testGBseq <- function(gb, db) {
    if (! file.exists("/tmp/seq")) {
        stop("Create /tmp/seq before running this script.")
    }
    res <- array(, dim=c(nrow(gb), 7), dimnames=list(NULL, c("seqNm", "sameLength",
                                           "seqLen1", "seqLen2",
                                           "nAmb1", "nAmb2",
                                           "distGenIs0")))
    lFiles <- character(nrow(gb))
    for (i in 1:nrow(gb)) {
        res[i, 1] <- gb$genbankNb[i]
        algPth <- "/tmp/seq"
        fileNm <- paste(gb$genbankNb[i], ".fas", sep="")
        algNm <- gsub("fas$", "afa", fileNm)
        tmpDB <- subset(db, GenBankSubmission == gb$genbankNb[i])
        if (nrow(tmpDB) == 0) {
            res[i, ] <- c(NA, NA, NA, NA, NA, NA, NA)
        }
        else {
            lFiles[i] <- algNm
            ## Sequence 1 - what's in the database
            ## Sequence 2 - what's in GenBank
            seqNm1 <- paste(">", tmpDB$GenBankSubmission, "_", tmpDB$Sample, sep="")
            seqNm2 <- paste(">", gb$genbankNb[i], "_", gb$vou[i], sep="")
            seq1 <- tmpDB$Sequence
            seq2 <- gb$vdb.seq[i]
            lSeq1 <- length(gregexpr("[actgACTG]", seq1)[[1]])
            lSeq2 <- length(gregexpr("[actgACTG]", seq2)[[1]])
            if (lSeq1 < 100 || lSeq2 < 100) {
                warning("sequence too short to be true")
                browser()
            }
            res[i, 2] <- lSeq1 == lSeq2
            res[i, 3] <- lSeq1
            res[i, 4] <- lSeq2
            cat(seqNm1, "\n", seq1, "\n", file=file.path(algPth, fileNm), append=FALSE, sep="")
            cat(seqNm2, "\n", seq2, "\n", file=file.path(algPth, fileNm), append=TRUE, sep="")
            res[i, 5] <- length(gregexpr("[^-]", seq1)[[1]]) - lSeq1
            res[i, 6] <- length(gregexpr("[^-]", seq2)[[1]]) - lSeq2
            system(paste("muscle -in", file.path(algPth, fileNm), "-out", file.path(algPth, algNm)))
            res[i, 7] <- dist.dna(read.dna(file=file.path(algPth, algNm), format="fasta"))
        }
    }
    oFile <- paste("/tmp/", format(Sys.time(), "%Y%m%d-%H%M%S"), "allseq.fas", sep="")
    mSeq <- mergeAlignment(lFiles[nzchar(lFiles)], output=oFile, seqFolder="/tmp/seq")
    res
}

if (FALSE) {
    ufgb <- read.csv(file="data/UF_genbankSequences.csv", stringsAsFactors=FALSE)

    compareSeqTmp <- testGBseq(gb=ufgb, db=allDB)

    compareSeq <- data.frame(compareSeqTmp, stringsAsFactors=FALSE)
    compareSeq$sameLength <- as.logical(compareSeq$sameLength)
    allGood <- compareSeq$sameLength & compareSeq$nAmb1 == compareSeq$nAmb2 & compareSeq$distGenIs0 == 0
    compareSeqPb <- compareSeq[!allGood, ]
    write.csv(compareSeqPb, file="/tmp/compareSeq.csv")
    }
