genLabel <- function(dbTmp) {
  seqNm <- paste(dbTmp$family, dbTmp$genusorhigher, dbTmp$modifier, dbTmp$species, dbTmp$Loc,
                 paste(dbTmp$"Collection.Code", dbTmp$Catalog_number, sep=""),
                 dbTmp$Sample, dbTmp$Type, sep="_")
  if ( dbTmp$Type != "")
      seqNm <- paste(seqNm, dbTmp$Sample, sep="_")
  seqNm <- gsub("_{2,}", "_", seqNm)
  seqNm <- gsub("_$", "", seqNm)
  seqNm <- gsub("\\s+$", "", seqNm)
  seqNm
}

genSp <- function(dbTmp) {
  seqNm <- paste(dbTmp$genusorhigher, dbTmp$modifier, dbTmp$species, sep="_")
  seqNm <- gsub("_{2,}", "_", seqNm)
  seqNm
}


genFasta <- function(db, out=file.path(tempdir(), paste(format(Sys.time(), "%Y%m%d-%H%M%S"), "seq.fas", sep="_"))) {
### db -- database in which the data is stored
### out -- file name of the fasta that will be generated
  for (i in 1:nrow(db)) {
    dbTmp <- db[i, ]
    seqNm <- genLabel(dbTmp)
    seqNm <- paste(">", seqNm, sep="")
    seqNm <- gsub("\\s+", "", seqNm)
    seqTmp <- dbTmp$Sequence
    cat(seqNm, "\n", dbTmp$Sequence, "\n", file=out,
        append=ifelse(i == 1, FALSE, TRUE), sep="")
  }
  TRUE
}

generate_unaligned_cuke_fasta <- function(cuke_seqs, out = file.path("data", "seq", "cuke_alg_unaligned.fas")) {
    unlink(out)
    seqs <- apply(cuke_seqs, 1, function(x) {
        paste0(">", x[1], "\n", x[2], collapse = "\n")
    })
    cat(seqs, sep = "\n", file = out, append = TRUE)
}

generate_aligned_cuke_fasta <- function(unaligned, aligned = file.path("data", "seq", "cuke_alg_aligned.fas")) {
    mafftCmd <- paste("mafft --auto --op 10 --thread -1",
                      unaligned, ">", aligned)

    system(mafftCmd)
}

generate_cleaned_cuke_fasta <- function(aligned, cleaned = file.path("data", "seq", "cuke_alg_cleaned.fas")) {
    ## Identify sequences with internal gaps
    seq_hol_c <- ape::read.dna(file = aligned, format = "fasta",
                               as.character = TRUE)
    seq_hol_c <- apply(seq_hol_c, 1,
                       function(x) paste(x, sep = "", collapse = ""))
    int_gap <- sapply(seq_hol_c, function(x)
        gregexpr("[actgnrmsykw]-+[actgnrmsykw]", x)[[1]][1] != -1)
    seq_with_gap <- names(int_gap[int_gap])

    ## Identify sequences with stop codons
    seq_hol <- ape::read.dna(file = aligned, format = "fasta")
    trans_seq <- foreach (i = 1:nrow(seq_hol)) %dopar% {
        seqinr::translate(as.character(seq_hol[i, ]), frame = 1,
                          numcode = 9)
    }
    seq_with_stop <- dimnames(seq_hol)[[1]][grep("\\*", trans_seq)]

    ## Remove sequences with internal gaps and stop codons
    to_rm <- union(seq_with_stop, seq_with_gap)

    message("Removing: ", to_rm)
    to_rm_ind <- match(to_rm, dimnames(seq_hol)[[1]])
    seq_hol <- seq_hol[-to_rm_ind, ]

    ## Write working copy of fasta file
    ape::write.dna(seq_hol, file = cleaned, format = "fasta", colsep = "")

}

load_cuke_alg <- function(alg) {

    ## identify sequences with ambiguities and rename them
    amb_seq <- chopper::checkAmbiguity(file = alg, quiet = TRUE)
    old_nm <- names(amb_seq)
    new_nm <- paste(old_nm, "_", sapply(amb_seq, length), "amb", sep = "")

    cuke_alg <- ape::read.dna(file = alg, format = "fasta")

    cuke_alg <- cuke_alg[-match(old_nm, dimnames(cuke_alg)[[1]]), ]
    cuke_alg <- chopper::cleanSeqLabels(cuke_alg, software = "RAxML")
    dimnames(cuke_alg)[[1]] <- gsub("\\\"", "", dimnames(cuke_alg)[[1]])
    invisible(cuke_alg)
}
