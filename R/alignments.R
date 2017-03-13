
generate_unaligned_cuke_fasta <- function(cuke_seqs, out = file.path("data", "seq", "cuke_alg_unaligned.fas")) {
    unlink(out)
    seqs <- apply(cuke_seqs, 1, function(x) {
        paste0(">", x[1], "\n", x[2], collapse = "\n")
    })
    cat(seqs, sep = "\n", file = out, append = TRUE)
    if (file.info(out)$size < 1)
        stop("Something is wrong... empty file for ", sQuote(out), ".")
}

generate_aligned_cuke_fasta <- function(unaligned, aligned = file.path("data", "seq", "cuke_alg_aligned.fas")) {
    mafftCmd <- paste("mafft --auto --op 10 --thread -1",
                      unaligned, ">", aligned)
    system(mafftCmd)
    if (file.info(aligned)$size < 1)
        stop("Something is wrong... empty file for ", sQuote(aligned), ". Make sure mafft is installed.")
}


guess_reading_frame <- function(alg, subsample = 25, verbose = TRUE) {
    alg_sub <- alg[sample(seq_len(dim(alg)[1]), size = subsample), ]
    n_nonsense <- integer(3)
    for (i in 1:3) {
        trans_alg <- mclapply(seq_len(nrow(alg_sub)), function(j) {
            seqinr::translate(as.character(alg_sub[j, ]), frame = i, numcode = 9)
        }, mc.preschedule = FALSE)
        n_stops_total <- vapply(trans_alg, function(x)
            sum(grepl("\\*", x)),
            numeric(1))
        n_nonsense[i] <- sum(n_stops_total > 0)
    }
    if (verbose) {
        message("Number of non sense codons for positions 1, 2 and 3 respectively: ",
                paste(n_nonsense, collapse = ", "), ".")
    }
    if (diff(sort(n_nonsense))[1] < subsample/5)
        warning("Not a lot of difference to be sure it's the correct frame.")
    if (min(n_nonsense) >  subsample * 0.05)
        warning("It seems like a very high error rate in the best frame.")
    which.min(n_nonsense)
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
    frame <- guess_reading_frame(seq_hol)
    trans_seq <- mclapply(seq_len(nrow(seq_hol)), function(i) {
        seqinr::translate(as.character(seq_hol[i, ]), frame = frame,
                          numcode = 9)
    }, mc.preschedule = FALSE)
    seq_with_stop <- dimnames(seq_hol)[[1]][grep("\\*", trans_seq)]

    ## Remove sequences with internal gaps and stop codons
    to_rm <- union(seq_with_stop, seq_with_gap)

    message("Removing: ", paste(to_rm, collapse = ", "))
    to_rm_ind <- match(to_rm, dimnames(seq_hol)[[1]])
    seq_hol <- seq_hol[-to_rm_ind, ]

    ## Write working copy of fasta file
    ape::write.dna(seq_hol, file = cleaned, format = "fasta", colsep = "")

}

load_cuke_alg <- function(alg) {

    ## identify sequences with ambiguities and rename them
    #amb_seq <- chopper::checkAmbiguity(file = alg, quiet = TRUE)
    #old_nm <- names(amb_seq)
    #new_nm <- paste(old_nm, "_", sapply(amb_seq, length), "amb", sep = "")

    cuke_alg <- ape::read.dna(file = alg, format = "fasta")

    #cuke_alg <- cuke_alg[-match(old_nm, dimnames(cuke_alg)[[1]]), ]
    #cuke_alg <- chopper::cleanSeqLabels(cuke_alg, software = "RAxML")
    #dimnames(cuke_alg)[[1]] <- gsub("\\\"", "", dimnames(cuke_alg)[[1]])
    invisible(cuke_alg)
}
