
build_raxml_tree <- function(alg_file, raxml_path = "raxmlHPC-PTHREADS-SSE3") {
    ## targets to create:
    ## - the directory where the magic takes place
    raxml_output <- "data/raxml"
    ## - the partition file (by default, each codon position will be its own partition)
    part_file <- "data/raxml/cukeAlg.partition"
    ## - the suffix for the RAxML analyses
    raxml_suffix <- "cukeBarcodes"


    if (!file.exists(raxml_output)) {
        dir.create(raxml_output)
    } else {
        ## if the folder already exists, delete previous files with
        ## the same prefix so RAxML doesn't complain about them
        lst_raxml_suffix <- list.files(path = raxml_output,
                                    pattern = paste0(raxml_suffix, "$"),
                                    full.names = TRUE)
        file.remove(lst_raxml_suffix)
    }

    alg_phy <- fas2phy(alg_file, format = "sequential")
    alg_phy <- names(alg_phy) # get the file name

    raxml_phy <- file.path(raxml_output, basename(alg_phy))

    ## move the generated file into the data/raxml folder
    file.rename(alg_phy, raxml_phy)

    ## create partition file
    raxmlPartitionCreate(raxml_phy, file.out = part_file, overwrite = TRUE)

    ## clean up the sequence names
    alg <- ape::read.dna(file = raxml_phy, format = "sequential")
    alg <- cleanSeqLabels(alg)

    ## save the correspondance between the replaced names and the
    ## original names for later if needed.
    table_names <- data.frame(
        original_names = attr(alg, "oldnames")[[1]],
        cleaned_names = dimnames(alg)[[1]],
        stringsAsFactors = FALSE
    )
    write.csv(table_names, file = "data/raxml/table_names.csv",
              row.names = FALSE)

    ## overwrite the alignment file with the cleaned up names
    ape::write.dna(alg, file = raxml_phy, format = "sequential", colw = 99999, colsep = "")

    raxmlCmd <- paste(raxml_path, "-s", raxml_phy, "-m GTRGAMMA",
                      "-q", part_file,
                      "-f a -p 10101 -x 10101 -# 500 -n ", raxml_suffix,
                      "-T8 -w", file.path(getwd(), raxml_output))
    system(raxmlCmd)
}
