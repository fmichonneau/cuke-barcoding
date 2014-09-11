source("R/packages.R")

build_raxml_tree <- function(algFile="data/cukeBarcodes-flagAmb.phy") {
    partFile <- "data/cukeBarcodes-partition"
    raxmlDir <- "data/raxml"
    if (!file.exists(raxmlDir)) dir.create(raxmlDir)
    raxmlPartitionCreate(algFile, file.out=partFile, overwrite=TRUE)
    raxmlCmd <- paste("raxmlHPC-PTHREADS-SSE3", "-s", algFile, "-m GTRGAMMA",
                      "-q", partFile,
                      "-f a -p 10101 -x 10101 -# 500 -n cukeBarcodes",
                      "-T8 -w", file.path(getwd(), raxmlDir))
    system(raxmlCmd)
}

build_PTP_tree <- function() {
    ptpDir <- "data/raxml_ptp"
    if (!file.exists(ptpDir)) dir.create(ptpDir)
    raxmlTree <- ape::read.tree(file="data/raxml/RAxML_bestTree.cukeBarcodes")
    raxmlAlg <- ape::read.dna(file="data/cukeBarcodes-flagAmb.phy.reduced",
                              format="sequential")
    toDrop <- raxmlTree$tip.label[! raxmlTree$tip.label %in% dimnames(raxmlAlg)[[1]]]
    tree <- drop.tip(raxmlTree, toDrop)
    ape::write.tree(tree, file=file.path(ptpDir, "RAxML_bestTree_reduced.cukeBarcodes"))
}

build_PTP_results <- function(pathPTP="~/sandbox/SpeciesCounting/") {
    ptpCmd <- paste("python", paste0(pathPTP, "bPTP.py"),
                    "-t data/raxml_ptp/RAxML_bestTree_reduced.cukeBarcodes",
                    "-s 10101 -o data/raxml_ptp/bPTPres -r",
                    "-i 50000 -n 500")
    system(ptpCmd)
}
