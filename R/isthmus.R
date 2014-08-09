
build_isthmus_alg <- function() {
    isthmusId <- read.table(file="data/isthmus-id.txt", header=FALSE)
    cukeAlg <- load_cukeAlg()
    mtchId <- lapply(isthmusId[, 1], function(x) {
        grep(paste0("_", x, "(_|$)"), dimnames(cukeAlg)[[1]])
    })
    lId <- sapply(mtchId, length)
    stopifnot(identical(as.character(isthmusId[which(lId > 1), 1]),
                        "UF9650"))
    stopifnot(all(lId != 0))
    mtchId <- sapply(mtchId, function(x) x[1])
    isthmusAlg <- cukeAlg[mtchId, ]
    ape::write.dna(isthmusAlg, file="data/cukeBarcodes-isthmus.phy",
                   format="sequential", colsep="", colw=dim(cukeAlg)[2])
    nexusPartitionCreate(alg="data/cukeBarcodes-isthmus.phy",
                         file.out="tmp/cukeBarcodes-isthmus.part",
                         overwrite=TRUE)
    alg2nex(file="data/cukeBarcodes-isthmus.phy", partition.file="tmp/cukeBarcodes-isthmus.part")
}
