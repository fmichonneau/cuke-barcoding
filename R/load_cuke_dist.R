load_cuke_dist <- function(alg, model) {
    model <- match.arg(model,  c("raw", "K80"))
    cukeDist <- ape::dist.dna(alg, as.matrix=TRUE, model=model,
                              pairwise.deletion=FALSE)
    cukeDist
}

load_cuke_dist_raw <- function(cuke_alg) {
   load_cuke_dist(cuke_alg, model = "raw")
}

load_cuke_dist_k2p <- function(cuke_alg) {
    load_cuke_dist(cuke_alg, model = "K80")
}
