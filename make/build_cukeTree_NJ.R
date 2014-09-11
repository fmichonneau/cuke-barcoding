
source("R/parseRscriptArgs.R")
args <- commandArgs(TRUE)
args <- parseRscriptArgs(args)

model <- args$model
Nrep <- args$Nrep

set.seed(10101)

source("R/load.R")

model <- match.arg(model, c("raw", "K80"))

if (identical(model, "raw")) {
    load_cukeTree_raw(overwrite=TRUE, Nrep=Nrep)
    load_cukeTree_raw_phylo4(overwrite=TRUE)
} else if (identical(model, "K80")) {
    load_cukeTree_k2p(overwrite=TRUE, Nrep=Nrep)
    load_cukeTree_k2p_phylo4(overwrite=TRUE)
}
