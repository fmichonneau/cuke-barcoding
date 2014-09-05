source("R/parseRscriptArgs.R")
args <- commandArgs(TRUE)
args <- parseRscriptArgs(args)

f <- args$file

fcont <- readLines(f)
fcont <- gsub("caption\\[.+\\]", "caption", fcont)
cat(fcont, sep="\n", file=gsub("\\.tex", "_noCapt.tex", f))
