
library(ape)

plot.pdf.tree <- function(file, treeRDS, ...) {
    tr <- readRDS(treeRDS)
    pdf(file=file, ...)
    plot(tr, no.margin=TRUE, cex=0.8)
    dev.off()
}
