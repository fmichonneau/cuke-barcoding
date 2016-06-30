
distFromTip <- function(tr, node, trPost) {
### trying to use Rcpp to do the intersect doesn't make it faster
    if (missing(trPost)) {
        trPost <- reorder(tr, order="postorder")
    }
    allDesc <- descendants(tr, node, type="all")
    lDesc <- allDesc[allDesc < nTips(tr)]
    tmpD <- mclapply(lDesc, function(x) {
        pth <- intersect(allDesc, ancestors(trPost, x, type="ALL"))
        sumEdgeLength(tr, pth)
    }, mc.preschedule = FALSE)
    unlist(tmpD[which.max(tmpD)])
}

findGroups <- function(tr, threshold=.015) {
### idea for optimization see if it makes a difference to supply the tree
### both in pre- and post-order when computing the distance to tip
    stopifnot(inherits(tr, "phylo4"))
    stopifnot(!hasDuplicatedLabels(tr))

    grp <- vector("list", nTips(tr))

    ## find the distance to the tips for each internal node and select the nodes
    ## below the threshold
    intNodes <- nodeId(tr, "internal")
    trPost <- reorder(tr, order = "postorder")
    lGrp <- mclapply(intNodes, function(x) distFromTip(tr, x, trPost),
                     mc.preschedule = FALSE)
    lGrp <- unlist(lGrp)
    if (length(lGrp) != length(intNodes)) {
        summary(lGrp)
        summary(intNodes)
    }
    names(lGrp) <- intNodes

    lGrp <- lGrp[lGrp <= threshold]

    ## find all the descendants for the nodes below the threshold
    descGrp <- mclapply(as.numeric(names(lGrp)), function(x) descendants(tr, x),
                        mc.preschedule = FALSE)

    ## remove overlapping sets
    snglGrp <- mclapply(descGrp, function(x) length(x) == 1)
    snglGrp <- unlist(snglGrp)
    edgeGrp <- do.call("rbind", mclapply(descGrp, function(x) {
                                    if(length(x) > 1) cbind(head(x, -1), tail(x, -1)) else NULL
                                }))
    if (!is.null(edgeGrp)) {
        graphGrp <- graph.data.frame(edgeGrp)
        descGrp <- c(split(V(graphGrp)$name, clusters(graphGrp)$membership),
                     descGrp[snglGrp])

        ## add singletons
        missingGrps <- setdiff(nodeId(tr, "tip"), as.numeric(unlist(descGrp)))
        descGrp <- c(descGrp, as.list(missingGrps))
    } else {
        descGrp <- nodeId(tr, "tip")
    }

    ## Return tip labels
    grp <- sapply(descGrp, function(x) tipLabels(tr)[as.numeric(x)])

    ## build a phylo4d object for the results
    dTip <- data.frame(Groups=rep(1:length(grp), sapply(grp, length)))
    rownames(dTip) <- unlist(grp)
    if (inherits(tr, "phylo4d")) {
        return(addData(tr, tip.data=dTip, rownamesAsLabels = TRUE))
    } else if (inherits(tr, "phylo4"))
        return(phylo4d(tr, tip.data=dTip, rownamesAsLabels = TRUE))
}


### Example code
## spGrp <- findGroups(trP4, threshold=.025)
## spGrpCopy <- spGrp
## tipLabels(spGrp) <- paste(tipData(spGrp)$Group, tipLabels(spGrp), sep="_")

## grpLbl <- paste("^", 1:max(tipData(spGrp)$Groups), "_", sep="")

## pdf(file="treeWithBars-025.pdf", height=100, width=10)
## par(mai=c(0.5,0,2,0), xpd=T)
## plot(as(spGrp, "phylo"), cex=.5, show.tip.label=T, no.margin=F, label.offset=0)
## barMonophyletic(grpLbl, as(spGrp, "phylo"), extra.space=0.01, cex.plot=.5, cex.text=.5,
##                 bar.at.tips=TRUE, include.tip.label=TRUE)
## add.scale.bar()
## dev.off()
