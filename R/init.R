initialize <- function() {
    cmptr_name <- system("uname -n", intern = TRUE)
    if (identical(cmptr_name, "francois-laptop") ||
        identical(cmptr_name, "francois-XPS-15-9560")) {
        options("mc.cores" = 8)
    } else if (identical(cmptr_name, "ryanlab.whitney.ufl.edu")) {
        options("mc.cores" = 12)
    } else {
        options("mc.cores" = 1)
        message("no parallelization")
    }
    message("number of cores ", getOption("mc.cores"))
    if (!dir.exists("tmp"))
        dir.create("tmp")
}
initialize()
