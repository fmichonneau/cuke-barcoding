
parseRscriptArgs <- function(args) {
    args <- paste0("list(", args, ")")
    args <- eval(parse(text=args))
    args
}
