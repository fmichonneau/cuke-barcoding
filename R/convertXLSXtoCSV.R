
######## Converts XLSX spreadsheet into CSV
### See here for more info https://github.com/dagwieers/unoconv
convertXLSXtoCSV <- function(file) {
    cmd <- paste("unoconv -f csv", file)
    if (file.exists(file)) {
        target <- gsub("\\.xlsx$", ".csv", file)
        if (file.exists(target))
            file.remove(target)
    }
    else {
        stop(file, "cannot be found.")
    }
    system("unoconv -l&") # start listener
    system("sleep 0.5;")
    system(cmd) # converts document
    system("pkill unoconv") # kill process
    if (file.exists(target)) {
        message("file conversion worked ", target, " was created")
    }
    else {
        stop("Something went wrong.")
    }
}
