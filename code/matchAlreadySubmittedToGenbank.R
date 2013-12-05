
### Sequences obtained with the following query in genbank:
### (holothuroidea[Organism]) AND (cytochrome OR COI[gene name]) AND (michonneau OR paulay OR kim OR werner OR Honey-Escandon)

### as of 2013-12-04 -- 605 sequences

ingb <- read.dna(file="/tmp/sequence.fasta", format="fasta", as.character=TRUE)

vposbtmp <- gregexpr("UF|MOL|FLMNH|WAM|NMV|GP|FRM|[ANG]([0-9]{2,3})\\s+|MRAC|GZG|BMOO|GUSTAV", names(ingb))
vposetmp <- gregexpr("cytochrome", names(ingb))

vposb <- lapply(vposbtmp, function(x) x[1])
vpose <- lapply(vposetmp, function(x) x[1])

vpos <- character(nrow(dtgb))
for (i in 1:nrow(dtgb)) {
    vpos[i] <- substr(names(ingb)[i], vposb[[i]], vpose[[i]]-2)
}

gbnb <- sapply(names(ingb), function(x) unlist(strsplit(x, "\\|"))[4])


dtgb <- data.frame(seqname=names(ingb), genbankNb=gbnb, voucherInfo=vpos,
                   seq=sapply(ingb, function(x) paste(x, collapse="")),
                   stringsAsFactors=FALSE)


write.csv(dtgb, file="/tmp/subgb.csv")

### After some manual edits to standardize voucher info
###   trying to match the (disparate) voucher info in
###   Genbank with the ones in the marbol spreadsheet

matchVoucher <- function(file="~/Dropbox/CukeBarcoding/GenbankMatch.csv",
                         seqdbfile="~/Documents/CukeBarcoding/data/MARBoL_Echinos_VIII_2013.csv") {

    vdb <- read.csv(file, stringsAsFactors=FALSE)
    seqdb <- read.csv(seqdbfile, stringsAsFactors=FALSE)

    ext <- paste(seqdb$Plate.Extraction.Genbank, seqdb$Cell, sep="_")
    ext <- gsub("_+$", "", ext)

    catnb <- paste(seqdb$Collection.Code, seqdb$Catalog_number, sep="")
    catnb <- sapply(catnb, function(x) gsub("\\([0-9]+\\)", "", x))

    samp <- gsub("_", "", seqdb$Sample)
     
    res <- vou <- character(nrow(vdb))

    for (i in 1:nrow(vdb)) {
        vinf <- vdb$voucherInfo_edited[i]
        if (vdb$typeOfVoucherInfo[i] == "Extraction") {            
            if (length(grep("FLMNH|GUSTAV", vinf)) > 0) {
                vinf <- gsub("^GUSTAV", "FLMNH", vinf)
                vinf <- unlist(strsplit(vinf, "_"))
                vinf[1] <- add.zeros(vinf[1], 3, prefix.length=5,
                                     quiet=TRUE)
                vinf <- paste(vinf, collapse="_")
                vinf <- gsub("^FLMNH", "FLMNH_", vinf)
                xx <- match(vinf, ext)
                if (is.na(xx)) message("no match for ", vinf)
                res[i] <- seqdb$Counter[xx]
                vou[i] <- vinf
            }
            else if (length(grep("^UF", vinf)) > 0) {
                vinf <- gsub("^UF\\-", "", vinf)
                xx <- match(vinf, ext)
                if (is.na(xx)) message("no match for ", vinf)
                res[i] <- seqdb$Counter[xx]
                vou[i] <- vinf
            }
            else {
                xx <- match(vinf, ext)
                if (is.na(xx)) message("no match for ", vinf)
                res[i] <- seqdb$Counter[xx]
                vou[i] <- vinf
            }
        }
        else if (vdb$typeOfVoucherInfo[i] == "Catalog") {
            vinf <- gsub("_|\\s+", "", vinf)
            vinf <- gsub("UFE", "UF", vinf)
            vinf <- gsub("UF0", "UF", vinf)
            xx <- match(vinf, catnb)
            if (is.na(xx)) message("no match for ", vinf)
            res[i] <- seqdb$Counter[xx]
            vou[i] <- vinf           
        }
        else if (vdb$typeOfVoucherInfo[i] == "Sample") {
            vinf <- gsub("AF", "AF0", vinf)
            xx <- match(vinf, samp)
            if (is.na(xx)) message("no match for ", vinf)
            res[i] <- seqdb$Counter[xx]
            vou[i] <- vinf
        }
        else message("no match for ", vinf)
    }
    data.frame(vou, res, vdb$seq)
}


tt <- matchVoucher()
write.csv(tt, file="toEdit.csv")

###

tt <- read.csv(file="toEdit.csv", stringsAsFactors=FALSE)

seqdb <- read.csv(file="data/MARBoL_Echinos_VIII_2013.csv", stringsAsFactors=FALSE)

tmpRes <- merge(seqdb, tt, by.x="Counter", by.y="res", all.x=TRUE)

write.csv(tmpRes, file="seqDB_withGenbankInfo.csv")
