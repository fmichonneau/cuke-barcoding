
RSCRIPT = Rscript

DATA_RAW = data/raw/MARBoL_Echinos_VIII_2013.xlsx

DATA_PREPROCESSED = data/raw/MARBoL_Echinos_VIII_2013.csv

DATA_PROCESSED_FASTA = data/cukeBarcodes.fas

DATA_PROCESSED = data/cukeBarcodes.rds

### Initialization
init: init-tmp

init-tmp:
	${RSCRIPT} -e "if (!file.exists('tmp/')) dir.create('tmp/')"

### Data preparation

data-raw: ${DATA_RAW}

data-preprocessed: ${DATA_PREPROCESSED}

data-processed: ${DATA_PROCESSED}

data/raw/MARBoL_Echinos_VIII_2013.csv: make/build_cukeBarcodesCSV.R ${DATA_RAW} 
	${RSCRIPT} $<

data/raw/cukeBarcodes.csv.rds: ${DATA_PREPROCESSED}
	${RSCRIPT} -e "saveRDS(read.csv(file='$<', stringsAsFactors=FALSE), file='data/raw/cukeBarcodes.csv.rds')"

data/cukeBarcodes-aligned.fas: make/build_alignedFasta.R data/raw/cukeBarcodes.csv.rds R/genFasta.R init-tmp
	${RSCRIPT} $<

data/cukeBarcodes-cleaned.fas: make/build_cleanedFasta.R data/cukeBarcodes-aligned.fas
	${RSCRIPT} $<

data/cukeBarcodes-flagAmb.rds: make/build_flagAmbiguities.R data/cukeBarcodes-cleaned.fas
	${RSCRIPT} $<

### clean

clean-tmp:
	-rm -f tmp/*.*

