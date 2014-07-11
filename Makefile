
RSCRIPT = Rscript

DATA_RAW = data/raw/MARBoL_Echinos_VIII_2013.xlsx

DATA_PREPROCESSED = data/raw/MARBoL_Echinos_VIII_2013.csv

DATA_PROCESSED_FASTA = data/cukeBarcodes.fas

DATA_PROCESSED = data/cukeBarcodes.rds

### Initialization
init: init-tmp

init-tmp:
	${RSCRIPT} -e "if (!file.exists('tmp/')) dir.create('tmp/'); if (!file.exists('figures/')) dir.create('figures/')"

### Data preparation

data-raw: ${DATA_RAW}

data-preprocessed: ${DATA_PREPROCESSED}

data-processed: ${DATA_PROCESSED}

data/raw/MARBoL_Echinos_VIII_2013.csv: make/build_cukeBarcodesCSV.R ${DATA_RAW} 
	${RSCRIPT} $<

data/raw/cukeBarcodes.csv.rds: ${DATA_PREPROCESSED}
	${RSCRIPT} -e "saveRDS(read.csv(file='$<', stringsAsFactors=FALSE), file='data/raw/cukeBarcodes.csv.rds')"

data/cukeBarcodes-aligned.fas: make/build_alignedFasta.R data/raw/cukeBarcodes.csv.rds R/genFasta.R
	${RSCRIPT} -e "if (!file.exists('tmp/')) dir.create('tmp/');"
	${RSCRIPT} $<

data/cukeBarcodes-cleaned.fas: make/build_cleanedFasta.R data/cukeBarcodes-aligned.fas
	${RSCRIPT} $<

data/cukeBarcodes-flagAmb.rds: make/build_flagAmbiguities.R data/cukeBarcodes-cleaned.fas
	${RSCRIPT} $<

data/cukeTree-k2p.rds: make/build_cukeTree.k2p.rds.R data/cukeBarcodes-flagAmb.rds
	${RSCRIPT} $<

data/cukeTree-raw.rds: make/build_cukeTree.raw.rds.R data/cukeBarcodes-flagAmb.rds
	${RSCRIPT} $<

### figures

figures/cukeTree-k2p.pdf: data/cukeTree-k2p.rds R/plot.pdf.tree.R
	${RSCRIPT} -e "if (!file.exists('figures/')) dir.create('figures/');"
	${RSCRIPT} -e "source('R/plot.pdf.tree.R'); plot.pdf.tree('figures/cukeTree-k2p.pdf', treeRDS='data/cukeTree-k2p.rds', height=350, width=50)"

figures/cukeTree-raw.pdf: data/cukeTree-raw.rds R/plot.pdf.tree.R
	${RSCRIPT} -e "if (!file.exists('figures/')) dir.create('figures/');"
	${RSCRIPT} -e "source('R/plot.pdf.tree.R'); plot.pdf.tree('figures/cukeTree-raw.pdf', treeRDS='data/cukeTree-raw.rds', height=350, width=50)"


### clean

clean-tmp:
	-rm -f tmp/*.*

