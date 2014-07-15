
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

#data/raw/MARBoL_Echinos_VIII_2013.csv: make/build_cukeBarcodesCSV.R ${DATA_RAW} 
#	${RSCRIPT} $<

data/raw/cukeBarcodes.csv.rds: ${DATA_PREPROCESSED}
	${RSCRIPT} -e "saveRDS(read.csv(file='$<', stringsAsFactors=FALSE), file='data/raw/cukeBarcodes.csv.rds')"

data/cukeBarcodes-aligned.fas: make/build_alignedFasta.R data/raw/cukeBarcodes.csv.rds R/genFasta.R
	${RSCRIPT} -e "if (!file.exists('tmp/')) dir.create('tmp/');"
	${RSCRIPT} $<

data/cukeBarcodes-cleaned.fas: make/build_cleanedFasta.R data/cukeBarcodes-aligned.fas
	${RSCRIPT} $<

data/cukeBarcodes-flagAmb.rds: make/build_flagAmbiguities.R data/cukeBarcodes-cleaned.fas
	${RSCRIPT} $<

data/cukeTree-k2p.rds: make/build_cukeTree_nj.R data/cukeBarcodes-flagAmb.rds
	${RSCRIPT} $< "file.in='data/cukeBarcodes-flagAmb.rds', model='K80', file.out='$@', Nrep=200"

data/cukeTree-raw.rds: make/build_cukeTree_nj.R data/cukeBarcodes-flagAmb.rds
	${RSCRIPT} $< "file.in='data/cukeBarcodes-flagAmb.rds', model='raw', file.out='$@', Nrep=200"

data/cukeBarcodes-flagAmb.phy: data/cukeBarcodes-flagAmb.rds
	${RSCRIPT} -e "library(ape); write.dna(readRDS('$<'), format='sequential', colw=1000, file='data/cukeBarcodes-flagAmb.phy')"

data/cukeBarcodes-raxml.tre: data/cukeBarcodes-flagAmb.phy
	${RSCRIPT} -e "library(seqManagement); raxmlPartitionCreate('$<', file.out='data/cukeBarcodes-partition', overwrite=TRUE)"
	raxmlHPC-PTHREADS-SSE3 -s $< -m GTRGAMMA -q data/cukeBarcodes-partition -T8 -f a -p 10101 -x 10101 -# 500 -n cukeBarcodes
	cp tmp/RAxML_bipartitions.cukeBarcodes $@

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

