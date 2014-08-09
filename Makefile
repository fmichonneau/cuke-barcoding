
RSCRIPT = Rscript

### Manuscript
all: cuke-barcoding.tex cuke-barcoding_nourl.bib clean-partial
	-xelatex -interaction=nonstopmode "\input" cuke-barcoding.tex
	-bibtex cuke-barcoding
	-xelatex -interaction=nonstopmode "\input" cuke-barcoding.tex
	xelatex -interaction=nonstopmode "\input" cuke-barcoding.tex

cuke-barcoding.tex: cuke-barcoding.Rnw code/cuke-barcoding.R
	${RSCRIPT} -e "library(knitr); knit('cuke-barcoding.Rnw')"

cuke-barcoding_nourl.bib: cuke-barcoding.Rnw
	-cp ~/Library/Barcoding.bib cuke-barcoding.bib
	${RSCRIPT} parseURLs.R

clean-partial:
	-rm *.bbl
	-rm *.blg
	-rm *.aux
	-rm *.log
	-rm *~

clean: clean-partial
	-rm impatiens_phylogeography.pdf
	-rm impatiens_phylogeography.tex

### Data preparation

#data/raw/MARBoL_Echinos_VIII_2013.csv: make/build_cukeBarcodesCSV.R data/raw/MARBoL_Echinos_VIII_2013.xlsx
#	${RSCRIPT} $<

data/raw/cukeBarcodes.csv.rds: data/raw/MARBoL_Echinos_VIII_2013.csv
	${RSCRIPT} -e "saveRDS(read.csv(file='$<', stringsAsFactors=FALSE), file='$@')"

data/cukeBarcodes-aligned.fas: make/build_alignedFasta.R data/raw/cukeBarcodes.csv.rds R/genFasta.R
	${RSCRIPT} -e "if (!file.exists('tmp/')) dir.create('tmp/');"
	${RSCRIPT} $<

data/cukeBarcodes-cleaned.fas: make/build_cleanedFasta.R data/cukeBarcodes-aligned.fas
	${RSCRIPT} $<

data/cukeBarcodes-flagAmb.rds: make/cukeBarcodes-flagAmb.rds.R data/cukeBarcodes-cleaned.fas
	${RSCRIPT} $<

data/cukeTree-k2p.rds: make/build_cukeTree_nj.R data/cukeBarcodes-flagAmb.rds
	${RSCRIPT} $< "file.in='data/cukeBarcodes-flagAmb.rds', model='K80', file.out='$@', Nrep=200"

data/cukeTree-raw.rds: make/build_cukeTree_nj.R data/cukeBarcodes-flagAmb.rds
	${RSCRIPT} $< "file.in='data/cukeBarcodes-flagAmb.rds', model='raw', file.out='$@', Nrep=200"

data/cukeBarcodes-flagAmb.phy: data/cukeBarcodes-flagAmb.rds
	${RSCRIPT} -e "library(ape); write.dna(readRDS('$<'), format='sequential', colw=1000, file='data/cukeBarcodes-flagAmb.phy')"

data/cukeBarcodes-raxml.tre: data/cukeBarcodes-flagAmb.phy ## not tested
	${RSCRIPT} -e "source('R/build.R'); build_raxml_tree('$<');"
	cp data/raxml/RAxML_bipartitions.cukeBarcodes $@

data/raxml_ptp/bPTPres: #data/cukeBarcodes-raxml.tre
	${RSCRIPT} -e "source('R/build.R'); build_PTP_tree(); build_PTP_results();"


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
