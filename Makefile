
RSCRIPT = Rscript

### Manuscript
all: cuke-barcoding.tex cuke-barcoding_nourl.bib clean-partial
	-xelatex -interaction=nonstopmode "\input" cuke-barcoding.tex
	-bibtex cuke-barcoding
	-xelatex -interaction=nonstopmode "\input" cuke-barcoding.tex
	xelatex -interaction=nonstopmode "\input" cuke-barcoding.tex

docx: cuke-barcoding.tex cuke-barcoding_nourl.bib
	${RSCRIPT} removeShortCaptions.R "file='$<'"
	pandoc --to=docx cuke-barcoding_noCapt.tex -o cuke-barcoding.docx --bibliography cuke-barcoding_nourl.bib --csl methods-in-ecology-and-evolution.csl
	-rm cuke-barcoding_noCapt.tex


cuke-barcoding.tex: cuke-barcoding.Rnw code/cuke-barcoding.R
	${RSCRIPT} -e "library(knitr); knit('cuke-barcoding.Rnw')"

cuke-barcoding_nourl.bib: cuke-barcoding.Rnw ~/Library/Barcoding.bib
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

## database
data/cukeDB_noLabels.rds: data/raw/MARBoL_Echinos_VIII_2013.csv R/load_cukeDB.R
	${RSCRIPT} -e "source('R/load.R'); load_cukeDB_noLabels(overwrite=TRUE)"

## alignments
data/cukeAlg-cleaned.fas: make/build_alignedFasta.R data/cukeDB_noLabels.rds R/genFasta.R
	${RSCRIPT} $<

data/cukeAlg-flagAmb.rds: R/load_cukeAlg.R data/cukeAlg-cleaned.fas
	${RSCRIPT} -e "source('R/load.R'); load_cukeAlg(overwrite=TRUE);"

## distance matrices
data/cukeDist-raw.rds: R/load_cukeDist.R data/cukeAlg-flagAmb.rds
	${RSCRIPT} -e "source('R/load.R'); load_cukeDist_raw(overwrite=TRUE)"

data/cukeDist-k2p.rds: R/load_cukeDist.R data/cukeAlg-flagAmb.rds
	${RSCRIPT} -e "source('R/load.R'); load_cukeDist_k2p(overwrite=TRUE)"

## NJ trees
data/cukeTree-k2p.rds: make/build_cukeTree_NJ.R R/load_cukeTree_NJ.R data/cukeAlg-flagAmb.rds data/cukeDist-k2p.rds
	${RSCRIPT} $< "model='K80', Nrep=1"

data/cukeTree-k2p-phylo4.rds: data/cukeTree-k2p.rds
	#dummy recipe

data/cukeTree-raw.rds: make/build_cukeTree_NJ.R R/load_cukeTree_NJ.R data/cukeAlg-flagAmb.rds data/cukeDist-raw.rds
	${RSCRIPT} $< "model='raw', Nrep=1"

data/cukeTree-raw-phylo4.rds: data/cukeTree-raw.rds
	# dummy recipe

## cukeDB with labels (the real deal)
data/cukeDB_withLabels.rds: data/cukeDB_noLabels.rds data/cukeTree-raw-phylo4.rds R/genFasta.R R/load_cukeDB.R
	${RSCRIPT} -e "source('R/load.R'); load_cukeDB(overwrite=TRUE)"

## taxonomyDf
data/taxonomyDf.rds: data/cukeDB_withLabels.rds R/load_taxonomyDf.R
	${RSCRIPT} -e "source('R/load.R'); load_taxonomyDf(overwrite=TRUE)"

## localGap
data/localGap-manualESUs.rds: data/raw/manualESUs.csv data/cukeDist-raw.rds data/cukeDB_withLabels.rds
	${RSCRIPT} -e "source('R/load.R'); load_localGap(overwrite=TRUE)"

## clustering/distances
## TODO --- better way to deal with the 100's of target file?
data/cukeTree-raw-all-002.rds: data/cukeDB_withLabels.rds data/taxonomyDf.rds data/cukeTree-k2p-phylo4.rds data/cukeTree-raw-phylo4.rds
	${RSCRIPT} -e "source('R/load.R'); build_cukeTree_clusterGrps(overwrite=TRUE)"

## Should I just keep the file in the raxml folder and create an rds file
##   for it in the data/ folder instead?
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
