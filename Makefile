RSCRIPT = Rscript

## Manuscript
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

### clean

clean-tmp:
	-rm -f tmp/*.*
