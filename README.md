
# File organization

The file dependency in the `R/` folder are loaded in this order (--> = is loaded
by):

packages.R --> build.R --> load.R --> build_*.R

# BEAST analyses

## 20140809.isthmus

* used nexus file generated by `build_isthmus_alg()`
* linked clocks and trees to have 1 clock and 1 tree
* used GTR+G+I for each codon positions
* log normal clock set at 1
