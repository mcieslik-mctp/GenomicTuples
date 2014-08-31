GenomicTuples
================================================================================

`GenomicTuples` defines general purpose containers for storing genomic tuples. 
It aims to provide functionality for tuples of genomic co-ordinates that are 
analogous to those available for genomic ranges in the [`GenomicRanges` Bioconductor package](http://www.bioconductor.org/packages/release/bioc/html/GenomicRanges.html).

__This package is in early development and requires the use of the development version of Bioconductor ([instructions here](http://bioconductor.org/developers/how-to/useDevel/)). __ It can be installed (provided it is passing the Travis CI build) with:

```R
source("http://bioconductor.org/biocLite.R")
useDevel()
biocLite(c('BiocGenerics', 'IRanges', 'Rcpp', 'Biobase', 'GenomicRanges', 'S4Vectors', 'GenomeInfoDb', 'testthat', 'knitr'))
devtools::install_github("PeteHaitch/GenomicTuples")
```

Alternatively, you may like to try using the [`packrat`](http://rstudio.github.io/packrat/) private library that is included in this git repository, which includes the necessary packages for using the `GenomicTuples` package.