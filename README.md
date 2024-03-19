# lsea - Lipid Structure Enrichment Analysis

An R package for perfoming enrichment analysis of lipid classes, chain lengths, and degrees of unsaturation.

## Dependencies
[compositions](https://cran.r-project.org/web/packages/compositions/index.html) <br />
[fgsea](https://bioconductor.org/packages/release/bioc/html/fgsea.html)

## Installation
```R
if(!require("devtools", quietly = TRUE))
    install.packages("devtools")

devtools::install_github("nbay13/lsea")
```
## Usage
See [vignette](https://nbay13.github.io/lsea/)
