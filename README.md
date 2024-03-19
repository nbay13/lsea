# lsea

`lsea` is an R package for perfoming Lipid Structure Enrichment Analysis (LSEA) of lipid classes, chain lengths, and degrees of unsaturation.

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
See vignette for a tutorial using example data to perform differential composition testing and Lipid Structure Enrichment Analysis (LSEA)

https://nbay13.github.io/lsea/
