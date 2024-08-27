## BiasDetect

Our goal is to use feature selection methods that can account for batch effects to flag SVGs that may be strongly influenced by technical factors. By examining the per-gene variance/deviance with and without a batch effect, we aim to identify features that could be interfering with the generation of clusters that correspond to known DLPFC spatial domains.

#### Installation

`BiasDetect` is a R package available in *Bioconductor* version 3.19 and later. You can install `BiasDetect` by using the following commands in R session from *Bioconductor*:

``` r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("BiasDetect")

## Check Bioconductor installation
BiocManager::valid()
```

Additionally, you can install development version from [GitHub](https://christinehou11.github.io/BiasDetect):

``` r
BiocManager::install("christinehou11/BiasDetect")
```

Install additional required packages before running package codes in vignettes.

``` r
pkgs <- c("dplyr", "tidyr")
required_pkgs <- pkgs[!pkgs %in% rownames(installed.packages())]
BiocManager::install(required_pkgs)
```

#### Tutorial

View the [tutorial](https://christinehou11.github.io/BiasDetect/articles/BiasDetec_vignettes.html) written by Christine Hou.

#### Helpful Link

View the additional documentation [Find Biased Features](https://jac-thom.github.io/findBiasedFeatures/) written by Jacqui Thompson.
