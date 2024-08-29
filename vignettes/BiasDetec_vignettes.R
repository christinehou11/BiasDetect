## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
    collapse = TRUE,
    comment = "#>"
)

## ----'install dev', eval = FALSE----------------------------------------------
#  remotes::install("christinehou11/BiasDetect")

## ----'pkg install', eval=FALSE------------------------------------------------
#  pkgs <- c("dplyr", "tidyr", "scran", "scry", "nnSVG")
#  required_pkgs <- pkgs[!pkgs %in% rownames(installed.packages())]
#  BiocManager::install(required_pkgs)

## ----'library', message=FALSE-------------------------------------------------
library(BiasDetect)

