#' Human Hippocampus .
#'
#' The data was generated from spatially-resolved transcriptomics (SRT) and
#' single-nucleus RNA-sequencing (snRNA-seq) data from adjacent tissue sections
#' of the anterior human hippocampus across ten adult neurotypical donors.
#' SRT data was generated using 10x Genomics Visium (n=36 capture areas) and
#' 10x Genomics Visium Spatial Proteogenomics (SPG) (n=8 capture areas).
#' snRNA-seq data was generated using 10x Genomics Chromium (n=26 total
#' snRNA-seq libraries). The original data was processed using `nnSVG()` feature
#' selection, and `spe` dataset is the `SpatialExperiment` object
#' containing the 2098 spatially variable genes (SVGs).
#'
#'
#' @docType data
#'
#' @usage data(spe)
#'
#' @format An SpatialExperiment object.
#' \describe{}
#'
#' @keywords datasets
#'
#' @source \href{https://github.com/LieberInstitute/spatial_hpc}{spatialHPC}
#'
#' @examples
#' data(spe)
"spe"
