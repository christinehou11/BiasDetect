#' @rdname biasDetect
#'
#' @name featureSelect
#'
#' @title Feature Selection
#'
#' @description Function to conduct feature selection on snRNA-seq and spatial
#'     transcriptomics data and calculate the difference of deviance and rank
#'     with and without the selected batch effect.
#'
#' @importFrom scry devianceFeatureSelection
#' @importFrom dplyr left_join filter
#' @importFrom SummarizedExperiment colData
#' @importFrom stats sd
#'
#' @param input  \code{SpatialExperiment} or \code{SingleCellExperiment}: Input
#'     data which can be either \code{SpatialExperiment} or
#'     \code{SingleCellExperiment}. It is assumed to have an \code{assays} slot
#'     containing \code{counts} assay for \code{devianceFeatureSelection()} to
#'     successfully operate and calculate the deviance and rank. The
#'     \code{logcounts} is strongly recommended to be included in \code{assays}.
#'     Also, the input data is assumed to set up \code{rownames(input)} as
#'     \code{genes} and \code{gene_name} in \code{rowData(input)}.
#'
#' @param batch_effect \code{character}: Either batch on \code{slide} or
#'     \code{subject} based on what metadata is available within input data.
#'     The name of \code{slide} or \code{subject} within each input object
#'     should be specified based on different scenarios.
#'
#' @param VGs \code{character}: Highly Variable Genes (HVGs) for
#'     \code{SingleCelleExperiment} object or Spatially Variable Genes (SVGs)
#'     for \code{SpatialExperiment} object. If it is a data frame, it is assumed
#'     to contain one column of identified variable genes with column name as
#'     "gene_name".
#'
#' @return If the input was provided as a \code{SpatialExperiment} or
#'     \code{SingleCellExperiment} object, the output values are returned as
#'     a data frame containing the deviance and rank with and without the bacth
#'     effect, the corresponding difference, the corresponding nSD, gene, gene
#'     name, and whether the gene is outlier defined by the chosen deviance and
#'     rank cutoff.
#'
#' @export
#'
#' @examples
#' # library(spatialLIBD)
#' # spe <- fetch_data(type = "spe")
#'
#' # featureSelect(spe, batch_effect = "brain", VGs = SVGs)
#'
featureSelect <- function(input, batch_effect = NULL, VGs = NULL) {

    stopifnot(.check_feature(input, batch_effect, VGs))

    if (!is.null(batch_effect)) {
        if (!batch_effect %in% names(colData(input))) {
            stop("The batch_effect is not a valid column")}
        batch_effect <- colData(input)[[batch_effect]]
    } else { stop("Please provide a valid batch_effect.")}

    message("Step 1: Running feature selection without batch...")
    bd <- devianceFeatureSelection(input, assay = "counts", fam = "binomial")
    bd_df <- .df_maker(bd)

    message("Step 2: Running feature selection with batch...")
    bd_batch <- devianceFeatureSelection(input, assay = "counts",
            fam = "binomial",batch = as.factor(batch_effect))
    bd_batch_df <- .df_maker(bd_batch)

    message("Step 3: Calculating deviance and rank difference...")
    batch_df <- left_join(bd_df, bd_batch_df, by=c("gene", "gene_name"),
                            suffix=c("_default","_batch"))
    batch_df <- filter(batch_df, .data$gene_name %in% VGs)

    batch_df$d_diff <- (batch_df$dev_default-batch_df$dev_batch)/
                        batch_df$dev_batch
    batch_df$r_diff <- batch_df$rank_batch-batch_df$rank_default

    batch_df
    }

#' @rdname biasDetect
#'
#' @name biasDetect
#'
#' @title Find Bias Genes
#'
#' @description Function to filter out the bias genes based on selected
#'    selected threshold of nSD of deviance and rank in data frame or plot.
#'
#' @importFrom dplyr filter
#' @importFrom tibble is_tibble
#' @importFrom rlang .data
#' @importFrom ggplot2 scale_x_log10 scale_y_reverse
#'            scale_y_log10 aes geom_abline
#' @importFrom gridExtra grid.arrange
#'
#' @param batch_df \code{data.frame} : Input data frame generated from
#'    `featureSelection()` function using \code{SpatialExperiment} or
#'    \code{SingleCellExperiment} object containing the raw data.
#'
#' @param nSD_dev \code{integer}: Number of standard deviation (nSD) on deviance
#'     difference. The default value is 5.
#'
#' @param nSD_rank \code{integer}: Number of standard deviation (nSD) on rank
#'     difference. The default value is 5.
#'
#' @param visual \code{logical}: Whether to display the detected bias genes by
#'     visualizations of \code{deviance} and \code{rank}. Default = FALSE to
#'     return the data frame. If it is TRUE, the returned format will be
#'     two parallel plots presenting bias genes based on both \code{nSD_dev} and
#'     \code{nSD_rank}.
#'
#' @return The output values will be the identified bias genes. The returned
#'     format can be either a data frame or a list with two elements, giving
#'     option to see either \code{deviance} or \code{rank} plot.
#'
#' @export
#'
#' @examples
#' # library(spatialLIBD)
#'
#' # batch_df <- featureSelect(spe, batch_effect = "brain", VGs = SVGs)
#' # biasDetect(batch_df, nSD_dev = 5, nSD_rank = 6)
#' # biasDetect(batch_df, nSD_dev = 5, nSD_rank = 6, visual = TRUE)
#'

biasDetect <- function(batch_df, nSD_dev = 5, nSD_rank = 5, visual = FALSE) {

    stopifnot(.check_bias(batch_df, nSD_dev, nSD_rank, visual))

    batch_df <- .sd_cutoff(batch_df, nSD_dev, "nSD_dev","dev")
    batch_df <- .sd_cutoff(batch_df, nSD_rank, "nSD_rank", "rank")

    biased.genes.df <- filter(batch_df,
                            .data$dev_outlier==TRUE | .data$rank_outlier==TRUE)

    if (visual == TRUE) {
        plot_dev <- .plot(batch_df,
                            "dev_default", "dev_batch", nSD_dev,"dev") +
            scale_x_log10() + scale_y_log10() +
            geom_abline(aes(slope = 1, intercept = 0), lty = 2)

        plot_rank <- .plot(batch_df,
                            "rank_default", "rank_batch", nSD_rank, "rank") +
            scale_y_reverse() +
            geom_abline(aes(slope = -1, intercept = 0), lty = 2)

        output <- list(deviance = plot_dev, rank = plot_rank)


    }
    else {
        output <- biased.genes.df
    }

    output
    }


