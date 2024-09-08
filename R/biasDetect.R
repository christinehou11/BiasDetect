#' @rdname biasDetect
#'
#' @name biasDetect
#'
#' @title Find Bias Genes
#'
#' @description Function to run 'BiasDetect' method for spatially variable
#'     genes (SVGs) identification on spatial-resolved transcriptomics (SRT)
#'     data with some types of batch effect.
#'
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
#'
#' @param batch_effect \code{character}: either batch on \code{slide} or
#'     \code{subject} based on what metadata is available within input data.
#'     The name of \code{slide} or \code{subject} within each input object
#'     should be specified based on different scenarios.
#'
#'
#' @param nSD_dev \code{integer}: number of standard deviation (nSD) on deviance
#'     difference. The default value is 5.
#'
#' @param nSD_rank \code{integer}: number of standard deviation (nSD) on rank
#'     difference. The default value is 5.
#'
#' @param VGs \code{character}: highly variable genes (HVGs) for
#'     \code{SingleCelleExperiment} object or spatially variable genes (SVGs)
#'     for \code{SpatialExperiment} object. If it is a data frame, it is assumed
#'     to contain one column of identified variable genes with column name as
#'     "gene_name".
#'
#' @param visual \code{logical}: Whether to display the detected bias genes by
#'     visualizations of \code{deviance} and \code{rank}. Default = FALSE to
#'     return the data frame. If it is TRUE, the returned format will be
#'     two parallel plots presenting bias genes based on both \code{nSD_dev} and
#'     \code{nSD_rank}.
#'
#' @return If the input was provided as a \code{SpatialExperiment} or
#'     \code{SingleCellExperiment} object, the output values are returned as
#'     either a data frame or two parallel plots to present the identified bias
#'     genes.
#'
#'
#' @importFrom scry devianceFeatureSelection
#' @importFrom dplyr left_join filter
#' @importFrom SummarizedExperiment rowData colData
#' @importFrom utils read.csv
#' @importFrom stats sd
#' @importFrom ggplot2 scale_x_log10 scale_y_reverse
#'            scale_y_log10 aes geom_abline
#' @importFrom gridExtra grid.arrange
#'
#' @export
#'
#' @examples
#' \donttest{
#' set.seed(123)
#' spe <- load("inst/data/spe_norm.Rdata")
#' libd <- read.csv("inst/data/libd-all_nnSVG_p-05-features-df.csv/")
#' SVGs <- libd$gene_name
#' biasDetect(spe, batch_effect = "brain", VGs = SVGs, nSD_dev = 5,
#'         nSD_rank = 6, visual = FALSE)
#'  biasDetect(spe, batch_effect = "brain", VGs = SVGs, nSD_dev = 5,
#'         nSD_rank = 6, visual = TRUE)
#' }
biasDetect <- function(input, batch_effect = NULL, VGs = NULL, nSD_dev = 5,
                        nSD_rank = 5, visual = FALSE) {

    stopifnot(.check(input, batch_effect, VGs, nSD_dev, nSD_rank))

    if (!is.null(batch_effect)) {
        if (!batch_effect %in% names(colData(input))) {
            stop("The batch_effect is not a valid column")}
        batch_effect <- colData(input)[[batch_effect]]
    } else { stop("Please provide a valid batch_effect.")}

    message("Step 1: Running feature selection without batch...")
    bd <- devianceFeatureSelection(input, assay = "counts", fam = "binomial")
    bd.df <- .df_maker(bd)

    message("Step 2: Running feature selection with batch...")
    bd.batch <- devianceFeatureSelection(input, assay = "counts",
            fam = "binomial",batch = as.factor(batch_effect))
    bd.batch.df <- .df_maker(bd.batch)

    message("Step 3: Calculating deviance and rank difference...")
    batch.df <- left_join(bd.df, bd.batch.df, by=c("gene", "gene_name"),
                            suffix=c("_default","_batch"))
    batch.df <- filter(batch.df, .data$gene_name %in% VGs)
    batch.df$d.diff <- (batch.df$dev_default-batch.df$dev_batch)/
                        batch.df$dev_batch
    batch.df$r.diff <- batch.df$rank_batch-batch.df$rank_default
    batch.df <- .sd_cutoff(batch.df, nSD_dev, "dev")
    batch.df <- .sd_cutoff(batch.df, nSD_rank, "rank")
    biased.genes.df <- filter(batch.df,
                    .data$dev_outlier==TRUE | .data$rank_outlier==TRUE)
    biased.gene <- biased.genes.df$gene
    names(biased.genes) <- rowData(input)[biased.genes,"gene_name"]

    if (visual == TRUE) {
        plot_dev <- .plot(batch.df, "dev_default", "dev_batch",
                        nSD_dev, "nSD_dev", "dev") +
                        scale_x_log10() + scale_y_log10() +
                        geom_abline(aes(slope = 1, intercept = 0), lty = 2)

        plot_rank <- .plot(batch.df, "rank_default", "rank_batch",
                        nSD_rank, "nSD_rank", "rank") +
                        scale_y_reverse() +
                        geom_abline(aes(slope = -1, intercept = 0), lty = 2)
        output <- grid.arrange(plot_dev, plot_rank, ncol=2)
    } else {
        output <- biased.genes.df }
    output
    }
