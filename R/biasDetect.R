#' biasDetect
#'
#' Function to run 'BiasDetect' method for spatially variable genes (SVGs)
#' identification on spatial-resolved transcriptomics (SRT) data with some types
#' of batch effect.
#'
#' `BiasDetect` is a feature-based Quality Control (QC) to identify SVGs on
#' spatial transcriptomics data with sometypes of batch effect. Regarding to the
#' spatial transcriptomics data experients, the batch can be defined as
#' "slide" or "subject".The `BiasDetect` method is based on binomial deviance
#' model (Townes et al, 2019) and applies cutoffs based on the number of
#' standard deviation (nSD) of deviance and rank difference as the data-driven
#' thresholding approach to detect the batch-biased features.
#'
#'
#' @param input Input data: assumed to be formatted as a
#'     \code{SpatialExperiment} object or \code{SingleCellExperiment} with an
#'     \code{assays} slot named \code{counts} containing raw expression counts.
#'
#'
#' @param batch_effect Batch effect: either \code{slide} or \code{subject} based
#'     on what metadata is available within input data. The name of \code{slide}
#'     or \code{subject} within each input object can be different.
#'
#' @param nSD_dev number of standard deviation (nSD) on deviance difference
#'
#' @param nSD_rank number of standard deviation (nSD) on rank difference
#'
#' @param visual visualization of nSD on deviance or rank difference. Default F
#'
#' @return If the input was provided as a \code{SpatialExperiment} object, the
#'   output values are returned as additional columns in the \code{rowData} slot
#'   of the input object. If the input was provided as a \code{numeric} matrix
#'   of values, the output is returned as a \code{numeric} matrix. The output
#'   values include spatial variance parameter estimates, likelihood ratio (LR)
#'   statistics, effect sizes (proportion of spatial variance), p-values, and
#'   multiple testing adjusted p-values.
#'
#'
#' @importFrom scry devianceFeatureSelection
#' @importFrom dplyr left_join filter
#' @importFrom SummarizedExperiment rowData
#' @importFrom ggplot2 ggplot geom_point scale_color_manual
#'            geom_abline labs theme theme_bw scale_x_log10 scale_y_reverse
#'            scale_y_log10 aes
#' @importFrom gridExtra grid.arrange
#' @importFrom RColorBrewer brewer.pal
#' @importFrom stats sd
#' @importFrom rlang .data
#'
#' @export
#'
#' @examples
biasDetect <-
    function(input, batch_effect, nSD_dev = 5, nSD_rank = 5, visual = FALSE) {

    stopifnot(

    )

    spe <- input

    message("Step 1: Running feature selection without batch...")
    bd <- devianceFeatureSelection(input, assay = "counts", fam = "binomial")
    bd.df <- .df_maker(bd)

    message("Step 2: Running feature selection with batch...")
    bd.batch <- devianceFeatureSelection(input, assay = "counts",
                                        fam = "binomial",
                                        batch = as.factor(batch_effect))
    bd.batch.df <- .df_maker(bd.batch)

    message("Step 3: Merging data frames...")
    batch.df <- left_join(bd.df, bd.batch.df,
                            by=c("gene", "gene_name"),
                            suffix=c("_default","_batch"))

    message("Step 4: Calculating deviance and rank difference...")
    batch.df$d.diff <- (batch.df$dev_default-batch.df$dev_batch)/
                        batch.df$dev_batch
    batch.df$r.diff <- batch.df$rank_batch-batch.df$rank_default

    message("Step 6: SD cutoff: deviance")
    mean_deviance <- mean(batch.df$d.diff)
    sd_deviance <- sd(batch.df$d.diff)
    batch.df$nSD_dev <- (batch.df$d.diff-mean_deviance)/sd_deviance
    batch.df$dev_outlier <- batch.df$nSD_dev >= nSD_dev

    message("Step 7: SD cutoff: rank")
    mean_rank <- mean(batch.df$r.diff)
    sd_rank <- sd(batch.df$r.diff)
    batch.df$nSD_rank <- (batch.df$r.diff-mean_rank)/sd_rank
    batch.df$rank_outlier <- batch.df$nSD_rank >= nSD_rank

    message("Step 8: Cluster results")
    biased.genes <- filter(batch.df,
                        .data$dev_outlier==TRUE | .data$rank_outlier==TRUE)$gene
    names(biased.genes) <- rowData(spe)[biased.genes,"gene_name"]

    message("Function execution complete. The bias genes are:")

    if (visual == TRUE) {

        col.pal <- brewer.pal(length(unique(batch.df$nSD.bin_dev)), "YlOrRd")
        col.pal[1] <- "grey"
        sd.interval <- nSD_dev
        batch.df$nSD.bin_dev <- cut(abs(batch.df$nSD_dev), right=FALSE,
                                breaks=seq(0,
                                            max(batch.df$nSD_dev)+sd.interval,
                                            by=sd.interval),
                                include.lowest=TRUE)
        dev <- ggplot(batch.df,
                        aes(x=.data$dev_default,
                            y=.data$dev_batch,
                            color=.data$nSD.bin_dev)) +
                geom_point()+
                scale_x_log10()+
                scale_y_log10()+
                scale_color_manual(values=col.pal)+
                geom_abline(aes(slope=1, intercept=0), lty=2)+
                labs(x="deviance (no batch)",
                        y="deviance (batch)")+
                theme_bw() +
                theme(legend.position="none")

        sd.interval <- nSD_rank
        batch.df$nSD.bin_rank <- cut(abs(batch.df$nSD_rank), right=FALSE,
                                    breaks=seq(0,
                                            max(batch.df$nSD_dev)+sd.interval,
                                            by=sd.interval),
                                    include.lowest=TRUE)
        col.pal2 <- brewer.pal(length(unique(batch.df$nSD.bin_rank)),
                                "YlOrRd")
        col.pal2[1] <- "grey"
        rank <- ggplot(batch.df,
                        aes(x=.data$rank_default,
                            y=.data$rank_batch,
                            color=.data$nSD.bin_rank))+
                geom_point()+
                scale_y_reverse()+
                scale_color_manual(values=col.pal2)+
                geom_abline(aes(slope=-1, intercept=0), lty=2)+
                labs(x="rank (no batch)",
                        y="rank (batch)")+
                theme_bw()+
                theme(legend.position="none")

        output <- grid.arrange(dev, rank, ncol=2)
    }
    else {
        output <- biased.genes
    }

    output

    }
