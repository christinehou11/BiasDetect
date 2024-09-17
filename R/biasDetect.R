#' @rdname biasDetect
#'
#' @name biasDetect
#'
#' @title Find Bias Genes
#'
#' @description Function to filter out the bias genes based on selected
#'    selected threshold of nSD of deviance and rank in data frame or plot.
#'
#' @importFrom dplyr filter select
#' @importFrom tibble is_tibble
#' @importFrom rlang .data
#' @importFrom ggplot2 ggplot geom_point scale_color_manual scale_y_log10
#'            geom_abline theme_bw aes labs scale_x_log10 scale_y_reverse
#' @importFrom ggrepel geom_text_repel
#' @importFrom RColorBrewer brewer.pal
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
#' data(spe)
#' spe_svgs <- spe
#'
#' # biasDetect()
#'

biasDetect <- function(batch_df, nSD_dev = 5, nSD_rank = 5, visual = FALSE) {

    tol <- .Machine$double.eps^0.5
    stopifnot(
        is.data.frame(batch_df) | is_tibble(batch_df), is.logical(visual),
        abs(nSD_dev - round(nSD_dev)) < tol,
        abs(nSD_rank - round(nSD_rank)) < tol
    )

    batch_df$dev_outlier <- batch_df$nSD_dev >= nSD_dev
    batch_df$nSD_bin_dev <- cut(abs(batch_df$nSD_dev), right=FALSE,
                      breaks=seq(0, max(batch_df$nSD_dev)+ nSD_dev,
                                 by= nSD_dev), include.lowest=TRUE)

    batch_df$rank_outlier <- batch_df$nSD_rank >= nSD_rank
    batch_df$nSD_bin_rank <- cut(abs(batch_df$nSD_rank), right=FALSE,
                      breaks=seq(0, max(batch_df$nSD_rank)+ nSD_rank,
                                 by= nSD_rank), include.lowest=TRUE)

    if (visual == TRUE) {
        options(ggrepel.max.overlaps = Inf)
        
        sd.interval <- nSD_dev
        col.pal <- brewer.pal(length(unique(
            batch_df[["nSD_bin_dev"]])),"YlOrRd")
        col.pal[1] <- "grey"
        plot_dev <- batch_df |> 
            ggplot(aes(x=.data[["dev_default"]], y=.data[["dev_batch"]], 
                color=.data[["nSD_bin_dev"]])) +
            geom_point() + scale_x_log10() + scale_y_log10() +
            scale_color_manual(values=col.pal) +
            geom_text_repel(
                aes(label = ifelse(.data[["nSD_dev"]] > sd.interval,
                .data[["gene_name"]], "")), size = 3) +
            labs(x= "dev (no batch)", y="dev (batch)",
                color = "nSD_bin") + 
            theme_bw() +
            geom_abline(aes(slope = 1, intercept = 0), lty = 2)
        
        sd.interval <- nSD_rank
        col.pal2 <- brewer.pal(length(unique(
            batch_df[["nSD_bin_rank"]])),"YlOrRd")
        col.pal2[1] <- "grey"
        plot_rank <- ggplot(batch_df, aes(x=.data[["rank_default"]],
                y=.data[["rank_batch"]], 
                color=.data[["nSD_bin_rank"]])) +
            geom_point() +
            scale_y_reverse() + scale_color_manual(values=col.pal2) +
            geom_abline(aes(slope = -1, intercept = 0), lty = 2) +
            labs(x= "rank (no batch)", y= "rank (batch)",
                color = "nSD_bin") + 
            theme_bw() +
            geom_text_repel(aes(label = ifelse(
                .data[["nSD_rank"]] > sd.interval,
                .data[["gene_name"]], "")), size = 3)

        output <- list(deviance = plot_dev, rank = plot_rank)
    }
    else {
        biased.genes.df <- filter(batch_df,
                    .data$dev_outlier==TRUE | .data$rank_outlier==TRUE)
        
        bias <- biased.genes.df$gene
        names(bias) <- biased.genes.df$gene_name
        output <- bias
    }

    output
    }
