# helper function

.is_wholenumber <- function(x, tol = .Machine$double.eps^0.5) {
    abs(x - round(x)) < tol
    }

#' @importFrom SummarizedExperiment rowData
.df_maker <- function(input) {
    df <- cbind.data.frame("gene"=rownames(input),
                "gene_name"=rowData(input)$gene_name,
                "dev"= rowData(input)$binomial_deviance,
                "rank"=(nrow(input)+1)-rank(rowData(input)$binomial_deviance))

    rownames(df) <- df$gene

    df
    }

.check_feature <- function(input, batch_effect, VGs) {
    inherits(input, c("SpatialExperiment", "SingleCellExperiment")) &
    is.character(batch_effect) & is.character(VGs)
    }

#' @importFrom tibble is_tibble
.check_bias <- function(batch_df, nSD_dev = 5, nSD_rank = 5, visual = FALSE) {
    is.data.frame(batch_df) | is_tibble(batch_df) &
    .is_wholenumber(nSD_dev) & .is_wholenumber(nSD_rank) &
    is.logical(visual)
    }

#' @importFrom stats sd
.sd_cutoff <- function(df, nSD_cutoff, nSD_var, type = c("dev", "rank")) {

    type <- match.arg(type)
    col <- ifelse(type == "dev", "d_diff", "r_diff")
    col2 <- ifelse(type == "dev", "nSD_dev", "nSD_rank")
    col3 <- ifelse(type == "dev", "dev_outlier", "rank_outlier")
    col4 <- ifelse(type == "dev", "nSD_bin_dev", "nSD_bin_rank")

    mean <- mean(df[[col]])
    sd <- sd(df[[col]])
    df[[col2]] <- (df[[col]] - mean) / sd
    df[[col3]] <- df[[col2]] >= nSD_cutoff

    df[[col4]] <- cut(abs(df[[nSD_var]]), right=FALSE,
                        breaks=seq(0, max(df[[nSD_var]])+ nSD_cutoff,
                                   by= nSD_cutoff),
                        include.lowest=TRUE)

    df
    }

#' @importFrom ggplot2 ggplot geom_point scale_color_manual
#'            geom_abline theme_bw aes labs
#' @importFrom ggrepel geom_text_repel
#' @importFrom RColorBrewer brewer.pal
#' @importFrom rlang .data
.plot <- function(data, x, y, nSD_cutoff, type = c("dev", "rank")) {

    type <- match.arg(type)
    col <- ifelse(type == "dev", "nSD_bin_dev", "nSD_bin_rank")
    col2 <- ifelse(type == "dev", "nSD_dev", "nSD_rank")

    col.pal <- brewer.pal(length(unique(data[[col]])), "YlOrRd")
    col.pal[1] <- "grey"

    plot <- ggplot(data, aes(x=.data[[x]], y=.data[[y]], color=.data[[col]])) +
        geom_point()+
        scale_color_manual(values=col.pal)+
        labs(x= paste(type, " (no batch)"), y= paste(type, " (batch)"))+
        geom_text_repel(aes(label = ifelse(.data[[col2]] >= nSD_cutoff,
                    .data[["gene_name"]], "")), size = 3, max.overlaps = 10) +
        theme_bw()

    plot
    }
