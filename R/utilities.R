# helper function

.is_wholenumber <- function(x, tol = .Machine$double.eps^0.5) {
    abs(x - round(x)) < tol
    }

#' @importFrom SummarizedExperiment assayNames assays rowData
.df_maker <- function(input) {
    df <- cbind.data.frame("gene"=rownames(input),
                "gene_name"=rowData(input)$gene_name,
                "dev"= rowData(input)$binomial_deviance,
                "rank"=(nrow(input)+1)-rank(rowData(input)$binomial_deviance))

    rownames(df) <- df$gene

    df
    }

.check <- function(input, batch_effect, VGs, nSD_dev, nSD_rank) {
    inherits(input, c("SpatialExperiment", "SingleCellExperiment")) &
        is.character(batch_effect) & is.character(VGs) &
        .is_wholenumber(nSD_dev) & .is_wholenumber(nSD_rank)
    }

.sd_cutoff <- function(df, nSD_cutoff, type = c("dev", "rank")) {

    type <- match.arg(type)
    col <- ifelse(type == "dev", "d.diff", "r.diff")
    col2 <- ifelse(type == "dev", "nSD_dev", "nSD_rank")
    col3 <- ifelse(type == "dev", "dev_outlier", "rank_outlier")

    mean <- mean(df[[col]])
    sd <- sd(df[[col]])
    df[[col2]] <- (df[[col]] - mean) / sd
    df[[col3]] <- df[[col2]] >= nSD_cutoff

    df
    }

#' @importFrom ggplot2 ggplot geom_point scale_color_manual
#'            geom_abline theme_bw aes labs
#' @importFrom RColorBrewer brewer.pal
#' @importFrom rlang .data
.plot <- function(data, x, y, nSD_cutoff, nSD_var, type = c("dev", "rank")) {

    type <- match.arg(type)
    col1 <- ifelse(type == "dev", "nSD.bin_dev", "nSD.bin_rank")

    data[[col1]] <- cut(abs(data[[nSD_var]]), right=FALSE,
                                    breaks=seq(0,
                                        max(data[[nSD_var]])+ nSD_cutoff,
                                        by= nSD_cutoff),
                                    include.lowest=TRUE)
    col.pal <- brewer.pal(length(unique(data[[col1]])), "YlOrRd")
    col.pal[1] <- "grey"

    plot <- ggplot(data, aes(x=.data[[x]], y=.data[[y]], color=.data[[col1]])) +
        geom_point()+
        scale_color_manual(values=col.pal)+
        labs(x= paste(type, " (no batch)"), y= paste(type, " (batch)"))+
        theme_bw()

    plot
    }
