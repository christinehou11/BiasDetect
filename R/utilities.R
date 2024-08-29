# helper function

#' @importFrom SummarizedExperiment assayNames assays rowData
.df_maker <- function(input) {
    df <- cbind.data.frame("gene"=rownames(input),
                "gene_name"=rowData(input)$gene_name,
                "dev"= rowData(input)$binomial_deviance,
                "rank"=(nrow(input)+1)-rank(rowData(input)$binomial_deviance))

    rownames(df) = df$gene

    df
    }
