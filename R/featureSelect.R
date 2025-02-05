#' @rdname featureSelect
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
#' @importFrom SummarizedExperiment colData rowData
#' @importFrom stats sd
#'
#' @param input  \code{SpatialExperiment} or \code{SingleCellExperiment}: Input
#'     data which can be either \code{SpatialExperiment} or
#'     \code{SingleCellExperiment}. It is assumed to have an \code{assays} slot
#'     containing \code{counts} assay for \code{devianceFeatureSelection()} to
#'     successfully operate and calculate the deviance and rank. The
#'     \code{logcounts} is strongly recommended to be included in \code{assays}.
#'     Also, the input data is assumed to set up \code{rownames(input)} as
#'     \code{genes} and \code{gene_name} in \code{rowData(input)}. The input can
#'     be either the raw and complete data object or the filtered data object
#'     containing only highly variable genes (HVGs) or spatially varaible genes
#'     (SVGs).
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
#'     gene name/ID, e.g. "ENSG00000002330".
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
#' library(dplyr)
#' library(SummarizedExperiment)
#' # retrieve raw spe input
#' ehub <- ExperimentHub::ExperimentHub()
#' spe <- ehub[["EH9605"]]
#' # subset 4 samples out
#' fix_order <- distinct(
#'     as.data.frame(colData(spe)), slide, array, 
#'     brnum, sample_id, position, sex) %>% 
#'     arrange(slide, array)
#' sub4 <- fix_order$sample_id[c(14, 16, 20, 21)]
#' spe_sub4 <- spe[,spe$sample_id %in% sub4]
#' # retrieve SVGs
#' res_ranks <- read.csv(
#'     system.file("extdata","res_ranks.csv",package = "BiasDetect"),
#'     row.names = 1, check.names = FALSE)
#' res_df_sub4 <- tidyr::pivot_longer(
#'     tibble::rownames_to_column(as.data.frame(res_ranks), var<-"gene_id"), 
#'     colnames(res_ranks), 
#'     names_to = "sample_id", 
#'     values_to = "rank", 
#'     values_drop_na = TRUE)
#' # subset to 4 samples (rank <= 2000, n > 1)
#' res_df2_sub4 <- filter(res_df_sub4, 
#'     sample_id %in%
#'     c("V11L05-333_B1","V11L05-333_D1","V11L05-335_D1","V11L05-336_A1"),
#'     rank <= 2000)
#' svgs_sub4 <- group_by(res_df2_sub4, gene_id) %>% 
#'     tally() %>% 
#'     filter(n > 1)
#' SVGs <- svgs_sub4$gene_id
#' # feature selection
#' batch_df <- featureSelect(spe_sub4, batch_effect = "sample_id", VGs = SVGs)
#' 
featureSelect <- function(input, batch_effect = NULL, VGs = NULL) {

    stopifnot(
        inherits(input, c("SpatialExperiment", "SingleCellExperiment")),
        is.character(batch_effect),
        is.character(VGs)
    )

    if (!is.null(batch_effect)) {
        if (!batch_effect %in% names(colData(input))) {
            stop("The batch_effect is not a valid column")}
        batch_effect <- colData(input)[[batch_effect]]
    } else { stop("Please provide a valid batch_effect.")}
    
    # Re-organize raw input
    input <- input[rowData(input)$gene_id %in% VGs, ]
    rownames(input) <- rowData(input)$gene_id
    
    message("Step 1: Running feature selection without batch...")
    bd <- devianceFeatureSelection(input, assay = "counts", fam = "binomial")
    bd_df <- cbind.data.frame("gene"=rownames(bd),
                            "gene_name"=rowData(bd)$gene_name,
                            "dev"= rowData(bd)$binomial_deviance,
                            "rank"=(nrow(bd)+1)
                                -rank(rowData(bd)$binomial_deviance))
    rownames(bd_df) <- bd_df$gene

    message("Step 2: Running feature selection with batch...")
    bd_batch <- devianceFeatureSelection(input, assay = "counts",
            fam = "binomial",batch = as.factor(batch_effect))
    bd_batch_df <- cbind.data.frame("gene"=rownames(bd_batch),
                              "gene_name"=rowData(bd_batch)$gene_name,
                              "dev"= rowData(bd_batch)$binomial_deviance,
                              "rank"=(nrow(bd_batch)+1)
                                  -rank(rowData(bd_batch)$binomial_deviance))
    rownames(bd_batch_df) <- bd_batch_df$gene

    message("Step 3: Calculating deviance and rank difference...")
    batch_df <- left_join(bd_df, bd_batch_df, by=c("gene", "gene_name"),
                            suffix=c("_default","_batch"))

    batch_df$d_diff <- (batch_df$dev_default-batch_df$dev_batch)/
                        batch_df$dev_batch
    mean_dev <- mean(batch_df$d_diff)
    sd_dev <- sd(batch_df$d_diff)
    batch_df$nSD_dev <- (batch_df$d_diff - mean_dev) / sd_dev

    batch_df$r_diff <- batch_df$rank_batch-batch_df$rank_default
    mean_rank <- mean(batch_df$r_diff)
    sd_rank <- sd(batch_df$r_diff)
    batch_df$nSD_rank <- (batch_df$r_diff - mean_rank) / sd_rank

    batch_df
    }


