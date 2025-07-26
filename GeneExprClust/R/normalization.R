#' Perform normalization on an expression matrix
#'
#' Supports "log2", "quantile", "zscore", and "deseq" normalization methods.
#'
#' @param expr Numeric matrix of expression values (genes x samples).
#' @param method Normalization method, one of "log2", "quantile", "zscore", or "deseq".
#' @return Normalized expression matrix.
#' @importFrom preprocessCore normalize.quantiles
#' @importFrom DESeq2 DESeqDataSetFromMatrix estimateSizeFactors counts
perform_normalization <- function(expr, method) {
  method <- match.arg(method, c("log2", "quantile", "zscore", "deseq"))

 if (method == "log2") {
    if (any(is.na(expr))) warning("NA values detected in expression matrix.")
    if (any(expr < 0, na.rm = TRUE)) stop("Negative values found: log2 normalization is not appropriate.")
    expr <- log2(expr + 1)
  } else if (method == "quantile") {
    if (!requireNamespace("preprocessCore", quietly = TRUE)) {
      stop("Package 'preprocessCore' is required for quantile normalization.")
    }
    if (any(is.na(expr))) warning("NA values detected; preprocessCore::normalize.quantiles may fail.")
    expr <- preprocessCore::normalize.quantiles(expr)
  } else if (method == "zscore") {
    if (any(is.na(expr))) warning("NA values detected; scale() will propagate NAs.")
    expr <- t(scale(t(expr)))
  } else if (method == "deseq") {
    if (!requireNamespace("DESeq2", quietly = TRUE)) {
      stop("Package 'DESeq2' is required for DESeq2 normalization.")
    }
    if (any(expr < 0, na.rm = TRUE)) stop("Negative values found: DESeq2 expects count data.")
    dds <- DESeq2::DESeqDataSetFromMatrix(countData = round(expr),
                                          colData = data.frame(row.names = colnames(expr)),
                                          design = ~1)
    dds <- DESeq2::estimateSizeFactors(dds)
    expr <- DESeq2::counts(dds, normalized = TRUE)
  }
  
  if (any(is.nan(expr))) warning("NaN values detected after normalization.")
  return(expr)
}

#' Perform batch correction on expression matrix using ComBat
#'
#' @param expr Numeric matrix of expression values (genes x samples).
#' @param batch Factor or vector indicating batch assignment per sample.
#' @return Batch-corrected expression matrix.
#' @importFrom sva ComBat
perform_batch_correction <- function(expr, batch) {
  if (!requireNamespace("sva", quietly = TRUE)) {
    stop("Package 'sva' is required for batch correction.")
  }
  if (length(batch) != ncol(expr)) stop("Length of 'batch' must equal number of samples (columns) in expression matrix.")
  expr <- sva::ComBat(expr, batch = batch)
  return(expr)
}

#' Filter low-expressed genes from expression matrix and ExpressionSet
#'
#' @param expr Numeric expression matrix (genes x samples).
#' @param eset ExpressionSet object.
#' @param filter_threshold Numeric, minimum expression threshold to consider a gene expressed.
#' @param filter_min_samples Integer, minimum number of samples that must meet the threshold.
#' @return List with filtered expression matrix (`expr`) and ExpressionSet (`eset`).
filter_low_expression <- function(expr, eset, filter_threshold = 1, filter_min_samples = 3) {
  keep <- rowSums(expr > filter_threshold) >= filter_min_samples
  
  if (sum(!keep) > 0) {
    message("Filtering out ", sum(!keep), " low-expressed genes.")
  }
  
  if (sum(keep) == 0) {
    warning("All genes filtered out! Returning empty results. Adjust filter thresholds.")
    return(list(
      expr = expr[integer(0), , drop = FALSE],  # Empty matrix with same structure
      eset = eset[integer(0), ]                 # Empty eset with same structure
    ))
  }
  
  list(
    expr = expr[keep, , drop = FALSE],
    eset = eset[rownames(expr)[keep], ]
  )
}

#' Normalize expression data in an ExpressionSet with optional batch correction and filtering
#'
#' Supports log2 transformation, quantile normalization, z-score standardization,
#' DESeq2 median-of-ratios normalization, optional batch correction (ComBat),
#' and low-expression filtering.
#'
#' @param eset ExpressionSet object to normalize.
#' @param method Normalization method: "log2", "quantile", "zscore", or "deseq".
#' @param batch Optional batch factor for ComBat batch correction (default NULL).
#' @param filter Logical, whether to filter low-expressed genes (default FALSE).
#' @param filter_threshold Numeric, minimum count threshold for filtering (default 10).
#' @param filter_min_samples Integer, minimum number of samples meeting threshold (default 3).
#' @return Normalized ExpressionSet.
#' @importFrom Biobase exprs sampleNames featureNames
#' @examples
#' data(eset)
#' # Log2 normalization
#' eset_log2 <- normalize_data(eset, method = "log2")
#' @export
normalize_data <- function(eset, method = c("log2", "quantile", "zscore", "deseq"),
                           batch = NULL, filter = FALSE,
                           filter_threshold = 10, filter_min_samples = 3) {
  if (!inherits(eset, "ExpressionSet")) stop("Input must be an ExpressionSet.")
  method <- match.arg(method)
  
  expr <- Biobase::exprs(eset)
  
  expr <- perform_normalization(expr, method)
  
  if (!is.null(batch)) {
    if (is.character(batch) && length(batch) == 1) {
      if (!(batch %in% colnames(Biobase::pData(eset)))) {
        stop(paste0("Batch variable '", batch, "' not found in phenotype data."))
      }
      batch_vec <- Biobase::pData(eset)[[batch]]
    } else {
      batch_vec <- batch
    }

    expr <- perform_batch_correction(expr, batch_vec)
  }
  
  if (filter) {
    res <- filter_low_expression(expr, eset, filter_threshold, filter_min_samples)
    expr <- res$expr
    eset <- res$eset
  }
  
  Biobase::exprs(eset) <- expr
  return(eset)
}