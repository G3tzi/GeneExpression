#' Detect and optionally remove outlier genes based on variability
#'
#' @param eset An ExpressionSet object
#' @param method Method for variability detection: "mad" (default) or "iqr"
#' @param threshold Threshold multiplier (e.g., 3 MADs above median)
#' @param plot Logical. If TRUE, shows a distribution plot of variability
#' @param remove Logical. If TRUE (default), removes outlier genes from the returned ExpressionSet.
#' @return A list with:
#'   \item{filtered}{ExpressionSet (original or filtered depending on `remove`)}
#'   \item{outliers}{Character vector of outlier gene names}
#' @importFrom matrixStats rowMads rowIQRs
#' @importFrom Biobase exprs featureNames
#' @importFrom graphics hist abline
#' @importFrom stats median mad quantile
#' @examples
#' data(eset)
#' # Detect and remove outlier genes using MAD method
#' eset_filtered <- detect_outliers(eset, plot = TRUE)
#' @export
detect_outliers <- function(eset, method = c("mad", "iqr"), threshold = 3, plot = FALSE, remove = TRUE) {
  if (!inherits(eset, "ExpressionSet")) stop("Input must be an ExpressionSet.")

  method <- match.arg(method)
  expr <- Biobase::exprs(eset)
  
  # Compute variability metric
  var_metric <- switch(method,
    mad = matrixStats::rowMads(expr),
    iqr = matrixStats::rowIQRs(expr)
  )
  
  # Determine outlier cutoff based on method
  if (method == "mad") {
    metric_median <- median(var_metric)
    mad_var <- mad(var_metric, constant = 1)
    cutoff <- metric_median + threshold * mad_var
  } else if (method == "iqr") {
    Q1 <- quantile(var_metric, 0.25)
    Q3 <- quantile(var_metric, 0.75)
    IQR_val <- Q3 - Q1
    cutoff <- Q3 + threshold * IQR_val
  }
  
  outliers <- which(var_metric > cutoff)
  outlier_genes <- Biobase::featureNames(eset)[outliers]

  if (plot) {
    hist(var_metric, breaks = 50, main = paste("Gene Variability (", method, ")"), 
         xlab = paste("Variability (", method, ")", sep = ""), col = "skyblue")
    abline(v = cutoff, col = "red", lwd = 2)
  }
  
  if (remove && length(outliers) > 0) {
    message(length(outliers), " outlier genes removed.")
    eset_filtered <- eset[-outliers, ]
    attr(eset_filtered, "outliers") <- outlier_genes
  } else {
    message(length(outliers), " outliers detected, but not removed.")
    eset_filtered <- eset
    attr(eset_filtered, "outliers") <- character(0)
  }
  
  return(eset_filtered) 
}

