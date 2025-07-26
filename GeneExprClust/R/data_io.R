#' Read gene expression data from file
#'
#' Supports CSV, TSV, or ExpressionSet from .RData file.
#'
#' @param path Path to data file (.csv, .tsv, or .RData)
#' @param type File type: "csv", "tsv", or "eset"
#' @param pheno_data Optional data frame or file path to sample metadata
#'
#' @return ExpressionSet object
#' @importFrom Biobase ExpressionSet AnnotatedDataFrame
#' @importFrom utils read.table
#' @examples 
#' # Load existing ExpressionSet from RData file
#' file_path <- system.file("data", "eset.rda", package = "GeneExprClust")
#' eset <- read_expression(file_path, type = "eset")
#' 
#' # Check the resulting ExpressionSet
#' print(eset)
#' dim(Biobase::exprs(eset))  # Get expression matrix dimensions
#' Biobase::pData(eset)       # Get phenotype data
#' @export
read_expression <- function(path, type = c("csv", "tsv", "eset"), pheno_data = NULL) {
  type <- match.arg(type)
  
  if (!file.exists(path)) stop("File not found: ", path)
  
  if (type %in% c("csv", "tsv")) {
    sep <- ifelse(type == "csv", ",", "\t")
    expr <- as.matrix(read.table(path, sep = sep, header = TRUE, row.names = 1, check.names = FALSE))

    pdata <- NULL
    if (!is.null(pheno_data)) {
      if (is.character(pheno_data) && file.exists(pheno_data)) {
        pdata_df <- read.table(pheno_data, sep = sep, header = TRUE, row.names = 1)
      } else if (is.data.frame(pheno_data)) {
        pdata_df <- pheno_data
      } else {
        stop("Invalid pheno_data input. Must be path or data frame.")
      }
      pdata <- AnnotatedDataFrame(pdata_df)
    } else {
      pdata_df <- data.frame(row.names = colnames(expr), condition = factor(rep("unknown", ncol(expr))))
      pdata <- AnnotatedDataFrame(pdata_df)
    }
    
    eset <- ExpressionSet(assayData = expr, phenoData = pdata)
    if (!inherits(eset, "ExpressionSet")) stop("Failed to create valid ExpressionSet.")
    return(eset)
  }

  if (type == "eset") {
    loaded <- load(path)
    eset_objs <- mget(loaded)
    eset <- Filter(function(x) inherits(x, "ExpressionSet"), eset_objs)
    if (length(eset) == 0) stop("No ExpressionSet found in loaded file.")
    return(eset[[1]])
  }
}
