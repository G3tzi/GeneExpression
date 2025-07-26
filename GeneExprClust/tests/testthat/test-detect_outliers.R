# Setup (shared test data)
setup_test_eset <- function() {
  test_exprs <- matrix(c(
    rnorm(90, mean = 10, sd = 0.1),  # Normal genes
    rnorm(10, mean = 10, sd = 2)     # Outlier genes
  ), nrow = 10, byrow = TRUE)
  rownames(test_exprs) <- paste0("gene_", 1:10)
  colnames(test_exprs) <- paste0("sample_", 1:10)
  Biobase::ExpressionSet(assayData = test_exprs)
}

test_that("detect_outliers returns ExpressionSet with outliers attribute", {
  eset <- setup_test_eset()
  result <- detect_outliers(eset)

  testthat::expect_s4_class(result, "ExpressionSet")
  testthat::expect_true("outliers" %in% names(attributes(result)))
  testthat::expect_type(attr(result, "outliers"), "character")
})

test_that("detect_outliers removes outliers when remove=TRUE", {
  eset <- setup_test_eset()
  result <- detect_outliers(eset, remove = TRUE)
  testthat::expect_true(nrow(result) < nrow(eset))
})

test_that("detect_outliers keeps outliers when remove=FALSE", {
  eset <- setup_test_eset()
  result <- detect_outliers(eset, remove = FALSE)
  testthat::expect_equal(nrow(result), nrow(eset))
})

test_that("detect_outliers works with method='iqr'", {
  eset <- setup_test_eset()
  result <- detect_outliers(eset, method = "iqr")
  outliers <- attr(result, "outliers")
  testthat::expect_type(outliers, "character")
  testthat::expect_gt(length(outliers), 0)
})

test_that("detect_outliers detects more outliers with lower threshold", {
  eset <- setup_test_eset()
  result_thresh1 <- detect_outliers(eset, threshold = 1)
  result_thresh5 <- detect_outliers(eset, threshold = 5)
  testthat::expect_gte(length(attr(result_thresh1, "outliers")), length(attr(result_thresh5, "outliers")))
})

test_that("detect_outliers throws error for non-ExpressionSet input", {
  testthat::expect_error(detect_outliers(matrix(1:10, nrow = 2)), "Input must be an ExpressionSet")
})

test_that("detect_outliers handles no outliers case", {
  no_outlier_exprs <- matrix(rnorm(100, mean = 10, sd = 0.1), nrow = 10)
  rownames(no_outlier_exprs) <- paste0("gene_", 1:10)
  no_outlier_eset <- Biobase::ExpressionSet(assayData = no_outlier_exprs)
  result <- detect_outliers(no_outlier_eset, threshold = 10)
  testthat::expect_length(attr(result, "outliers"), 0)
})
