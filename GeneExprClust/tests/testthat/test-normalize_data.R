create_test_eset <- function() {
  # Create matrix where some genes meet thresholds and some don't
  expr <- matrix(
    c(
      # 3 genes that should be kept (meet threshold in >=5 samples)
      rep(c(10, 10, 10, 10, 10, 0, 0), 3),
      # 3 genes that should be filtered (don't meet threshold)
      rep(c(1, 1, 0, 0, 0, 0, 0), 3)
    ),
    nrow = 6,
    byrow = TRUE,
    dimnames = list(paste0("gene", 1:6), paste0("sample", 1:7))
  )

  Biobase::ExpressionSet(assayData = expr)
}

testthat::test_that("normalize_data validates input correctly", {
  testthat::expect_error(normalize_data(matrix(1:10, nrow=2)), "Input must be an ExpressionSet")
  
  eset <- create_test_eset()
  testthat::expect_error(normalize_data(eset, method="invalid"), "should be one of")
  
  neg_eset <- eset
  Biobase::exprs(neg_eset)[1,1] <- -1
  testthat::expect_error(normalize_data(neg_eset, method="log2"), "Negative values found")
})

testthat::test_that("log2 normalization works correctly", {
  eset <- create_test_eset()
  norm_eset <- normalize_data(eset, method="log2")
  
  testthat::expect_s4_class(norm_eset, "ExpressionSet")
  
  original <- Biobase::exprs(eset)
  transformed <- Biobase::exprs(norm_eset)
  testthat::expect_equal(transformed, log2(original + 1))
  
  testthat::expect_equal(rownames(transformed), rownames(original))
  testthat::expect_equal(colnames(transformed), colnames(original))
})

testthat::test_that("quantile normalization works correctly", {
  eset <- create_test_eset()
  original_expr <- Biobase::exprs(eset)
  
  # Create mock
  mock_norm <- function(x) {
    out <- apply(x, 2, rank)
    dimnames(out) <- dimnames(x)
    return(out)
  }
  
  # Mock using the correct namespace
  testthat::with_mocked_bindings(
    normalize.quantiles = mock_norm,
    .package = "preprocessCore",
    {
      norm_eset <- normalize_data(eset, method = "quantile")
    }
  )
  
  testthat::expect_s4_class(norm_eset, "ExpressionSet")
  testthat::expect_equal(dim(norm_eset), dim(eset))
  testthat::expect_equal(exprs(norm_eset), mock_norm(original_expr))
})

testthat::test_that("zscore normalization works correctly", {
  eset <- create_test_eset()
  norm_eset <- normalize_data(eset, method="zscore")
  
  # Check row means are ~0 and sd ~1
  row_stats <- t(apply(Biobase::exprs(norm_eset), 1, function(x) c(mean(x), sd(x))))
  testthat::expect_true(all(abs(row_stats[,1]) < 1e-10))
  testthat::expect_true(all(abs(row_stats[,2] - 1) < 1e-10))
})

testthat::test_that("deseq normalization works correctly", {
  # Create test data with row/col names
  gene_names <- paste0("gene", 1:10)
  sample_names <- paste0("sample", 1:10)
  
  expr <- matrix(
    rpois(100, lambda = 10),
    nrow = 10,
    dimnames = list(gene_names, sample_names)
  )
  eset <- Biobase::ExpressionSet(assayData = expr)
  
  # Create mock counts WITH DIMNAMES
  mock_counts <- matrix(
    rpois(100, lambda = 10),
    nrow = 10,
    dimnames = list(gene_names, sample_names)  # Same as input
  )
  
  # Set up mocks (unchanged)
  mockery::stub(normalize_data, "perform_normalization", function(expr, method) mock_counts)
  
  # Test and compare
  norm_eset <- normalize_data(eset, method = "deseq")
  testthat::expect_equal(unname(Biobase::exprs(norm_eset)), unname(mock_counts))  # Compare values only
  testthat::expect_equal(dimnames(Biobase::exprs(norm_eset)), dimnames(mock_counts))  # Compare names separately
})

testthat::test_that("batch correction works when requested", {
  eset <- create_test_eset()
  original_expr <- exprs(eset)
  
  # Explicitly match batch vector to 7 samples
  batch <- c(1, 1, 1, 1, 1, 2, 2)
  
  # Define mock ComBat
  mock_combat <- function(x, batch) {
    x + 1  # Simple mock transformation
  }

  # Use mock
  testthat::with_mocked_bindings(
    ComBat = mock_combat,
    .package = "sva",
    {
      norm_eset <- normalize_data(eset, method = "log2", batch = batch)
    }
  )

  expected <- log2(original_expr + 1) + 1
  testthat::expect_equal(exprs(norm_eset), expected)
})

testthat::test_that("filtering works correctly", {
  expr <- matrix(
    c(
      # Genes 1–5: pass the threshold in 3 samples (will be kept)
      rep(c(10, 10, 10), 5),
      # Genes 6–10: only 2 values > threshold (will be filtered out)
      rep(c(10, 10, 1), 5)
    ),
    nrow = 10,
    byrow = TRUE,
    dimnames = list(paste0("gene", 1:10), paste0("sample", 1:3))
  )
  eset <- Biobase::ExpressionSet(assayData = expr)
  
  norm_eset <- normalize_data(
    eset,
    method = "log2",  # Or your expected normalization method
    filter = TRUE,
    filter_threshold = 2,
    filter_min_samples = 3
  )
  
  testthat::expect_lt(nrow(norm_eset), nrow(eset))
  
  expr_after <- Biobase::exprs(norm_eset)
  expr_before <- Biobase::exprs(normalize_data(eset, method = "log2"))  # Same normalization, no filtering
  
  kept_genes <- rownames(expr_after)
  removed_genes <- setdiff(rownames(expr_before), kept_genes)
  
  if (length(kept_genes) > 0) {
    kept_expr <- expr_before[kept_genes, ]
    testthat::expect_true(all(rowSums(kept_expr > 2) >= 3),
                info = "All kept genes should meet threshold")
  }
  
  if (length(removed_genes) > 0) {
    removed_expr <- expr_before[removed_genes, ]
    testthat::expect_true(all(rowSums(removed_expr > 2) < 3),
                info = "All removed genes should fail threshold")
  }
})


testthat::test_that("normalization preserves sample and feature metadata", {
  eset <- create_test_eset()
  Biobase::pData(eset)$group <- rep(c("A","B"), length.out = 7)
  Biobase::fData(eset)$symbol <- paste0("GENE_", 1:6)
  
  norm_eset <- normalize_data(eset, method="log2")
  
  testthat::expect_equal(Biobase::pData(norm_eset), Biobase::pData(eset))
  testthat::expect_equal(Biobase::fData(norm_eset), Biobase::fData(eset))
  testthat::expect_equal(sampleNames(norm_eset), sampleNames(eset))
  testthat::expect_equal(featureNames(norm_eset), featureNames(eset))
})