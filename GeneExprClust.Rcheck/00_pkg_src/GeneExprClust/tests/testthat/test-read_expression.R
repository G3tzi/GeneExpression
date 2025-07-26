library(testthat)
library(Biobase)

test_that("read_expression reads CSV correctly", {
  expr_data <- matrix(1:9, nrow = 3, dimnames = list(c("gene1", "gene2", "gene3"), c("sample1", "sample2", "sample3")))
  tmpfile <- tempfile(fileext = ".csv")
  write.csv(expr_data, tmpfile, row.names = TRUE)
  
  eset <- read_expression(tmpfile, type = "csv")
  expect_s4_class(eset, "ExpressionSet")
  expect_equal(exprs(eset), expr_data)
})

test_that("read_expression reads TSV correctly", {
  expr_data <- matrix(1:9, nrow = 3, dimnames = list(c("gene1", "gene2", "gene3"), c("sample1", "sample2", "sample3")))
  tmpfile <- tempfile(fileext = ".tsv")
  write.table(expr_data, tmpfile, sep = "\t", row.names = TRUE, col.names = NA)
  
  eset <- read_expression(tmpfile, type = "tsv")
  expect_s4_class(eset, "ExpressionSet")
  expect_equal(exprs(eset), expr_data)
})

test_that("read_expression reads ExpressionSet from RData", {
  expr_data <- matrix(rnorm(9), nrow = 3,
                      dimnames = list(paste0("gene", 1:3), paste0("sample", 1:3)))
  pheno <- data.frame(row.names = paste0("sample", 1:3), group = c("A", "B", "A"))
  eset <- ExpressionSet(assayData = expr_data, phenoData = AnnotatedDataFrame(pheno))
  
  tmpfile <- tempfile(fileext = ".RData")
  save(eset, file = tmpfile)
  
  eset_loaded <- read_expression(tmpfile, type = "eset")
  expect_s4_class(eset_loaded, "ExpressionSet")
  expect_equal(exprs(eset_loaded), expr_data)  # Now rownames will match
})

test_that("read_expression handles pheno_data input", {
  expr_data <- matrix(1:4, nrow = 2, dimnames = list(c("gene1", "gene2"), c("sample1", "sample2")))
  expr_file <- tempfile(fileext = ".csv")
  write.csv(expr_data, expr_file, row.names = TRUE)
  
  pheno_df <- data.frame(row.names = c("sample1", "sample2"), condition = c("A", "B"))
  eset <- read_expression(expr_file, type = "csv", pheno_data = pheno_df)
  expect_equal(pData(eset)$condition, c("A", "B"))
  
  pheno_file <- tempfile(fileext = ".csv")
  write.csv(pheno_df, pheno_file, row.names = TRUE)
  eset2 <- read_expression(expr_file, type = "csv", pheno_data = pheno_file)
  expect_equal(pData(eset2)$condition, c("A", "B"))
})

test_that("read_expression handles missing files", {
  expect_error(read_expression("nonexistent.csv", type = "csv"), "File not found")
})

test_that("read_expression handles invalid pheno_data", {
  expr_data <- matrix(1:4, nrow = 2, dimnames = list(c("gene1", "gene2"), c("sample1", "sample2")))
  expr_file <- tempfile(fileext = ".csv")
  write.csv(expr_data, expr_file, row.names = TRUE)
  
  expect_error(read_expression(expr_file, type = "csv", pheno_data = 42), "Invalid pheno_data input")
})
