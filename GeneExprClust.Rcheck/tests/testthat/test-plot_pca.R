library(testthat)
library(Biobase)
library(ggplot2)

test_that("plot_pca produces expected output", {
  expr <- matrix(rnorm(100), nrow=10, 
                 dimnames=list(paste0("G",1:10), paste0("S",1:10)))
  pdata <- data.frame(Condition=rep(c("A","B"), each=5),
                      row.names=paste0("S",1:10))
  eset <- ExpressionSet(assayData=expr, phenoData=AnnotatedDataFrame(pdata))
  
  p <- plot_pca(eset, "Condition")
  expect_s3_class(p, "ggplot")
  expect_match(p$labels$x, "^PC1 \\([0-9.]+%\\)$")
  expect_match(p$labels$y, "^PC2 \\([0-9.]+%\\)$")
  expect_identical(p$labels$colour, "Group")
})

test_that("plot_pca handles ellipse parameter", {
  expr <- matrix(rnorm(100), nrow=10, 
                 dimnames=list(paste0("G",1:10), paste0("S",1:10)))
  pdata <- data.frame(Condition=rep(c("A","B"), each=5),
                      row.names=paste0("S",1:10))
  eset <- ExpressionSet(assayData=expr, phenoData=AnnotatedDataFrame(pdata))
  
  p_ellipse <- plot_pca(eset, "Condition", ellipse=TRUE)
  expect_true(any(sapply(p_ellipse$layers, function(x) inherits(x$geom, "GeomPolygon"))))
  
  p_no_ellipse <- plot_pca(eset, "Condition", ellipse=FALSE)
  expect_false(any(sapply(p_no_ellipse$layers, function(x) inherits(x$geom, "GeomPolygon"))))
})

test_that("plot_pca validates inputs", {
  # Add colnames to expression matrix to match pdata row names
  expr <- matrix(rnorm(100), nrow=10, ncol=10,
                 dimnames=list(paste0("G",1:10), paste0("S",1:10)))

  # This should error because input is not ExpressionSet
  expect_error(plot_pca(expr, "Condition"), "must be an ExpressionSet")
  
  # pdata row names must match colnames(expr)
  pdata <- data.frame(Wrong=rep(1,10), row.names=paste0("S",1:10))
  
  # Now this ExpressionSet is valid (sampleNames match)
  eset2 <- ExpressionSet(assayData=expr, phenoData=AnnotatedDataFrame(pdata))
  
  # Test should expect the exact error message from the function
  expect_error(plot_pca(eset2, "Condition"), 
               sprintf("Column '%s' not found in phenotype data.", "Condition"))
})

test_that("plot_pca handles edge cases", {
  # Case with only one sample - no phenotype data provided
  eset1 <- ExpressionSet(assayData=matrix(1:10, ncol=1, 
                                          dimnames=list(paste0("G",1:10), "S1")))
  expect_error(plot_pca(eset1, "Condition"))  # Will error because "Condition" column missing
  
  # Proper matching names between expr and pdata
  expr2 <- matrix(rnorm(20), nrow=10, dimnames=list(paste0("G",1:10), c("S1","S2")))
  pdata2 <- data.frame(Condition=c("A","B"), row.names=c("S1","S2"))  # rownames match colnames(expr2)
  eset2 <- ExpressionSet(assayData=expr2, phenoData=AnnotatedDataFrame(pdata2))  # Valid ExpressionSet
  
  p <- plot_pca(eset2, "Condition")
  expect_s3_class(p, "ggplot")
})
