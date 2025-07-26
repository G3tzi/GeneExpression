library(testthat)
library(Biobase)
library(ggplot2)

test_that("plot_elbow returns a ggplot object", {
  mat <- matrix(rnorm(100), nrow = 10)
  eset <- ExpressionSet(assayData = mat)
  
  p <- plot_elbow(eset, max_k = 5, target = "gene")
  expect_s3_class(p, "ggplot")
})

test_that("plot_elbow produces decreasing WSS values", {
  mat <- matrix(rnorm(100), nrow = 10)
  eset <- ExpressionSet(assayData = mat)
  
  p <- plot_elbow(eset, max_k = 5, target = "sample")
  
  df <- ggplot2::ggplot_build(p)$data[[1]]
  expect_equal(df$x, 1:5)
  expect_true(all(diff(df$y) <= 0 | diff(df$y) < 0))
})
