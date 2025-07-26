library(testthat)
library(stats)
library(mockery)

test_that("plot_dendrogram produces correct plot output", {
  hc <- hclust(dist(matrix(rnorm(100), ncol=10)))
  
  expect_silent(plot_dendrogram(hc))
  expect_silent(plot_dendrogram(hc, "Custom Title"))
  
  mock_plot <- mockery::mock(function(x, ...) {
    args <- list(...)
    expect_equal(args$main, "Test Title")
    expect_equal(args$xlab, "")
    expect_equal(args$sub, "")
    expect_equal(args$cex, 0.6)
  })
  
  mockery::stub(plot_dendrogram, "plot", mock_plot)
  plot_dendrogram(hc, "Test Title")
})
    
test_that("plot_dendrogram maintains hclust structure", {
  hc <- hclust(dist(matrix(rnorm(100), ncol=10)))
  
  # Store original structure
  original_order <- hc$order
  original_labels <- hc$labels
  original_method <- hc$method
  
  # Call plotting function
  plot_dendrogram(hc)
  
  # Verify structure unchanged
  expect_identical(hc$order, original_order)
  expect_identical(hc$labels, original_labels)
  expect_identical(hc$method, original_method)
})