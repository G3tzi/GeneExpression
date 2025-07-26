library(testthat)
library(Biobase)

test_that("plot_heatmap produces expected output", {
  expr <- matrix(rnorm(100), nrow=10, dimnames=list(paste0("G",1:10), paste0("S",1:10)))
  eset <- ExpressionSet(assayData=expr)
  
  hm <- plot_heatmap(eset)
  expect_s4_class(hm, "Heatmap")
  expect_equal(hm@matrix, expr)

  clusters <- setNames(rep(1:2, each=5), rownames(expr))
  hm_clust <- plot_heatmap(eset, clusters=clusters)
  expect_s4_class(hm_clust@right_annotation, "HeatmapAnnotation")
  
  hm_rowscale <- plot_heatmap(eset, scale="row")
  expect_equal(hm_rowscale@matrix, t(scale(t(expr))))
})

test_that("plot_heatmap handles edge cases", {

  eset1 <- ExpressionSet(assayData=matrix(1:5, nrow=1, dimnames=list("G1", paste0("S",1:5))))
  expect_s4_class(plot_heatmap(eset1), "Heatmap")
  
  eset2 <- ExpressionSet(assayData=matrix(1:5, ncol=1, dimnames=list(paste0("G",1:5), "S1")))
  expect_s4_class(plot_heatmap(eset2), "Heatmap")
})

test_that("plot_heatmap validates inputs", {
  expr <- matrix(rnorm(100), nrow=10)
  eset <- ExpressionSet(assayData=expr)
  
  expect_error(plot_heatmap(expr), "must be an ExpressionSet")
  
  expect_error(plot_heatmap(eset, scale="invalid"), "should be one of")
})