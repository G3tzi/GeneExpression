## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.align = "center",
  warning = FALSE,
  message = FALSE
)

## ----libraries----------------------------------------------------------------
library(GeneExprClust)
library(ALL)              # Bioconductor dataset
library(Biobase)          # ExpressionSet manipulation
library(ggplot2)          # Enhanced plotting
library(dplyr)            # Data manipulation
library(knitr)            # Table formatting

# Load the ALL dataset
data(ALL)

# Examine the dataset structure
message("Dataset dimensions: ", dim(ALL))
message("Number of genes: ", nrow(ALL)) 
message("Number of samples: ", ncol(ALL))

# Look at sample information
head(pData(ALL)[, c("age", "sex", "BT", "mol.biol")], 10)

## ----data_prep----------------------------------------------------------------
# Filter samples - focus on B-cell samples with molecular biology information
selected_samples <- which(!is.na(pData(ALL)$mol.biol) & 
                         pData(ALL)$BT %in% c("B1", "B2", "B3"))

message("Selected ", length(selected_samples), " samples for analysis\n")

# Subset the data
ALL_subset <- ALL[, selected_samples]

# Filter genes - keep most variable genes to focus on interesting patterns
gene_vars <- apply(exprs(ALL_subset), 1, var, na.rm = TRUE)
top_variable_genes <- order(gene_vars, decreasing = TRUE)[1:2000]

ALL_filtered <- ALL_subset[top_variable_genes, ]

message("Final dataset dimensions: ", dim(ALL_filtered))

# Examine the molecular biology subtypes in our subset
table(pData(ALL_filtered)$mol.biol)

## ----basic_pipeline, results='hide'-------------------------------------------
# Create output directory
output_dir <- "ALL_analysis_basic"

# Run basic pipeline
basic_results <- run_pipeline(
  expression_set = ALL_filtered,
  output_dir = output_dir,
  
  # Use standard normalization and filtering
  norm_method = "quantile",
  filter_low_expr = FALSE,
  filter_threshold = 100,  # Higher threshold for microarray data
  
  # Try multiple clustering methods
  clustering_methods = c("hierarchical", "kmeans"),
  k_clusters = 6,  # Expecting several molecular subtypes
  target = "gene",
  
  # Enable validation and visualization
  perform_validation = TRUE,
  create_plots = TRUE,
  
  # Use online gene sets for enrichment
  perform_enrichment = TRUE,
  gene_sets_source = "online",
  gene_sets_collection = "H",  # Hallmark pathways
  
  verbose = TRUE
)

## ----basic_results------------------------------------------------------------
# Summary of the analysis
message(basic_results$summary)

# Clustering results summary
for (method in names(basic_results$clustering_results)) {
  clusters <- basic_results$clustering_results[[method]]$clusters
  message("\n", stringr::str_to_title(method), " Clustering Results:\n")
  message("Number of clusters: ", length(unique(clusters)))
  message("Cluster sizes: ", table(clusters))
}

# Validation metrics
if (length(basic_results$validation_results) > 0) {
  message("\nClustering Validation Metrics:\n")
  for (method in names(basic_results$validation_results)) {
    val <- basic_results$validation_results[[method]]
    if (!is.null(val$avg_silhouette)) {
      message(method, " - Average Silhouette Score: ", round(val$avg_silhouette, 3))
    }
  }
}

## ----advanced_pipeline, results='hide'----------------------------------------
# Advanced pipeline with multiple methods and parameters
advanced_results <- run_pipeline(
  expression_set = ALL_filtered,
  output_dir = "ALL_analysis_advanced",
  
  # Advanced normalization
  norm_method = "quantile",
  filter_low_expr = FALSE,
  filter_threshold = 100,
  
  # Outlier detection
  detect_outliers = TRUE,
  outlier_method = "mad",
  outlier_threshold = 3,
  remove_outliers = TRUE,
  
  # Multiple clustering methods
  clustering_methods = c("hierarchical", "kmeans", "dbscan"),
  k_clusters = 6,
  max_k = 12,  # For elbow plot
  dbscan_eps = 0.8,  # DBSCAN parameters
  dbscan_minPts = 5,
  
  # Enhanced validation
  perform_validation = TRUE,
  gap_statistic = TRUE,
  gap_max_k = 10,
  gap_B = 50,  # More bootstrap samples
  
  # Rich visualization
  create_plots = TRUE,
  pca_color_by = "mol.biol",  # Color PCA by molecular subtype
  heatmap_scale = "row",
  
  # Comprehensive enrichment analysis
  perform_enrichment = TRUE,
  gene_sets_source = "online",
  gene_sets_collection = "C2",  # Curated gene sets (larger collection)
  enrichment_p_threshold = 0.05,
  min_gene_set_size = 10,
  max_gene_set_size = 200,
  
  verbose = TRUE
)

## ----clustering_comparison----------------------------------------------------
# Compare clustering methods
methods <- names(advanced_results$clustering_results)
comparison_df <- data.frame(
  Method = character(),
  N_Clusters = integer(),
  Avg_Silhouette = numeric(),
  WSS = numeric(),
  stringsAsFactors = FALSE
)

for (method in methods) {
  clusters <- advanced_results$clustering_results[[method]]$clusters
  
  # Handle DBSCAN noise points
  if (method == "dbscan") {
    n_clusters <- length(unique(clusters[clusters != 0]))
    n_noise <- sum(clusters == 0)
  } else {
    n_clusters <- length(unique(clusters))
    n_noise <- 0
  }
  
  # Get validation metrics
  val <- advanced_results$validation_results[[method]]
  avg_sil <- if(!is.null(val$avg_silhouette)) round(val$avg_silhouette, 3) else NA
  wss <- if(!is.null(val$wss)) round(val$wss, 0) else NA
  
  comparison_df <- rbind(comparison_df, data.frame(
    Method = stringr::str_to_title(method),
    N_Clusters = n_clusters,
    N_Noise = if(method == "dbscan") n_noise else NA,
    Avg_Silhouette = avg_sil,
    WSS = wss
  ))
}

kable(comparison_df, caption = "Clustering Methods Comparison")

## ----biological_interpretation------------------------------------------------
# Get the molecular biology data for the samples that were actually clustered
mol_biol <- pData(ALL_filtered)$mol.biol

for (method in methods) {
  clusters <- advanced_results$clustering_results[[method]]$clusters
  
  message("\n", stringr::str_to_title(method), " Clustering vs Molecular Biology Subtypes:\n")
  message("Length check - mol_biol: ", length(mol_biol), " clusters: ", length(clusters))
  
  # Ensure both vectors have the same length
  if (length(mol_biol) != length(clusters)) {
    message("Warning: Length mismatch detected. Attempting to align data...\n")
    
    # Option 1: If clusters is shorter (due to outlier removal), subset mol_biol
    if (length(clusters) < length(mol_biol)) {
      # Get the sample names or indices that were actually used in clustering
      if (!is.null(names(clusters))) {
        # If clusters has names, match by names
        sample_names <- names(clusters)
        mol_biol_aligned <- mol_biol[match(sample_names, colnames(ALL_filtered))]
      } else {
        # If no names, assume first n samples were used
        mol_biol_aligned <- mol_biol[1:length(clusters)]
      }
    } else {
      # Option 2: If mol_biol is shorter, subset clusters
      clusters_aligned <- clusters[1:length(mol_biol)]
      mol_biol_aligned <- mol_biol
      clusters <- clusters_aligned
    }
  } else {
    mol_biol_aligned <- mol_biol
  }
  
  # Verify alignment
  if (length(mol_biol_aligned) != length(clusters)) {
    message("Error: Could not align data lengths. Skipping ", method)
    next
  }
  
  # Create contingency table
  if (method == "dbscan") {
    # Remove noise points (cluster 0) for comparison
    non_noise <- clusters != 0
    if (sum(non_noise) == 0) {
      message("All points classified as noise in DBSCAN. Skipping comparison.\n")
      next
    }
    cont_table <- table(mol_biol_aligned[non_noise], clusters[non_noise])
  } else {
    cont_table <- table(mol_biol_aligned, clusters)
  }
  
  print(cont_table)
  
  # Calculate adjusted rand index for clustering quality
  if (requireNamespace("mclust", quietly = TRUE)) {
    tryCatch({
      if (method == "dbscan") {
        ari <- mclust::adjustedRandIndex(mol_biol_aligned[non_noise], clusters[non_noise])
      } else {
        ari <- mclust::adjustedRandIndex(mol_biol_aligned, clusters)
      }
      message("Adjusted Rand Index: ", round(ari, 3))
    }, error = function(e) {
      message("Could not calculate ARI: ", e$message)
    })
  }
}

## ----enrichment_analysis------------------------------------------------------
# Examine enrichment results
if (length(advanced_results$enrichment_results) > 0) {
  message("Pathway Enrichment Analysis Results:\n")
  
  for (method in names(advanced_results$enrichment_results)) {
    enrich <- advanced_results$enrichment_results[[method]]
    
    if (nrow(enrich$summary) > 0) {
      message("\n", stringr::str_to_title(method), " - Top Enriched Pathways:\n")
      
      # Show top 5 pathways per cluster
      top_pathways <- enrich$summary %>%
        group_by(cluster) %>%
        slice_head(n = 3) %>%
        select(cluster, pathway, p_adjusted, enrichment_ratio, genes_overlap)
      
      print(kable(top_pathways, digits = 4))
    } else {
      message("\n", stringr::str_to_title(method), " - No significant pathways found\n")
    }
  }
}

## ----visualization_analysis---------------------------------------------------
# Display key plots if they were generated
plots_available <- names(advanced_results$plots)
message("Generated plots: ", paste(plots_available, collapse = ", "))

# PCA plot colored by molecular subtype
if ("pca" %in% plots_available) {
  print(advanced_results$plots$pca)
}

# Elbow plot for optimal k selection
if ("elbow" %in% plots_available) {
  print(advanced_results$plots$elbow)
}

## ----sessionInfo--------------------------------------------------------------
sessionInfo()

