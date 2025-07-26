#' Complete Gene Expression Analysis Pipeline
#'
#' This function provides a comprehensive workflow for gene expression analysis,
#' including data loading, preprocessing, normalization, exploratory analysis,
#' clustering, visualization, and pathway enrichment analysis.
#'
#' @param data_path Character, path to expression data file or ExpressionSet .RData file
#' @param expression_set The ExpressionSet object to be analyzed
#' @param data_type Character, type of input data: "csv", "tsv", or "eset"
#' @param pheno_data Optional data frame or file path to sample metadata
#' @param output_dir Character, directory to save results and plots (default: "pipeline_results")
#' 
#' # Normalization parameters
#' @param normalize Logical, whether to perform normalization (default: TRUE)
#' @param norm_method Character, normalization method: "log2", "quantile", "zscore", "deseq"
#' @param batch_correction Character or vector, batch variable for correction (default: NULL)
#' @param filter_low_expr Logical, whether to filter low-expressed genes (default: TRUE)
#' @param filter_threshold Numeric, minimum expression threshold (default: 1)
#' @param filter_min_samples Integer, minimum samples meeting threshold (default: 3)
#' 
#' # Exploratory analysis parameters
#' @param detect_outliers Logical, whether to detect and remove outlier genes (default: TRUE)
#' @param outlier_method Character, method for outlier detection: "mad" or "iqr"
#' @param outlier_threshold Numeric, threshold multiplier for outlier detection (default: 3)
#' @param remove_outliers Logical, whether to remove detected outliers (default: TRUE)
#' 
#' # Clustering parameters
#' @param clustering_methods Character vector, methods to use: "hierarchical", "kmeans", "dbscan"
#' @param target Character, cluster "gene" (default) or "sample"
#' @param k_clusters Integer, number of clusters for hierarchical/kmeans (default: 5)
#' @param max_k Integer, maximum k for elbow plot (default: 15)
#' @param hclust_method Character, linkage method for hierarchical clustering (default: "ward.D2")
#' @param distance Character, distance metric (default: "euclidean")
#' @param kmeans_nstart Integer, number of random starts for k-means (default: 25)
#' @param dbscan_eps Numeric, DBSCAN epsilon parameter (required if using DBSCAN)
#' @param dbscan_minPts Integer, DBSCAN minimum points parameter (default: 5)
#' 
#' # Validation parameters
#' @param perform_validation Logical, whether to perform clustering validation (default: TRUE)
#' @param gap_statistic Logical, whether to compute gap statistic (default: TRUE)
#' @param gap_max_k Integer, maximum k for gap statistic (default: 10)
#' @param gap_B Integer, bootstrap samples for gap statistic (default: 20)
#' 
#' # Visualization parameters
#' @param create_plots Logical, whether to create visualization plots (default: TRUE)
#' @param heatmap_scale Character, scaling for heatmap: "none", "row", "column" (default: "row")
#' @param pca_color_by Character, phenotype variable to color PCA by (default: NULL)
#' @param save_plots Logical, whether to save plots to files (default: TRUE)
#' 
#' # Enrichment analysis parameters
#' @param perform_enrichment Logical, whether to perform pathway enrichment (default: TRUE)
#' @param gene_sets_source Character, source for gene sets: "online", "gmt", "example", or "custom"
#' @param gene_sets_collection Character, MSigDB collection for online source (default: "H")
#' @param gene_sets_species Character, species for online gene sets (default: "Homo sapiens")
#' @param gene_sets_file Character, path to GMT file if using "gmt" source
#' @param custom_gene_sets List, custom gene sets if using "custom" source
#' @param enrichment_p_threshold Numeric, p-value threshold for enrichment (default: 0.05)
#' @param min_gene_set_size Integer, minimum gene set size (default: 5)
#' @param max_gene_set_size Integer, maximum gene set size (default: 500)
#' 
#' # General parameters
#' @param verbose Logical, whether to print progress messages (default: TRUE)
#' 
#' @return A list containing all analysis results:
#' \itemize{
#'   \item original_data - Original ExpressionSet
#'   \item processed_data - Processed ExpressionSet after normalization/filtering
#'   \item outlier_info - Information about detected outliers
#'   \item clustering_results - Results from all clustering methods
#'   \item validation_results - Clustering validation metrics
#'   \item plots - List of generated plots
#'   \item enrichment_results - Pathway enrichment analysis results
#'   \item summary - Summary of the analysis
#' }
#' 
#' @examples
#' data(eset)
#' results <- run_pipeline(
#'   expression_set = eset,
#'   clustering_methods = c("hierarchical", "kmeans"),
#'   k_clusters = 4
#' )
#' @importFrom stringr str_to_title
#' @importFrom grDevices dev.off png
#' @importFrom utils write.csv
#' @export
run_pipeline <- function(
  # Data input parameters
  data_path = NULL,
  expression_set = NULL,
  data_type = c("csv", "tsv", "eset"),
  pheno_data = NULL,
  output_dir = "pipeline_results",
  
  # Normalization parameters
  normalize = TRUE,
  norm_method = c("log2", "quantile", "zscore", "deseq"),
  batch_correction = NULL,
  filter_low_expr = TRUE,
  filter_threshold = 1,
  filter_min_samples = 3,
  
  # Exploratory analysis parameters
  detect_outliers = TRUE,
  outlier_method = c("mad", "iqr"),
  outlier_threshold = 3,
  remove_outliers = TRUE,
  
  # Clustering parameters
  clustering_methods = c("hierarchical", "kmeans"),
  target = c("gene", "sample"),
  k_clusters = 5,
  max_k = 15,
  hclust_method = "ward.D2",
  distance = "euclidean",
  kmeans_nstart = 25,
  dbscan_eps = NULL,
  dbscan_minPts = 5,
  
  # Validation parameters
  perform_validation = TRUE,
  gap_statistic = TRUE,
  gap_max_k = 10,
  gap_B = 20,
  
  # Visualization parameters
  create_plots = TRUE,
  heatmap_scale = c("row", "none", "column"),
  pca_color_by = NULL,
  save_plots = TRUE,
  
  # Enrichment analysis parameters
  perform_enrichment = TRUE,
  gene_sets_source = c("online", "gmt", "example", "custom"),
  gene_sets_collection = "H",
  gene_sets_species = "Homo sapiens",
  gene_sets_file = NULL,
  custom_gene_sets = NULL,
  enrichment_p_threshold = 0.05,
  min_gene_set_size = 5,
  max_gene_set_size = 500,
  
  # General parameters
  verbose = TRUE
) {
  
  data_type <- match.arg(data_type)
  norm_method <- match.arg(norm_method)
  outlier_method <- match.arg(outlier_method)
  target <- match.arg(target)
  heatmap_scale <- match.arg(heatmap_scale)
  gene_sets_source <- match.arg(gene_sets_source)
  
  if (is.null(data_path) && is.null(expression_set)) {
    stop("Either 'data_path' or 'expression_set' must be provided")
  }

  if (!is.null(data_path) && !is.null(expression_set)) {
    stop("Provide either 'data_path' OR 'expression_set', not both")
  }
  
  if (!is.null(expression_set) && !inherits(expression_set, "ExpressionSet")) {
    stop("'expression_set' must be an ExpressionSet object")
  }

  valid_methods <- c("hierarchical", "kmeans", "dbscan")
  if (!all(clustering_methods %in% valid_methods)) {
    stop("Invalid clustering methods. Must be one or more of: ", paste(valid_methods, collapse = ", "))
  }
  
  if ("dbscan" %in% clustering_methods && is.null(dbscan_eps)) {
    stop("dbscan_eps parameter is required when using DBSCAN clustering")
  }
  
  if (save_plots && !dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    if (verbose) message("Created output directory: ", output_dir)
  }
  
  results <- list()
  
  # ============================================================================
  # STEP 1: DATA LOADING AND INITIAL PROCESSING
  # ============================================================================
  
  if (verbose) message("\n=== STEP 1: Loading and Initial Processing ===\n")
  
  if (!is.null(expression_set)) {
    if (verbose) message("Using provided ExpressionSet object\n")
    original_eset <- expression_set
    data_source <- "provided ExpressionSet"
  } else {
    if (verbose) message("Loading expression data from: ", data_path)
    original_eset <- read_expression(data_path, type = data_type, pheno_data = pheno_data)
    data_source <- data_path
  }
  
  if (verbose) {
    message("Loaded expression data:\n")
    message("  Genes: ", nrow(original_eset))
    message("  Samples: ", ncol(original_eset))
  }
  
  results$original_data <- original_eset
  processed_eset <- original_eset
  
  # ============================================================================
  # STEP 2: NORMALIZATION AND FILTERING
  # ============================================================================
  
  if (normalize) {
    if (verbose) message("\n=== STEP 2: Normalization and Filtering ===\n")
    
    if (verbose) message("Applying ", norm_method, " normalization...\n")
    
    processed_eset <- normalize_data(
      processed_eset,
      method = norm_method,
      batch = batch_correction,
      filter = filter_low_expr,
      filter_threshold = filter_threshold,
      filter_min_samples = filter_min_samples
    )
    
    if (verbose) {
      message("After normalization and filtering:\n")
      message("  Genes: ", nrow(processed_eset))
      message("  Samples: ", ncol(processed_eset))
    }
  }
  
  results$processed_data <- processed_eset
  
  # ============================================================================
  # STEP 3: EXPLORATORY ANALYSIS AND OUTLIER DETECTION
  # ============================================================================
  
  outlier_info <- list()
  if (detect_outliers) {
    if (verbose) message("\n=== STEP 3: Exploratory Analysis and Outlier Detection ===\n")
    
    processed_eset <- detect_outliers(
      processed_eset,
      method = outlier_method,
      threshold = outlier_threshold,
      plot = create_plots && save_plots,
      remove = remove_outliers
    )
    
    outlier_genes <- attr(processed_eset, "outliers")
    outlier_info$method <- outlier_method
    outlier_info$threshold <- outlier_threshold
    outlier_info$outlier_genes <- outlier_genes
    outlier_info$n_outliers <- length(outlier_genes)
    
    if (verbose) {
      message("Outlier detection using ", outlier_method, " method:")
      message("  Outliers detected: ", length(outlier_genes))
      if (remove_outliers) {
        message("  Final dataset - Genes: ", nrow(processed_eset), " Samples: ", ncol(processed_eset))
      }
    }
  }
  
  results$outlier_info <- outlier_info
  results$processed_data <- processed_eset
  
  # ============================================================================
  # STEP 4: CLUSTERING ANALYSIS
  # ============================================================================
  
  if (verbose) message("\n=== STEP 4: Clustering Analysis ===\n")
  
  clustering_results <- list()
  
  for (method in clustering_methods) {
    if (verbose) message("Running ", method, " clustering...")
    
    if (method == "hierarchical") {
      result <- cluster_expression(
        processed_eset,
        method = "hierarchical",
        k = k_clusters,
        target = target,
        hclust_method = hclust_method,
        distance = distance
      )
    } else if (method == "kmeans") {
      result <- cluster_expression(
        processed_eset,
        method = "kmeans",
        k = k_clusters,
        target = target,
        nstart = kmeans_nstart
      )
    } else if (method == "dbscan") {
      result <- cluster_expression(
        processed_eset,
        method = "dbscan",
        target = target,
        eps = dbscan_eps,
        minPts = dbscan_minPts,
        distance = distance
      )
    }
    
    clustering_results[[method]] <- result
    
    if (verbose) {
      n_clusters <- length(unique(result$clusters))
      if (method == "dbscan") {
        n_noise <- sum(result$clusters == 0)
        n_clusters <- n_clusters - (n_noise > 0)  # Don't count noise as cluster
        message("  Found ", n_clusters, " clusters and ", n_noise, " noise points")
      } else {
        message("  Created ", n_clusters, " clusters")
      }
    }
  }
  
  results$clustering_results <- clustering_results
  
  # ============================================================================
  # STEP 5: CLUSTERING VALIDATION
  # ============================================================================
  
  validation_results <- list()
  if (perform_validation && length(clustering_results) > 0) {
    if (verbose) message("\n=== STEP 5: Clustering Validation ===\n")
    
    expr_data <- if (target == "gene") {
      t(exprs(processed_eset))
    } else {
      exprs(processed_eset)
    }
    
    # Compute validation metrics for each clustering method
    for (method_name in names(clustering_results)) {
      if (verbose) message("Validating ", method_name, " clustering...")
      
      clusters <- clustering_results[[method_name]]$clusters
      method_validation <- list()
      
      # WSS score
      tryCatch({
        method_validation$wss <- wss_score(expr_data, clusters)
      }, error = function(e) {
        if (verbose) message("  Warning: Could not compute WSS - ", e$message)
      })
      
      # Silhouette analysis
      if (length(unique(clusters)) >= 2) {
        tryCatch({
          sil_result <- silhouette_score(expr_data, clusters, distance)
          method_validation$silhouette <- sil_result
          method_validation$avg_silhouette <- mean(sil_result[, "sil_width"])
        }, error = function(e) {
          if (verbose) message("  Warning: Could not compute silhouette - ", e$message)
        })
      }
      
      validation_results[[method_name]] <- method_validation
    }
    
    # Gap statistic for k-means
    if (gap_statistic && target == "gene") {
      if (verbose) message("Computing gap statistic...\n")
      tryCatch({
        gap_result <- gap_statistic(t(exprs(processed_eset)), max_k = gap_max_k, B = gap_B)
        validation_results$gap_statistic <- gap_result
      }, error = function(e) {
        if (verbose) message("Warning: Could not compute gap statistic - ", e$message)
      })
    }
  }
  
  results$validation_results <- validation_results
  
  # ============================================================================
  # STEP 6: VISUALIZATION
  # ============================================================================
  
  plots <- list()
  if (create_plots) {
    if (verbose) message("\n=== STEP 6: Creating Visualizations ===\n")
    
    # PCA plot
    if (verbose) message("Creating PCA plot...\n")
    tryCatch({
      plots$pca <- plot_pca(processed_eset, color_by = pca_color_by, main = "PCA Analysis")
      if (save_plots) {
        ggsave(file.path(output_dir, "pca_plot.png"), plots$pca, width = 8, height = 6)
      }
    }, error = function(e) {
      if (verbose) message("Warning: Could not create PCA plot - ", e$message)
    })
    
    # Elbow plot for k-means
    if ("kmeans" %in% clustering_methods) {
      if (verbose) message("Creating elbow plot...\n")
      tryCatch({
        plots$elbow <- plot_elbow(processed_eset, max_k = max_k, target = target)
        if (save_plots) {
          ggsave(file.path(output_dir, "elbow_plot.png"), plots$elbow, width = 8, height = 6)
        }
      }, error = function(e) {
        if (verbose) message("Warning: Could not create elbow plot - ", e$message)
      })
    }
    
    # Clustering-specific plots
    for (method_name in names(clustering_results)) {
      method_result <- clustering_results[[method_name]]
      
      # Heatmap with cluster annotation
      if (verbose) message("Creating ", method_name, " heatmap...")
      tryCatch({
        plot_title <- paste(stringr::str_to_title(method_name), "Clustering Heatmap")
        plots[[paste0(method_name, "_heatmap")]] <- plot_heatmap(
          processed_eset,
          clusters = method_result$clusters,
          main = plot_title,
          scale = heatmap_scale
        )
        if (save_plots) {
          png(file.path(output_dir, paste0(method_name, "_heatmap.png")), 
              width = 1200, height = 800, res = 150)
          print(plots[[paste0(method_name, "_heatmap")]])
          dev.off()
        }
      }, error = function(e) {
        if (verbose) message("Warning: Could not create ", method_name, " heatmap - ", e$message)
      })
      
      # Dendrogram for hierarchical clustering
      if (method_name == "hierarchical") {
        if (verbose) message("Creating dendrogram...\n")
        tryCatch({
          if (save_plots) {
            plots[[method_name]] <- plot_dendrogram(method_result$model,
                                                    main = "Hierarchical Clustering Dendrogram")
            
            ggsave(file.path(output_dir, paste0(method_name, ".png")),
                   plots[[method_name]], width = 10, height = 6)
          }
        }, error = function(e) {
          if (verbose) message("Warning: Could not create dendrogram - ", e$message)
        })
      }
      
      # Cluster profiles
      if (target == "gene") {
        if (verbose) message("Creating ", method_name, " cluster profiles...")
        tryCatch({
          plots[[paste0(method_name, "_profiles")]] <- plot_cluster_profiles(
            processed_eset, 
            method_result$clusters, 
            scale = TRUE
          )
          if (save_plots) {
            ggsave(file.path(output_dir, paste0(method_name, "_profiles.png")), 
                   plots[[paste0(method_name, "_profiles")]], width = 10, height = 6)
          }
        }, error = function(e) {
          if (verbose) message("Warning: Could not create cluster profiles - ", e$message)
        })
      }
      
      # DBSCAN-specific visualization
      if (method_name == "dbscan") {
        if (verbose) message("Creating DBSCAN results plot...\n")
        tryCatch({
          plots$dbscan_pca <- plot_dbscan_results(processed_eset, method_result$clusters, target)
          if (save_plots) {
            ggsave(file.path(output_dir, "dbscan_pca.png"), plots$dbscan_pca, width = 8, height = 6)
          }
        }, error = function(e) {
          if (verbose) message("Warning: Could not create DBSCAN plot - ", e$message)
        })
      }
    }
  }
  
  results$plots <- plots
  
  # ============================================================================
  # STEP 7: PATHWAY ENRICHMENT ANALYSIS
  # ============================================================================
  
  enrichment_results <- list()
  if (perform_enrichment && target == "gene" && length(clustering_results) > 0) {
    if (verbose) message("\n=== STEP 7: Pathway Enrichment Analysis ===\n")
    
    # Load gene sets based on source
    gene_sets <- NULL
    if (gene_sets_source == "online") {
      if (verbose) message("Loading gene sets from MSigDB...\n")
      tryCatch({
        gene_sets <- load_online_gene_sets(
          species = gene_sets_species,
          collection = gene_sets_collection,
          min_size = min_gene_set_size,
          max_size = max_gene_set_size
        )
      }, error = function(e) {
        if (verbose) message("Warning: Could not load online gene sets - ", e$message)
        if (verbose) message("Falling back to example gene sets...\n")
        gene_sets <- create_example_gene_sets(
          gene_universe = featureNames(processed_eset),
          n_pathways = 50,
          pathway_size = 30
        )
      })
    } else if (gene_sets_source == "gmt") {
      if (is.null(gene_sets_file)) {
        stop("gene_sets_file must be provided when using 'gmt' source")
      }
      if (verbose) message("Loading gene sets from GMT file: ", gene_sets_file)
      gene_sets <- load_gmt_gene_sets(gene_sets_file, min_gene_set_size, max_gene_set_size)
    } else if (gene_sets_source == "example") {
      if (verbose) message("Creating example gene sets...\n")
      gene_sets <- create_example_gene_sets(
        gene_universe = featureNames(processed_eset),
        n_pathways = 50,
        pathway_size = 30
      )
    } else if (gene_sets_source == "custom") {
      if (is.null(custom_gene_sets)) {
        stop("custom_gene_sets must be provided when using 'custom' source")
      }
      gene_sets <- custom_gene_sets
      if (verbose) message("Using provided custom gene sets...\n")
    }
    
    # Perform enrichment analysis for each clustering method
    if (!is.null(gene_sets)) {
      for (method_name in names(clustering_results)) {
        if (verbose) message("Performing enrichment analysis for ", method_name, " clusters...\n")
        
        clusters <- clustering_results[[method_name]]$clusters
        
        # Remove noise points for DBSCAN
        if (method_name == "dbscan") {
          clusters <- clusters[clusters != 0]
        }
        
        if (length(clusters) > 0) {
          tryCatch({
            ora_result <- perform_ora(
              clusters = clusters,
              gene_sets = gene_sets,
              universe = featureNames(processed_eset),
              p_threshold = enrichment_p_threshold,
              min_set_size = min_gene_set_size,
              max_set_size = max_gene_set_size
            )
            
            # Create summary
            summary_result <- summarize_ora_results(
              ora_result,
              top_n = 10,
              p_threshold = enrichment_p_threshold
            )
            
            enrichment_results[[method_name]] <- list(
              detailed_results = ora_result,
              summary = summary_result
            )
            
            # Save enrichment results
            if (save_plots && nrow(summary_result) > 0) {
              write.csv(summary_result, 
                       file.path(output_dir, paste0(method_name, "_enrichment_summary.csv")),
                       row.names = FALSE)
            }
            
          }, error = function(e) {
            if (verbose) message("Warning: Enrichment analysis failed for ", method_name, " - ", e$message)
          })
        }
      }
    }
  }
  
  results$enrichment_results <- enrichment_results
  
  # ============================================================================
  # STEP 8: GENERATE SUMMARY
  # ============================================================================
  
  if (verbose) message("\n=== STEP 8: Generating Summary ===\n")
  
  summary_info <- list(
    input_data = list(
      source = data_source,
      type = if(!is.null(data_path)) data_type else "ExpressionSet",
      original_genes = nrow(original_eset),
      original_samples = ncol(original_eset),
      final_genes = nrow(processed_eset),
      final_samples = ncol(processed_eset)
    ),
    processing = list(
      normalization = ifelse(normalize, norm_method, "none"),
      batch_correction = !is.null(batch_correction),
      filtering = filter_low_expr,
      outlier_detection = detect_outliers,
      outliers_removed = if(detect_outliers) length(outlier_info$outlier_genes) else 0
    ),
    clustering = list(
      methods = clustering_methods,
      target = target,
      k_clusters = k_clusters
    ),
    validation = list(
      performed = perform_validation,
      gap_statistic = gap_statistic
    ),
    enrichment = list(
      performed = perform_enrichment,
      gene_sets_source = if(perform_enrichment) gene_sets_source else "none",
      methods_analyzed = names(enrichment_results)
    )
  )
  
  results$summary <- summary_info
  
  # Save complete summary
  if (save_plots) {
    saveRDS(results, file.path(output_dir, "complete_results.rds"))
    
    # Create text summary
    summary_text <- paste0(
      "Gene Expression Analysis Pipeline Summary\n",
      "========================================\n\n",
      "Input Data:\n",
      "  - Source: ", data_source,
    "  - Type: ", if(!is.null(data_path)) data_type else "ExpressionSet",
      "  - Original dimensions: ", nrow(original_eset), " genes x ", ncol(original_eset), " samples\n",
      "  - Final dimensions: ", nrow(processed_eset), " genes x ", ncol(processed_eset), " samples\n\n",
      "Processing:\n",
      "  - Normalization: ", ifelse(normalize, norm_method, "none"),
      "  - Batch correction: ", ifelse(!is.null(batch_correction), "yes", "no"),
      "  - Low expression filtering: ", ifelse(filter_low_expr, "yes", "no"),
      "  - Outlier detection: ", ifelse(detect_outliers, paste(outlier_method, "method"), "no"),
      "  - Outliers removed: ", ifelse(detect_outliers, length(outlier_info$outlier_genes), 0), "\n\n",
      "Clustering:\n",
      "  - Methods: ", paste(clustering_methods, collapse = ", "),
      "  - Target: ", target,
      "  - Number of clusters: ", k_clusters, "\n\n",
      "Validation: ", ifelse(perform_validation, "performed", "not performed"),
      "Enrichment Analysis: ", ifelse(perform_enrichment, paste("performed using", gene_sets_source, "gene sets"), "not performed"), "\n\n",
      "Output saved to: ", output_dir
    )
    
    writeLines(summary_text, file.path(output_dir, "analysis_summary.txt"))
  }
  
  if (verbose) {
    message("\n=== PIPELINE COMPLETED SUCCESSFULLY ===\n")
    message("Results saved to: ", output_dir)
    message("Total clustering methods analyzed: ", length(clustering_results))
    if (perform_enrichment) {
      message("Enrichment analysis completed for: ", length(enrichment_results), " methods\n")
    }
  }
  
  return(results)
}