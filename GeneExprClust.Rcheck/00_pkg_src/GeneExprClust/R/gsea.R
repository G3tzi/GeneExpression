#' Perform Over-Representation Analysis (ORA) for Gene Clusters
#'
#' Performs hypergeometric test for pathway enrichment in gene clusters.
#' This is the main enrichment analysis function for clustered gene expression data.
#'
#' @param clusters Named vector of cluster assignments (gene names as names, cluster IDs as values)
#' @param gene_sets Named list of gene sets/pathways (each element is a vector of gene names)
#' @param universe Vector of background genes. If NULL, uses all genes from clusters and gene_sets
#' @param p_adjust_method Method to adjust p-values for multiple testing (default: "BH")
#' @param min_set_size Minimum number of genes required in a gene set to be tested
#' @param max_set_size Maximum number of genes allowed in a gene set to be tested
#' @param p_threshold P-value threshold for significance (used in summary functions)
#'
#' @return List of ORA results per cluster, each containing a data frame with enrichment statistics
#' @importFrom stats phyper p.adjust
#' @examples
#' # Basic ORA analysis
#' clusters <- c(1, 1, 2, 2, 3, 3, 1, 2)
#' names(clusters) <- c("GENE1", "GENE2", "GENE3", "GENE4", "GENE5", 
#'                      "GENE6", "GENE7", "GENE8")
#' gene_sets <- list(
#'   Pathway1 = c("GENE1", "GENE2", "GENE7", "GENE9", "GENE10"),
#'   Pathway2 = c("GENE3", "GENE4", "GENE8", "GENE11", "GENE12")
#' )
#' ora_results <- perform_ora(clusters, gene_sets)
#' @export
perform_ora <- function(clusters, gene_sets, universe = NULL,
                        p_adjust_method = "BH", min_set_size = 5, 
                        max_set_size = 500, p_threshold = 0.05) {

  if (!is.vector(clusters) || is.null(names(clusters))) {
    stop("'clusters' must be a named vector with gene names as names and cluster IDs as values")
  }
  
  if (!is.list(gene_sets) || is.null(names(gene_sets))) {
    stop("'gene_sets' must be a named list of gene vectors")
  }
  
  # Set universe if not provided
  if (is.null(universe)) {
    universe <- unique(c(names(clusters), unlist(gene_sets)))
    message("Universe not provided. Using ", length(universe), " unique genes from clusters and gene sets.")
  }
  
  # Filter gene sets by size
  original_n_sets <- length(gene_sets)
  gene_sets <- gene_sets[sapply(gene_sets, function(x) {
    valid_genes <- intersect(x, universe)
    length(valid_genes) >= min_set_size && length(valid_genes) <= max_set_size
  })]
  
  message("Filtered gene sets: ", original_n_sets, " -> ", length(gene_sets), 
      " (min size: ", min_set_size, ", max size: ", max_set_size, ")\n")
  
  if (length(gene_sets) == 0) {
    stop("No gene sets passed the size filters")
  }
  
  # Get unique clusters
  unique_clusters <- unique(clusters)
  message("Analyzing ", length(unique_clusters), " clusters...\n")
  
  results <- list()
  
  # Perform ORA for each cluster
  for (cluster_id in unique_clusters) {
    message("Processing cluster ", cluster_id, "...")
    
    # Get genes in this cluster that are also in universe
    cluster_genes <- intersect(names(clusters)[clusters == cluster_id], universe)
    
    if (length(cluster_genes) == 0) {
      message(" No valid genes found. Skipping.\n")
      next
    }
    
    ora_results <- data.frame()
    
    # Test enrichment for each pathway
    for (pathway_name in names(gene_sets)) {
      pathway_genes <- intersect(gene_sets[[pathway_name]], universe)
      genes_in_both <- intersect(cluster_genes, pathway_genes)
      
      # Skip if no overlap
      if (length(genes_in_both) == 0) next
      
      # Hypergeometric test parameters
      k <- length(genes_in_both)        # genes in both cluster and pathway
      N <- length(universe)             # total genes in universe
      K <- length(pathway_genes)        # genes in pathway
      n <- length(cluster_genes)        # genes in cluster
      
      # Hypergeometric test (one-tailed, testing for over-representation)
      pval <- phyper(k - 1, K, N - K, n, lower.tail = FALSE)
      
      # Calculate enrichment metrics
      expected_overlap <- (n * K) / N
      enrichment_ratio <- k / expected_overlap
      
      # Store results
      ora_results <- rbind(ora_results, data.frame(
        pathway = pathway_name,
        p_value = pval,
        genes_in_pathway = K,
        genes_in_cluster = n,
        genes_overlap = k,
        expected_overlap = round(expected_overlap, 2),
        enrichment_ratio = round(enrichment_ratio, 3),
        overlap_genes = paste(genes_in_both, collapse = ";"),
        stringsAsFactors = FALSE
      ))
    }
    
    # Adjust p-values if we have results
    if (nrow(ora_results) > 0) {
      ora_results$p_adjusted <- p.adjust(ora_results$p_value, method = p_adjust_method)
      ora_results <- ora_results[order(ora_results$p_value), ]
      
      n_significant <- sum(ora_results$p_adjusted < p_threshold)
      message(" Found ", nrow(ora_results), " pathways tested, ", n_significant, " significant.\n")
    } else {
      message(" No pathways with overlapping genes found.\n")
    }
    
    results[[paste0("cluster_", cluster_id)]] <- ora_results
  }
  
  message("ORA analysis completed.\n")
  return(results)
}

#' Extract Gene Names from ExpressionSet
#'
#' Helper function to extract gene names from Bioconductor ExpressionSet objects.
#'
#' @param eset ExpressionSet object
#' @param id_type Type of gene identifiers to extract ("featureNames", "SYMBOL", "ENTREZID", etc.)
#'
#' @return Character vector of gene names
#' @examples
#' # Create a sample ExpressionSet
#' library(Biobase)
#' 
#' # Create expression matrix
#' expr_matrix <- matrix(rnorm(100), nrow = 10, ncol = 10)
#' rownames(expr_matrix) <- paste0("probe_", 1:10)
#' colnames(expr_matrix) <- paste0("sample_", 1:10)
#' 
#' # Create feature data with gene symbols
#' feature_data <- data.frame(
#'   SYMBOL = c("GENE1", "GENE2", "GENE3", "GENE4", "GENE5", 
#'              "GENE6", "GENE7", "GENE8", "GENE9", "GENE10"),
#'   ENTREZID = paste0("ID_", 1:10),
#'   stringsAsFactors = FALSE
#' )
#' rownames(feature_data) <- rownames(expr_matrix)
#' 
#' # Create phenotype data
#' pheno_data <- data.frame(
#'   condition = rep(c("control", "treatment"), each = 5),
#'   stringsAsFactors = FALSE
#' )
#' rownames(pheno_data) <- colnames(expr_matrix)
#' 
#' # Create ExpressionSet
#' eset <- ExpressionSet(
#'   assayData = expr_matrix,
#'   phenoData = AnnotatedDataFrame(pheno_data),
#'   featureData = AnnotatedDataFrame(feature_data)
#' )
#' 
#' # Extract gene names using different methods
#' probe_ids <- extract_gene_names(eset)  # Default: featureNames
#' gene_symbols <- extract_gene_names(eset, id_type = "SYMBOL")
#' entrez_ids <- extract_gene_names(eset, id_type = "ENTREZID")
#' @export
extract_gene_names <- function(eset, id_type = "featureNames") {
  
  if (!requireNamespace("Biobase", quietly = TRUE)) {
    stop("Biobase package is required. Install with: BiocManager::install('Biobase')")
  }
  
  if (!inherits(eset, "ExpressionSet")) {
    stop("Input must be an ExpressionSet object")
  }
  
  if (id_type == "featureNames") {
    return(Biobase::featureNames(eset))
  } else {
    fdata <- Biobase::fData(eset)
    if (id_type %in% colnames(fdata)) {
      gene_names <- fdata[[id_type]]
      gene_names <- gene_names[!is.na(gene_names)]
      return(as.character(gene_names))
    } else {
      stop("ID type '", id_type, "' not found in feature data. Available columns: ", 
           paste(colnames(fdata), collapse = ", "))
    }
  }
}

#' Perform ORA Analysis from ExpressionSet and Clusters
#'
#' Convenience wrapper function that takes an ExpressionSet and cluster results
#' to perform ORA analysis.
#'
#' @param eset ExpressionSet object containing gene expression data
#' @param clusters Named vector of cluster assignments
#' @param gene_sets Named list of gene sets for enrichment testing
#' @param gene_id_type Type of gene IDs to use from ExpressionSet
#' @param ... Additional arguments passed to perform_ora()
#'
#' @return List of ORA results per cluster
#' @examples
#' \dontrun{
#' # ORA analysis from ExpressionSet
#' ora_results <- ora_from_expression_set(eset, clusters, gene_sets)
#' }
#' @export
ora_from_expression_set <- function(eset, clusters, gene_sets, 
                                   gene_id_type = "featureNames", ...) {
  
  universe_genes <- extract_gene_names(eset, gene_id_type)
  
  if (!all(names(clusters) %in% universe_genes)) {
    warning("Some cluster gene names don't match ExpressionSet gene names. ",
            "Check gene ID consistency.")
  }
  
  perform_ora(clusters = clusters, 
              gene_sets = gene_sets, 
              universe = universe_genes, 
              ...)
}

#' Summarize ORA Results
#'
#' Create a summary table of top enriched pathways across all clusters.
#'
#' @param ora_results Output from perform_ora()
#' @param top_n Number of top results per cluster to include
#' @param p_threshold Adjusted p-value threshold for significance
#' @param min_enrichment_ratio Minimum enrichment ratio threshold
#'
#' @return Data frame summarizing top enrichment results
#' @importFrom utils head
#' @examples
#' \dontrun{
#' # Summarize ORA results
#' summary_table <- summarize_ora_results(ora_results, top_n = 5)
#' }
#' @export
summarize_ora_results <- function(ora_results, top_n = 10, p_threshold = 0.05, 
                                 min_enrichment_ratio = 1.5) {
  
  if (length(ora_results) == 0) {
    return(data.frame())
  }
  
  summary_list <- list()
  
  for (cluster_name in names(ora_results)) {
    df <- ora_results[[cluster_name]]
    
    if (nrow(df) == 0) next
    
    # Filter by significance and enrichment ratio
    significant <- df[df$p_adjusted < p_threshold & 
                     df$enrichment_ratio >= min_enrichment_ratio, ]
    
    if (nrow(significant) > 0) {
      # Take top N results
      top_results <- head(significant, top_n)
      top_results$cluster <- cluster_name
      
      # Select relevant columns for summary
      summary_df <- top_results[, c("cluster", "pathway", "p_adjusted", 
                                   "enrichment_ratio", "genes_overlap", 
                                   "genes_in_pathway", "genes_in_cluster")]
      
      summary_list[[cluster_name]] <- summary_df
    }
  }
  
  if (length(summary_list) == 0) {
    message("No significant enrichments found with current thresholds.\n")
    return(data.frame())
  }
  
  # Combine all results
  combined_summary <- do.call(rbind, summary_list)
  rownames(combined_summary) <- NULL
  
  combined_summary <- combined_summary[order(combined_summary$p_adjusted), ]
  
  return(combined_summary)
}

#' Get Detailed Results for Specific Cluster
#'
#' Extract detailed ORA results for a specific cluster.
#'
#' @param ora_results Output from perform_ora()
#' @param cluster_id Cluster ID to extract results for
#' @param p_threshold P-value threshold for filtering
#' @param include_genes Whether to include the list of overlapping genes
#'
#' @return Data frame of detailed results for the specified cluster
#' @examples
#' \dontrun{
#' # Get results for specific cluster
#' cluster1_results <- get_cluster_results(ora_results, cluster_id = 1)
#' }
#' @export
get_cluster_results <- function(ora_results, cluster_id, p_threshold = 0.05, 
                               include_genes = TRUE) {
  
  cluster_key <- paste0("cluster_", cluster_id)
  
  if (!cluster_key %in% names(ora_results)) {
    stop("Cluster ", cluster_id, " not found in results. Available clusters: ",
         paste(gsub("cluster_", "", names(ora_results)), collapse = ", "))
  }
  
  df <- ora_results[[cluster_key]]
  
  if (nrow(df) == 0) {
    message("No results found for cluster ", cluster_id)
    return(data.frame())
  }
  
  significant <- df[df$p_adjusted < p_threshold, ]
  
  if (!include_genes) {
    significant$overlap_genes <- NULL
  }
  
  return(significant)
}

#' Load Gene Sets from Online Databases
#'
#' Fetch gene sets directly from online databases using msigdbr package.
#' This connects to servers automatically without needing file downloads.
#'
#' @param species Species name (e.g., "Homo sapiens", "Mus musculus")
#' @param collection MSigDB category (e.g., "H" for Hallmark, "C2" for curated, "C5" for GO)
#' @param subcategory Optional subcategory (e.g., "GO:BP" for GO Biological Process)
#' @param min_size Minimum size of gene sets to keep
#' @param max_size Maximum size of gene sets to keep
#'
#' @return Named list of gene sets
#' @examples
#' # Get Hallmark pathways for human
#' hallmark_sets <- load_online_gene_sets("Homo sapiens", "H")
#' @importFrom msigdbr msigdbr
#' @examples
#' # Load Hallmark pathways for human
#' hallmark_sets <- load_online_gene_sets("Homo sapiens", "H")
#' @export
load_online_gene_sets <- function(species = "Homo sapiens", collection = "H", 
                                 subcategory = NULL, min_size = 5, max_size = 500) {
  
  if (!requireNamespace("msigdbr", quietly = TRUE)) {
    stop("msigdbr package is required. Install with: install.packages('msigdbr')")
  }
  
  message("Fetching gene sets from MSigDB for ", species, " collection ", collection, "...\n")
  
  # Fetch from MSigDB
  if (is.null(subcategory)) {
    gene_sets_df <- msigdbr::msigdbr(species = species, collection = collection)
  } else {
    gene_sets_df <- msigdbr::msigdbr(species = species, collection = collection, subcategory = subcategory)
  }
  
  if (nrow(gene_sets_df) == 0) {
    stop("No gene sets found for the specified criteria")
  }
  
  # Convert to named list format
  gene_sets <- split(gene_sets_df$gene_symbol, gene_sets_df$gs_name)
  
  # Filter by size
  original_count <- length(gene_sets)
  gene_sets <- gene_sets[sapply(gene_sets, function(x) length(x) >= min_size & length(x) <= max_size)]
  
  message("Loaded ", length(gene_sets), " gene sets (filtered from ", original_count, ")")
  return(gene_sets)
}

#' Load Gene Sets from GMT File
#'
#' Helper function to load gene sets from GMT (Gene Matrix Transposed) format files.
#' This is a common format for pathway databases like MSigDB.
#'
#' @param gmt_file Path to GMT file
#' @param min_size Minimum size of gene sets to keep
#' @param max_size Maximum size of gene sets to keep
#'
#' @return Named list of gene
#' @examples
#' \dontrun{
#' # Load gene sets from GMT file
#' gene_sets <- load_gmt_gene_sets("pathways.gmt")
#' }
#' @export
load_gmt_gene_sets <- function(gmt_file, min_size = 5, max_size = 500) {
  
  if (!file.exists(gmt_file)) {
    stop("GMT file not found: ", gmt_file)
  }
  
  # Read GMT file
  lines <- readLines(gmt_file)
  gene_sets <- list()
  
  for (line in lines) {
    parts <- strsplit(line, "\t")[[1]]
    
    if (length(parts) < 3) next
    
    pathway_name <- parts[1]
    description <- parts[2]
    genes <- parts[3:length(parts)]
    
    # Remove empty genes
    genes <- genes[genes != "" & !is.na(genes)]
    
    # Filter by size
    if (length(genes) >= min_size && length(genes) <= max_size) {
      gene_sets[[pathway_name]] <- genes
    }
  }
  
  message("Loaded ", length(gene_sets), " gene sets from ", gmt_file)
  return(gene_sets)
}

#' Create Simple Example Gene Sets
#'
#' Creates example gene sets for testing purposes.
#'
#' @param gene_universe Vector of all available gene names
#' @param n_pathways Number of pathways to create
#' @param pathway_size Average size of each pathway
#'
#' @return Named list of example gene sets
#' @importFrom stats rnorm
#' @examples
#' data(eset)
#' # Create example gene sets
#' example_sets <- create_example_gene_sets(rownames(eset), n_pathways = 20)
#' @export
create_example_gene_sets <- function(gene_universe, n_pathways = 50, pathway_size = 30) {
  
  gene_sets <- list()
  
  for (i in 1:n_pathways) {
    n_genes <- max(5, round(rnorm(1, pathway_size, pathway_size * 0.3)))
    n_genes <- min(n_genes, length(gene_universe))
    
    pathway_genes <- sample(gene_universe, n_genes)
    pathway_name <- paste0("Pathway_", sprintf("%03d", i))
    
    gene_sets[[pathway_name]] <- pathway_genes
  }
  
  return(gene_sets)
}

#' Parallel Over-Representation Analysis (ORA) for Gene Clusters
#' @param clusters Named vector of cluster assignments
#' @param gene_sets Named list of gene sets/pathways
#' @param universe Vector of background genes
#' @param p_adjust_method Method to adjust p-values (default: "BH")
#' @param min_set_size Minimum gene set size
#' @param max_set_size Maximum gene set size
#' @param p_threshold P-value threshold for significance
#' @param BPPARAM BiocParallel backend
#' @return List of ORA results per cluster
#' @importFrom stats phyper p.adjust
#' @examples
#' \dontrun{
#' # Parallel ORA analysis
#' ora_results <- perform_ora_parallel(clusters, gene_sets)
#' }
#' @export
perform_ora_parallel <- function(clusters, gene_sets, universe = NULL,
                                p_adjust_method = "BH", min_set_size = 5, 
                                max_set_size = 500, p_threshold = 0.05,
                                BPPARAM = BiocParallel::bpparam()) {
  
  if (!requireNamespace("BiocParallel", quietly = TRUE)) {
    stop("BiocParallel package is required for parallel ORA")
  }
  
  if (!is.vector(clusters) || is.null(names(clusters))) {
    stop("'clusters' must be a named vector")
  }
  
  if (!is.list(gene_sets) || is.null(names(gene_sets))) {
    stop("'gene_sets' must be a named list of gene vectors")
  }
  
  if (is.null(universe)) {
    universe <- unique(c(names(clusters), unlist(gene_sets)))
    message("Universe not provided. Using ", length(universe), " unique genes.")
  }
  
  gene_sets <- gene_sets[sapply(gene_sets, function(x) {
    valid_genes <- intersect(x, universe)
    length(valid_genes) >= min_set_size && length(valid_genes) <= max_set_size
  })]
  
  if (length(gene_sets) == 0) {
    stop("No gene sets passed the size filters")
  }
  
  unique_clusters <- unique(clusters)
  message("Analyzing ", length(unique_clusters), " clusters in parallel...")
  
  # Parallel execution across clusters
  results <- BiocParallel::bplapply(unique_clusters, function(cluster_id) {
    cluster_genes <- intersect(names(clusters)[clusters == cluster_id], universe)
    
    if (length(cluster_genes) == 0) {
      return(data.frame())
    }
    
    # Perform ORA for this cluster
    ora_results <- data.frame()
    
    for (pathway_name in names(gene_sets)) {
      pathway_genes <- intersect(gene_sets[[pathway_name]], universe)
      genes_in_both <- intersect(cluster_genes, pathway_genes)
      
      if (length(genes_in_both) == 0) next
      
      # Hypergeometric test
      k <- length(genes_in_both)
      N <- length(universe)
      K <- length(pathway_genes)
      n <- length(cluster_genes)
      
      pval <- phyper(k - 1, K, N - K, n, lower.tail = FALSE)
      expected_overlap <- (n * K) / N
      enrichment_ratio <- k / expected_overlap
      
      ora_results <- rbind(ora_results, data.frame(
        pathway = pathway_name,
        p_value = pval,
        genes_in_pathway = K,
        genes_in_cluster = n,
        genes_overlap = k,
        expected_overlap = round(expected_overlap, 2),
        enrichment_ratio = round(enrichment_ratio, 3),
        overlap_genes = paste(genes_in_both, collapse = ";"),
        stringsAsFactors = FALSE
      ))
    }
    
    # Adjust p-values
    if (nrow(ora_results) > 0) {
      ora_results$p_adjusted <- p.adjust(ora_results$p_value, method = p_adjust_method)
      ora_results <- ora_results[order(ora_results$p_value), ]
    }
    
    return(ora_results)
  }, BPPARAM = BPPARAM)
  
  names(results) <- paste0("cluster_", unique_clusters)
  
  message("Parallel ORA analysis completed.\n")
  return(results)
}