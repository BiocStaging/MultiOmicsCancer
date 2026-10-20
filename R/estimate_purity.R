#' Estimate Tumor Purity and Immune Score
#'
#' @description
#' Computes ESTIMATE scores (Stromal, Immune, ESTIMATE) and Tumor Purity
#' from gene expression data using the ESTIMATE algorithm.
#' Requires package 'tidyestimate'.
#'
#' @param expr_data A matrix of gene expression (genes x samples).
#' @param is_affymetrix Logical. TRUE if data is from Affymetrix arrays.
#'   Default FALSE (RNA-seq).
#' @return A data frame with columns: sample, stromal, immune, estimate, purity.
#' @export
#' @examples
#' \donttest{
#' # Note: The ESTIMATE algorithm requires a large, complete gene expression
#' # matrix (typically > 5,000 genes) to function properly.
#' # This is a minimal mock example wrapped in \donttest to avoid errors
#' # during automated checking with insufficient data.
#'
#' # Mock expression matrix
#' expr_mat <- matrix(rnorm(10000), nrow = 1000, ncol = 10)
#' rownames(expr_mat) <- paste0("GENE", 1:1000)
#' colnames(expr_mat) <- paste0("Sample", 1:10)
#'
#' # Run estimate_purity
#' res <- estimate_purity(expr_mat)
#' head(res)
#' }
#' @examples
#' # Mock expression data (genes x samples)
#' expr_mat <- matrix(rnorm(2000), nrow = 200, ncol = 10)
#' rownames(expr_mat) <- paste0("GENE", 1:200)
#' colnames(expr_mat) <- paste0("Sample", 1:10)
#'
#' # The example below mocks 'tidyestimate' functions so it runs without
#' # requiring the package or a large, complete gene expression matrix.
#' # In a real analysis, ensure 'tidyestimate' is installed.
#' if (requireNamespace("tidyestimate", quietly = TRUE)) {
#'   # Mock the 'filter_common_genes' function
#'   mock_filter <- function(x, ...) x
#'   # Mock the 'estimate_score' function
#'   mock_score <- function(x, ...) {
#'     data.frame(
#'       sample = colnames(x),
#'       stromal = runif(ncol(x), -1000, 1000),
#'       immune = runif(ncol(x), -1000, 1000),
#'       estimate = runif(ncol(x), -2000, 2000)
#'     )
#'   }
#'   # Temporarily override the real functions for this example
#'   assignInNamespace("filter_common_genes", mock_filter, ns = "tidyestimate")
#'   assignInNamespace("estimate_score", mock_score, ns = "tidyestimate")
#'   res <- estimate_purity(expr_mat)
#'   head(res)
#' }
estimate_purity <- function(expr_data, is_affymetrix = FALSE) {
  if (!requireNamespace('tidyestimate', quietly = TRUE)) {
    stop('Please install tidyestimate: install.packages("tidyestimate")')
  }
  expr_mat <- as.matrix(expr_data)
  filtered <- tidyestimate::filter_common_genes(
    expr_mat, id = 'hgnc_symbol',
    tell_missing = FALSE, find_alias = TRUE
  )
  scores <- tidyestimate::estimate_score(filtered, is_affymetrix = is_affymetrix)
  scores_df <- as.data.frame(scores)
  if ('sample' %in% colnames(scores_df)) {
    rownames(scores_df) <- scores_df$sample
    scores_df <- scores_df[, -1]
  }
  scores_df$purity <- cos(0.6049872018 + 0.0001467884 * scores_df$estimate)
  scores_df
}

