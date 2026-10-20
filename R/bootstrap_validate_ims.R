#' Nested Bootstrap Validation for IMS
#'
#' @param meth_data Matrix of methylation (CpGs x samples).
#' @param surv_data Data frame with columns 'time' and 'status'.
#' @param immune_genes Character vector of immune gene symbols.
#' @param B Integer. Number of bootstrap replicates. Default 100.
#' @param p_threshold Numeric. P-value cutoff. Default 0.05.
#' @param min_cpgs Integer. Minimum CpGs required. Default 3.
#' @param seed Integer. Random seed. Default 123.
#' @return A list with apparent C-index, optimism, corrected C-index,
#'   and the bootstrap distribution.
#' @export
#' @examples
#' # Mock methylation data (5 CpGs, 15 samples)
#' # CpG names must match the immune gene list format (e.g., Gene_CpG)
#' meth_data <- matrix(runif(75), nrow = 5, ncol = 15)
#' rownames(meth_data) <- c("CD274_cg1", "CD274_cg2", "CD274_cg3", "PDCD1_cg1", "PDCD1_cg2")
#' colnames(meth_data) <- paste0("Sample", 1:15)
#'
#' # Mock survival data
#' surv_data <- data.frame(
#'   time = runif(15, 10, 100),
#'   status = sample(c(0, 1), 15, replace = TRUE)
#' )
#' rownames(surv_data) <- colnames(meth_data)
#'
#' # Run bootstrap (using B=5 and p_threshold=1 for a fast, error-free example)
#' res <- bootstrap_validate_ims(
#'   meth_data = meth_data,
#'   surv_data = surv_data,
#'   immune_genes = c("CD274", "PDCD1"),
#'   B = 5,
#'   p_threshold = 1.0,
#'   min_cpgs = 1
#' )
bootstrap_validate_ims <- function(meth_data, surv_data, immune_genes = NULL,
                                   B = 100L, p_threshold = 0.05,
                                   min_cpgs = 3L, seed = 123) {
  if (is.null(immune_genes)) {
    immune_genes <- c(
      'CD274','PDCD1','CTLA4','LAG3','TIGIT','HAVCR2','PDCD1LG2','VSIR','BTLA',
      'IFNG','TNF','IL6','IL10','TGFB1','CXCL9','CXCL10',
      'HLA-A','HLA-B','HLA-C','B2M','TAP1','TAP2',
      'CD8A','CD8B','CD4','FOXP3','GZMA','GZMB','PRF1',
      'STAT1','IRF1','JAK1','JAK2','IDO1')
  }
  meth_data <- as.matrix(meth_data)
  surv_data <- as.data.frame(surv_data)
  common <- intersect(colnames(meth_data), rownames(surv_data))
  meth_data <- meth_data[, common, drop = FALSE]
  surv_data <- surv_data[common, , drop = FALSE]
  pattern <- paste0('^(', paste(immune_genes, collapse = '|'), ')_')
  immune_cpgs <- rownames(meth_data)[grepl(pattern, rownames(meth_data))]
  meth_sub <- meth_data[immune_cpgs, , drop = FALSE]
  if (anyNA(meth_sub)) {
    for (i in seq_len(nrow(meth_sub))) {
      na_idx <- is.na(meth_sub[i, ])
      if (any(na_idx)) meth_sub[i, na_idx] <- mean(meth_sub[i, ], na.rm = TRUE)
    }
  }
  meth_z <- t(scale(t(meth_sub)))
  surv_data <- surv_data[colnames(meth_z), , drop = FALSE]

  fit_pipeline <- function(idx) {
    sel_cpgs <- character(0); sel_coefs <- numeric(0)
    for (j in seq_len(nrow(meth_z))) {
      df <- data.frame(time = surv_data$time[idx], status = surv_data$status[idx],
                       meth = as.numeric(meth_z[j, idx]))
      fit <- tryCatch(survival::coxph(survival::Surv(time, status) ~ meth, data = df),
                      error = function(e) NULL)
      if (is.null(fit)) next
      s <- summary(fit)$coefficients
      if (s[1, 'Pr(>|z|)'] < p_threshold) {
        hr <- exp(s[1, 'coef'])
        if (hr > 0.1 && hr < 10) {
          sel_cpgs <- c(sel_cpgs, rownames(meth_z)[j])
          sel_coefs <- c(sel_coefs, s[1, 'coef'])
        }
      }
    }
    if (length(sel_cpgs) < min_cpgs) return(NULL)
    list(cpgs = sel_cpgs, coefs = sel_coefs)
  }

  compute_cindex <- function(sel, idx) {
    if (is.null(sel)) return(NA)
    ims <- as.numeric(t(meth_z[sel$cpgs, idx, drop = FALSE]) %*% matrix(sel$coefs, ncol = 1))
    df <- data.frame(time = surv_data$time[idx], status = surv_data$status[idx], IMS = ims)
    fit <- tryCatch(survival::coxph(survival::Surv(time, status) ~ IMS, data = df),
                    error = function(e) NULL)
    if (is.null(fit)) return(NA)
    summary(fit)$concordance[1]
  }

  n <- nrow(surv_data)
  full_sel <- fit_pipeline(seq_len(n))
  apparent <- compute_cindex(full_sel, seq_len(n))

  optimism <- numeric(B)
  for (b in seq_len(B)) {
    boot_idx <- sample(seq_len(n), replace = TRUE)
    sel <- fit_pipeline(boot_idx)
    c_boot <- compute_cindex(sel, boot_idx)
    c_orig <- compute_cindex(sel, seq_len(n))
    optimism[b] <- c_boot - c_orig
  }
  optimism_mean <- mean(optimism, na.rm = TRUE)
  corrected <- apparent - optimism_mean
  list(apparent = apparent, optimism = optimism_mean,
       corrected = corrected, boot_distribution = optimism)
}

