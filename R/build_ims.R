#' Build an Immune Methylation Score (IMS)
#'
#' @param meth_data A matrix of methylation beta values (CpGs x samples).
#' @param surv_data A data frame with columns 'time' and 'status'.
#' @param immune_genes Character vector of immune gene symbols. If NULL,
#'   uses a default 35-gene immune panel.
#' @param p_threshold Numeric. P-value threshold for CpG selection. Default 0.05.
#' @param hr_range Numeric vector length 2. Allowed HR range. Default c(0.1, 10).
#' @param min_cpgs Integer. Minimum CpGs required. Default 3.
#' @return A list with IMS, selected_cpgs, coefs, model.
#' @examples
#' # Source: Mock data simulated for demonstration purposes.
#' set.seed(123)
#' meth_data <- matrix(runif(100), nrow = 10, ncol = 10)
#' rownames(meth_data) <- paste0("CD274_cg", 1:10)
#' colnames(meth_data) <- paste0("Sample", 1:10)
#' surv_data <- data.frame(
#'   time = runif(10, 10, 100),
#'   status = sample(c(0, 1), 10, replace = TRUE)
#' )
#' rownames(surv_data) <- colnames(meth_data)
#' res <- build_ims(
#'   meth_data = meth_data, surv_data = surv_data,
#'   immune_genes = c("CD274"), p_threshold = 1,
#'   hr_range = c(0, 100), min_cpgs = 2
#' )
#' @export
build_ims <- function(meth_data, surv_data, immune_genes = NULL,
                      p_threshold = 0.05, hr_range = c(0.1, 10),
                      min_cpgs = 3L) {
  if (is.null(immune_genes)) {
    immune_genes <- c('CD274','PDCD1','CTLA4','LAG3','TIGIT','HAVCR2',
      'PDCD1LG2','VSIR','BTLA','IFNG','TNF','IL6','IL10','TGFB1',
      'CXCL9','CXCL10','HLA-A','HLA-B','HLA-C','B2M','TAP1','TAP2',
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
  cox_res <- do.call(rbind, lapply(seq_len(nrow(meth_z)), function(i) {
    df <- data.frame(time = surv_data[['time']], status = surv_data[['status']],
                     meth = as.numeric(meth_z[i, ]))
    fit <- tryCatch(survival::coxph(survival::Surv(time, status) ~ meth, data = df),
                    error = function(e) NULL)
    if (is.null(fit)) return(NULL)
    s <- summary(fit)$coefficients
    data.frame(CpG = rownames(meth_z)[i], coef = s[1, 'coef'],
               HR = exp(s[1, 'coef']), p_value = s[1, 'Pr(>|z|)'],
               stringsAsFactors = FALSE)
  }))
  sel <- cox_res[cox_res$p_value < p_threshold &
                 cox_res$HR > hr_range[1] & cox_res$HR < hr_range[2], , drop = FALSE]
  if (nrow(sel) < min_cpgs) stop('Not enough CpGs passed filter.')
  coefs <- stats::setNames(sel$coef, sel$CpG)
  ims_matrix <- t(meth_z[sel$CpG, , drop = FALSE])
  IMS <- as.numeric(ims_matrix %*% coefs)
  names(IMS) <- colnames(meth_z)
  list(IMS = IMS, selected_cpgs = sel, coefs = coefs, model = cox_res)
}

