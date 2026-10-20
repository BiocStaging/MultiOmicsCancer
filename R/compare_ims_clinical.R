#' Compare IMS with Clinical Variables
#'
#' @description
#' Compares prognostic performance of IMS vs a clinical variable
#' (e.g., Stage) using Cox models and likelihood ratio test.
#'
#' @param ims_values Named numeric vector of IMS values (from build_ims).
#' @param surv_data Data frame with columns 'time' and 'status'.
#' @param clinical_var Named vector or factor of clinical variable
#'   (e.g., Stage_binary) aligned with surv_data rownames.
#' @param var_name Character. Name of clinical variable for output labels.
#' @return A list with C-index for IMS alone, clinical alone, both,
#'   Delta C-index, and LRT p-value.
#' @examples
#' # Source: Mock data simulated for demonstration purposes.
#' set.seed(123)
#' n <- 60
#' ims_values <- setNames(rnorm(n), paste0("Sample", 1:n))
#' surv_data <- data.frame(
#'   time = runif(n, 10, 100),
#'   status = sample(c(0, 1), n, replace = TRUE),
#'   row.names = names(ims_values)
#' )
#' clinical_var <- setNames(
#'   sample(c("Low", "High"), n, replace = TRUE),
#'   names(ims_values)
#' )
#' res <- compare_ims_clinical(
#'   ims_values = ims_values, surv_data = surv_data,
#'   clinical_var = clinical_var, var_name = "Stage"
#' )
#' @export
compare_ims_clinical <- function(ims_values, surv_data, clinical_var,
                                 var_name = 'Clinical') {
  common <- intersect(names(ims_values), rownames(surv_data))
  df <- data.frame(
    time = surv_data[common, 'time'],
    status = surv_data[common, 'status'],
    IMS = ims_values[common],
    Clinical = clinical_var[common],
    stringsAsFactors = FALSE
  )
  df <- df[stats::complete.cases(df), ]
  if (nrow(df) < 50) stop('Not enough complete samples (< 50).')
  if (length(unique(df$Clinical)) < 2) {
    stop('Clinical variable must have at least 2 levels.')
  }

  m_ims <- survival::coxph(survival::Surv(time, status) ~ IMS, data = df)
  m_clin <- survival::coxph(survival::Surv(time, status) ~ Clinical, data = df)
  m_both <- survival::coxph(survival::Surv(time, status) ~ IMS + Clinical, data = df)

  c_ims <- summary(m_ims)$concordance[1]
  c_clin <- summary(m_clin)$concordance[1]
  c_both <- summary(m_both)$concordance[1]
  delta_c <- c_both - c_clin

  lrt <- tryCatch({
    a <- stats::anova(m_clin, m_both, test = 'LRT')
    a$`Pr(>|Chi|)`[2]
  }, error = function(e) NA)

  list(
    var_name = var_name,
    cindex_ims = c_ims,
    cindex_clinical = c_clin,
    cindex_both = c_both,
    delta_c = delta_c,
    lrt_pvalue = lrt,
    n = nrow(df),
    cox_ims = summary(m_ims)$coefficients,
    cox_clinical = summary(m_clin)$coefficients,
    cox_both = summary(m_both)$coefficients
  )
}

