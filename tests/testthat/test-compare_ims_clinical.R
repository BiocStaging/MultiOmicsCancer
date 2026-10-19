test_that('compare_ims_clinical returns expected values', {
  set.seed(1)
  n <- 100
  samples <- paste0('S', 1:n)
  ims <- rnorm(n)
  names(ims) <- samples
  stage <- sample(c('Early', 'Late'), n, replace = TRUE, prob = c(0.7, 0.3))
  names(stage) <- samples
  risk <- 0.5 * ims + ifelse(stage == 'Late', 0.8, 0)
  surv_time <- -log(runif(n)) / (0.1 * exp(risk))
  surv <- data.frame(time = surv_time, status = rbinom(n, 1, 0.7),
                     row.names = samples)
  res <- compare_ims_clinical(ims, surv, stage, var_name = 'Stage')
  expect_type(res, 'list')
  expect_true('delta_c' %in% names(res))
  expect_true('lrt_pvalue' %in% names(res))
})

