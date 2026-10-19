test_that('bootstrap_validate_ims returns expected list', {
  set.seed(1)
  n <- 100
  samples <- paste0('S', 1:n)
  cpgs <- c('CD8A_cg1', 'PDCD1_cg2', 'TNF_cg3', 'HLA-A_cg4', 'GZMB_cg5')
  meth <- matrix(runif(length(cpgs) * n), nrow = length(cpgs),
                 dimnames = list(cpgs, samples))
  risk <- 2 * meth[1, ] - 1.5 * meth[2, ]
  surv_time <- -log(runif(n)) / (0.1 * exp(risk))
  surv <- data.frame(time = surv_time, status = rbinom(n, 1, 0.7),
                     row.names = samples)
  res <- bootstrap_validate_ims(meth, surv, B = 5, min_cpgs = 2L)
  expect_type(res, 'list')
  expect_true('apparent' %in% names(res))
  expect_true('corrected' %in% names(res))
})

