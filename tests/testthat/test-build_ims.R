test_that('build_ims works with valid data', {
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
  res <- build_ims(meth, surv, min_cpgs = 2L)
  expect_type(res, 'list')
  expect_true('IMS' %in% names(res))
  expect_true(length(res$IMS) == n)
})

test_that('build_ims errors with no immune CpGs', {
  meth <- matrix(1:10, nrow = 2, dimnames = list(c('X_cg1', 'Y_cg2'), paste0('S', 1:5)))
  surv <- data.frame(time = 1:5, status = c(1,0,1,0,1), row.names = paste0('S', 1:5))
  expect_error(build_ims(meth, surv, min_cpgs = 2L))
})

