context("binomialfreq")
test_that("the binomial frequentist RAR output is", {
  set.seed(100211)
  expect_equal(max(binomialfreq(p_control = 0.7, p_treatment = 0.7, N_total = 200,
                           block_number = 2, simulation = 10)$N_enrolled), 200)
  expect_equal(binomialfreq(p_control = 0.1, p_treatment = 0.5, N_total = 200,
                               block_number = 3, simulation = 10)$power, 1)
  expect_equal(binomialfreq(p_control = 0.1, p_treatment = 0.5, N_total = 100,
                            block_number = 2, simulation = 10, replace = FALSE)$power, 1)
  expect_equal(binomialfreq(p_control = 0.1, p_treatment = 0.8, N_total = 200,
                            block_number = 100, simulation = 10)$power, 1)
  expect_equal(binomialfreq(p_control = 0.1, p_treatment = 0.6, N_total = 200,
                           block_number = 3, simulation = 10,
                           alternative = "greater")$power, 1)
  expect_equal(binomialfreq(p_control = 0.1, p_treatment = 0.6, N_total = 200,
                           block_number = 3, simulation = 10,
                           alternative = "greater", replace = TRUE)$power, 1)
  expect_equal(binomialfreq(p_control = 0.99, p_treatment = 0.01, 120,
                           simulation = 20)$power, 0)
  expect_equal(binomialfreq(p_control = 0.1, p_treatment = 0.8, N_total = 200,
                            block_number = 100, simulation = 10, early_stop = TRUE)$early_stop,
                rep(1, 10))
  expect_equal(binomialfreq(p_control = 0.1, p_treatment = 0.1, N_total = 200,
                           block_number = 3, simulation = 10,
                           alternative = "less", replace = TRUE)$power, 0)
  expect_equal(binomialfreq(p_control = 0.3, p_treatment = 0.3, N_total = 100,
                            block_number = 3, simulation = 10,
                            alternative = "less")$power, 0.1)
  expect_equal(binomialfreq(p_control = 0.3, p_treatment = 0.3, N_total = 400,
                            block_number = 1, simulation = 10,
                            alternative = "less")$power, 0.1)
  expect_equal(binomialfreq(p_control = 0.1, p_treatment = 0.2, N_total = 200,
                           block_number = 4, simulation = 10,
                           alternative = "less", replace = TRUE)$power, 0)
  expect_equal(binomialfreq(p_control = 0.1, p_treatment = 0.2, N_total = 100,
                            block_number = 60, replace = "FALSE", simulation = 1)$N_enrolled, 100)
  expect_equal(min(binomialfreq(p_control = 0.1, p_treatment = 0.01, N_total = 200,
                           block_number = 2, simulation = 10)$p_treatment_estimate), 0)
  expect_equal(min(binomialfreq(p_control = 0.01, p_treatment = 0.2, N_total = 200,
                               block_number = 2, simulation = 10)$p_control_estimate), 0)
  expect_equal(binomialfreq(p_control = 0.1, p_treatment = 0.2, N_total = 1000,
                            block_number = 1, simulation = 1)$N_enrolled, 1000)
  expect_error(binomialfreq(p_control = 1.1, p_treatment = 0.5, N_total = 200,
                           block_number = 3, simulation = 10))
  expect_error(binomialfreq(p_control = 0.1, p_treatment = 1.2, N_total = 200,
                           block_number = 3, simulation = 10))
  expect_error(binomialfreq(p_control = 0.1, p_treatment = 0.2, N_total = -100,
                           block_number = 3, simulation = 10))
  expect_error(binomialfreq(p_control = 0.1, p_treatment = 0.2, N_total = 100,
                           block_number = 3, simulation = 10.2))
  expect_error(binomialfreq(p_control = 0.1, p_treatment = 0.2, N_total = 100,
                           block_number = 3, drift = -0.15, simulation = 10))
  expect_error(binomialfreq(p_control = 0.1, p_treatment = 0.2, N_total = 100,
                           block_number = 3, replace = "YES"))
  expect_error(binomialfreq(p_control = 0.1, p_treatment = 0.2, N_total = 120,
                           conf_int = 1.2))
  expect_error(binomialfreq(p_control = 0.1, p_treatment = 0.2, N_total = 100,
                           block_number = 1.2))
  expect_error(binomialfreq(p_control = 0.1, p_treatment = 0.2, N_total = 100,
                           alternative = "two-sided"))
  expect_error(binomialfreq(p_control = 0.1, p_treatment = 0.2, N_total = 100,
                           alternative = "two-sided"))
  expect_error(binomialfreq(p_control = 0.1, p_treatment = 0.2, N_total = 100,
                           zvalue = c(-1, 1)))
  expect_error(binomialfreq(p_control = 0.1, p_treatment = 0.2, N_total = 100,
                           rand_ratio = c(-1, 1, 2)))
  expect_error(binomialfreq(p_control = 0.1, p_treatment = 0.2, N_total = 100,
                           block_number = 120))
  expect_error(binomialfreq(p_control = 0.1, p_treatment = 0.2, N_total = 100,
                            early_stop = "NO"))
})
