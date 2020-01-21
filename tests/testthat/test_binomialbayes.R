context("binomialbayes")
test_that("the binomial Bayesian RAR output is", {
  set.seed(100211)
  expect_equal(max(binomialbayes(p_control = 0.7, p_treatment = 0.7, N_total = 200,
                                 block_number = 2, simulation = 10)$N_enrolled), 200)
  expect_equal(binomialbayes(p_control = 0.1, p_treatment = 0.5, N_total = 200,
                            block_number = 3, simulation = 10)$power, 1)
  expect_equal(binomialbayes(p_control = 0.1, p_treatment = 0.5, N_total = 200,
                             block_number = 100, simulation = 10)$early_success, rep(1, 10))
  expect_equal(binomialbayes(p_control = 0.1, p_treatment = 0.5, N_total = 200,
                             block_number = 1, simulation = 10)$power, 1)
  expect_equal(binomialbayes(p_control = 0.1, p_treatment = 0.6, N_total = 200,
                            block_number = 3, simulation = 10,
                            alternative = "greater")$power, 1)
  expect_equal(binomialbayes(p_control = 0.99, p_treatment = 0.01, 120,
                            simulation = 20)$power, 0)
  expect_equal(binomialbayes(p_control = 0.1, p_treatment = 0.2, N_total = 200,
                            block_number = 4, simulation = 10,
                            alternative = "less")$power, 0)
  expect_equal(binomialbayes(p_control = 0.1, p_treatment = 0.2, N_total = 200,
                             block_number = 200, simulation = 10,
                             alternative = "less")$power, 0)
  expect_equal(binomialbayes(p_control = 0.1, p_treatment = 0.1, N_total = 1,
                             block_number = 1, simulation = 10)$power, 0)
  expect_equal(min(binomialbayes(p_control = 0.01, p_treatment = 0.2, N_total = 200,
                                block_number = 2, simulation = 10)$N_enrolled[1]), 100)
  expect_error(binomialbayes(p_control = 1.1, p_treatment = 0.5, N_total = 200,
                            block_number = 3, simulation = 10))
  expect_error(binomialbayes(p_control = 0.1, p_treatment = 1.2, N_total = 200,
                            block_number = 3, simulation = 10))
  expect_error(binomialbayes(p_control = 0.1, p_treatment = 0.2, N_total = -100,
                            block_number = 3, simulation = 10))
  expect_error(binomialbayes(p_control = 0.1, p_treatment = 0.2, N_total = 100,
                            block_number = 3, simulation = 10.2))
  expect_error(binomialbayes(p_control = 0.1, p_treatment = 0.2, N_total = 100,
                            block_number = 3, drift = -0.15, simulation = 10))
  expect_error(binomialbayes(p_control = 0.1, p_treatment = 0.2, N_total = 100,
                            block_number = 3, replace = "YES"))
  expect_error(binomialbayes(p_control = 0.1, p_treatment = 0.2, N_total = 120,
                            conf_int = 1.2))
  expect_error(binomialbayes(p_control = 0.1, p_treatment = 0.2, N_total = 100,
                            block_number = 1.2))
  expect_error(binomialbayes(p_control = 0.1, p_treatment = 0.2, N_total = 100,
                            alternative = "two-sided"))
  expect_error(binomialbayes(p_control = 0.1, p_treatment = 0.2, N_total = 100,
                            alternative = "two-sided"))
  expect_error(binomialbayes(p_control = 0.1, p_treatment = 0.2, N_total = 100,
                            block_number = 120, simulation = 2))
  expect_error(binomialbayes(p_control = 0.1, p_treatment = 0.2, N_total = 100,
                             block_number = 10, simulation = 2, prob_accept_ha = 1.2))
})
