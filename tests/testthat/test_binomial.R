context("binomial RAR")
test_that("the binomial RAR output is", {
  set.seed(100211)
  expect_equal(max(binomialRAR(p_control = 0.7, p_treatment = 0.7, N_total = 200,
                           block_number = 2, simulation = 10)$N_enrolled), 200)
  expect_equal(binomialRAR(p_control = 0.1, p_treatment = 0.5, N_total = 200,
                               block_number = 3, simulation = 10)$power, 0)
  expect_equal(binomialRAR(p_control = 0.1, p_treatment = 0.6, N_total = 200,
                           block_number = 3, simulation = 10,
                           alternative = "greater")$power, 1)
  expect_error(binomialRAR(p_control = 1.1, p_treatment = 0.5, N_total = 200,
                           block_number = 3, simulation = 10))
  expect_error(binomialRAR(p_control = 0.1, p_treatment = 1.2, N_total = 200,
                           block_number = 3, simulation = 10))
  expect_error(binomialRAR(p_control = 0.1, p_treatment = 1.2, N_total = -100,
                           block_number = 3, simulation = 10))
  expect_error(binomialRAR(p_control = 0.1, p_treatment = 0.2, N_total = -100,
                           block_number = 1.2))
  expect_error(binomialRAR(p_control = 0.1, p_treatment = 0.2, N_total = 100,
                           alternative = "two-sided"))
  expect_error(binomialRAR(p_control = 0.1, p_treatment = 0.2, N_total = 100,
                           alternative = "two-sided"))
  expect_error(binomialRAR(p_control = 0.1, p_treatment = 0.2, N_total = 100,
                           zvalue = c(-1, 1, 2)))
  expect_error(binomialRAR(p_control = 0.1, p_treatment = 0.2, N_total = 100,
                           rand_ratio = c(-1, 1, 2)))
  expect_error(binomialRAR(p_control = 0.1, p_treatment = 0.2, N_total = 100,
                           block_number = 51))

})
