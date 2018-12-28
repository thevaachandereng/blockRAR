context("binomial RAR")
test_that("the binomial RAR output is", {
  set.seed(100211)
  expect_equal(max(binomialRAR(p_control = 0.7, p_treatment = 0.7, N_total = 200,
                           block_number = 2, simulation = 10)$N_enrolled), 200)
  expect_equal(binomialRAR(p_control = 0.1, p_treatment = 0.5, N_total = 200,
                               block_number = 3, simulation = 10)$power, 0)

})
