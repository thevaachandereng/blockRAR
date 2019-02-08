#' @title Block Design for Response-Adaptive Randomization for Binomial Data
#'
#' @description Simulation for binomial counts for block design for
#'    response-adaptive randomization with time as a confounding
#'
#' @inheritParams binomialfreq
#' @param number_mcmc scalar. Number of Monte Carlo Markov Chain draws in
#'   sampling posterior.
#' @param futility_prob scalar. Probability of stopping early for futility.
#' @param early_success_prob scalar. Probability of stopping early for success.
#' @param prob_accept_ha scalar. Probability of accepting
#'   alternative hypothesis.
#'
#' @return a list with details on the simulation.
#' \describe{
#'   \item{\code{power}}{
#'     scalar. The power of the trial, ie. the proportion of success over the
#'     number of simulation ran.}
#'   \item{\code{p_control_estimate}}{
#'     scalar. The estimated proportion of events under the control group.}
#'   \item{\code{p_treatment_estimate}}{
#'     scalar. The estimated proportion of events under the treatment group.}
#'   \item{\code{N_enrolled}}{
#'     vector. The number of patients enrolled in the trial (sum of control
#'     and experimental group for each simulation. )}
#'   \item{\code{N_control}}{
#'     vector. The number of patients enrolled in the control group for
#'     each simulation.}
#'   \item{\code{N_control}}{
#'     vector. The number of patients enrolled in the experimental group for
#'     each simulation.}
#' }
#'
#' @importFrom stats rbinom binomial coef
#' @importFrom arm bayesglm sim
#' @importFrom bayesDP bdpbinomial
#' @importFrom dplyr mutate group_by summarize
#' @importFrom tibble as.tibble
#'
#' @export binomialbayes
#'

binomialbayes <- function(
  p_control,
  p_treatment,
  N_total,
  block_number       = 4,
  drift              = 0,
  simulation         = 10000,
  a0                 = 0.5,
  b0                 = 0.5,
  number_mcmc        = 10000,
  prob_accept_ha     = 0.95,
  early_success_prob = 0.99,
  futility_prob      = 0.10,
  alternative        = "greater"
  ){
  # stop if proportion of control is not between 0 and 1
  if((p_control <= 0 | p_control >= 1)){
    stop("The proportion of event for the control group needs to between 0 and 1!")
  }

  # stop if proportion of treatment is not between 0 and 1
  if((p_treatment <= 0 | p_treatment >= 1)){
    stop("The proportion of event for the treatment group needs to between 0 and 1!")
  }

  ## make sure sample size is an integer!
  if((N_total <= 0 | N_total %% 1 != 0)){
    stop("The sample size needs to be a positive integer!")
  }

  #
  if((block_number <= 0 | block_number %% 1 != 0)){
    stop("The number of blocks needs to be a positve integer!")
  }

  if((N_total / block_number <= 2) & replace == FALSE){
    warning("The sampling is done with replacement and replace input is ignored!")
  }

  if((simulation <= 0 | simulation %% 1 != 0)){
    stop("The number of simulation needs to be a positve integer!")
  }

  if((prob_accept_ha <= 0 |prob_accept_ha >= 1)){
    stop("The confidence interval needs to between 0 and 1!")
  }

  if((alternative != "less" & alternative != "greater")){
    stop("The alternative can only be less or greater!")
  }


  if(drift + p_control >= 1 | drift + p_control <= 0 |
     drift + p_treatment >= 1 | drift + p_treatment <= 0){
    stop("The drift value is too high causing the proportion of event to exceed 1
         in either the control or treatment group, pick a lower value for drift!")
  }

  group <- rep(floor(N_total / block_number), block_number)
  if((N_total - sum(group)) > 0){
    index <- sample(1:block_number, N_total - sum(group))
    group[index] <- group[index] + 1
  }

  power              <- 0
  N_control          <- NULL
  N_treatment        <- NULL
  sample_size        <- NULL
  prop_diff_estimate <- NULL

  for(k in 1:simulation){
    data_total            <- NULL
    test_stat             <- 0
    index                 <- block_number
    stop_success          <- 0
    stop_futility         <- 0
    for(i in 1:block_number){

      if(is.null(data_total)){
        yt <- 0
        Nt <- 0
        yc <- 0
        Nc <- 0

        est_interim <- bdpbinomial(y_t         = yt,
                                   N_t         = Nt,
                                   y_c         = yc,
                                   N_c         = Nc,
                                   a0          = a0,
                                   b0          = b0,
                                   number_mcmc = number_mcmc)

        if(alternative == "greater"){
          rr <- mean(est_interim$posterior_treatment$posterior -
                       est_interim$posterior_control$posterior > 0)
        }
        else{
          rr <- mean(est_interim$posterior_treatment$posterior -
                       est_interim$posterior_control$posterior < 0)
        }

      }

      data <- data.frame(
        treatment = sample(0:1, replace = T, group[i], prob = c(1 - rr, rr)),
        outcome   = rep(NA, group[i]))

      data$outcome <- rbinom(dim(data)[1], 1, prob = data$treatment * p_treatment +
                               (1 - data$treatment) * p_control +
                               drift * sum(group[1:i]) / N_total)

      data_total <- rbind(data_total, data)

      yt <- sum(data_total$outcome[data_total$treatment == 1])
      Nt <- length(data_total$outcome[data_total$treatment == 1])
      yc <- sum(data_total$outcome[data_total$treatment == 0])
      Nc <- length(data_total$outcome[data_total$treatment == 0])

      est_interim <- bdpbinomial(y_t         = yt,
                                 N_t         = Nt,
                                 y_c         = yc,
                                 N_c         = Nc,
                                 a0          = a0,
                                 b0          = b0,
                                 number_mcmc = number_mcmc)


      if(alternative == "greater"){
        rr <- mean(est_interim$posterior_treatment$posterior -
                     est_interim$posterior_control$posterior > 0)
      }
      else{
        rr <- mean(est_interim$posterior_treatment$posterior -
                     est_interim$posterior_control$posterior < 0)
      }


      if(rr > early_success_prob){
        index        <- i
        stop_success <- 1
        break
      }

      if(rr < futility_prob){
        index         <- i
        stop_futility <- 1
        break
      }

    }

    data_total <- data_total %>%
      mutate(time = factor(rep(1:index, group[1:index])))

    if(N_total / block_number > 2 | block_number < 2){
      fit0 <- bayesglm(formula = outcome ~ as.factor(treatment) + as.factor(time),
                       family  = binomial(link="logit"),
                       data    = data_total)
      post_trt <- coef(sim(fit0, n.sims = number_mcmc))[, 2]

      diff_est <- mean(post_trt)

      if(alternative == "greater"){
        prob_ha <- mean(post_trt > 0)
      }
      else{
        prob_ha <- mean(post_trt < 0)
      }
    }
    else{
      yt <- sum(data_total$outcome[data_total$treatment == 1])
      Nt <- length(data_total$outcome[data_total$treatment == 1])
      yc <- sum(data_total$outcome[data_total$treatment == 0])
      Nc <- length(data_total$outcome[data_total$treatment == 0])

      est_final <- bdpbinomial(y_t         = yt,
                               N_t         = Nt,
                               y_c         = yc,
                               N_c         = Nc,
                               a0          = a0,
                               b0          = b0,
                               number_mcmc = number_mcmc)
      diff_est <- mean(est_final$posterior_treatment$posterior -
                       est_final$posterior_control$posterior)

      if(alternative == "greater"){
        prob_ha <- mean(est_final$posterior_treatment$posterior -
                        est_final$posterior_control$posterior > 0)
      }
      else{
        prob_ha <- mean(est_final$posterior_treatment$posterior -
                        est_final$posterior_control$posterior < 0)
      }
    }

    N_control          <- c(N_control, sum(data_total$treatment == 0))
    N_treatment        <- c(N_treatment, sum(data_total$treatment == 1))
    sample_size        <- c(sample_size, dim(data_total)[1])
    prop_diff_estimate <- c(prop_diff_estimate, diff_est)

    if(prob_ha > (prob_accept_ha)){
      power <- power + 1
    }


  }

  output <- list(
    power                 = power / simulation,
    prop_diff_estimate    = prop_diff_estimate,
    N_enrolled            = sample_size,
    N_control             = N_control,
    N_treatment           = N_treatment
  )

  return(output)

}

## quiets concerns of R CMD check re: the .'s that appear in pipelines
if(getRversion() >= "2.15.1")  utils::globalVariables(c("treatment", "outcome"))


