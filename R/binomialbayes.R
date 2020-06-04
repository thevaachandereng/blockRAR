#' @title Block Design for Response-Adaptive Randomization for Binomial Data
#'
#' @description Simulation for binomial counts for block design for
#'    response-adaptive randomization with time as a confounding
#'
#' @inheritParams binomialfreq
#' @param a0 scalar. Prior value for the beta rate \code{Beta(a0, b0)}.
#'  Default is 0.5.
#' @param b0 scalar. Prior value for the beta rate \code{Beta(a0, b0)}.
#'  Default is 0.5.
#' @param p scalar. Power for randomization ratio.
#' @param number_mcmc scalar. Number of Monte Carlo Markov Chain draws in
#'   sampling posterior.
#' @param prob_accept_ha scalar. Probability of accepting
#'   alternative hypothesis.
#' @param early_success_prob scalar. Probability of stopping early for success.
#' @param futility_prob scalar. Probability of stopping early for futility.
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
#'   \item{\code{randomization_ratio}}{
#'     matrix. The randomization ratio allocated for each block.}
#' }
#'
#' @importFrom stats rbinom binomial coef quasi rbeta
#' @importFrom arm bayesglm sim
#' @importFrom dplyr mutate group_by summarize
#' @importFrom tibble as.tibble
#'
#' @export binomialbayes
#'
#' @examples
#' binomialbayes(p_control = 0.20, p_treatment = 0.30, N_total = 100, simulation = 3)
#' binomialbayes(p_control = 0.50, p_treatment = 0.30, N_total = 100, simulation = 3)
#'

binomialbayes <- function(
  p_control,
  p_treatment,
  N_total,
  block_number              = 4,
  drift                     = 0,
  simulation                = 10000,
  a0                        = 0.5,
  b0                        = 0.5,
  p                         = "n/2N",
  number_mcmc               = 10000,
  prob_accept_ha            = 0.95,
  early_success_prob        = 0.99,
  futility_prob             = 0.01,
  alternative               = "greater",
  size_equal_randomization  = 20,
  min_patient_earlystop     = 20,
  max_prob                  = 0.8
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

  # make sure the block size is a positive integer!
  if((block_number <= 0 | block_number %% 1 != 0)){
    stop("The number of blocks needs to be a positve integer!")
  }

  # make sure the simulation is a positive integer
  if((simulation <= 0 | simulation %% 1 != 0)){
    stop("The number of simulation needs to be a positve integer!")
  }

  # confidence interval between 0 and 1
  if((prob_accept_ha <= 0 | prob_accept_ha >= 1)){
    stop("The confidence interval needs to between 0 and 1!")
  }

  # the alternative is either less or greater than
  if((alternative != "less" & alternative != "greater")){
    stop("The alternative can only be less or greater!")
  }

  # making sure the drift didnt make the prop of control/treatment > 1 / < 0
  if(drift + p_control >= 1 | drift + p_control <= 0 |
     drift + p_treatment >= 1 | drift + p_treatment <= 0){
    stop("The drift value is too high causing the proportion of event to exceed 1
         in either the control or treatment group, pick a lower value for drift!")
  }

  # making sure the drift didnt make the prop of control/treatment > 1 / < 0
  if(N_total < block_number){
    stop("The number of blocks can't exceed the number of patients!")
  }

  # computing the group size if its symmetric or not
  group <- rep(floor(N_total / block_number), block_number)

  if((N_total - sum(group)) > 0){
    index        <- sample(1:block_number, N_total - sum(group))
    group[index] <- group[index] + 1
  }

  # storing values all important variables as well as the drift
  power              <- 0
  N_control          <- NULL
  N_treatment        <- NULL
  sample_size        <- NULL
  prop_diff_estimate <- NULL
  early_success      <- NULL
  early_futility     <- NULL
  drift_p            <- seq(drift / N_total, drift,  length.out = N_total)
  randomization      <- array(NA, c(simulation, N_total))

  # going through all the simulations
  for(k in 1:simulation){
    #storing each simulation values
    data_total            <- data.frame()
    test_stat             <- 0
    time                  <- rep(1:block_number, group[1:block_number])
    stop_success          <- 0
    stop_futility         <- 0

    for(i in 1:N_total){

      if(any((i - 1) == cumsum(group)) | i == 1){

        # if data_total is null, set all the outcome to 0
        if(nrow(data_total) < size_equal_randomization){
          yt <- 0
          Nt <- 0
          yc <- 0
          Nc <- 0
        }

        # estimating the interim value estimates
        est_interim <- bdpbinomial(y_t         = yt,
                                   N_t         = Nt,
                                   y_c         = yc,
                                   N_c         = Nc,
                                   a0          = a0,
                                   b0          = b0,
                                   number_mcmc = number_mcmc)

        if(p == "n/2N"){
          pi <- nrow(data_total) / (2 * N_total)
        }
        else if(p >= 0 & p <= 1){
          pi <- p
        }
        else{
          pi <- 0.5
        }
        # altering the randomization ratio based on Thall and Wathen's paper
        if(alternative == "greater"){
          diff <- est_interim$posterior_treatment$posterior -
            est_interim$posterior_control$posterior

          rr <- mean(diff > 0)^pi / (mean(diff > 0)^pi + mean(diff < 0)^pi)
        }
        else{
          diff <- est_interim$posterior_treatment$posterior -
            est_interim$posterior_control$posterior
          rr <- mean(diff < 0)^pi / (mean(diff > 0)^pi + mean(diff < 0)^pi)
        }

        # maximum probability assigning to the treatment group is 0.8
        if(rr > max_prob){
          rr <- max_prob
        }

        # maximum probability assigning to the control group is 0.8
        else if(rr < (1 - max_prob)){
          rr <- 1 - max_prob
        }
      }

      randomization[k, i] <- rr

      # creating the dataset for each patient
      data <- data.frame(
        treatment = sample(0:1, size = 1, replace = T, prob = c(1 - rr, rr)),
        outcome   = rep(NA, 1))

      # adding the outcome with time trends (linear time trend)
      data$outcome <- rbinom(nrow(data), 1, prob = data$treatment * p_treatment +
                               (1 - data$treatment) * p_control +
                               drift_p[i])

      # joining the dataset using rbind
      data_total <- rbind(data_total, data)

      data_interim <- data_total %>%
        mutate(time = time[1:i])


      if(any(N_total / block_number < 2 | block_number == 1 | all(data_interim$time == 1))){
        # accumulating information on the events
        yt <- sum(data_interim$outcome[data_interim$treatment == 1])
        Nt <- length(data_interim$outcome[data_interim$treatment == 1])
        yc <- sum(data_interim$outcome[data_interim$treatment == 0])
        Nc <- length(data_interim$outcome[data_interim$treatment == 0])

        # performing interim analysis to allow for stopping early
        est_interim <- bdpbinomial(y_t         = yt,
                                   N_t         = Nt,
                                   y_c         = yc,
                                   N_c         = Nc,
                                   a0          = a0,
                                   b0          = b0,
                                   number_mcmc = number_mcmc)

        # calculating the mean posterior difference for treatment estimate
        if(alternative == "greater"){
          int_analysis <- mean(est_interim$posterior_treatment$posterior -
                                 est_interim$posterior_control$posterior > 0)
        }
        else{
          int_analysis <- mean(est_interim$posterior_treatment$posterior -
                                 est_interim$posterior_control$posterior < 0)
        }
      }
      else{
        fit_int <- tryCatch(expr = bayesglm(formula = outcome ~ as.factor(treatment) + as.factor(time),
                                            data    = data_interim,
                                            family  = quasi(link = "identity", variance = "mu(1-mu)"),
                                            start   = rep(0.1, 1 + length(levels(data_interim$time)))),
                            error   = function(data = data_interim){
                              rbeta(number_mcmc,
                                    sum(data$outcome[data$time == data$time[1] & data$treatment == 1]) + a0,
                                    length(data$outcome[data$time == data$time[1] & data$treatment == 1]) + b0) -
                                rbeta(number_mcmc,
                                      sum(data$outcome[data$time == data$time[1] & data$treatment == 0]) + a0,
                                      length(data$outcome[data$time == data$time[1] & data$treatment == 0])) + b0})
        if(length(fit_int) == number_mcmc){
          int_trt <- fit_int
        }
        else{
          int_trt <- coef(sim(fit_int, n.sims = number_mcmc))[, 2]
        }

        if(alternative == "greater"){
          int_analysis <- mean(int_trt > 0)
        }
        else{
          int_analysis <- mean(int_trt < 0)
        }
      }

      # check for early stopping for success
      if(int_analysis > early_success_prob & nrow(data_total) >= min_patient_earlystop){
        time         <- time[1:i]
        stop_success <- 1
        if(i < N_total){
          randomization[k, (i+1):N_total] <- 0
        }
        break
      }

      # check for early stopping for futility
      if(int_analysis < futility_prob & nrow(data_total) >= min_patient_earlystop){
        time          <- time[1:i]
        stop_futility <- 1
        if(i < N_total){
          randomization[k, (i+1):N_total] <- 0
        }
        break
      }

    }

    # mutate time factor for time column
    data_total <- data_total %>%
      mutate(time = as.factor(time))

    # setting data final same as data_total
    data_final <- data_total

    # if block size is less than 2 or number of time block is less than 2
    # or the number o patient in each block is less than 2, do not do stratified
    # analysis
    if(N_total / block_number < 2 | all(data_total$time == 1) | block_number < 2){
      # accumulating details for analysis
      yt <- sum(data_total$outcome[data_total$treatment == 1])
      Nt <- length(data_total$outcome[data_total$treatment == 1])
      yc <- sum(data_total$outcome[data_total$treatment == 0])
      Nc <- length(data_total$outcome[data_total$treatment == 0])

      ## non-stratified analysis
      est_final <- bdpbinomial(y_t         = yt,
                               N_t         = Nt,
                               y_c         = yc,
                               N_c         = Nc,
                               a0          = a0,
                               b0          = b0,
                               number_mcmc = number_mcmc)

      # estimating posterior treatment difference
      diff_est <- mean(est_final$posterior_treatment$posterior -
                         est_final$posterior_control$posterior)

      # for different condition, computing the probability of accepting
      # alternative hypothesis
      if(alternative == "greater"){
        prob_ha <- mean(est_final$posterior_treatment$posterior -
                          est_final$posterior_control$posterior > 0)
      }
      else{
        prob_ha <- mean(est_final$posterior_treatment$posterior -
                          est_final$posterior_control$posterior < 0)
      }
    }

    # stratified analysis for factor time
    else{
      # perform bayes generalized linear models with quasi family and identity link
      fit0 <- tryCatch(expr = bayesglm(formula = outcome ~ as.factor(treatment) + as.factor(time),
                                       data    = data_total,
                                       family  = quasi(link = "identity", variance = "mu(1-mu)"),
                                       start   = rep(0.1, 1 + length(levels(data_total$time)))),
                       error   = function(data = data_total){
                         rbeta(number_mcmc,
                               sum(data$outcome[data$time == data$time[1] & data$treatment == 1]) + a0,
                               length(data$outcome[data$time == data$time[1] & data$treatment == 1]) + b0) -
                           rbeta(number_mcmc,
                                 sum(data$outcome[data$time == data$time[1] & data$treatment == 0]) + a0,
                                 length(data$outcome[data$time == data$time[1] & data$treatment == 0])) + b0})




      # estimating the treatment effect
      if(length(fit0) == number_mcmc){
        post_trt <- fit0
      }
      else{
        post_trt <- coef(sim(fit0, n.sims = number_mcmc))[, 2]
      }

      # computing mean posterior treatment effect
      diff_est <- mean(post_trt)

      # computing the probability of accepting alternative hypothesis
      if(alternative == "greater"){
        prob_ha <- mean(post_trt > 0)
      }
      else{
        prob_ha <- mean(post_trt < 0)
      }

    }

    # storing all the control, treatment information for each trial simulation
    N_control          <- c(N_control, sum(data_final$treatment == 0))
    N_treatment        <- c(N_treatment, sum(data_final$treatment == 1))
    sample_size        <- c(sample_size, nrow(data_final))
    prop_diff_estimate <- c(prop_diff_estimate, diff_est)
    early_success      <- c(early_success, stop_success)
    early_futility     <- c(early_futility, stop_futility)

    # computing the power for all the trial simulation
    if(prob_ha > (prob_accept_ha)){
      power <- power + 1
    }


  }

  # storing the information for outputs
  output <- list(
    power                 = power / simulation,
    prop_diff_estimate    = prop_diff_estimate,
    N_enrolled            = sample_size,
    N_control             = N_control,
    N_treatment           = N_treatment,
    early_success         = early_success,
    early_futilty         = early_futility,
    prob_trt              = randomization
  )

  # return output
  return(output)

}

## quiets concerns of R CMD check re: the .'s that appear in pipelines
if(getRversion() >= "2.15.1")  utils::globalVariables(c("treatment", "outcome"))


