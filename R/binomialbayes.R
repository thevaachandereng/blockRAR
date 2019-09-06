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
#' }
#'
#' @importFrom stats rbinom binomial coef quasi rbeta
#' @importFrom arm bayesglm sim
#' @importFrom bayesDP bdpbinomial
#' @importFrom dplyr mutate group_by summarize
#' @importFrom tibble as.tibble
#'
#' @export binomialbayes
#'
#' @examples
#' binomialbayes(p_control = 0.20, p_treatment = 0.30, N_total = 100, simulation = 10)
#' binomialbayes(p_control = 0.50, p_treatment = 0.30, N_total = 100, simulation = 5)
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
  p                  = 0.5,
  number_mcmc        = 10000,
  prob_accept_ha     = 0.95,
  early_success_prob = 0.99,
  futility_prob      = 0.01,
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

  # computing the group size if its symmetric or not
  group <- rep(floor(N_total / block_number), block_number)
  if((N_total - sum(group)) > 0){
    index <- sample(1:block_number, N_total - sum(group))
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

  # going through all the simulations
  for(k in 1:simulation){
    #storing each simulation values
    data_total            <- data.frame()
    test_stat             <- 0
    index                 <- block_number
    stop_success          <- 0
    stop_futility         <- 0

    for(i in 1:block_number){

      # if data_total is null, set all the outcome to 0
      if(dim(data_total)[1] == 0){
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

        # altering the randomization ratio based on Thall and Wathen's paper
        if(alternative == "greater"){
          diff <- est_interim$posterior_treatment$posterior -
                  est_interim$posterior_control$posterior

          rr <- mean(diff > 0)^p / (mean(diff > 0)^p + mean(diff < 0)^p)
        }
        else{
          diff <- est_interim$posterior_treatment$posterior -
                  est_interim$posterior_control$posterior
          rr <- mean(diff < 0)^p / (mean(diff > 0)^p + mean(diff < 0)^p)
        }

      # creating the dataset for each block
      data <- data.frame(
        treatment = sample(0:1, replace = T, group[i], prob = c(1 - rr, rr)),
        outcome   = rep(NA, group[i]))

      # adding the outcome with time trends (linear time trend)
      data$outcome <- rbinom(dim(data)[1], 1, prob = data$treatment * p_treatment +
                               (1 - data$treatment) * p_control +
                               drift_p[((1:dim(data)[1]) + dim(data_total)[1])])

      # joining the dataset using rbind
      data_total <- rbind(data_total, data)

      # accumulating information on the events
      yt <- sum(data_total$outcome[data_total$treatment == 1])
      Nt <- length(data_total$outcome[data_total$treatment == 1])
      yc <- sum(data_total$outcome[data_total$treatment == 0])
      Nc <- length(data_total$outcome[data_total$treatment == 0])

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
        rr <- mean(est_interim$posterior_treatment$posterior -
                     est_interim$posterior_control$posterior > 0)
      }
      else{
        rr <- mean(est_interim$posterior_treatment$posterior -
                     est_interim$posterior_control$posterior < 0)
      }

      # check for early stopping for success
      if(rr > early_success_prob){
        index        <- i
        stop_success <- 1
        break
      }

      # check for early stopping for futility
      if(rr < futility_prob){
        index         <- i
        stop_futility <- 1
        break
      }

    }

    # mutate time factor for time column
    data_total <- data_total %>%
      mutate(time = factor(rep(1:index, group[1:index])))

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

      for(i in levels(data_total$time)){

        # selecting the treatment and control group
        trt_grp <- data_total$outcome[data_total$treatment == 1 & data_total$time == i]
        ctr_grp <- data_total$outcome[data_total$treatment == 0 & data_total$time == i]

        # only if there is at least one person in the
        # treatment group and one control group,  perform the analysis
        if(length(trt_grp) == 0 | length(ctr_grp) == 0){
          data_total <- data_total[-which(data_total$time == i), ]
        }
      }

      # drop levels for all the timepoints dropped due to missing data
      data_total <- droplevels.data.frame(data_total)

      ## if there is more than 1 factor level, fit time as a factor level
      if(length(levels(data_total$time)) > 1){
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
      }

      # else fit just the treatment effect model
      else if (length(levels(data_total$time)) == 1){
        # perform bayes generalized linear models with quasi family and identity link
        fit0 <- bayesglm(formula = outcome ~ as.factor(treatment),
                         data    = data_total,
                         family  = quasi(link = "identity", variance = "mu(1-mu)"),
                         start   = rep(0.1, 1 + length(levels(data_total$time))))
        # estimating the treatment effect
        post_trt <- coef(sim(fit0, n.sims = number_mcmc))[, 2]
      }

      # if all the columns are dropped, fit equal to 0
      else{
        fit0 <- 0
        post_trt <- 0
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
    sample_size        <- c(sample_size, dim(data_final)[1])
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
    early_futilty         = early_futility
  )

  # return output
  return(output)

}

## quiets concerns of R CMD check re: the .'s that appear in pipelines
if(getRversion() >= "2.15.1")  utils::globalVariables(c("treatment", "outcome"))


