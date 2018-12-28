#' @title Block Design for Response-Adaptive Randomization
#'
#' @description Simulation for binomial counts for block design for
#'    response-adaptive randomization with time as a confounding
#' @param p_treatment scalar. Proportion of events under the treatment arm.
#' @param p_control scalar. Proportion of events under the control arm.
#' @param N_total scalar. Total sample size.
#' @param block_number scalar. Number of blocks or time levels. The default is set to 4.
#'   If \code{block_number} is set to 1. This is a traditional RCT design.
#' @param simulation scalar. Number of simulation to be ran. The default is set to 10000.
#' @param zvalue vector. The z-value cutoff for corrected chi-square test statistics.
#' @param rand_ratio vector. The randomization ratio is set based on the corrected
#'   chi-square test statistic. The length of rand_ratio should be the same length of
#'   zvalue.
#' @param conf.int scalar. Confidence level of the interval.
#' @param alternative character. A string specifying the alternative hypothesis,
#'    must be one of "less" (default) or "greater".
#' @param replace character. should sampling be with replacement? If replace is set to
#'    FALSE (default), the 0 for control, 1 for treatment is replicated to the closest
#'    integer and this vector is sampled with no replacement. If replace is set to TRUE,
#'    the sampling is done based on randomization ratio provided with replacement.
#'
#' @return a list with details on the simulation.
#' \describe{
#'   \item{\code{power}}{
#'     scalar. The power of the trial, ie. the proportion of success over the
#'     number of simulation ran.}
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
#' @importFrom stats rbinom mantelhaen.test chisq.test
#' @importFrom ldbounds bounds
#' @importFrom dplyr mutate
#'
#' @export binomialRAR
#'
#' @examples
#' binomialRAR(p_control = 0.7, p_treatment = 0.65, N_total = 200,
#'             block_number = 2, simulation = 100)

binomialRAR <- function(
  p_control,
  p_treatment,
  N_total,
  block_number    = 4,
  simulation      = 10000,
  zvalue          = NULL,
  rand_ratio      = NULL,
  conf.int        = 0.95,
  alternative     = "less",
  replace         = FALSE
){

  if(is.null(zvalue)){
    zvalue <- c(1, 1.5, 2, 2.5)
  }

  if(is.null(rand_ratio)){
    rand_ratio <- c(1, 1.5, 2, 2.5)
  }

  if((p_control <= 0 | p_control >= 1)){
    stop("The proportion of event for the control group needs to between 0 and 1!")
  }

  if((p_treatment <= 0 | p_treatment >= 1)){
    stop("The proportion of event for the treatment group needs to between 0 and 1!")
  }

  if((N_total < 0 | N_total %% 1 != 0)){
    stop("The sample size needs to be a positive integer!")
  }

  if((block_number < 0 | block_number %% 1 != 0)){
    stop("The number of blocks needs to be a positve integer!")
  }

  if((simulation < 0 | simulation %% 1 != 0)){
    stop("The number of simulation needs to be a positve integer!")
  }

  group <- rep(floor(N_total / block_number), block_number)
  if((N_total - sum(group)) > 0){
    index <- sample(1:block_number, N_total - sum(group))
    group[index] <- group[index] + 1
  }


  time <- seq(1 / block_number, 1, 1 / block_number)
  bounds <- bounds(time, iuse = c(1, 1), alpha = c(1 - conf.int, 1 - conf.int))$upper.bounds

  power <- 0

  N_control   <- NULL
  N_treatment <- NULL
  sample_size <- NULL

  for(k in 1:simulation){
    data_total <- NULL
    test_stat  <- 0
    index      <- block_number
    for(i in 1:block_number){

      rr <- rand_ratio[which.max(rand_ratio >= test_stat)]

      data <- data.frame(
        treatment =
          if(replace == TRUE){
            if(alternative == "less"){
              sample(0:1, replace = T, group[i], prob = c(1, rr))
            }
            else{
              sample(0:1, replace = T, group[i], prob = c(rr, 1))
            }
          }
        else{
          if(alternative == "less"){
            sum_ratio <- rr + 1
            sampling <- rep(c(0, 1), round(c(group[i] * 1 / sum_ratio - 0.0001,
                                           group[i] * rr / sum_ratio + 0.0001)))
            sample(sampling, length(sampling))
          }
          else{
            sum_ratio <- rr + 1
            sampling <- rep(c(0, 1), round(c(group[i] * rr / sum_ratio - 0.0001,
                                             group[i] * 1 / sum_ratio + 0.0001)))
            sample(sampling, length(sampling))
          }

        },
        outcome = rep(NA, group[i]))

      data$outcome <- rbinom(dim(data)[1], 1, prob = data$treatment * p_treatment +
                               (1 - data$treatment) * p_control)

      data_total <- rbind(data_total, data)

      data_total$treatment <- as.factor(data_total$treatment)
      data_total$outcome <- as.factor(data_total$outcome)

      if(is.na(match("1", levels(data_total$outcome)))){
        return(factor(data_total$outcome, levels=c(levels(data_total$outcome), "1")))
      }

      if(is.na(match("0", levels(data_total$outcome)))){
        return(factor(data_total$outcome, levels=c(levels(data_total$outcome), "0")))
      }

      test_stat <- sqrt(as.numeric(chisq.test(data_total$treatment,
                                              data_total$outcome)$statistic))

      if(test_stat > bounds[i]){
        index <- i
        break
      }

    }

    data_total <- data_total %>%
      mutate(time = factor(rep(1:index, group[1:index])))

    if(all(data_total$time == 1)){
      p.val <- chisq.test(data_total$treatment, data_total$outcome)$p.value
    }
    else{
      p.val <- mantelhaen.test(table(data_total), alternative = alternative)$p.val
    }

    N_control <- c(N_control, sum(data_total$treatment == 0))
    N_treatment <- c(N_treatment, sum(data_total$treatment == 1))
    sample_size <- c(sample_size, dim(data_total)[1])

    if(p.val < (1 - conf.int)){
      power <- power + 1
    }


  }

  output <- list(
    power        = power / simulation,
    N_enrolled   = sample_size,
    N_control    = N_control,
    N_treatment  = N_treatment
    )

  return(output)

}






