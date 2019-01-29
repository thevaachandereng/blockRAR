#' @title Block Design for Response-Adaptive Randomization for Binomial Data
#'
#' @description Simulation for binomial counts for block design for
#'    response-adaptive randomization with time as a confounding
#' @param p_treatment scalar. Proportion of events under the treatment arm.
#' @param p_control scalar. Proportion of events under the control arm.
#' @param N_total scalar. Total sample size.
#' @param block_number scalar. Number of blocks or time levels. The default is set to 4.
#'   If \code{block_number} is set to 1. This is a traditional RCT design.
#' @param drift scalar. The increase or decrease in proportion of event over time.
#'   In this case, the proportion of failure changes in each block by the number of
#'   patient accured over the total sample size. The full drift effect is seen in the
#'   final block.
#' @param simulation scalar. Number of simulation to be ran. The default is set to 10000.
#' @param zvalue vector. The z-value cutoff for corrected chi-square test statistics.
#' @param rand_ratio vector. The randomization ratio is set based on the corrected
#'   chi-square test statistic. The length of rand_ratio should be the same length of
#'   zvalue.
#' @param conf_int scalar. Confidence level of the interval.
#' @param alternative character. A string specifying the alternative hypothesis,
#'    must be one of "less" (default) or "greater".
#' @param correct logical. a logical indicating whether to apply continuity correction
#'    when computing the test statistic: one half is subtracted from all |O - E|
#'    differences; however, the correction will not be bigger than the differences themselves.
#' @param replace logical. should sampling be with replacement? If replace is set to
#'    FALSE (default), the 0 for control, 1 for treatment is replicated to the closest
#'    integer and this vector is sampled with no replacement. If replace is set to TRUE,
#'    the sampling is done based on randomization ratio provided with replacement.
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
#' @importFrom stats rbinom mantelhaen.test chisq.test
#' @importFrom ldbounds bounds
#' @importFrom dplyr mutate group_by summarize
#' @importFrom tibble as.tibble
#'
#' @export binomialRAR
#'
#' @examples
#' binomialRAR(p_control = 0.7, p_treatment = 0.65, N_total = 200,
#'             block_number = 2, simulation = 100)
#' binomialRAR(p_control = 0.5, p_treatment = 0.40, N_total = 200,
#'             block_number = 2, simulation = 100, drift = -0.15)

binomialRAR <- function(
  p_control,
  p_treatment,
  N_total,
  block_number    = 4,
  drift           = 0,
  simulation      = 10000,
  zvalue          = c(1, 1.5, 2),
  rand_ratio      = c(1, 1.5, 2, 2.5),
  conf_int        = 0.95,
  alternative     = "less",
  correct         = FALSE,
  replace         = TRUE
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

  if((conf_int <= 0 | conf_int >= 1)){
    stop("The confidence interval needs to between 0 and 1!")
  }

  if((alternative != "less" & alternative != "greater")){
    stop("The alternative can only be less or greater!")
  }

  if((replace != "TRUE" & replace != "FALSE")){
    stop("The replacement for sampling is either TRUE or FALSE!")
  }

  if(any(zvalue <= 0)){
    stop("The zvalue can only be greater than 0!")
  }

  if(any(rand_ratio < 1) | (length(rand_ratio) - 1) != length(zvalue)){
    stop("The randomization ratio needs to be greater than or equal to 1 and the length of
         randomization ratio needs to be greater than the length of zvalue by one!")
  }

  if(drift + p_control >= 1 | drift + p_control <= 0 |
     drift + p_treatment >= 1 | drift + p_treatment <= 0){
    stop("The drift value is too high causing the proportion of event to exceed 1
         in either the control or treatment group, pick a lower value for drift!")
  }

  if(replace == FALSE & N_total / block_number <= 2){
    replace <- TRUE
  }

  group <- rep(floor(N_total / block_number), block_number)
  if((N_total - sum(group)) > 0){
    index <- sample(1:block_number, N_total - sum(group))
    group[index] <- group[index] + 1
  }


  time <- seq(1 / block_number, 1, 1 / block_number)
  bounds <- bounds(time, iuse = c(1, 1), alpha = c(1 - conf_int, 1 - conf_int))$upper.bounds

  power <- 0

  N_control   <- NULL
  N_treatment <- NULL
  sample_size <- NULL
  p_control_estimate    <- NULL
  p_treatment_estimate  <- NULL

  for(k in 1:simulation){
    data_total            <- NULL
    test_stat             <- 0
    index                 <- block_number
    for(i in 1:block_number){
      if(length(rand_ratio) == 1){
        rr <- rand_ratio
      }
      else{
        zval <- c(-0.1, zvalue)
        rr <- rand_ratio[max(which(test_stat >= zval))]
      }

      if(!is.null(data_total) & length(levels(factor(data_total$treatment))) == 2){
        data_summary <- data_total %>%
          group_by(treatment) %>%
          summarize(prop = mean(as.numeric(outcome) - 1))
      }
      else{
        data_summary <- as.tibble(data.frame(treatment = c(0, 1), prop = c(0, 0)))
      }

      data <- data.frame(
        treatment =
          if(replace == TRUE){
            if((alternative == "less" &
               (data_summary$prop[1] >= data_summary$prop[2])) |
               (alternative == "greater" &
                (data_summary$prop[2] > data_summary$prop[1]))){
              sample(0:1, replace = T, group[i], prob = c(1, rr))
            }
            else{
              sample(0:1, replace = T, group[i], prob = c(rr, 1))
            }
          }
        else{
          if((alternative == "less" &
              (data_summary$prop[1] >= data_summary$prop[2])) |
             (alternative == "greater" &
              (data_summary$prop[2] > data_summary$prop[1]))){
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
                               (1 - data$treatment) * p_control +
                               drift * sum(group[1:i]) / N_total)

      data_total <- rbind(data_total, data)

      data_total$treatment <- as.factor(data_total$treatment)
      data_total$outcome <- as.factor(data_total$outcome)

      if(is.na(match("1", levels(data_total$outcome)))){
        data_total$outcome <- factor(data_total$outcome,
                                     levels=c(levels(data_total$outcome), "1"))
      }

      if(is.na(match("0", levels(data_total$outcome)))){
        data_total$outcome <- factor(data_total$outcome,
                                     levels=c(levels(data_total$outcome), "0"))
      }

      if(is.na(match("1", levels(data_total$treatment)))){
        data_total$treatment <- factor(data_total$treatment,
                                     levels=c(levels(data_total$treatment), "1"))
      }

      if(is.na(match("0", levels(data_total$treatment)))){
        data_total$treatment <- factor(data_total$treatment,
                                     levels=c(levels(data_total$treatment), "0"))
      }

      summ_data <-  data_total %>%
        group_by(treatment) %>%
        summarize(prop = mean(as.numeric(outcome) - 1))

      if(all(data_total$outcome == 1) | all(data_total$outcome == 0) |
         all(data_total$treatment == 1) | all(data_total$treatment == 0)){
        test_stat <- 0
      }
      else if(((summ_data$prop[2] - summ_data$prop[1] > 0) & alternative == "less") |
              ((summ_data$prop[1] - summ_data$prop[2] > 0) & alternative == "greater")){
        test_stat <- 0

      }
      else{
        if(i == 1 | N_total / block_number <= 2 |
           all(data_total[data_total$outcome == 1, 1] == 1) |
           all(data_total[data_total$outcome == 0, 1] == 1) |
           all(data_total[data_total$outcome == 1, 1] == 0) |
           all(data_total[data_total$outcome == 0, 1] == 0)){
          test_stat <- sqrt(as.numeric(chisq.test(data_total$treatment,
                                                  data_total$outcome,
                                                  correct = correct)$statistic))
        }
        else{
          temp_data <- data_total %>%
                          mutate(time = factor(rep(1:i, group[1:i])))
          test_stat <- sqrt(as.numeric(mantelhaen.test(table(temp_data),
                                       alternative = alternative,
                                       correct = correct)$statistic))

        }
      }

      if(test_stat > bounds[i]){
        index <- i
        break
      }

    }

    data_total <- data_total %>%
      mutate(time = factor(rep(1:index, group[1:index])))

    summary_data <-  data_total %>%
      group_by(treatment) %>%
      summarize(prop = mean(as.numeric(outcome) - 1))

    if(all(data_total$time == 1) | N_total / block_number < 2){
      if(((summary_data$prop[1] - summary_data$prop[2] > 0) & alternative == "less") |
         ((summary_data$prop[2] - summary_data$prop[1] > 0) & alternative == "greater")){
        p.val <- chisq.test(data_total$treatment, data_total$outcome,
                            correct = correct)$p.value
      }
      else{
        p.val <- 1
      }
    }
    else{
      p.val <- mantelhaen.test(table(data_total), alternative = alternative,
                               correct = correct)$p.val
    }

    N_control <- c(N_control, sum(data_total$treatment == 0))
    N_treatment <- c(N_treatment, sum(data_total$treatment == 1))
    sample_size <- c(sample_size, dim(data_total)[1])
    p_control_estimate <- c(p_control_estimate, summary_data$prop[1])
    p_treatment_estimate <- c(p_treatment_estimate, summary_data$prop[2])

    if(p.val < (1 - conf_int)){
      power <- power + 1
    }


  }

  output <- list(
    power                 = power / simulation,
    p_control_estimate    = p_control_estimate ,
    p_treatment_estimate  = p_treatment_estimate,
    N_enrolled            = sample_size,
    N_control             = N_control,
    N_treatment           = N_treatment
    )

  return(output)

}

## quiets concerns of R CMD check re: the .'s that appear in pipelines
if(getRversion() >= "2.15.1")  utils::globalVariables(c("treatment", "outcome"))




