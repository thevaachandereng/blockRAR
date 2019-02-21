#' @title Stratified posterior mean for Binomial Data
#'
#' @description Computing the posterior mean of treatment difference for
#'   stratified data with Bayesian method. The stratification is done over time.
#'
#' @inheritParams prop_strata
#' @inheritParams binomialbayes
#'
#' @return the weighted mean of posterior mean of difference (treatment - control).
#'
#' @importFrom bayesDP bdpbinomial
#' @export prop_stratabayes
#'
#' @examples
#' set.seed(20999)
#' prop_stratabayes(c(0, 1, 1, 1, 0, 1, 0, 0, 1, 1, 0, 1, 0, 1, 0),
#'                  c(0, 0, 1, 1, 0, 1, 1, 1, 0, 1, 1, 0, 1, 0, 1),
#'                  as.factor(rep(1:3, each = 5)))
prop_stratabayes <- function(treatment,
                             outcome,
                             block,
                             a0          = 0.5,
                             b0          = 0.5,
                             number_mcmc = 10000){
  # making sure the outcome is numeric
  outcome <- as.numeric(as.character(outcome))

  # stop if the all the variables length are not the same
  stopifnot(length(treatment) == length(outcome),
            length(treatment) == length(block))


  # if there is one block or one patient per block, compute simple mean
  if(all(block == 1) | length(treatment) == length(unique(block))){
    yt <- sum(outcome[treatment == 1])
    Nt <- length(outcome[treatment == 1])
    yc <- sum(outcome[treatment == 0])
    Nc <- length(outcome[treatment == 0])

    est_final <- bdpbinomial(y_t         = yt,
                             N_t         = Nt,
                             y_c         = yc,
                             N_c         = Nc,
                             a0          = a0,
                             b0          = b0,
                             number_mcmc = number_mcmc)
    diff_est <- mean(est_final$posterior_treatment$posterior -
                     est_final$posterior_control$posterior)

  }
  # else compute the stratified mean
  else{
    # assigning prop_diff
    prop_diff <- 0
    weight    <- 0
    # go through every block and compute the weighted prop difference
    for(i in levels(block)){
      # selecting the
      trt_grp <- outcome[treatment == 1 & block == i]
      ctr_grp <- outcome[treatment == 0 & block == i]

      # only if there is at least one treatment group and one control group, perform the analysis
      if(length(trt_grp) > 0 & length(ctr_grp) > 0){

        yt <- sum(trt_grp)
        Nt <- length(trt_grp)
        yc <- sum(ctr_grp)
        Nc <- length(ctr_grp)

        est_final <- bdpbinomial(y_t         = yt,
                                 N_t         = Nt,
                                 y_c         = yc,
                                 N_c         = Nc,
                                 a0          = a0 / length(levels(block)) ,
                                 b0          = b0 / length(levels(block)),
                                 number_mcmc = number_mcmc)


        # computing the weighted prop difference and the weight
        prop_diff <- prop_diff + 1 / (1 / length(ctr_grp) + 1 / length(trt_grp)) *
                    (mean(est_final$posterior_treatment$posterior -
                          est_final$posterior_control$posterior))
        weight    <- weight + 1 / (1 / length(ctr_grp) + 1 / length(trt_grp))
      }
    }
    # dividing it by the overall weight
    prop_diff <- prop_diff / weight
  }
  # return the stratified prop difference
  return(prop_diff)
}
