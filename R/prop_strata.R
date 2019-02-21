#' @title Stratified Proportion Estimate for Binomial Data
#'
#' @description Computing the proportion of treatment difference for
#'   stratified data. The stratification is done over time.
#' @param treatment vector. The vector with treatment assignment, 0 for control
#'  and 1 for treatment group.
#' @param outcome vector. The vector with outcome, 0 for failure
#'  and 1 for success. Must be the same length as treatment variable.
#' @param block vector. The vector with factor level of the block.
#'  Must be same lenhth as treatment variable.
#'
#' @return the weighted mean of proportion difference (treatment - control).
#'
#' @export prop_strata
#'
#' @examples
#' set.seed(20999)
#' prop_strata(c(0, 1, 1, 1, 0, 1, 0, 0, 1, 1, 0, 1, 0, 1, 0),
#'             c(0, 0, 1, 1, 0, 1, 1, 1, 0, 1, 1, 0, 1, 0, 1),
#'             as.factor(rep(1:3, each = 5)))


prop_strata <- function(treatment, outcome, block){
  # making sure the outcome is numeric
  outcome <- as.numeric(as.character(outcome))

  # stop if the all the variables length are not the same
  stopifnot(length(treatment) == length(outcome),
            length(treatment) == length(block))

  # if there is one block or one patient per block, compute simple mean
  if(all(block == 1) | length(treatment) == length(unique(block))){
    prop_diff <- mean(outcome[treatment == 1]) - mean(outcome[treatment == 0])
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
        # computing the weighted prop difference and the weight
        prop_diff <- prop_diff + 1 / (1 / length(ctr_grp) + 1 / length(trt_grp)) * (mean(trt_grp) - mean(ctr_grp))
        weight    <- weight + 1 / (1 / length(ctr_grp) + 1 / length(trt_grp))
      }
    }
    # dividing it by the overall weight
    prop_diff <- prop_diff / weight
  }
  # return the stratified prop difference
  return(prop_diff)
}
