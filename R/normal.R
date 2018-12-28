#' @title Block Design for Response-Adaptive Randomization
#'
#' @description Simulation for normal mean for block design for
#'    response-adaptive randomization with time as a confounding
#' @param mu_control scalar. Mean outcome in the control arm.
#' @param mu_treatment scalar. Mean outcome in the treatment arm.
#' @param sd_control scalar. Standard deviation of outcome in the control arm.
#' @param sd_treatment scalar. Standard deviation of outcome in the treatment
#' @inheritParams binomialRAR
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
#' @importFrom stats rnorm anova lm
#' @importFrom ldbounds bounds
#' @importFrom dplyr mutate
#'
#' @export normalRAR
#'
#' @examples
#' normalRAR(mu_control = 5.7, mu_treatment = 0.65,
#'           sd_control = 1.1, sd_treatment = 1.2,
#'           N_total = 200, block_number = 2, simulation = 100)
#'

normalRAR <- function(
  mu_control,
  mu_treatment,
  sd_control,
  sd_treatment,
  N_total,
  block_number    = 4,
  simulation      = 10000,
  zvalue          = NULL,
  rand_ratio      = NULL,
  conf_int        = 0.95,
  alternative     = "less",
  replace         = FALSE
){

  if(is.null(zvalue)){
    zvalue <- c(1, 1.5, 2, 2.5)
  }

  if(is.null(rand_ratio)){
    rand_ratio <- c(1, 1.5, 2, 2.5)
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




  output <- list(
    power        = power / simulation,
    N_enrolled   = sample_size,
    N_control    = N_control,
    N_treatment  = N_treatment
  )

  return(output)

}
