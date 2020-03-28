## ----setup, echo=FALSE, results="hide"-----------------------------------
knitr::opts_chunk$set(
  echo = FALSE,
  collapse = TRUE,
  comment = "#>"
)

## ---- echo = FALSE, fig=TRUE, fig.width = 6, fig.height = 6--------------
library(ggplot2)
suppressWarnings(RNGversion("3.5.0"))
set.seed(20999)


dat <- data.frame(
prob_response <- 0.5 + cumsum(c(0, runif(49, 0, 0.01))) + c(0, rnorm(49, 0, 0.01)),
time <- 0:49
)
dat$rand.ratio <- dat$prob_response


p <- ggplot(dat, aes(x = time))
p <- p + geom_line(aes(y = prob_response), color = "blue")
p <- p + scale_y_continuous(sec.axis = sec_axis(~., name = "Randomization Fraction for Treatment"))
p <- p + scale_colour_manual(values = c("blue"))
p <- p + labs(y = "Probability of Response (Clinical)",
              x = "Time (months)")
p <- p + theme(axis.text=element_text(size=12),
               axis.title=element_text(size=16,face="bold"))
p

## ---- eval = FALSE, echo=TRUE--------------------------------------------
#  install.packages("blockRAR")

## ---- eval = FALSE, echo=TRUE--------------------------------------------
#  devtools::install_github("thevaachandereng/blockRAR@vx.xx.x")
#  # or
#  devtools::install_version("blockRAR", version = "x.x.x", repos = "http://cran.us.r-project.org")

## ---- cache=FALSE, warning=FALSE, comment=FALSE, eval = FALSE, echo=TRUE, results="hide"----
#  devtools::install_github("thevaachandereng/blockRAR")

## ----lib, results="asis", eval=TRUE, echo=TRUE---------------------------
library(blockRAR)

## ---- warning=FALSE------------------------------------------------------
binomialfreq(p_control    = 0.25, 
             p_treatment  = 0.45, 
             N_total      = 200, 
             block_number = 5,  
             drift        = 0, 
             simulation   = 10, 
             conf_int     = 0.95,
             alternative  = "greater",
             early_stop   = FALSE)

## ------------------------------------------------------------------------
binomialbayes(p_control          = 0.35, 
              p_treatment        = 0.35, 
              N_total            = 150, 
              block_number       = 2,
              drift              = 0.10,
              simulation         = 10,
              a0                 = 0.5,
              b0                 = 0.5, 
              number_mcmc        = 10000, 
              prob_accept_ha     = 0.95,
              early_success_prob = 1,
              futility_prob      = 0,
              alternative        = "greater")

## ------------------------------------------------------------------------
sessionInfo()

