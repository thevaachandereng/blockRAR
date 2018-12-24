## ----setup, echo=FALSE, results="hide"-----------------------------------
knitr::opts_chunk$set(
  echo = FALSE,
  collapse = TRUE,
  comment = "#>"
)
set.seed(43232)

## ---- echo = FALSE, fig=TRUE, fig.width = 6, fig.height = 6--------------
library(ggplot2)
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
#  
#  # or
#  
#  devtools::install_version("blockRAR", version = "x.x.x", repos = "http://cran.us.r-project.org")

## ---- cache=FALSE, warning=FALSE, comment=FALSE, eval = TRUE, echo=TRUE, results="hide"----
devtools::install_github("thevaachandereng/blockRAR")

## ----lib, results="asis", eval=TRUE, echo=TRUE---------------------------
library(blockRAR)

