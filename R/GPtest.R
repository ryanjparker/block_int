# code to analyze sequential design setting
library(fields)
library(lhs)
library(MASS)
library(multicore)
library(sensitivity)
library(lattice)

source("R/cov.R")
source("R/estimate.R")

# load test data
test.dat <- read.csv("data/GPtest.csv")

# fit model to data
data <- list(
		Yobs=test.dat$y,
		X=cbind(test.dat$S1,test.dat$S2,test.dat$S3)
	)

n <- length(test.dat$y)
p <- 3
starts <- log( c(4,4,4) )

fit.full <- hd.estimate(data, 1, rep(1,n), matrix(1,nrow=1,ncol=1), p, starts, rep(FALSE,p), ce_cov, ce_partial, FALSE)
