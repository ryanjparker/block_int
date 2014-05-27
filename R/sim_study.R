# code to run simulation study
library(fields)
library(lhs)
library(MASS)
library(multicore)
library(sensitivity)

#source("R/create_blocks.R")
source("R/cov.R")
source("R/estimate.R")

# function to execte the simulation study based on given factors
"sim_exp" <- function(design, factors, which.exp) {

	#res <- mclapply(1:design$Nreps, function(i) {
	#res <- mclapply(1:3, function(i) {
	#res <- lapply(1:design$Nreps, function(i) {
	res <- lapply(1:1, function(i) { #design$Nreps, function(i) {
		seed <- 1983 + i + design$Nreps*(which.exp-1)
		set.seed(seed)  # set a seed for reproducibility

		# generate data
		data <- generate_data(design, factors)

		# get results for models...

			# .. oracle
			res.orac <- eval.orac(design, factors, data)

			# ... full
			res.full <- eval.full(design, factors, data)

			# 250 obs per block
			factors$Nblock_obs_ind <- 250
				res.indr <- eval.indr(design, factors, data)
				res.indc <- eval.indc(design, factors, data, res.indr$theta)
				rir_250 <- res.indr
				ric_250 <- res.indc

			# 100 obs per block
			factors$Nblock_obs_ind <- 100
				res.indr <- eval.indr(design, factors, data)
				res.indc <- eval.indc(design, factors, data, res.indr$theta)
				rir_100 <- res.indr
				ric_100 <- res.indc

			# 50 obs per block
			factors$Nblock_obs_ind <- 50
				res.indr <- eval.indr(design, factors, data)
				res.indc <- eval.indc(design, factors, data, res.indr$theta)
				rir_50 <- res.indr
				ric_50 <- res.indc

			# 25 obs per block
			factors$Nblock_obs_ind <- 25
				res.indr <- eval.indr(design, factors, data)
				res.indc <- eval.indc(design, factors, data, res.indr$theta)
				rir_25 <- res.indr
				ric_25 <- res.indc

			# ... dependent blocks/random assignment
			#res.depr <- eval.depr(design, factors, data)

			# ... dependent blocks/cluster assignment
			#res.depc <- eval.depc(design, factors, data, res.depr$theta)

		# return results
		r <- list(seed=seed, p=factors$p, tau=factors$tau,
			# oracle
			orac.status=res.orac$status, orac.time=res.orac$time, orac.rmse_p=res.orac$rmse_p,
			# full
			full.status=res.full$status, full.time=res.full$time, full.rmse_t=res.full$rmse_t, full.rmse_p=res.full$rmse_p,
			full.rmse_s_S=res.full$rmse_s_S, full.rmse_s_T=res.full$rmse_s_T,
			# ind random/cluster
			ir250.status=rir_250$status, ir250.time=rir_250$time, ir250.rmse_t=rir_250$rmse_t, ir250.rmse_full_p=rir_250$rmse_full_p,
			ir250.rmse_s_S=rir_250$rmse_s_T, ir250.rmse_s_T=rir_250$rmse_s_T,
			ic250.status=ric_250$status, ic250.time=ric_250$time, ic250.rmse_t=ric_250$rmse_t, ic250.rmse_full_p=ric_250$rmse_full_p,
			ic250.rmse_s_S=ric_250$rmse_s_T, ic250.rmse_s_T=ric_250$rmse_s_T,
			ir100.status=rir_100$status, ir100.time=rir_100$time, ir100.rmse_t=rir_100$rmse_t, ir100.rmse_full_p=rir_100$rmse_full_p,
			ir100.rmse_s_S=rir_100$rmse_s_T, ir100.rmse_s_T=rir_100$rmse_s_T,
			ic100.status=ric_100$status, ic100.time=ric_100$time, ic100.rmse_t=ric_100$rmse_t, ic100.rmse_full_p=ric_100$rmse_full_p,
			ic100.rmse_s_S=ric_100$rmse_s_T, ic100.rmse_s_T=ric_100$rmse_s_T,
			ir50.status=rir_50$status, ir50.time=rir_50$time, ir50.rmse_t=rir_50$rmse_t, ir50.rmse_full_p=rir_50$rmse_full_p,
			ir50.rmse_s_S=rir_50$rmse_s_T, ir50.rmse_s_T=rir_50$rmse_s_T,
			ic50.status=ric_50$status, ic50.time=ric_50$time, ic50.rmse_t=ric_50$rmse_t, ic50.rmse_full_p=ric_50$rmse_full_p,
			ic50.rmse_s_S=ric_50$rmse_s_T, ic50.rmse_s_T=ric_50$rmse_s_T,
			ir25.status=rir_25$status, ir25.time=rir_25$time, ir25.rmse_t=rir_25$rmse_t, ir25.rmse_full_p=rir_25$rmse_full_p,
			ir25.rmse_s_S=rir_25$rmse_s_T, ir25.rmse_s_T=rir_25$rmse_s_T,
			ic25.status=ric_25$status, ic25.time=ric_25$time, ic25.rmse_t=ric_25$rmse_t, ic25.rmse_full_p=ric_25$rmse_full_p,
			ic25.rmse_s_S=ric_25$rmse_s_T, ic25.rmse_s_T=ric_25$rmse_s_T,

			# sensitivity indices
			orac.Si_S=data$Si_S, orac.Si_T=data$Si_T,
			full.Si_S=res.full$Si_S, full.Si_T=res.full$Si_T,
			ir250.Si_S=rir_250$Si_S, ir250.Si_T=rir_250$Si_T,
			ic250.Si_S=ric_250$Si_S, ic250.Si_T=ric_250$Si_T,
			ir100.Si_S=rir_100$Si_S, ir100.Si_T=rir_100$Si_T,
			ic100.Si_S=ric_100$Si_S, ic100.Si_T=ric_100$Si_T,
			ir50.Si_S=rir_50$Si_S, ir50.Si_T=rir_50$Si_T,
			ic50.Si_S=ric_50$Si_S, ic50.Si_T=ric_50$Si_T,
			ir25.Si_S=rir_25$Si_S, ir25.Si_T=rir_25$Si_T,
			ic25.Si_S=ric_25$Si_S, ic25.Si_T=ric_25$Si_T
#			# ind/random
#			indr.status=res.indr$status, indr.time=res.indr$time, indr.rmse_t=res.indr$rmse_t, indr.rmse_s=res.indr$rmse_s,
#				indr.rmse_full_p=res.indr$rmse_full_p, #indr.rmse_block_p=res.indr$rmse_block_p, indr.rmse_local_p=res.indr$rmse_local_p,
#			# ind/cluster
#			indc.status=res.indc$status, indc.time=res.indc$time, indc.rmse_t=res.indc$rmse_t, indc.rmse_s=res.indc$rmse_s,
#				indc.rmse_full_p=res.indc$rmse_full_p#, indc.rmse_block_p=res.indc$rmse_block_p, indc.rmse_local_p=res.indc$rmse_local_p,
#			# dep/random
#			depr.status=res.depr$status, depr.time=res.depr$time, depr.rmse_t=res.depr$rmse_t, depr.rmse_s=res.depr$rmse_s,
#				depr.rmse_full_p=res.depr$rmse_full_p, depr.rmse_block_p=res.depr$rmse_block_p, depr.rmse_local_p=res.depr$rmse_local_p,
#			# dep/cluster
#			depc.status=res.depc$status, depc.time=res.depc$time, depc.rmse_t=res.depc$rmse_t, depc.rmse_s=res.depc$rmse_s,
#				depc.rmse_full_p=res.depc$rmse_full_p, depc.rmse_block_p=res.depc$rmse_block_p, depc.rmse_local_p=res.depc$rmse_local_p

    )
print(round(unlist(r),3)[1:20])

		r
	})


	# return results
	res.df <- as.data.frame(do.call("rbind",res))
print(colnames(res.df))
	start_Si <- which(colnames(res.df) == "orac.Si_S")
	for (i in 1:(start_Si-1)) { res.df[,i] <- unlist(res.df[,i]) }   # unlist the columns

	res.df
}

# generate training, and test data sets
"generate_data" <- function(design, factors) {
	data <- list()

	data$n <- factors$p*factors$Nper_p

	# locations to sample from
	#data$Xobs <- matrix(seq(0,1,length=data$n), nrow=data$n, ncol=1)
	data$Xobs <- randomLHS(data$n, factors$p)

	# locations to predict at
	#data$Xpred <- randomLHS(factors$Npred, factors$p)
	data$Xpred <- matrix(runif(factors$Npred*factors$p),nrow=factors$Npred,ncol=factors$p)

	data$X <- rbind(data$Xobs, data$Xpred)

	# locations for sensitivity analysis
	#data$X1_s <- data.frame(matrix(runif(factors$p*design$Nsens),nrow=design$Nsens,ncol=factors$p))
	#data$X2_s <- data.frame(matrix(runif(factors$p*design$Nsens),nrow=design$Nsens,ncol=factors$p))
	data$X1_s <- data.frame(randomLHS(design$Nsens, factors$p))
	data$X2_s <- data.frame(randomLHS(design$Nsens, factors$p))

	# set covariance parameters
	#b <- 3
	jvec <- 1:factors$p
	#data$theta <- 400
	data$theta <- factors$tau * ( (1 - (jvec-1)/factors$p)^factors$b - (1 - jvec/factors$p)^factors$b )
	#data$theta[factors$Fzero] <- 0.01    # set some to have small effect

#	# construct distance matrices
#	data$D <- array(NA, dim=c(factors$p,data$n+factors$Npred,data$n+factors$Npred))
#	lapply(1:factors$p, function(k) {
#		data$D[k,,] <<- rdist(c(data$Xobs[,k],data$Xpred[,k]))^2
#		diag(data$D[k,,]) <<- 0
#	})

	# construct covariance
#	data$Sigma <- factors$sigma2 * ce_cov(data$theta, data$D)
	data$Sigma <- factors$sigma2 * ce_cov(data$theta, data$X)

	# generate response
	y <- t(chol(data$Sigma)) %*% rnorm(data$n+factors$Npred)

	data$Yobs  <- y[1:data$n]
	data$Ypred <- y[data$n+1:factors$Npred]

	# compute sensitivity indices
	iy <- chol2inv(chol(data$Sigma[1:data$n,1:data$n])) %*% data$Yobs
	Si <- sobol2002(model = ce_full_pred_X, data$X1_s, data$X2_s, nboot = 0, Xobs=data$Xobs, iy=iy, theta=data$theta)
	data$Si_S <- Si$S[,1]
	data$Si_T <- Si$T[,1]

	data
}

# fit/evaluate models

# oracle
"eval.orac" <- function(design, factors, data) {
	status <- FALSE
	mse_p  <- NA

	t1 <- proc.time()
	try({
		# MSE of predictions
		preds <- ce_full_pred(data$Yobs, data$n, factors$Npred, data$Sigma)

		mse_p <- mean( (preds-data$Ypred)^2 )

		status <- TRUE
	})
	t2 <- proc.time()-t1

	list(
		status=status, time=as.vector(t2[3]), rmse_p=sqrt(mse_p)
	)
}

# full
"eval.full" <- function(design, factors, data) {
	status <- FALSE
	mse_t  <- NA
	mse_p  <- NA
	mse_s_S<- NA
	mse_s_T<- NA
	Si_S   <- NA
	Si_T   <- NA

	t1 <- proc.time()
	try({
		# fit full model
		fit <- hd.estimate(data, 1, rep(1,data$n), matrix(1,nrow=1,ncol=1), factors$p, log(data$theta), rep(FALSE,factors$p), ce_cov, ce_partial, TRUE)

		if (fit$convergence) {
			mse_t <- mean( (fit$theta-data$theta)^2 )

			# compute covariance for fitted model
#			fitSigma <- factors$sigma2 * ce_cov(fit$theta, data$D)
			fitSigma <- factors$sigma2 * ce_cov(fit$theta, data$X)

			preds <- ce_full_pred(data$Yobs, data$n, factors$Npred, fitSigma)
			mse_p <- mean( (preds-data$Ypred)^2 )

			# compute sensitivity indices
			iy <- chol2inv(chol(fitSigma[1:data$n,1:data$n])) %*% data$Yobs
			Si <- sobol2002(model = ce_full_pred_X, data$X1_s, data$X2_s, nboot = 0, Xobs=data$Xobs, iy=iy, theta=fit$theta)
			Si_S <- Si$S[,1]
			Si_T <- Si$T[,1]
			mse_s_S <- mean( (data$Si_S-Si_S)^2 )
			mse_s_T <- mean( (data$Si_T-Si_T)^2 )

			status <- TRUE
		}
	})
	t2 <- proc.time()-t1

	list(
		status=status, time=as.vector(t2[3]),
		rmse_t=sqrt(mse_t), rmse_p=sqrt(mse_p), rmse_s_S=sqrt(mse_s_S), rmse_s_T=sqrt(mse_s_T),
		Si_S=Si_S, Si_T=Si_T
	)
}

# independent blocks/random assignment
"eval.indr" <- function(design, factors, data) {
	status      <- FALSE
	mse_t       <- NA
	mse_full_p  <- NA
	mse_block_p <- NA
	mse_local_p <- NA
	mse_s_S     <- NA
	mse_s_T     <- NA
	theta       <- NA
	Si_S        <- NA
	Si_T        <- NA

	t1 <- proc.time()
	try({
		# assign to blocks
		sy <- sort(data$Yobs,index.return=TRUE)$ix
		nb <- ceiling(data$n/factors$Nblock_obs_ind)
		oy <- unlist(lapply(1:factors$Nblock_obs_ind, function(b) seq(1:nb)[sample(nb)]))[1:data$n]
		B  <- oy[sy]
		newB <- sample(nb,factors$Npred,replace=TRUE)

		# fit model
		fit <- hd.estimate(data, nb, B, cbind(1:nb,1:nb), factors$p, log(data$theta), rep(FALSE,factors$p), ce_cov, ce_partial, FALSE)
	})
	t2 <- proc.time()-t1

	try({
		if (fit$convergence) {
			theta <- fit$theta

			mse_t <- mean( (fit$theta-data$theta)^2 )

			# compute covariance for fitted model
#			fitSigma <- factors$sigma2 * ce_cov(fit$theta, data$D)
			fitSigma <- factors$sigma2 * ce_cov(fit$theta, data$X)

			preds <- ce_full_pred(data$Yobs, data$n, factors$Npred, fitSigma)
			mse_full_p <- mean( (preds-data$Ypred)^2 )

			#preds <- ce_block_pred(fit$theta, data$Yobs, data$n, factors$Npred, B, newB, data$D, cbind(1:nb,1:nb))
			#mse_block_p <- mean( (preds-data$Ypred)^2 )
			#preds <- ce_local_pred(data$Yobs, data$n, factors$Npred, fitSigma)
			#mse_local_p <- mean( (preds-data$Ypred)^2 )

			# compute sensitivity indices
			iy <- chol2inv(chol(fitSigma[1:data$n,1:data$n])) %*% data$Yobs
			Si <- sobol2002(model = ce_full_pred_X, data$X1_s, data$X2_s, nboot = 0, Xobs=data$Xobs, iy=iy, theta=fit$theta)
			Si_S <- Si$S[,1]
			Si_T <- Si$T[,1]
			mse_s_S <- mean( (data$Si_S-Si_S)^2 )
			mse_s_T <- mean( (data$Si_T-Si_T)^2 )

			status <- TRUE
		}
	})

	list(
		status=status, time=as.vector(t2[3]), theta=theta,
		rmse_t=sqrt(mse_t), rmse_s_S=sqrt(mse_s_S), rmse_s_T=sqrt(mse_s_T),
		rmse_full_p=sqrt(mse_full_p), #rmse_block_p=sqrt(mse_block_p), rmse_local_p=sqrt(mse_local_p),
		Si_S=Si_S, Si_T=Si_T
	)
}

# independent blocks/cluster assignment
"eval.indc" <- function(design, factors, data, init.theta) {
	status      <- FALSE
	mse_t       <- NA
	mse_full_p  <- NA
	mse_block_p <- NA
	mse_local_p <- NA
	mse_s_S     <- NA
	mse_s_T     <- NA
	Si_S        <- NA
	Si_T        <- NA

	t1 <- proc.time()
	try({
		# assign to blocks
#		initSigma <- factors$sigma2 * ce_cov(init.theta, data$D)
		initSigma <- factors$sigma2 * ce_cov(init.theta, data$X)
		hc <- hclust( as.dist( 1-initSigma[1:data$n,1:data$n] ) )
		nb <- ceiling(data$n/factors$Nblock_obs_ind)
		B  <- cutree(hc, k=nb)
		#newB <- B[data$n+1:factors$Npred]
		#B    <- B[1:data$n]

		# fit model
		fit <- hd.estimate(data, nb, B, cbind(1:nb,1:nb), factors$p, log(init.theta), rep(FALSE,factors$p), ce_cov, ce_partial, FALSE)
	})
	t2 <- proc.time()-t1

	try({
		if (fit$convergence) {
			mse_t <- mean( (fit$theta-data$theta)^2 )

			# compute covariance for fitted model
#			fitSigma <- factors$sigma2 * ce_cov(fit$theta, data$D)
			fitSigma <- factors$sigma2 * ce_cov(fit$theta, data$X)

			preds <- ce_full_pred(data$Yobs, data$n, factors$Npred, fitSigma)
			mse_full_p <- mean( (preds-data$Ypred)^2 )
			#preds <- ce_block_pred(fit$theta, data$Yobs, data$n, factors$Npred, B, newB, data$D, cbind(1:nb,1:nb))
			#mse_block_p <- mean( (preds-data$Ypred)^2 )
			#preds <- ce_local_pred(data$Yobs, data$n, factors$Npred, fitSigma)
			#mse_local_p <- mean( (preds-data$Ypred)^2 )

			# compute sensitivity indices
			iy <- chol2inv(chol(fitSigma[1:data$n,1:data$n])) %*% data$Yobs
			Si <- sobol2002(model = ce_full_pred_X, data$X1_s, data$X2_s, nboot = 0, Xobs=data$Xobs, iy=iy, theta=fit$theta)
			mse_s <- mean( (data$Si$T[,1]-Si$T[,1])^2 )
			Si_S <- Si$S[,1]
			Si_T <- Si$T[,1]
			mse_s_S <- mean( (data$Si_S-Si_S)^2 )
			mse_s_T <- mean( (data$Si_T-Si_T)^2 )

			status <- TRUE
		}
	})

	list(
		status=status, time=as.vector(t2[3]),
		rmse_t=sqrt(mse_t), rmse_s_S=sqrt(mse_s_S), rmse_s_T=sqrt(mse_s_T),
		rmse_full_p=sqrt(mse_full_p), #rmse_block_p=sqrt(mse_block_p), rmse_local_p=sqrt(mse_local_p)
		Si_S=Si_S, Si_T=Si_T
	)
}

# dependent blocks/random assignment
"eval.depr" <- function(design, factors, data) {
	status <- FALSE
	mse_t  <- NA
	mse_full_p  <- NA
	mse_block_p  <- NA
	mse_local_p  <- NA
	mse_s  <- NA

	t1 <- proc.time()
	try({
		# assign to blocks
		sy <- sort(data$Yobs,index.return=TRUE)$ix
		nb <- ceiling(data$n/factors$Nblock_obs_dep)
		oy <- unlist(lapply(1:factors$Nblock_obs_dep, function(b) seq(1:nb)[sample(nb)]))[1:data$n]
		B  <- oy[sy]
		newB <- sample(nb,factors$Npred,replace=TRUE)

		# get neighbors
		block_order <- sample(1:nb, nb)
		lag <- 3
		nmat <- c()
		sapply(1:nb, function(i) {
			for (j in i+1:lag) {
				if (j > nb) break;
				nmat <<- rbind( nmat, c(block_order[i],block_order[j]) )
			}
		})

		# fit model
		fit <- hd.estimate(data, nb, B, nmat, factors$p, log(data$theta), rep(FALSE,factors$p), ce_cov, ce_partial, FALSE)

		if (fit$convergence) {
			mse_t <- mean( (fit$theta-data$theta)^2 )

			# compute covariance for fitted model
#			fitSigma <- factors$sigma2 * ce_cov(fit$theta, data$D)
			fitSigma <- factors$sigma2 * ce_cov(fit$theta, data$X)

			preds <- ce_full_pred(data$Yobs, data$n, factors$Npred, fitSigma)
			mse_full_p <- mean( (preds-data$Ypred)^2 )
#			preds <- ce_block_pred(fit$theta, data$Yobs, data$n, factors$Npred, B, newB, data$D, nmat)
			mse_block_p <- mean( (preds-data$Ypred)^2 )
			preds <- ce_local_pred(data$Yobs, data$n, factors$Npred, fitSigma)
			mse_local_p <- mean( (preds-data$Ypred)^2 )

			status <- TRUE
		}
	})
	t2 <- proc.time()-t1

	list(
		status=status, time=as.vector(t2[3]),
		rmse_t=sqrt(mse_t), rmse_s=sqrt(mse_s),
		rmse_full_p=sqrt(mse_full_p), rmse_block_p=sqrt(mse_block_p), rmse_local_p=sqrt(mse_local_p)
	)
}

# dependent blocks/cluster assignment
"eval.depc" <- function(design, factors, data, init.theta) {
	status <- FALSE
	mse_t  <- NA
	mse_full_p  <- NA
	mse_block_p  <- NA
	mse_local_p  <- NA
	mse_s  <- NA

	t1 <- proc.time()
	try({
		# assign to blocks
		sy <- sort(data$Yobs,index.return=TRUE)$ix
		nb <- ceiling(data$n/factors$Nblock_obs_dep)
		oy <- unlist(lapply(1:factors$Nblock_obs_dep, function(b) seq(1:nb)[sample(nb)]))[1:data$n]
		B  <- oy[sy]
		newB <- sample(nb,factors$Npred,replace=TRUE)

		# get neighbors
		block_order <- unique(B[hc$order])
		lag <- 3
		nmat <- c()
		sapply(1:nb, function(i) {
			for (j in i+1:lag) {
				if (j > nb) break;
				nmat <<- rbind( nmat, c(block_order[i],block_order[j]) )
			}
		})

		# fit model
		fit <- hd.estimate(data, nb, B, nmat, factors$p, log(data$theta), rep(FALSE,factors$p), ce_cov, ce_partial, FALSE)

		if (fit$convergence) {
			mse_t <- mean( (fit$theta-data$theta)^2 )

			# compute covariance for fitted model
#			fitSigma <- factors$sigma2 * ce_cov(fit$theta, data$D)
			fitSigma <- factors$sigma2 * ce_cov(fit$theta, data$X)

			preds <- ce_full_pred(data$Yobs, data$n, factors$Npred, fitSigma)
			mse_full_p <- mean( (preds-data$Ypred)^2 )
#			preds <- ce_block_pred(fit$theta, data$Yobs, data$n, factors$Npred, B, newB, data$D, nmat)
			mse_block_p <- mean( (preds-data$Ypred)^2 )
			preds <- ce_local_pred(data$Yobs, data$n, factors$Npred, fitSigma)
			mse_local_p <- mean( (preds-data$Ypred)^2 )

			status <- TRUE
		}
	})
	t2 <- proc.time()-t1

	list(
		status=status, time=as.vector(t2[3]),
		rmse_t=sqrt(mse_t), rmse_s=sqrt(mse_s),
		rmse_full_p=sqrt(mse_full_p), rmse_block_p=sqrt(mse_block_p), rmse_local_p=sqrt(mse_local_p)
	)
}

if (FALSE) { # test sim design

# fixed design parameters
sim.design <- list(
	# number of replications
	Nreps=4
)

sim.factors <- expand.grid(
	# generate data for this number of inputs
	#p=c(5,10,15),
	p=c(10),
	#p=c(1),
	# b: sparsity
	b=c(1,3,9),
	# tau: difficulty of problem
	tau=c(3),
	# number of observations per input
	Nper_p=50,
	# number of locations to predict at
	Npred=100,
	# number of observations per block
	Nblock_obs_ind=100, Nblock_obs_dep=50,
	# which inputs are not important?
	Fzero=1,
	# observation variance?
	sigma2=1
)

}

if (FALSE) {
	set.seed(311)
	sf <- sim.factors[1,]
	dat <- generate_data(sim.design, sf)

	# number of obs per block
	npb <- 100

	# number of blocks
	nb <- ceiling(dat$n/npb)

	done
	#sx <- sort(dat$Xobs[,1],index.return=TRUE)$ix
	#pdf("pdf/data1.pdf");plot(dat$Xobs[sx,1], dat$Yobs[sx],type="l");graphics.off()
	#pdf("pdf/data2.pdf");plot(dat$Xobs[,2], dat$Yobs,type="l");graphics.off()
	#ef <- eval_full(log(dat$theta), sim.design, sim.factors[1,], dat)
	#opts <- list("algorithm"="NLOPT_LD_LBFGS", "ftol_rel"=1.0e-8)
	#res <- nloptr(x0=log(dat$theta), eval_f=eval_full, opts=opts, design=sim.design, factors=sim.factors[1,], data=dat)
	#res <- nlm(f=eval_full, p=log(dat$theta), ndigit=3, design=sim.design, factors=sim.factors[1,], data=dat, check.analyticals=FALSE)
	#print( res )
}

if (FALSE) {
	# full model
	t1 <- proc.time()
	#fit.full <- hd.estimate(dat, 1, rep(1,dat$n), matrix(1,nrow=1,ncol=1), sf$p, log(dat$theta), rep(FALSE,sf$p), ce_cov, ce_partial, TRUE)
#	pred.full <- ce_full_pred(dat$Yobs, dat$n, sf$Npred, ce_cov(fit.full$theta, dat$D))
	pred.full <- ce_full_pred(dat$Yobs, dat$n, sf$Npred, ce_cov(fit.full$theta, dat$X))
	time.full <- proc.time()-t1
	print(time.full)
}

if (FALSE) {
	# ind blocks/random assignment

	# randomly place observations into blocks
	sy <- sort(dat$Yobs,index.return=TRUE)$ix
	nb <- ceiling(dat$n/npb)
	oy <- unlist(lapply(1:npb, function(b) seq(1:nb)[sample(nb)]))[1:dat$n]
	B  <- oy[sy]
	newB <- sample(nb,sf$Npred,replace=TRUE)
	t1 <- proc.time()
	#fit.indr <- hd.estimate(dat, nb, B, cbind(1:nb,1:nb), sf$p, log(dat$theta), rep(FALSE,sf$p), ce_cov, ce_partial, TRUE)
#	pred.indr <- ce_block_pred(fit.indr$theta, dat$Yobs, dat$n, sf$Npred, B, newB, dat$D, cbind(1:nb,1:nb))
	pred.indr <- ce_block_pred(fit.indr$theta, dat$Yobs, dat$n, sf$Npred, B, newB, dat$X, cbind(1:nb,1:nb))
	time.indr <- proc.time()-t1
	print(time.indr)
}

if (FALSE) {
	# ind blocks/cluster assignment

	# cluster observations into blocks

	# cluster 1-Sigma
	hc <- hclust( as.dist( 1-dat$Sigma[1:dat$n,1:dat$n] ) )
	B  <- cutree(hc, k=nb)

	t1 <- proc.time()
	#fit.indc <- hd.estimate(dat, nb, B, cbind(1:nb,1:nb), sf$p, log(dat$theta), rep(FALSE,sf$p), ce_cov, ce_partial, TRUE)
#	pred.indc <- ce_block_pred(fit.indc$theta, dat$Yobs, dat$n, sf$Npred, B, newB, dat$D, cbind(1:nb,1:nb))
	pred.indc <- ce_block_pred(fit.indc$theta, dat$Yobs, dat$n, sf$Npred, B, newB, dat$X, cbind(1:nb,1:nb))
	time.indc <- proc.time()-t1
	print(time.indc)
}

#print(pred.full); print(pred.indr); print(pred.indc)

if (FALSE) {  # test run sims
	options(cores=4)

	# run the experiment for each combination of factors
	#res <- lapply(1:nrow(sim.factors), function(i) {
	res <- lapply(1:1, function(i) {
	  print(sim.factors[i,])
	  exp_res <- sim_exp(sim.design, sim.factors[i,], i)
		save(exp_res, file=paste0("output/exp_",i,".RData"))

print(head(exp_res))

	  exp_res
	})

	mres <- do.call(rbind, res)
}
