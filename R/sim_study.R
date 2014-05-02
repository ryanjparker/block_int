# run simulation study
library(fields)
library(lhs)
library(MASS)
library(multicore)

#source("R/create_blocks.R")
source("R/cov.R")
source("R/estimate.R")

# function to execte the simulation study based on given factors
"sim_exp" <- function(design, factors, which.exp) {

	res <- mclapply(1:design$Nreps, function(i) {
	#res <- lapply(1:design$Nreps, function(i) {
		seed <- 1983 + i + design$Nreps*(which.exp-1)
    set.seed(seed)  # set a seed for reproducibility

		# generate data
		data <- generate_data(design, factors)

		# get results for models...

			# .. oracle
			res.orac <- eval.orac(design, factors, data)

			# ... full
			#res.full <- eval.full(design, factors, data)

			# ... independent blocks/random assignment
			res.indr <- eval.indr(design, factors, data)

			# ... independent blocks/cluster assignment
			res.indc <- eval.indc(design, factors, data, res.indr$theta)

#			# ... dependent blocks/random assignment
#			res.depr <- eval.depr(design, factors, data)
#
#			# ... dependent blocks/cluster assignment
#			res.depc <- eval.depc(design, factors, data)

		# return results
		r <- list(seed=seed, p=factors$p, tau=factors$tau,
			# oracle
			orac.status=res.orac$status, orac.time=res.orac$time, orac.rmse_t=res.orac$rmse_t, orac.rmse_p=res.orac$rmse_p, orac.rmse_s=res.orac$rmse_s,
			# full
			#full.status=res.full$status, full.time=res.full$time, full.rmse_t=res.full$rmse_t, full.rmse_p=res.full$rmse_p, full.rmse_s=res.full$rmse_s,
			# ind/random
			indr.status=res.indr$status, indr.time=res.indr$time, indr.rmse_t=res.indr$rmse_t, indr.rmse_s=res.indr$rmse_s,
				indr.rmse_full_p=res.indr$rmse_full_p, #indr.rmse_block_p=res.indr$rmse_block_p, indr.rmse_local_p=res.indr$rmse_local_p,
			# ind/cluster
			indc.status=res.indc$status, indc.time=res.indc$time, indc.rmse_t=res.indc$rmse_t, indc.rmse_s=res.indc$rmse_s,
				indc.rmse_full_p=res.indc$rmse_full_p#, indc.rmse_block_p=res.indc$rmse_block_p, indc.rmse_local_p=res.indc$rmse_local_p,
#			# dep/random
#			depr.status=res.depr$status, depr.time=res.depr$time, depr.rmse_t=res.depr$rmse_t, depr.rmse_s=res.depr$rmse_s,
#				depr.rmse_full_p=res.depr$rmse_full_p, depr.rmse_block_p=res.depr$rmse_block_p, depr.rmse_local_p=res.depr$rmse_local_p,
#			# dep/cluster
#			depc.status=res.depc$status, depc.time=res.depc$time, depc.rmse_t=res.depc$rmse_t, depc.rmse_s=res.depc$rmse_s,
#				depc.rmse_full_p=res.depc$rmse_full_p, depc.rmse_block_p=res.depc$rmse_block_p, depc.rmse_local_p=res.depc$rmse_local_p

    )
print(round(unlist(r),3))

		r
	})

	# return results
	res.df <- as.data.frame(do.call("rbind",res))
	for (i in 1:ncol(res.df)) { res.df[,i] <- unlist(res.df[,i]) }   # unlist the columns

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

	# set covariance parameters
	#b <- 3
	jvec <- 1:factors$p
	#data$theta <- 400
	data$theta <- factors$tau * ( (1 - (jvec-1)/factors$p)^factors$b - (1 - jvec/factors$p)^factors$b )
	#data$theta[factors$Fzero] <- 0.01    # set some to have small effect

	# construct distance matrices
	data$D <- array(NA, dim=c(factors$p,data$n+factors$Npred,data$n+factors$Npred))
	lapply(1:factors$p, function(k) {
		data$D[k,,] <<- rdist(c(data$Xobs[,k],data$Xpred[,k]))^2
		diag(data$D[k,,]) <<- 0
	})

	# construct covariance
	data$Sigma <- factors$sigma2 * ce_cov(data$theta, data$D)

	# generate response
	y <- t(chol(data$Sigma)) %*% rnorm(data$n+factors$Npred)

	data$Yobs  <- y[1:data$n]
	data$Ypred <- y[data$n+1:factors$Npred]

	data
}

# fit/evaluate models

# oracle
"eval.orac" <- function(design, factors, data) {
	status <- FALSE
	mse_t  <- NA
	mse_p  <- NA
	mse_s  <- NA

	t1 <- proc.time()
	try({
		# MSE of predictions
		preds <- ce_full_pred(data$Yobs, data$n, factors$Npred, data$Sigma)

		mse_t <- 0
		mse_p <- mean( (preds-data$Ypred)^2 )
		mse_s <- 0

		status <- TRUE
	})
	t2 <- proc.time()-t1

	list(
		status=status, time=as.vector(t2[3]),
		rmse_t=sqrt(0), rmse_p=sqrt(mse_p), rmse_s=sqrt(0)
	)
}

# full
"eval.full" <- function(design, factors, data) {
	status <- FALSE
	mse_t  <- NA
	mse_p  <- NA
	mse_s  <- NA

	t1 <- proc.time()
	try({
		# fit full model
		fit <- hd.estimate(data, 1, rep(1,data$n), matrix(1,nrow=1,ncol=1), factors$p, log(data$theta), rep(FALSE,factors$p), ce_cov, ce_partial, FALSE)

		if (fit$convergence) {
			mse_t <- mean( (fit$theta-data$theta)^2 )

			# compute covariance for fitted model
			fitSigma <- factors$sigma2 * ce_cov(fit$theta, data$D)

			preds <- ce_full_pred(data$Yobs, data$n, factors$Npred, fitSigma)
			mse_p <- mean( (preds-data$Ypred)^2 )

			status <- TRUE
		}
	})
	t2 <- proc.time()-t1

	list(
		status=status, time=as.vector(t2[3]),
		rmse_t=sqrt(mse_t), rmse_p=sqrt(mse_p), rmse_s=sqrt(mse_s)
	)
}

# independent blocks/random assignment
"eval.indr" <- function(design, factors, data) {
	status <- FALSE
	mse_t  <- NA
	mse_full_p  <- NA
	mse_block_p  <- NA
	mse_local_p  <- NA
	mse_s  <- NA
	theta <- NA

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
			fitSigma <- factors$sigma2 * ce_cov(fit$theta, data$D)

			preds <- ce_full_pred(data$Yobs, data$n, factors$Npred, fitSigma)
			mse_full_p <- mean( (preds-data$Ypred)^2 )

			#preds <- ce_block_pred(fit$theta, data$Yobs, data$n, factors$Npred, B, newB, data$D, cbind(1:nb,1:nb))
			#mse_block_p <- mean( (preds-data$Ypred)^2 )
			#preds <- ce_local_pred(data$Yobs, data$n, factors$Npred, fitSigma)
			#mse_local_p <- mean( (preds-data$Ypred)^2 )

			status <- TRUE
		}
	})

	list(
		status=status, time=as.vector(t2[3]), theta=theta,
		rmse_t=sqrt(mse_t), rmse_s=sqrt(mse_s),
		rmse_full_p=sqrt(mse_full_p) #, rmse_block_p=sqrt(mse_block_p), rmse_local_p=sqrt(mse_local_p)
	)
}

# independent blocks/cluster assignment
"eval.indc" <- function(design, factors, data, init.theta) {
	status <- FALSE
	mse_t  <- NA
	mse_full_p  <- NA
	mse_block_p  <- NA
	mse_local_p  <- NA
	mse_s  <- NA

	t1 <- proc.time()
	try({
		# assign to blocks
		initSigma <- factors$sigma2 * ce_cov(init.theta, data$D)
		hc <- hclust( as.dist( 1-initSigma[1:data$n,1:data$n] ) )
		nb <- ceiling(data$n/factors$Nblock_obs_ind)
		B  <- cutree(hc, k=nb)
		#newB <- B[data$n+1:factors$Npred]
		#B    <- B[1:data$n]

		# fit model
		fit <- hd.estimate(data, nb, B, cbind(1:nb,1:nb), factors$p, log(init.theta), rep(FALSE,factors$p), ce_cov, ce_partial, TRUE)
	})
	t2 <- proc.time()-t1

	try({
		if (fit$convergence) {
			mse_t <- mean( (fit$theta-data$theta)^2 )

			# compute covariance for fitted model
			fitSigma <- factors$sigma2 * ce_cov(fit$theta, data$D)

			preds <- ce_full_pred(data$Yobs, data$n, factors$Npred, fitSigma)
			mse_full_p <- mean( (preds-data$Ypred)^2 )
			#preds <- ce_block_pred(fit$theta, data$Yobs, data$n, factors$Npred, B, newB, data$D, cbind(1:nb,1:nb))
			#mse_block_p <- mean( (preds-data$Ypred)^2 )
			#preds <- ce_local_pred(data$Yobs, data$n, factors$Npred, fitSigma)
			#mse_local_p <- mean( (preds-data$Ypred)^2 )

			status <- TRUE
		}
	})

	list(
		status=status, time=as.vector(t2[3]),
		rmse_t=sqrt(mse_t), rmse_s=sqrt(mse_s),
		rmse_full_p=sqrt(mse_full_p) #, rmse_block_p=sqrt(mse_block_p), rmse_local_p=sqrt(mse_local_p)
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
			fitSigma <- factors$sigma2 * ce_cov(fit$theta, data$D)

			preds <- ce_full_pred(data$Yobs, data$n, factors$Npred, fitSigma)
			mse_full_p <- mean( (preds-data$Ypred)^2 )
			preds <- ce_block_pred(fit$theta, data$Yobs, data$n, factors$Npred, B, newB, data$D, nmat)
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
"eval.depc" <- function(design, factors, data) {
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
			fitSigma <- factors$sigma2 * ce_cov(fit$theta, data$D)

			preds <- ce_full_pred(data$Yobs, data$n, factors$Npred, fitSigma)
			mse_full_p <- mean( (preds-data$Ypred)^2 )
			preds <- ce_block_pred(fit$theta, data$Yobs, data$n, factors$Npred, B, newB, data$D, nmat)
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
	pred.full <- ce_full_pred(dat$Yobs, dat$n, sf$Npred, ce_cov(fit.full$theta, dat$D))
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
	pred.indr <- ce_block_pred(fit.indr$theta, dat$Yobs, dat$n, sf$Npred, B, newB, dat$D, cbind(1:nb,1:nb))
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
	pred.indc <- ce_block_pred(fit.indc$theta, dat$Yobs, dat$n, sf$Npred, B, newB, dat$D, cbind(1:nb,1:nb))
	time.indc <- proc.time()-t1
	print(time.indc)
}

#print(pred.full); print(pred.indr); print(pred.indc)

if (TRUE) {
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
