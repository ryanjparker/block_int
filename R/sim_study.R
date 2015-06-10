# code to run simulation study
library(fields)
library(lhs)
library(MASS)
library(multicore)
library(sensitivity)
#library(spacious)

#source("R/create_blocks.R")
source("R/cov.R")
source("R/estimate.R")

# function to execte the estimation simulation study based on given factors
"sim_exp_est" <- function(design, factors, which.exp, which.part) {

	exp.step  <- design$Nreps/20
	exp.start <- ((which.part-1)*exp.step+1)
	exp.end   <- exp.start+exp.step-1

	#res <- mclapply(1:design$Nreps, function(i) {
	#res <- mclapply(1:3, function(i) {
	#res <- lapply(1:design$Nreps, function(i) {
	#res <- lapply(1:1, function(i) { #design$Nreps, function(i) {
	res <- lapply(exp.start:exp.end, function(i) {

		seed <- 1983 + i + design$Nreps*(which.exp-1)
		set.seed(seed)  # set a seed for reproducibility

		# generate data
		data <- generate_data(design, factors)

		# get results for models...

			# .. oracle
			res.orac <- eval.orac(design, factors, data)

			# ... full
			res.full <- eval.full(design, factors, data)

			# ... independent blocks

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

			# ... dependent blocks

			# 50 obs per block
			factors$Nblock_obs_dep <- 50
				res.depr <- eval.depr(design, factors, data)
				res.depc <- eval.depc(design, factors, data, res.depr$theta)
				rdr_50 <- res.depr
				rdc_50 <- res.depc

			# 25 obs per block
			factors$Nblock_obs_dep <- 25
				res.depr <- eval.depr(design, factors, data)
				res.depc <- eval.depc(design, factors, data, res.depr$theta)
				rdr_25 <- res.depr
				rdc_25 <- res.depc

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
			# dep random/cluster
			dr50.status=rdr_50$status, dr50.time=rdr_50$time, dr50.rmse_t=rdr_50$rmse_t, dr50.rmse_full_p=rdr_50$rmse_full_p,
			dr50.rmse_s_S=rdr_50$rmse_s_T, dr50.rmse_s_T=rdr_50$rmse_s_T,
			dc50.status=rdc_50$status, dc50.time=rdc_50$time, dc50.rmse_t=rdc_50$rmse_t, dc50.rmse_full_p=rdc_50$rmse_full_p,
			dc50.rmse_s_S=rdc_50$rmse_s_T, dc50.rmse_s_T=rdc_50$rmse_s_T,
			dr25.status=rdr_25$status, dr25.time=rdr_25$time, dr25.rmse_t=rdr_25$rmse_t, dr25.rmse_full_p=rdr_25$rmse_full_p,
			dr25.rmse_s_S=rdr_25$rmse_s_T, dr25.rmse_s_T=rdr_25$rmse_s_T,
			dc25.status=rdc_25$status, dc25.time=rdc_25$time, dc25.rmse_t=rdc_25$rmse_t, dc25.rmse_full_p=rdc_25$rmse_full_p,
			dc25.rmse_s_S=rdc_25$rmse_s_T, dc25.rmse_s_T=rdc_25$rmse_s_T,

			# sensitivity indices
			orac.Si_S=data$Si_S, orac.Si_T=data$Si_T,
			full.Si_S=res.full$Si_S, full.Si_T=res.full$Si_T,
			# ind
			ir250.Si_S=rir_250$Si_S, ir250.Si_T=rir_250$Si_T,
			ic250.Si_S=ric_250$Si_S, ic250.Si_T=ric_250$Si_T,
			ir100.Si_S=rir_100$Si_S, ir100.Si_T=rir_100$Si_T,
			ic100.Si_S=ric_100$Si_S, ic100.Si_T=ric_100$Si_T,
			ir50.Si_S=rir_50$Si_S, ir50.Si_T=rir_50$Si_T,
			ic50.Si_S=ric_50$Si_S, ic50.Si_T=ric_50$Si_T,
			ir25.Si_S=rir_25$Si_S, ir25.Si_T=rir_25$Si_T,
			ic25.Si_S=ric_25$Si_S, ic25.Si_T=ric_25$Si_T,
			# dep
			dr50.Si_S=rdr_50$Si_S, dr50.Si_T=rdr_50$Si_T,
			dc50.Si_S=rdc_50$Si_S, dc50.Si_T=rdc_50$Si_T,
			dr25.Si_S=rdr_25$Si_S, dr25.Si_T=rdr_25$Si_T,
			dc25.Si_S=rdc_25$Si_S, dc25.Si_T=rdc_25$Si_T

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

# function to execte the prediction simulation study based on given factors
"sim_exp_pred" <- function(design, factors, which.exp, which.part) {

	exp.step  <- design$Nreps/20
	exp.start <- ((which.part-1)*exp.step+1)
	exp.end   <- exp.start+exp.step-1

	#res <- mclapply(1:design$Nreps, function(i) {
	#res <- mclapply(1:3, function(i) {
	#res <- lapply(1:design$Nreps, function(i) {
	#res <- lapply(1:1, function(i) { #design$Nreps, function(i) {
	res <- lapply(exp.start:exp.end, function(i) {
		seed <- 1983 + i + design$Nreps*(which.exp-1)
		set.seed(seed)  # set a seed for reproducibility

		# generate data
		data <- generate_data(design, factors)
#str(data)
done

		r <- list()

		# compute prediction error...

if (TRUE) {
			# all data
				res <- pred.full(design, factors, data); r <- c(r, full=res)
print(r)
}

if (TRUE) {
			# ... local kriging
				# 25
				res <- pred.local(design, factors, data, 25);  r <- c(r, l25=res)
				# 50
				res <- pred.local(design, factors, data, 50);  r <- c(r, l50=res)
				# 100
				res <- pred.local(design, factors, data, 100); r <- c(r, l100=res)
				# 250
				#res <- pred.local(design, factors, data, 250); r <- c(r, l250=res)
				# 500
				#res <- pred.local(design, factors, data, 500); r <- c(r, l500=res)
				# 1000
				#res <- pred.local(design, factors, data, 1000); r <- c(r, l1000=res)
}

			# ... block kriging
if (FALSE) {
				# independent ordered
					# 25
					res <- pred.order(design, factors, data, 25, FALSE);  r <- c(r, o25_i=res)
					# 50
					res <- pred.order(design, factors, data, 50, FALSE);  r <- c(r, o50_i=res)
					# 100
					res <- pred.order(design, factors, data, 100, FALSE);  r <- c(r, o100_i=res)
					# 250
					res <- pred.order(design, factors, data, 250, FALSE);  r <- c(r, o250_i=res)
}

if (FALSE) {
				# dependent ordered
					# 25
					res <- pred.order(design, factors, data, 25, TRUE);  r <- c(r, o25_d=res)
					# 50
					res <- pred.order(design, factors, data, 50, TRUE);  r <- c(r, o50_d=res)
					# 100
					res <- pred.order(design, factors, data, 100, TRUE);  r <- c(r, o100_d=res)
					# 250
					res <- pred.order(design, factors, data, 250, TRUE);  r <- c(r, o250_d=res)
}

if (FALSE) {
if (FALSE) {
				# perform the clustering once
				t1 <- proc.time()
				try({
					hc <- hclust( as.dist( 1-data$Sigma[1:data$n,1:data$n] ) )
				})
				t2 <- proc.time()-t1
				r <- c(r, clust.time=t2[3])
print(r)
}

if (FALSE) {
str(data$Xobs)
				t1 <- proc.time()
				try({
str(km)
				})
				t2 <- proc.time()-t1
				r <- c(r, clust.time=t2[3])
print(r)
done
}

if (FALSE) {
				# independent clustering
					# 25
					res <- pred.clust(design, factors, data, 25, FALSE);  r <- c(r, c25_i=res)
					# 50
					res <- pred.clust(design, factors, data, 50, FALSE);  r <- c(r, c50_i=res)
					# 100
					res <- pred.clust(design, factors, data, 100, FALSE);  r <- c(r, c100_i=res)
					# 250
					res <- pred.clust(design, factors, data, 250, FALSE);  r <- c(r, c250_i=res)
					# 500
					res <- pred.clust(design, factors, data, 500, FALSE);  r <- c(r, c500_i=res)
}

if (FALSE) {
				# dependent clustering
					# 25
#					res <- pred.clust(design, factors, data, 25, TRUE);  r <- c(r, c25_d=res)
					# 50
#					res <- pred.clust(design, factors, data, 50, TRUE);  r <- c(r, c50_d=res)
					# 100
					res <- pred.clust(design, factors, data, 100, TRUE);  r <- c(r, c100_d=res)
					# 250
					res <- pred.clust(design, factors, data, 250, TRUE);  r <- c(r, c250_d=res)
					# 500
					res <- pred.clust(design, factors, data, 500, TRUE);  r <- c(r, c500_d=res)
					# 1000
#					res <- pred.clust(design, factors, data, 1000, TRUE);  r <- c(r, c1000_d=res)
}
}

if (TRUE) {
			# ... best subset
				# 25
				res <- pred.sub(design, factors, data, 25);  r <- c(r, s25=res)
				# 50
				res <- pred.sub(design, factors, data, 50);  r <- c(r, s50=res)
				# 100
				res <- pred.sub(design, factors, data, 100);  r <- c(r, s100=res)
				# 250
				res <- pred.sub(design, factors, data, 250);  r <- c(r, s250=res)
				# 500
				res <- pred.sub(design, factors, data, 500);  r <- c(r, s500=res)
				# 1000
				res <- pred.sub(design, factors, data, 1000);  r <- c(r, s1000=res)
				# 2000
#				res <- pred.sub(design, factors, data, 2000);  r <- c(r, s2000=res)
				# 3000
#				res <- pred.sub(design, factors, data, 3000);  r <- c(r, s3000=res)
print(r)
}

		# return results
		r <- c(seed=seed, p=factors$p, r)
print(round(unlist(r),3))

		unlist(r)
	})

	# return results
	res.df <- as.data.frame(do.call("rbind",res))
#print(colnames(res.df))
	for (i in 1:ncol(res.df)) { res.df[,i] <- unlist(res.df[,i]) }   # unlist the columns

	res.df
}

# generate training, and test data sets
"generate_data" <- function(design, factors) {
	data <- list()

	isp <- factors$pred

	if (isp) {
		# prediction study
		data$n <- factors$Nper_p
	} else {
		# estimation study
		data$n <- factors$p*factors$Nper_p
	}

	# locations to sample from
	if (isp & exists("predXobs")) {
		data$Xobs <- predXobs
	} else {
		#data$Xobs <- matrix(seq(0,1,length=data$n), nrow=data$n, ncol=1)
		data$Xobs <- randomLHS(data$n, factors$p+factors$Fzero)
		#data$Xobs <- maximinLHS(data$n, factors$p)
		if (isp) predXobs  <<- data$Xobs
	}

	# locations to predict at
	if (isp & exists("predXpred")) {
		data$Xpred <- predXpred
	} else {
		#data$Xpred <- randomLHS(factors$Npred, factors$p)
		data$Xpred <- matrix(runif(factors$Npred*(factors$p+factors$Fzero)),nrow=factors$Npred,ncol=factors$p+factors$Fzero)
		if (isp) predXpred  <<- data$Xpred
	}

	data$X <- rbind(data$Xobs, data$Xpred)

	# set covariance parameters
	#b <- 3
	jvec <- 1:factors$p
	#data$theta <- 400
	data$theta <- factors$tau * ( (1 - (jvec-1)/factors$p)^factors$b - (1 - jvec/factors$p)^factors$b )
	#data$theta[factors$Fzero] <- 0.01    # set some to have small effect
#print( round(data$theta, 4)); done

if (FALSE) {
print(sum( data$theta^2 ))

fx <- function(x) { sum( (factors$tau * ( (1 - (jvec-1)/factors$p)^x - (1 - jvec/factors$p)^x ))^2 ) - 0.5 }
print(uniroot(fx, c(1,100)))

print(sum( data$theta ))
print(sum( data$theta^2 ))
done
}

if (FALSE) {
cat("Constructing D\n")
	# construct distance matrices
	data$D <- array(NA, dim=c(factors$p,data$n+factors$Npred,data$n+factors$Npred))
	lapply(1:factors$p, function(k) {
		data$D[k,,] <<- rdist(c(data$Xobs[,k],data$Xpred[,k]))^2
		diag(data$D[k,,]) <<- 0
	})
}

	if (isp & exists("predSigma")) {
		data$Sigma     <- predSigma
		data$cholSigma <- predCholSigma
	} else {
#cat("Computing Sigma\n")
	# construct covariance
#	data$Sigma <- factors$sigma2 * ce_cov(data$theta, data$D)

		data$Sigma     <- factors$sigma2 * ce_cov(data$theta, data$X)

		if (isp) {
			#data$cholSigma <- gpuChol(data$Sigma)
			data$cholSigma <- chol(data$Sigma)
			predSigma     <<- data$Sigma
			predCholSigma <<- data$cholSigma
		} else {
			data$cholSigma <- chol(data$Sigma)
		}
	}

	if (isp & !exists("predInvSigma")) {
		#predInvSigma <<- gpuChol2Inv(data$Sigma[1:data$n,1:data$n])
		predInvSigma <<- chol2inv(data$Sigma[1:data$n,1:data$n])
	}

#cat("Generating y\n")
	# generate response
	#y <- t(chol(data$Sigma)) %*% rnorm(data$n+factors$Npred)
	#cholSigma <- chol(data$Sigma) #[1:2000,1:2000])
	z <- matrix(rnorm(data$n+factors$Npred), ncol=1)
#t1 <- proc.time()
#	y <- as.vector( gpuMM( t(data$cholSigma), z ) )
#print(proc.time()-t1)
#print( summary(y) )
#t1 <- proc.time()
	y <- as.vector( t(data$cholSigma) %*% z )
print(sd(y))
#print(proc.time()-t1)
#print( summary(y) )
#done

	data$Yobs  <- y[1:data$n]
	data$Ypred <- y[data$n+1:factors$Npred]

	if (!isp) {
		# locations for sensitivity analysis
		#data$X1_s <- data.frame(matrix(runif(factors$p*design$Nsens),nrow=design$Nsens,ncol=factors$p))
		#data$X2_s <- data.frame(matrix(runif(factors$p*design$Nsens),nrow=design$Nsens,ncol=factors$p))
		data$X1_s <- data.frame(randomLHS(design$Nsens, factors$p+factors$Fzero))
		data$X2_s <- data.frame(randomLHS(design$Nsens, factors$p+factors$Fzero))

		# compute sensitivity indices
		iy <- chol2inv(chol(data$Sigma[1:data$n,1:data$n])) %*% data$Yobs
		Si <- sobol2002(model = ce_full_pred_X, data$X1_s, data$X2_s, nboot = 0, Xobs=data$Xobs, iy=iy, theta=data$theta)
		data$Si_S <- Si$S[,1]
		data$Si_T <- Si$T[,1]
	}

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
		fit <- hd.estimate(data, 1, rep(1,data$n), matrix(1,nrow=1,ncol=1), factors$p, log(data$theta), rep(FALSE,factors$p), ce_cov, ce_partial, FALSE)

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
		nb <- ceiling(data$n/factors$Nblock_obs_ind)

if (FALSE) {
		sy <- sort(data$Yobs,index.return=TRUE)$ix
		oy <- unlist(lapply(1:factors$Nblock_obs_ind, function(b) seq(1:nb)[sample(nb)]))[1:data$n]
		#B  <- oy[sy]
		B  <- rep(NA, length(sy))
		B[sy]  <- oy
		newB <- sample(nb,factors$Npred,replace=TRUE)
} else { # no randomness
		sy <- sort(data$Yobs,index.return=TRUE)$ix
		oy <- unlist(lapply(1:factors$Nblock_obs_ind, function(b) seq(1:nb)))[1:data$n]
		B  <- rep(NA, length(sy))
		B[sy]  <- oy
}

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
		nb <- ceiling(data$n/factors$Nblock_obs_ind)

if (FALSE) { # use estimated cov matrix
#		initSigma <- factors$sigma2 * ce_cov(init.theta, data$D)
		initSigma <- factors$sigma2 * ce_cov(init.theta, data$X)
		hc <- hclust( as.dist( 1-initSigma[1:data$n,1:data$n] ) )
		B  <- cutree(hc, k=nb)
		#newB <- B[data$n+1:factors$Npred]
		#B    <- B[1:data$n]
} else { # use observation locations
		km <- kmeans(data$Xobs, nb, iter.max=100)
		B  <- km$cluster
}

		# fit model
		fit <- hd.estimate(data, nb, B, cbind(1:nb,1:nb), factors$p, log(data$theta), rep(FALSE,factors$p), ce_cov, ce_partial, FALSE)
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
		nb <- ceiling(data$n/factors$Nblock_obs_dep)

if (FALSE) {
		sy <- sort(data$Yobs,index.return=TRUE)$ix
		oy <- unlist(lapply(1:factors$Nblock_obs_dep, function(b) seq(1:nb)[sample(nb)]))[1:data$n]
		#B  <- oy[sy]
		B  <- rep(NA, length(sy))
		B[sy]  <- oy
		newB <- sample(nb,factors$Npred,replace=TRUE)

		# get neighbors
		block_order <- sample(1:nb, nb)
		lag <- 2
		nmat <- c()
		sapply(1:nb, function(i) {
			for (j in i+1:lag) {
				if (j > nb) break;
				nmat <<- rbind( nmat, c(block_order[i],block_order[j]) )
			}
		})
} else {
		sy <- sort(data$Yobs,index.return=TRUE)$ix
		oy <- unlist(lapply(1:factors$Nblock_obs_dep, function(b) seq(1:nb)))[1:data$n]
		B  <- rep(NA, length(sy))
		B[sy]  <- oy

		centerB <- do.call("cbind", lapply(1:nb, function(b) { colMeans(data$Xobs[B==b,]) }))

		# order blocks based on (1) distance from 0, then (2) distance from each other
		to0 <- rdist(matrix(0, nrow=1, ncol=factors$p), centerB)
		blockD <- rdist(centerB)
		diag(blockD) <- Inf

		Border <- rep(NA, nb)
		Border[1] <- which.min(to0)

		for (b in 2:nb) {
			Border[b] <- which.min(blockD[Border[b-1],])
			blockD[Border[1:(b-1)],] <- Inf
			blockD[,Border[1:(b-1)]] <- Inf
		}

		# get neighbors
		block_order <- Border
		lag <- 1
		nmat <- c()
		sapply(1:nb, function(i) {
			for (j in i+1:lag) {
				if (j > nb) break;
				nmat <<- rbind( nmat, c(block_order[i],block_order[j]) )
			}
		})
}

		# fit model
		fit <- hd.estimate(data, nb, B, nmat, factors$p, log(data$theta), rep(FALSE,factors$p), ce_cov, ce_partial, FALSE)
	})
	t2 <- proc.time()-t1

	try({
		if (fit$convergence) {
			theta <- fit$theta

			mse_t <- mean( (fit$theta-data$theta)^2 )

			# compute covariance for fitted model
			fitSigma <- factors$sigma2 * ce_cov(fit$theta, data$X)

			preds <- ce_full_pred(data$Yobs, data$n, factors$Npred, fitSigma)
			mse_full_p <- mean( (preds-data$Ypred)^2 )

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
		rmse_full_p=sqrt(mse_full_p),
		Si_S=Si_S, Si_T=Si_T
	)
}

# dependent blocks/cluster assignment
"eval.depc" <- function(design, factors, data, init.theta) {
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
		nb <- ceiling(data$n/factors$Nblock_obs_dep)

if (FALSE) { # estimated cov matrix
		initSigma <- factors$sigma2 * ce_cov(init.theta, data$X)
		hc <- hclust( as.dist( 1-initSigma[1:data$n,1:data$n] ) )
		B  <- cutree(hc, k=nb)

		# get neighbors
		block_order <- unique(B[hc$order])
		lag <- 2
		nmat <- c()
		sapply(1:nb, function(i) {
			for (j in i+1:lag) {
				if (j > nb) break;
				nmat <<- rbind( nmat, c(block_order[i],block_order[j]) )
			}
		})
} else { # k-means on observation locations
		km <- kmeans(data$Xobs, nb, iter.max=100)
		B  <- km$cluster

		# order blocks based on (1) distance from 0, then (2) distance from each other
		to0 <- rdist(matrix(0, nrow=1, ncol=factors$p), km$centers)
		blockD <- rdist(km$centers)
		diag(blockD) <- Inf

		Border <- rep(NA, nb)
		Border[1] <- which.min(to0)

		for (b in 2:nb) {
			Border[b] <- which.min(blockD[Border[b-1],])
			blockD[Border[1:(b-1)],] <- Inf
			blockD[,Border[1:(b-1)]] <- Inf
		}

		# get neighbors
		block_order <- Border
		lag <- 1
		nmat <- c()
		sapply(1:nb, function(i) {
			for (j in i+1:lag) {
				if (j > nb) break;
				nmat <<- rbind( nmat, c(block_order[i],block_order[j]) )
			}
		})
}

		# fit model
		fit <- hd.estimate(data, nb, B, nmat, factors$p, log(data$theta), rep(FALSE,factors$p), ce_cov, ce_partial, FALSE)

	})
	t2 <- proc.time()-t1

	try({
		if (fit$convergence) {
			mse_t <- mean( (fit$theta-data$theta)^2 )

			# compute covariance for fitted model
			fitSigma <- factors$sigma2 * ce_cov(fit$theta, data$X)

			preds <- ce_full_pred(data$Yobs, data$n, factors$Npred, fitSigma)
			mse_full_p <- mean( (preds-data$Ypred)^2 )

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
		rmse_full_p=sqrt(mse_full_p),
		Si_S=Si_S, Si_T=Si_T
	)
}

# prediction

# full
"pred.full" <- function(design, factors, data) {
	mse  <- NA

	t1 <- proc.time()
	try({
		#preds <- ce_full_pred(data$Yobs, data$n, factors$Npred, data$Sigma)
#as.vector( Sigma[Nfit+1:Npred,1:Nfit] %*% chol2inv(chol(Sigma[1:Nfit,1:Nfit])) %*% y )

		#preds <- gpuMM(data$Sigma[data$n+1:factors$Npred,1:data$n], predInvSigma) %*% data$Yobs
		preds <- data$Sigma[data$n+1:factors$Npred,1:data$n] %*% (predInvSigma %*% data$Yobs)
		mse <- mean( (preds-data$Ypred)^2 )
	})
	t2 <- proc.time()-t1

	list( time=as.vector(t2[3]), rmse=sqrt(mse) )
}

# local
"pred.local" <- function(design, factors, data, Nlocal) {
	mse  <- NA

	t1 <- proc.time()
	try({
#cat("Predicting for N=",factors$Npred," sites\n", sep="")
#cat("Theta:\n"); print(theta)

		# determine Nlocal closest observations
		#with(data, print(D[1,1,n+1:factors$Npred]) )
if (FALSE) {
		# closest by distance
		closest <- lapply(1:factors$Npred, function(px) {
			#dists <- apply(data$D[,data$n+px,1:data$n], 2, function(col) { sqrt( sum(col^2) ) })
			#dists <- apply(data$D[,data$n+px,1:data$n], 2, function(col) { sum(data$theta*col) })
			dists <- colSums( data$theta * sqrt( ( t(data$Xobs) - data$Xpred[px,])^2 ) )
			locs  <- sort(dists,index.return=TRUE)$ix[1:Nlocal]
			#print(round(dists,3)); print(locs)

			locs
		})
} else {
		# closest by covariance
		closest <- lapply(1:factors$Npred, function(px) {
			#dists <- apply(data$D[,data$n+px,1:data$n], 2, function(col) { sqrt( sum(col^2) ) })
			#dists <- apply(data$D[,data$n+px,1:data$n], 2, function(col) { sum(data$theta*col) })
			dists <- 1-data$Sigma[data$n+px,1:data$n]
			locs  <- sort(dists,index.return=TRUE)$ix[1:Nlocal]
			#print(round(dists,3)); print(locs)

			locs
		})
#print(closest)
}

		# get predictions
		#if (Nlocal >= 500) gpu <- TRUE
		#else               gpu <- FALSE
		preds <- sapply(1:factors$Npred, function(px) {
			ce_full_pred(data$Yobs[closest[[px]]], Nlocal, 1, data$Sigma[c(closest[[px]],data$n+px),c(closest[[px]],data$n+px)]) #,gpu=gpu)
		})

		mse <- mean( (preds-data$Ypred)^2 )
	})
	t2 <- proc.time()-t1

	list( time=as.vector(t2[3]), rmse=sqrt(mse) )
}

"do.pred_blocks" <- function(design, factors, data, NperB, Nclose, centerB, nb, B) {
	newB <- apply(data$Xpred, 1, function(row) {
		#print( ( data$theta*sqrt( (row-centerB)^2 ) ) )
		which.min( colSums( data$theta*sqrt( (row-centerB)^2 ) ) )
	})
#print(newB)

	closestB <- apply(data$Xpred, 1, function(row) {
		#print( ( data$theta*sqrt( (row-centerB)^2 ) ) )
		sorted <- sort( colSums( data$theta*sqrt( (row-centerB)^2 ) ), index.return=TRUE )

		list(vals=sorted$x[1:Nclose], id=sorted$ix[1:Nclose])
	})
#print(closestB[[1]])
#print( (1/closestB[[1]]$vals)/sum(1/closestB[[1]]$vals))

	# get inverses
	invs <- lapply(1:nb, function(b) {
		idb1 <- which(B==b)

		N1 <- length(idb1)

		if (N1 > 0) {
			chol2inv(chol(data$Sigma[idb1,idb1]))
		} else {
			NA
		}
	})

	# get predictions
	preds <- rep(NA, factors$Npred)

if (FALSE) { # closest block
	sapply(1:nb, function(b) {
		idb1 <- which(B==b)
		idb2 <- which(newB==b)

		N1 <- length(idb1)
		N2 <- length(idb2)

		if (N1 > 0 & N2 > 0) {
			preds[idb2] <<- ce_full_pred(data$Yobs[idb1], N1, N2, data$Sigma[c(idb1,data$n+idb2),c(idb1,data$n+idb2)])
		}
	})
} else {  # Nclose blocks
#as.vector( Sigma[Nfit+1:Npred,1:Nfit] %*% chol2inv(chol(Sigma[1:Nfit,1:Nfit])) %*% y )
	preds <- sapply(1:factors$Npred, function(i) {
		w <- (1/closestB[[i]]$vals)/sum(1/closestB[[i]]$vals)
		y <- sapply(1:Nclose, function(b) {
			idb1 <- which(B==closestB[[i]]$id[b])
			data$Sigma[data$n+i,idb1] %*% invs[[ closestB[[i]]$id[b] ]] %*% data$Yobs[idb1]
		})

		sum(w*y)
#print(y); print(sum(w*y)); print(data$Ypred[i]); done
	})
}
#print(preds)

	mse <- mean( (preds-data$Ypred)^2 )
#print(mse)

	mse
}

# do independent block prediction
"do.pred_blocks_ind" <- function(design, factors, data, centerB, nb, B, invy) {
	newB <- apply(data$Xpred, 1, function(row) {
		#print( ( data$theta*sqrt( (row-centerB)^2 ) ) )
		which.min( colSums( data$theta*sqrt( (row-centerB)^2 ) ) )
	})
#print(newB)

	# get predictions
	preds <- rep(NA, factors$Npred)

	sapply(1:nb, function(b) {
		idb1 <- which(B==b)
		idb2 <- which(newB==b)

		N1 <- length(idb1)
		N2 <- length(idb2)

		if (N1 > 0 & N2 > 0) {
			preds[idb2] <<- data$Sigma[data$n+idb2,idb1] %*% invy[[b]]
		}
	})
#print(preds)

	mse <- mean( (preds-data$Ypred)^2 )
#print(mse)

	mse
}

# do dependent block prediction
"do.pred_blocks_dep" <- function(design, factors, data, centerB, nb, B, nmat) {
	newB <- apply(data$Xpred, 1, function(row) {
		#print( ( data$theta*sqrt( (row-centerB)^2 ) ) )
		which.min( colSums( data$theta*sqrt( (row-centerB)^2 ) ) )
	})
#print(newB)

	# get predictions
	y_0 <- rep(NA, factors$Npred)

	unewB <- sort(unique(newB))
#print(unewB)

	sapply(unewB, function(b) {
		# neighbors of this block
		neighbors <- as.vector(nmat[which( rowSums( nmat==b ) == 1 ),])
		neighbors <- sort(neighbors[neighbors != b])

		# what new points are in this block?
		in.new   <- which(newB == b)+data$n
		n.in.new <- length(in.new)

		# what observed points are in this block?
		in.obs   <- which(B == b)
		n.in.obs <- length(in.obs)

		# new and observed in this block
		in.b   <- c(in.new, in.obs)
		n.in.b <- n.in.new + n.in.obs

		A_0 <- matrix(0, nrow=n.in.new, ncol=n.in.new)
		b_0 <- rep(0, n.in.new)

		sapply(neighbors, function(n1) {
			in.n1   <- which(B == n1)
			n.in.n1 <- length(in.n1)

			# take points so that we have Y = (Y_{b,new}, Y_{b,obs}, Y_{n1})
			in.pair.n1 <- c(in.b, in.n1)

			# covariance for b and n1
			Sigma.n1    <- data$Sigma[in.pair.n1,in.pair.n1]
			invSigma.n1 <- chol2inv(chol(Sigma.n1))

			A_0 <<- A_0 + invSigma.n1[1:n.in.new,1:n.in.new]
			b_0 <<- b_0 +
				matrix(invSigma.n1[1:n.in.new,n.in.new+1:n.in.obs],nrow=n.in.new,ncol=n.in.obs) %*% data$Yobs[in.obs] +
				matrix(invSigma.n1[1:n.in.new,n.in.new+n.in.obs+1:n.in.n1],nrow=n.in.new,ncol=n.in.n1) %*% data$Yobs[in.n1]
		})

		# complete predictions
		b_0 <- -b_0
		y_0[in.new-data$n] <<- chol2inv(chol(A_0)) %*% b_0
	})

	mse <- mean( (y_0-data$Ypred)^2 )
#print(sqrt(mse))

	mse
}

"pred.order" <- function(design, factors, data, NperB, dep) {
	mse  <- NA

	t1 <- proc.time()
	try({
		# assign to blocks
		sy <- sort(data$Yobs,index.return=TRUE)$ix
		nb <- ceiling(data$n/NperB)
		oy <- unlist(lapply(1:NperB, function(b) seq(1:nb)))[1:data$n]
		B  <- rep(NA, length(sy))
		B[sy]  <- oy

		# get prediction blocks using block centroids
		centerB <- do.call("cbind", lapply(1:nb, function(b) { colMeans(data$Xobs[B==b,]) }))
#print(centerB)

		if (dep) {
			# create neighbor structure

			nmat <- c()

if (FALSE) {
			block_order <- sample(1:nb, nb)
			lag <- 2
			sapply(1:nb, function(i) {
				for (j in i+1:lag) {
					if (j > nb) break;
					nmat <<- rbind( nmat, c(block_order[i],block_order[j]) )
				}
			})
}

if (TRUE) {
			sapply(1:nb, function(b) {
				# get closest neighbors in each direction
				closest <- apply(abs(centerB[,b] - centerB), 1, function(row) { sort(row, index.return=TRUE)$ix[2] })
				# get p closest neighbors
#				dists   <- sqrt( colSums( (centerB[,b] - centerB)^2 ) )
#				sorted  <- sort(dists, index.return=TRUE)
#				closest <- sorted$ix[1 + 1:factors$p]
#print(closest)

				larger <- unique(closest[which(closest > b)])
#print(larger)

				if (length(larger) > 0) {
					nmat <<- rbind(nmat, cbind(b, larger))
#print(larger); print(cbind(b,larger))
				}
			})
}

		} else {
			# get inverses x y
			invy <- lapply(1:nb, function(b) {
				idb1 <- which(B==b)

				N1 <- length(idb1)

				if (N1 > 0) {
					chol2inv(chol(data$Sigma[idb1,idb1])) %*% data$Yobs[idb1]
				} else {
					NA
				}
			})
		}
	})
	t2 <- proc.time()-t1
	ptime <- as.vector(t2[3])

	t1 <- proc.time()
	try({
		if (dep) {
			mse <- do.pred_blocks_dep(design, factors, data, centerB, nb, B, nmat)
		} else {
			mse <- do.pred_blocks_ind(design, factors, data, centerB, nb, B, invy)
		}
	})
	t2 <- proc.time()-t1

	list(ptime=ptime, time=as.vector(t2[3]), rmse=sqrt(mse) )
}

"pred.clust" <- function(design, factors, data, NperB, dep) {
	mse  <- NA

	t1 <- proc.time()
	try({
		# assign to blocks
		nb <- ceiling(data$n/NperB)
		#B  <- cutree(hc, k=nb)
		km <- kmeans(data$Xobs, nb, iter.max=100)
		B  <- km$cluster

		# get prediction blocks using block centroids
#		centerB <- do.call("cbind", lapply(1:nb, function(b) { colMeans(data$Xobs[B==b,]) }))
		centerB <- t(km$centers)
#print(centerB)

		if (dep) {
			# create neighbor structure

			nmat <- c()

if (FALSE) {
			block_order <- sample(1:nb, nb)
			lag <- 2
			sapply(1:nb, function(i) {
				for (j in i+1:lag) {
					if (j > nb) break;
					nmat <<- rbind( nmat, c(block_order[i],block_order[j]) )
				}
			})
}

			sapply(1:nb, function(b) {
				# get closest neighbors in each direction
				closest <- apply(abs(centerB[,b] - centerB), 1, function(row) {
					sort(row, index.return=TRUE)$ix[2]
				})
#print(closest)

				larger <- unique(closest[which(closest > b)])
#print(larger)

				if (length(larger) > 0) {
					nmat <<- rbind(nmat, cbind(b, larger))
#print(larger); print(cbind(b,larger))
				}
			})

		} else {
			# get inverses x y
			invy <- lapply(1:nb, function(b) {
				idb1 <- which(B==b)

				N1 <- length(idb1)

				if (N1 > 0) {
					chol2inv(chol(data$Sigma[idb1,idb1])) %*% data$Yobs[idb1]
				} else {
					NA
				}
			})
		}
	})
	t2 <- proc.time()-t1
	ptime <- as.vector(t2[3])

	t1 <- proc.time()
	try({
		if (dep) {
			mse <- do.pred_blocks_dep(design, factors, data, centerB, nb, B, nmat)
		} else {
			mse <- do.pred_blocks_ind(design, factors, data, centerB, nb, B, invy)
		}
	})
	t2 <- proc.time()-t1

	list(ptime=ptime, time=as.vector(t2[3]), rmse=sqrt(mse) )
}

# predict on subset
"pred.sub" <- function(design, factors, data, Nsub) {
	# create subset
	t1 <- proc.time()
	try({
#		p <- seq(0,1,length=Nsub)
#		s <- sample.int(data$n, Nsub)
		km <- kmeans(data$Xobs, Nsub, iter.max=100)
		s <- sapply(1:Nsub, function(b) {
			which(km$cluster == b)[1]
		})

		invy <- chol2inv(chol(data$Sigma[s,s])) %*% data$Yobs[s]
	})
	t2 <- proc.time()-t1
	ptime <- as.vector(t2[3])

	t1 <- proc.time()
	try({
		#preds <- gpuMM(data$Sigma[data$n+1:factors$Npred,1:data$n], predInvSigma) %*% data$Yobs
		preds <- data$Sigma[data$n+1:factors$Npred,s] %*% invy
		mse <- mean( (preds-data$Ypred)^2 )
	})
	t2 <- proc.time()-t1

	list(ptime=ptime, time=as.vector(t2[3]), rmse=sqrt(mse) )
}

if (FALSE) { # test sim design

# fixed design parameters
sim.design <- list(
	# number of replications
	Nreps=4
)

sim.factors <- expand.grid(
	# generate data for this number of inputs
	#p=c(5,10,15), p=c(10),
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
	#res <- lapply(1:1, function(i) {
	res <- lapply(which_exp, function(i) {
	  print(sim.factors[i,])
	  exp_res <- sim_exp_pred(sim.design, sim.factors[i,], i, which_part)
		save(exp_res, file=paste0("output/exp_pred_",i,".RData"))

print(head(exp_res))

	  exp_res
	})

	mres <- do.call(rbind, res)
}
