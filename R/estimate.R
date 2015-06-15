library(parallel)

# function to run algorithm for fitting a high domensiontal block composite model
"hd.estimate" <- function(
	data,
	nblocks, B, neighbors,
	R,
	theta, theta_fixed,
	f.cc, f.cp,
	verbose, tol=1e-8, maxIter=50, parallel=FALSE
) {
	# f.cc: function to compute covariance
	# f.cp: function to compute partial derivatives
	# nblocks: number of blocks
	# B: block memberships
	# neighbors: 2 column matrix listing unique neighbors of blocks
	# R: number of covariance function parameters
	# theta: initial values for covariance parameters
	# theta_fixed: TRUE/FALSE indicating if parameter is fixed

	# verbose: print messages?
	# tol: error tolerance for identifying convergence
	# maxIter: maximum number of Fisher scoring iterations
	# parallel: run operations in parallel?

	Ntheta <- R

	nH <- R*(R+1)/2

	seq.R <- 1:R
	seq.R2 <- 1:(R^2)
	seq.RH <- 1:nH

	# how many parameters are not fixed?
	Nnot_fixed <- length(theta_fixed)-sum(theta_fixed)
	which.fixed <- which(theta_fixed==TRUE)
	which.not_fixed <- which(theta_fixed==FALSE)

	# transformation functions
	t_theta <- function(theta) { exp(theta) }

	# function to update theta with fisher scoring
	u <- rep(0, R)
	W <- vector("list", R)
	H <- rep(0, R*(R+1)/2)
	FI <- matrix(0, nrow=R, ncol=R)

	if (parallel) {
		par_lapply <- mclapply
	} else {
		par_lapply <- lapply
	}

	"update_theta" <- function(theta) {
		if (Nnot_fixed == 0) {   # don't do any updating
			return(theta)
		}

		u[seq.R]   <<- 0
		H[seq.RH]  <<- 0
		FI[seq.R2] <<- 0

		#apply(neighbors, 1, function(row) {
		bres <- par_lapply(1:nrow(neighbors), function(irow) {
			u <- rep(0, R)
			W <- vector("list", R)
			H <- rep(0, R*(R+1)/2)

			row <- neighbors[irow,]

			in.pair <- which(B==row[1] | B==row[2])

			#Sigma <- f.cc(t_theta(theta), data$D[,in.pair,in.pair])
			Sigma <- f.cc(t_theta(theta), data$X[in.pair,])
			invSigma <- chol2inv(chol(Sigma))

			q <- invSigma %*% data$Yobs[in.pair]

			# compute the Ws
			for (r in 1:R) {
				if (theta_fixed[r]) { next; }

				#partial <- f.cp(r, theta, Sigma, data$D[,in.pair,in.pair])
				partial <- f.cp(r, theta, Sigma, data$X[in.pair,])
				W[[r]] <- invSigma %*% partial
				u[r]   <- u[r] -0.5 * sum( diag(W[[r]]) ) + 0.5 * t(q) %*% partial %*% q
			}

			# compute the hessian H
			# this is stored in a compact form in H and then expanded into FI (fisher info)
			index <- 1
			sapply(seq.R, function(r) {
				sapply(r:R, function(s) {
					if (!theta_fixed[r] & !theta_fixed[s]) {
						#H[index] <<- H[index] + 0.5 * sum(diag( W[[r]] %*% W[[s]] ))
						H[index] <<- H[index] + 0.5 * sum_diag_mm(W[[r]], W[[s]]) #sum(diag( W[[r]] %*% W[[s]] ))
					}
					index <<- index+1
				})
			})

			list(u=u, H=H)
		})

		# gather results
		sapply(1:length(bres), function(i) {
			u <<- u+bres[[i]]$u

			index <- 1
			sapply(1:R, function(r) {
				sapply(r:R, function(s) {
					H[index] <<- H[index] + bres[[i]]$H[index]

					index <<- index+1
				})
			})
		})

		index <- 1
		sapply(seq.R, function(r) {
			sapply(r:R, function(s) {
				if (!theta_fixed[r] & !theta_fixed[s]) {
					FI[r,s] <<- H[index]
					if (r != s) {
						FI[s,r] <<- H[index]
					}
				}
				index <<- index+1
			})
		})

		tryCatch({
			# fix params if they go to 0
			cholFI <- chol(FI[which.not_fixed,which.not_fixed])
		}, error = function(e) {
			cat("Unable to invert FI; trying to fix smallest param\n")
			smallest <- which.min(theta[which.not_fixed])
			which.fixed <<- c(which.fixed, which.not_fixed[smallest])
			which.not_fixed <<- which.not_fixed[-smallest]
			cholFI <<- chol(FI[which.not_fixed,which.not_fixed])
		}) #, finally = print("Hello"))

		# cap change at 1 since we're on log scale
		change <- chol2inv(cholFI) %*% u[which.not_fixed]
		change <- ifelse(abs(change) >= 1, sign(change)*1, change)

		theta[which.not_fixed] <- theta[which.not_fixed] + change

		theta
	}

	# compute log likelihood
	"loglik" <- function(theta) {
		ll <- sum( apply(neighbors, 1, function(row) {
			in.pair <- which(B==row[1] | B==row[2])

			#Sigma     <- f.cc(t_theta(theta), data$D[,in.pair,in.pair])
			Sigma     <- f.cc(t_theta(theta), data$X[in.pair,])
			cholSigma <- chol(Sigma)
			invSigma  <- chol2inv(cholSigma)

			-sum(log(diag(cholSigma))) -0.5 * t(data$Yobs[in.pair]) %*% invSigma %*% data$Yobs[in.pair]
		}) )
	}

	# estimate params
	names.show <- c("theta", rep("", Ntheta-1))
	names.show <- c(names.show, "log lik")

	ll <- loglik(theta)

	# save theta and log lik at each iteration
	iters_theta <- t_theta(theta)
	iters_ll    <- ll

	for (iter in 1:maxIter) {
		prev.theta <- theta

		# update theta
		theta <- update_theta(prev.theta)

		# get log likelihood
		ll <- loglik(theta)

		# save values at each iteration
		iters_theta <- rbind( iters_theta, t_theta(theta) )
		iters_ll    <- c( iters_ll, ll )

		if (verbose) {
			show <- round(c(t_theta(theta), ll),2)
			names(show) <- names.show
			cat("iter ",iter,":\n",sep="")
			print( show )
		}

		# have we converged?
		if ( abs(ll - iters_ll[iter])/abs(0.1 + ll) <= tol ) {
			if (verbose) {
				cat("Converged at iteration",iter,"\n")
			}

			break
		}
	}

	convergence <- TRUE
	if (iter == maxIter) {
		warning("Possible issues with convergence: maximum number of iterations reached.")
		convergence <- FALSE
	}

	# transform theta
	theta <- t_theta(theta)

	# return estimates and standard errors
	list(
		convergence=convergence, nIter=iter,
		iters_theta=iters_theta, iters_ll=iters_ll,
		theta=theta, ll=ll
	)
}



if (FALSE) {
	# compute log likelihood
	"loglik" <- function(beta, theta) {
		ll <- sum( apply(neighbors, 1, function(row) {
			in.pair <- B==row[1] | B==row[2]
			n.pair <- sum(in.pair)

			Sigma <- compute_cov(cov, t_theta(theta), D[in.pair,in.pair])
			cholSigma <- chol(Sigma)
			invSigma <- chol2inv(cholSigma)

			Xb <- X[in.pair,] %*% beta
			ymXb <- y[in.pair]-Xb

			-sum(log(diag(cholSigma))) -0.5 * t(ymXb) %*% invSigma %*% ymXb
		}) )
	}
}

# full model
"eval_full" <- function(x, design, factors, data) {
cat("This is broken.\n")
return(NA)
	# compute covariance
	#Sigma <- factors$sigma2 * ce_cov(exp(x), data$D[,1:data$n,1:data$n])
	Sigma <- factors$sigma2 * ce_cov(exp(x), data$X[1:data$n,])
	cholSigma <- chol(Sigma)
	invSigma <- chol2inv(cholSigma)

	nll <- sum(log(diag(cholSigma))) +0.5 * t(data$Yobs) %*% invSigma %*% data$Yobs
print(nll)
	grad <- rep(0, factors$p)
	q <- invSigma %*% data$Yobs
	W <- vector("list", factors$p)
	lapply(1:factors$p, function(t) {
#ACK
		partial <- -exp(x[t]) * data$D[t,1:data$n,1:data$n] * Sigma
		W[[t]] <<- invSigma %*% partial
		grad[t] <<- 0.5*sum(diag( W[[t]] )) -0.5 * t(q) %*% partial %*% q
	})

	hess <- matrix(0, nrow=factors$p, ncol=factors$p)
	sapply(1:factors$p, function(p1) {
		sapply(p1:factors$p, function(p2) {
			hess[p1,p2] <<- 0.5*sum(diag( W[[p1]] %*% W[[p2]] ))
			if (p1 != p2) hess[p2,p1] <<- hess[p1,p2]
		})
	})

	#list( "objective"=nll, "gradient"=grad)

	attr(nll, "gradient") <- grad
	#attr(nll, "hessian")  <- chol2inv(chol(hess))
	attr(nll, "hessian")  <- hess

	nll
}
