# function to compute anisotropic squared exponential covariance
library(compiler)
library(Rcpp)
sourceCpp("src/cov.cpp")
#print(ce_full_pred_X_cov)
#print(str(ce_full_pred_X_cov(matrix(runif(100*2),nrow=100), matrix(runif(100*2),nrow=100), c(0.3,0.3))))

"sum_diag_mm" <- function(A, B) {
	sum_diag_mm_Rcpp(as.matrix(A), as.matrix(B))
}

#"ce_cov" <- cmpfun( function(theta, X) {
#"ce_cov" <- function(theta, D) {
"ce_cov" <- function(theta, X) {
#	Sigma <- Reduce('+',
#		lapply(1:length(theta), function(k) { theta[k]*D[k,,] })
#	)
#	exp(-Sigma)
	ce_cov_Rcpp(as.vector(theta), as.matrix(X)) + diag(nrow(X))*0.001
}#)

"ce_partial" <- function(e, theta, Sigma, X) {
#	-D[e,,] * exp(theta[e]) * Sigma
	ce_partial_Rcpp(as.integer(e), as.vector(theta), as.matrix(Sigma), as.matrix(X))
}

"ce_full_pred" <- function(y, Nfit, Npred, Sigma) {
	as.vector( Sigma[Nfit+1:Npred,1:Nfit] %*% chol2inv(chol(Sigma[1:Nfit,1:Nfit])) %*% y )
}

"ce_full_pred_X" <- function(X, Xobs, iy, theta) {
	# compute covariance between prediction and observation locations

#	# compute covariance between X and Xobs
#	Sigma <- t(exp(-apply(X, 1, function(row) {
#		apply(Xobs, 1, function(obs) {
#			theta %*% (row-obs)^2
#		})
#	})))

	Sigma <- ce_full_pred_X_cov(as.matrix(X), as.matrix(Xobs), as.vector(theta))

	as.vector(Sigma %*% iy)
}

"ce_local_pred" <- function(y, Nfit, Npred, Sigma, Nlocal=100) {
	y_0 <- rep(0, Npred)

	for (i in 1:Npred) {
		# get closest
		close <- sort(Sigma[Nfit+i,1:Nfit], index.return=TRUE, decreasing=TRUE)$ix[1:Nlocal]
#sorted <- sort(Sigma[Nfit+i,1:Nfit], index.return=TRUE, decreasing=TRUE)
#print(head(sorted$x)); print(tail(sorted$x))
#print(summary(Sigma[Nfit+i,1:Nfit][close]))
#print(summary(Sigma[Nfit+i,1:Nfit][-close]))
		#y_0[i] <- Sigma[Nfit+i,close] %*% chol2inv(chol(Sigma[close,close])) %*% y[close]
		y_0[i] <- Sigma[Nfit+i,1:Nfit][close] %*% chol2inv(chol(Sigma[1:Nfit,1:Nfit][close,close])) %*% y[1:Nfit][close]
	}

	y_0
}

"ce_block_pred" <- function(theta, y, Nfit, Npred, B, newB, D, neighbors) {
	# predict when we have blocks

	# predict for each unique block
	y_0 <- rep(0, Npred)

	# unique blocks where we have points to predict
	uB <- unique(newB)

	for (b in uB) {
		# neighbors of this block
		close <- unique( c(neighbors[neighbors[,1]==b,2], neighbors[neighbors[,2]==b,1]) )

		# what new points are in this block?
		in.new   <- which(newB == b)
		n.in.new <- length(in.new)

		A_0 <- matrix(0, nrow=n.in.new, ncol=n.in.new)
		b_0 <- rep(0, n.in.new)

		for (n1 in close) {
			# observations in this block pair
			in.obs <- which(B==b | B==n1)
			n.in.obs <- length(in.obs)

			# take points so that we have Y = (Y_{new}, Y_{obs})
			in.pair <- c(in.new+Nfit, in.obs)

			# covariance for observations and new
			Sigma.pair    <- ce_cov(theta, D[,in.pair,in.pair])
			invSigma.pair <- chol2inv(chol(Sigma.pair))

			A_0 <- A_0 + invSigma.pair[1:n.in.new,1:n.in.new]
			b_0 <- b_0 + invSigma.pair[1:n.in.new,n.in.new+1:n.in.obs] %*% y[in.obs]
		}

		# complete predictions
		b_0 <- -b_0
		y_0[in.new] <- chol2inv(chol(A_0)) %*% b_0
	}

	y_0
}

