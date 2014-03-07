# function to compute anisotropic squared exponential covariance
library(compiler)

#"ce_cov" <- cmpfun( function(theta, X) {
"ce_cov" <- function(theta, D) {
	Sigma <- Reduce('+',
		lapply(1:length(theta), function(k) { theta[k]*D[k,,] })
	)
	exp(-Sigma)
}#)

"ce_partial" <- function(e, theta, Sigma, D) {
	-D[e,,] * exp(theta[e]) * Sigma
}

"ce_full_pred" <- function(y, Nfit, Npred, Sigma) {
	as.vector( Sigma[Nfit+1:Npred,1:Nfit] %*% chol2inv(chol(Sigma[1:Nfit,1:Nfit])) %*% y )
}
