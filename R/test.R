# test estimation
library(fields)
library(MASS)

# generate data
s <- seq(0, 1, length=20)
S <- as.matrix( expand.grid(s, s) )
#S <- as.matrix( expand.grid(s, s, s, s) )

#n <- 200; S <- S[sample.int(nrow(S),n),]
n <- nrow(S)

D <- vector("list", length=ncol(S))
D2 <- vector("list", length=ncol(S))

sapply(1:ncol(S), function(column) {
	D[[column]]       <<- rdist(S[,column])
	diag(D[[column]]) <- 0
	D2[[column]]      <<- (D[[column]])^2
})

"compute_cov" <- function(Ntheta, theta, n, D, D2) {
	Sigma <- matrix(0, nrow=n, ncol=n)
	for (i in 1:Ntheta) {
		Sigma <- Sigma + -theta[i] * D2[[i]]
	}
	exp(Sigma)
}

Ntheta <- ncol(S)
theta <- rep(1/0.1, Ntheta)
cat("Computing Sigma...\n")
Sigma <- compute_cov(Ntheta, theta, n, D, D2)

# simulate data
set.seed(311)
y <- mvrnorm(1, rep(0,n), Sigma)

if (TRUE) {
	# plot data
	pdf("pdf/sim_data.pdf")
		image.plot(matrix(y,length(s),length(s)))
	graphics.off()
}

# cluster 1-Sigma
hc  <- hclust( as.dist( 1-Sigma ) )
hcm <- cutree(hc, k=round(nrow(Sigma)/20) )
