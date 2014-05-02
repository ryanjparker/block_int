# generate random function with b-spline
library(fda)

x <- seq(0,1,len=100)
n <- length(x)
nbasis <- 4

Bspline.basis <- create.bspline.basis(c(0,1),norder=4,nbasis=nbasis)
weights <- getbasismatrix(x, Bspline.basis, nderiv=0)
nweights <- nrow(weights)

# 1D
set.seed(1983)
coefs <- rnorm(nbasis)
y <- apply(weights, 1, function(w) { crossprod(coefs,w) })
pdf("pdf/1d.pdf")
	plot(x, y, main="1D", type="l")
graphics.off()

# 2D
set.seed(1983)
coefs <- matrix(rnorm(nbasis^2), nrow=nbasis, ncol=nbasis)
y <- matrix(NA, nrow=nweights, ncol=nweights)
for (i in 1:nweights) {
	for (j in 1:nweights) {
		y[i,j] <- sapply(1:nbasis, function(b1) {
			sapply(1:nbasis, function(b2) {
			})
		})
	}
}
done
pdf("pdf/2d.pdf")
	#plot(x, y, main="2D", type="l")
graphics.off()
