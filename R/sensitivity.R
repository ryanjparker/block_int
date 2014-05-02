# test sensitivity analysis

library(cubature)

fx <- function(x) { exp(-theta %*% (x-xo)^2) }
fx2 <- function(x) { fx(x)^2 }
fxe <- function() {
	prod(sapply(1:length(theta), function(k) {
		sqrt(pi)*
		(
			erf(sqrt(theta[k])*xo[k])-
			erf(sqrt(theta[k])*(xo[k]-1))
		)/(2*sqrt(theta[k]))
	}))
}

fxe2 <- function() {
	prod(sapply(1:length(theta), function(k) {
		sqrt(pi/2)*
		(
			erf(sqrt(2*theta[k])*xo[k])-
			erf(sqrt(2*theta[k])*(xo[k]-1))
		)/(2*sqrt(theta[k]))
	}))
}

theta <- rep(2, 2)
xo <- c(0.5,0.5)

t1<-proc.time()
X <- matrix(runif(2*10000),ncol=2)
mean(apply(X, 1, fx))
mean(apply(X, 1, fx2))
t2<-proc.time()-t1;print(t2)

t1<-proc.time()
adaptIntegrate(fx, c(0,0), c(1,1))
adaptIntegrate(fx2, c(0,0), c(1,1))
t2<-proc.time()-t1;print(t2)

t1<-proc.time()
fxe()
fxe2()
t2<-proc.time()-t1;print(t2)

