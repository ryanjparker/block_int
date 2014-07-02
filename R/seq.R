# code to analyze sequential design setting
library(fields)
library(lhs)
library(MASS)
library(multicore)
library(sensitivity)
library(lattice)

source("R/cov.R")
source("R/estimate.R")

# find low error predictor
"find_predictor" <- function(f, p, Niter=500) {
	# get initial observations
	n0 <- 100
	X0 <- as.matrix(maximinLHS(n0,p))
	y0 <- as.vector(branin(X0))
	mu_y <- mean(y0)
	sd_y <- sd(y0)

	starts <- c(log(11.5),log(3.0))

	# fit model to initial design
	data <- list(Yobs=(y0-mu_y)/sd_y, X=X0)
	fit0 <- hd.estimate(data, 1, rep(1,n0), matrix(1,nrow=1,ncol=1), p, starts, rep(FALSE,p), ce_cov, ce_partial, FALSE)
	#print(fit0)
	theta <- fit0$theta

	ncur <- length(y0)
	ycur <- y0
	Xcur <- X0

	if (Niter == 0) return( list(y=ycur, X=Xcur, theta=theta) )

	Nc    <- 100
	for (iter in 1:Niter) {
		# get candidate points
		Xc <- matrix(runif(Nc*p), nrow=Nc, ncol=p)

		# compute prediction variance at candidate points
		Sigma <- ce_cov(as.vector(theta), rbind(Xcur,Xc))
		sigmas <- sqrt( diag( 1 - Sigma[ncur+1:Nc,1:ncur] %*% chol2inv(chol(Sigma[1:ncur,1:ncur])) %*% Sigma[1:ncur,ncur+1:Nc] ) )
		to_add <- which.max(sigmas)
		if (sigmas[to_add] < 0.01) break

		cat("Adding",Xc[to_add,],"(",sigmas[to_add],")\n")

		ncur <- ncur+1
		ycur <- c(ycur,f(Xc[which.max(sigmas),]))
		Xcur <- rbind(Xcur, Xc[which.max(sigmas),])
		mu_y <- mean(ycur)
		sd_y <- sd(ycur)

		data <- list(Yobs=(ycur-mu_y)/sd_y, X=Xcur)
		fit <- hd.estimate(data, 1, rep(1,ncur), matrix(1,nrow=1,ncol=1), p, log(theta), rep(FALSE,p), ce_cov, ce_partial, FALSE)
		theta <- fit$theta
	}

	list(y=ycur, X=Xcur, theta=theta)
}

branin <- function(x) {
	x <- matrix(x,ncol=2)
	x1 <- -5 + x[,1]*15
	x2 <- x[,2]*15
#print(x1); print(x2)

	a <- 1
	b <- 5.1/(4*pi^2)
	c <- 5/pi
	r <- 6
	s <- 10
	t <- 1/(8*pi)

	a*(x2 - b*x1^2 + c*x1 - r)^2 + s*(1-t)*cos(x1)+s
}

"eval.pred" <- function(pred, x) {
	x <- matrix(x, ncol=length(pred$theta))
	Sigma <- ce_cov(as.vector(pred$theta), pred$X)
	iy <- chol2inv(chol(Sigma)) %*% pred$y
	ynew <- ce_full_pred_X(x, pred$X, iy, pred$theta);
}

#branin.pred.0 <- find_predictor(branin, 2, 0)
#branin.pred.500 <- find_predictor(branin, 2, 500)

if (TRUE) {
	# plot branin on (0,1) x (0,1) with prediction surface
	x1<-x2<-seq(0,1,len=100)
	z <- outer(x1,x2,function(x,y)branin(cbind(x,y)))

	z2 <- outer(x1,x2,function(x,y)eval.pred(branin.pred.0,cbind(x,y)))
	z3 <- outer(x1,x2,function(x,y)eval.pred(branin.pred.500,cbind(x,y)))

	e1 <- abs(z2-z)
	e2 <- abs(z3-z)

	# colors for z
	nrz <- nrow(z)
	ncz <- ncol(z)
	jet.colors <- colorRampPalette( c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000") )
	nbcol <- 100
	color <- jet.colors(nbcol)
	zfacet <- z[-1, -1] + z[-1, -ncz] + z[-nrz, -1] + z[-nrz, -ncz]
	facetcol <- cut(zfacet, nbcol)

	# colors for z2
	zfacet2 <- z2[-1, -1] + z2[-1, -ncz] + z2[-nrz, -1] + z2[-nrz, -ncz]
	facetcol2 <- cut(zfacet2, nbcol)

	# colors for error
	vals <- seq(min(c(e1,e2)), max(c(e1,e2)), length=nbcol)
	e1.facet <- e1[-1, -1] + e1[-1, -ncz] + e1[-nrz, -1] + e1[-nrz, -ncz]
	e1.facetcol <- cut(e1.facet, nbcol)
	#e1.facetcol <- sapply(1:length(e1.facet), function(i) { which.min( abs(e1[i]-vals) ) })
	#e1.facetcol <- apply(e1.facet, 2, function(row) { sapply(row, function(e) { which.min( abs(e-vals) ) }) })
	e2.facet <- e2[-1, -1] + e2[-1, -ncz] + e2[-nrz, -1] + e2[-nrz, -ncz]
	e2.facetcol <- cut(e2.facet, nbcol)
	#e2.facetcol <- sapply(1:length(e2.facet), function(i) { which.min( abs(e2[i]-vals) ) })

if (FALSE) {
	pdf("pdf/branin_0.pdf", height=9, width=3)
		par(mfrow=c(3,1))
		par(mar=c(1,1,1,1))
		#contour(x1, x2, z, vfont = c("sans serif", "plain"))
		#persp(x1,x2,z,theta=35,phi=30,expand=0.5,col="lightblue")
		#persp(x1,x2,z,theta=45,expand=0.5,col="lightblue",shade=1)

		persp(x1,x2,z,theta=55,phi=30,expand=1,col=color[facetcol], zlim=c(0,400), lwd=0.1, border=NA)
		persp(x1,x2,z2,theta=55,phi=30,expand=1,col=color[facetcol2], zlim=c(0,400), lwd=0.1, border=NA)
		persp(x1,x2,abs(z2-z),theta=55,phi=30,expand=1,zlim=c(0,10), col="darkgreen", lwd=0.1, border="darkgray")
	graphics.off()
}

	lab.x <- expression(x[1])
	lab.y <- expression(x^2)

	# actual function
	pdf("z.pdf")
		#par(mar=c(1,1,1,1))
		#persp(x1,x2,z,theta=55,phi=30,expand=1,col=color[facetcol], zlim=c(0,400), lwd=0.1, border=NA, xlab=lab.x, ylab=lab.y, zlab="f(x,y)")
		wireframe(z)
	graphics.off()

	# initial surface error
	elim <- c(0, max(c(e1,e2)) )
	pdf("pdf/branin_error_initial.pdf")
		par(mar=c(1,1,1,1))
		persp(x1,x2,e1,theta=55,phi=30,expand=1,col=color[e1.facetcol], zlim=elim, lwd=0.1, border=NA, xlab=lab.x, ylab=lab.y, zlab="Absolute Error")
	graphics.off()

	# initial surface error
	pdf("pdf/branin_error_seq.pdf")
		par(mar=c(1,1,1,1))
		persp(x1,x2,e2,theta=55,phi=30,expand=1,col=color[e2.facetcol], zlim=elim, lwd=0.1, border=NA, xlab=lab.x, ylab=lab.y, zlab="Absolute Error")
	graphics.off()

}
