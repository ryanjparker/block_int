source("R/cov.R")
source("R/estimate.R")
library(fields)
library(lhs)
library(spacious)
library(sensitivity)

options(cores=4); options(mc.cores=4)
set.seed(1983)

if (!exists("d.train")) {
	# read in sarcos data
	load("data/sarcos.RData")

	d.all <- sarcos$sarcos.inv

	d.train <- sarcos$sarcos.inv
	#d.train <- d.train[sample(nrow(d.train), 10000),] # sample the data
	d.test  <- sarcos$sarcos.inv.test

if (TRUE) {
	cat("Removing test samples from training\n")

	in.both <- unlist( apply(d.test, 1, function(row) {
		which(colSums(t(d.train) - row)==0)
	}) )

	d.train <- d.train[-in.both,]
} else {
	cat("Test samples *not* removed from training!\n")
}

	nX <- 21

	col.X <- 1:nX
	col.Y <- 22

	# scale inputs to be in [0,1]
	X.max <- apply(d.train[,col.X], 2, max)
	X.min <- apply(d.train[,col.X], 2, min)

	sapply(col.X, function(col) {
		d.train[,col] <<- (d.train[,col] - X.min[col]) / (X.max[col]-X.min[col])
		d.test[,col] <<- (d.test[,col] - X.min[col]) / (X.max[col]-X.min[col])
	})

if (FALSE) {
	# scale Y to have mean 0 and SD 1
	Y.mean <- mean(d.train[,col.Y])
	Y.sd   <- sd(d.train[,col.Y])
	d.train[,col.Y] <- (d.train[,col.Y] - Y.mean)/Y.sd
	d.test[,col.Y]  <- (d.test[,col.Y] - Y.mean)/Y.sd
} else {
	cat("Removing quadratic trend.\n")
	Xmat <- as.data.frame( cbind(d.train[,col.Y], d.train[,col.X], d.train[,col.X]^2) ); names(Xmat) <- paste0("v", 1:ncol(Xmat))
	fit.lm <- lm(v1 ~ ., data=Xmat)
	d.train[,col.Y] <- rstandard(fit.lm)  # get standardized residuals
	Y.sd <- summary(fit.lm)$sigma

	Xmat <- as.data.frame( cbind(d.test[,col.Y], d.test[,col.X], d.test[,col.X]^2) ); names(Xmat) <- paste0("v", 1:ncol(Xmat))
	d.test[,col.Y] <- (d.test[,col.Y]-predict(fit.lm, newdata=Xmat))/Y.sd
}

	n.train <- nrow(d.train)
	n.test <- nrow(d.test)

} else {
	# data loaded
	cat("Data previously loaded and cleaned.\n")
}

if (FALSE) { # compare threads vs no threads

# assign to blocks
NperB <- 100
nb <- ceiling(n.train/NperB)
cat("# per block = ",NperB," => # of blocks = ",nb,"\n",sep="")

# cluster on locations
km <- kmeans(d.train[,col.X], nb, iter.max=100)
B  <- km$cluster

# fit model
t1 <- proc.time()
fit1 <- hd.estimate(data=list(X=d.train[,col.X], Yobs=d.train[,col.Y]), nb, B, cbind(1:nb,1:nb), nX, rep(log(1),nX), rep(FALSE,nX), ce_cov, ce_partial, verbose=TRUE)
ftime <- proc.time()-t1
print(ftime)

t1 <- proc.time()
fit2 <- hd.estimate(data=list(X=d.train[,col.X], Yobs=d.train[,col.Y]), nb, B, cbind(1:nb,1:nb), nX, rep(log(1),nX), rep(FALSE,nX), ce_cov, ce_partial, verbose=TRUE, parallel=TRUE)
ftime <- proc.time()-t1
print(ftime)

} # end threads vs no threads

if (FALSE) { # fit models

if (TRUE) {
# IND C / 100
	# assign to blocks
	NperB <- 100
	nb    <- ceiling(n.train/NperB)
	cat("# per block = ",NperB," => # of blocks = ",nb,"\n",sep="")

	# cluster on locations
	km <- kmeans(d.train[,col.X], nb, iter.max=100)
	B  <- km$cluster

	# fit the model
	t1 <- proc.time()
	fit.IC100 <- hd.estimate(data=list(X=d.train[,col.X], Yobs=d.train[,col.Y]), nb, B, cbind(1:nb,1:nb), nX, rep(log(1),nX), rep(FALSE,nX), ce_cov, ce_partial, verbose=TRUE, parallel=TRUE)
	time.IC100 <- proc.time()-t1
	print(time.IC100)
}

if (TRUE) {
# IND C / 50
	# assign to blocks
	NperB <- 50
	nb    <- ceiling(n.train/NperB)
	cat("# per block = ",NperB," => # of blocks = ",nb,"\n",sep="")

	# cluster on locations
	km <- kmeans(d.train[,col.X], nb, iter.max=100)
	B  <- km$cluster

	# fit the model
	t1 <- proc.time()
	fit.IC50 <- hd.estimate(data=list(X=d.train[,col.X], Yobs=d.train[,col.Y]), nb, B, cbind(1:nb,1:nb), nX, rep(log(1),nX), rep(FALSE,nX), ce_cov, ce_partial, verbose=TRUE, parallel=TRUE)
	time.IC50 <- proc.time()-t1
	print(time.IC50)
}

if (TRUE) {
# DEP C / 50
	# assign to blocks
	NperB <- 50
	nb    <- ceiling(n.train/NperB)
	cat("# per block = ",NperB," => # of blocks = ",nb,"\n",sep="")

#r <- sort(tapply(B,B,length)); print(summary(r));done
	# cluster on locations
	km <- kmeans(d.train[,col.X], nb, iter.max=100)
	B  <- km$cluster

	# order blocks based on (1) distance from 0, then (2) distance from each other
	to0 <- rdist(matrix(0, nrow=1, ncol=nX), km$centers)
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

	# fit the model
	t1 <- proc.time()
	fit.DC50 <- hd.estimate(data=list(X=d.train[,col.X], Yobs=d.train[,col.Y]), nb, B, nmat, nX, rep(log(1),nX), rep(FALSE,nX), ce_cov, ce_partial, verbose=TRUE, parallel=TRUE)
	time.DC50 <- proc.time()-t1
	print(time.DC50)
}

	save(fit.DC50,time.DC50,fit.IC50,time.IC50,fit.IC100,time.IC100,file="output/ex/sarcos.RData")

} # end fit models

"compute_preds" <- function(fit, train, test) {
	res <- list(local=list(), subset=list())

	Ntrain <- nrow(train)
	Ntest  <- nrow(test)

if (TRUE) {
	# local
	for (Nlocal in c(10,25,100)) { #,500,1000)) {
		# closest by covariance
		#preds <- sapply(1:Ntest, function(px) {
		preds <- unlist( mclapply(1:Ntest, function(px) {
			localSigma <- ce_cov_single(fit$theta, train[,col.X], matrix(test[px,col.X],nrow=1))
#if (max(localSigma) < 1) { print( c(px, max(localSigma), which.max(localSigma)) ); done }
#print(round( c(sum(localSigma > 0.75), sum(localSigma < 0.75 & localSigma > 0.5),sum(localSigma < 0.5 & localSigma > 0.25),sum(localSigma < 0.01)), 4))

			dists <- 1-localSigma
			locs  <- sort(dists,decreasing=FALSE,index.return=TRUE)$ix[1:Nlocal]
#print(head(round(as.vector(dists[locs]),3),20)); print(head(locs)); print(head(localSigma[locs])); #done

			predSigma <- ce_cov(fit$theta, rbind(train[locs,col.X], matrix(test[px,col.X],nrow=1))) - diag(Nlocal + 1)*0.001
			#invy <- gpuChol2Inv(predSigma[1:Nlocal,1:Nlocal]) %*% train[s,col.Y]
			invy <- chol2inv(chol(predSigma[1:Nlocal,1:Nlocal])) %*% train[locs,col.Y]

			#preds <- gpuMM(predSigma[Nlocal+1:Ntest,1:Nlocal], invy)
			preds <- predSigma[Nlocal+1,1:Nlocal] %*% invy

#print(c(preds,test[px,col.Y]))
			preds
		}, mc.cores=4) )

#print(c(preds,d.test[3,col.Y]))

    rmse <- sqrt(mean( (preds-test[,col.Y])^2 ))
print(rmse)
		res[["local"]][[paste(Nlocal)]] <- rmse
	}
}

	# subset
	for (Nsub in c(1000,5000,10000)) {
		# create subset
		t1 <- proc.time()
		try({
			km <- kmeans(train[,col.X], Nsub, iter.max=100)
			s <- sapply(1:Nsub, function(b) { which(km$cluster == b)[1] })

			# compute covariance
			predSigma <- ce_cov(fit$theta, rbind(train[s,col.X], test[,col.X])) - diag(Nsub + Ntest)*0.001

			invy <- gpuChol2Inv(predSigma[1:Nsub,1:Nsub]) %*% train[s,col.Y]
			#invy <- chol2inv(chol(predSigma[1:Nsub,1:Nsub])) %*% d.train[s,col.Y]
		})
		t2 <- proc.time()-t1
		ptime <- as.vector(t2[3])

		t1 <- proc.time()
		try({
			preds <- gpuMM(predSigma[Nsub+1:Ntest,1:Nsub], invy)
			#preds <- predSigma[Nsub+1:n.test,1:Nsub] %*% invy
			rmse <- sqrt(mean( (preds-test[,col.Y])^2 ))
		})
		t2 <- proc.time()-t1

		print(rmse)
		res[["subset"]][[paste(Nsub)]] <- rmse
	}

	res
} # end compute_preds()

if (FALSE) { # plot theta
	max.theta <- 1 #<- max(c(fit.DC50$theta, fit.IC50$theta, fit.IC100$theta))
	pdf("pdf/ex/theta.pdf", height=2.5)
		par(mar = c(3,3.1,0.1,0.1)) # b,l,t,r
		plot(1:21, (fit.DC50$theta)/max.theta, pch=1, xlab="", ylab="", xaxt="n")
		points(1:21, (fit.IC50$theta)/max.theta, pch=2)
		points(1:21, (fit.IC100$theta)/max.theta, pch=3)

		abline(v=7.5, lty=2)
		abline(v=14.5, lty=2)

		axis(side = 1, at = c(2:7,9:14,16:21), labels=NA)
		axis(side = 1, at = c(1,8,15))

		legend("topright", c("100 obs/block (Ind)", "50 obs/block (Ind)", "50 obs/block (Dep)"), pch=c(3,2,1), cex=0.75, inset=0.02)
		mtext(expression(widehat(theta)), side=2, line=2)
		mtext("Input", side=1, line=1.5)

	graphics.off()
}

if (FALSE) { # make predictions
	load("output/ex/sarcos.RData")

#compute_preds(fit.DC50, d.train, d.test); done
preds.DC50 <- compute_preds(fit.DC50, d.train, d.test)
preds.IC50 <- compute_preds(fit.IC50, d.train, d.test)
preds.IC100 <- compute_preds(fit.IC100, d.train, d.test)

save(preds.DC50,preds.IC50,preds.IC100,file="output/ex/sarcos_preds.RData")
done

} # end predictions


"do_sub_preds" <- function(x, Xobs, iy, theta) {
#str(x); str(Xobs); str(iy); str(theta)
	Xmat <- as.data.frame( cbind(x, x^2) ); names(Xmat) <- paste0("v", 1+1:ncol(Xmat))
	pred <- predict(fit.lm, newdata=Xmat) + Y.sd*ce_full_pred_X(x, Xobs, iy, theta)
}

"compute_sens" <- function(fit) {
set.seed(1983)
	Nsens <- 1000
	Nsub  <- 10000

	# compute sensitivity indices
	km <- kmeans(d.train[,col.X], Nsub, iter.max=100)
	s <- sapply(1:Nsub, function(b) { which(km$cluster == b)[1] })

	# compute covariance
	predSigma <- ce_cov(fit$theta, d.train[s,col.X]) - diag(Nsub)*0.001

	iy <- gpuChol2Inv(predSigma[1:Nsub,1:Nsub]) %*% d.train[s,col.Y]
	#iy <- chol2inv(chol(predSigma[1:Nsub,1:Nsub])) %*% d.train[s,col.Y]

	Nr <- 10
	res <- lapply(1:Nr, function(r) {
		X1_s <- data.frame(randomLHS(Nsens, nX))
		X2_s <- data.frame(randomLHS(Nsens, nX))

		#Si <- sobol2002(model = ce_full_pred_X, X1_s, X2_s, nboot = 0, Xobs=d.train[s,col.X], iy=iy, theta=fit$theta)
		Si <- sobol2002(model = do_sub_preds, X1_s, X2_s, nboot = 0, Xobs=d.train[s,col.X], iy=iy, theta=fit$theta)
print(r)
		list(S=Si$S[,1], T=Si$T[,1])
	})

	S <- colMeans( do.call("rbind", lapply(1:Nr, function(r) { res[[r]]$S })) )
	T <- colMeans( do.call("rbind", lapply(1:Nr, function(r) { res[[r]]$T })) )
print( round(T,4) )

	#list(S=Si$S[,1], T=Si$T[,1])
	list(S=S, T=T)
}

if (FALSE) { # perform sensitivity analysis
	load("output/ex/sarcos.RData")

sens.DC50 <- compute_sens(fit.DC50)
sens.IC50 <- compute_sens(fit.IC50)
sens.IC100 <- compute_sens(fit.IC100)
save(sens.DC50,sens.IC50,sens.IC100,file="output/ex/sarcos_sens.RData")
done

}

if (FALSE) { # plot total sensitivity
	max.s <- max(c(sens.DC50$T, sens.IC50$T, sens.IC100$T))
	pdf("pdf/ex/sens.pdf", height=2.5)
		par(mar = c(3,3.1,0.1,0.1)) # b,l,t,r
		s <- sens.DC50$T; s[s < 0] <- 0; s[s > 1] <- 1
		plot(1:21, s, pch=1, xlab="", ylab="", xaxt="n", ylim=c(0,1))
		s <- sens.IC50$T; s[s < 0] <- 0; s[s > 1] <- 1
		points(1:21, s, pch=2)
		s <- sens.IC100$T; s[s < 0] <- 0; s[s > 1] <- 1
		points(1:21, s, pch=3)

		abline(v=7.5, lty=2)
		abline(v=14.5, lty=2)

		axis(side = 1, at = c(2:7,9:14,16:21), labels=NA)
		axis(side = 1, at = c(1,8,15))

		#legend("top", c("100 obs/block (Ind)", "50 obs/block (Ind)", "50 obs/block (Dep)"), pch=c(3,2,1), cex=0.5, inset=c(0,-0.2))
		mtext("Total Sensitivity", side=2, line=2)
		mtext("Input", side=1, line=1.5)

	graphics.off()
}

if (TRUE) { # plot most important indices

if (TRUE) {
	y   <- seq(min(d.all[,col.Y]), max(d.all[,col.Y])+0.01, length=20)
	x15 <- seq(min(d.all[,col.X[15]]), max(d.all[,col.X[15]])+0.01, length=20)
	x18 <- seq(min(d.all[,col.X[18]]), max(d.all[,col.X[18]])+0.01, length=20)
	hmat1 <- matrix(NA, nrow=20, ncol=20)
	hmat2 <- matrix(NA, nrow=20, ncol=20)
	for (i2 in 2:length(y)) {
		sub1 <- d.all[,col.Y] >= y[i2-1] & d.all[,col.Y] < y[i2]

		for (i1 in 2:length(x15)) {
			sub2 <- d.all[,col.X[15]] >= x15[i1-1] & d.all[,col.X[15]] < x15[i1]
			if (sum(sub1&sub2) > 0) hmat1[i1,i2] <- sum(sub1&sub2) #length(d.all[ sub1&sub2, col.Y ])
		}

		for (i1 in 2:length(x18)) {
			sub2 <- d.all[,col.X[18]] >= x18[i1-1] & d.all[,col.X[18]] < x18[i1]
			if (sum(sub1&sub2) > 0) hmat2[i1,i2] <- sum(sub1&sub2) #length(d.all[ sub1&sub2, col.Y ])
		}
	}

	hmat1 <- log(hmat1)
	hmat2 <- log(hmat2)

}

if (FALSE) {
	x15 <- seq(min(d.all[,col.X[15]]), max(d.all[,col.X[15]])+0.01, length=20)
	x18 <- seq(min(d.all[,col.X[18]]), max(d.all[,col.X[18]])+0.01, length=20)
	tmat <- matrix(NA, nrow=20, ncol=20)
	for (i1 in 2:length(x15)) {
		for (i2 in 2:length(x18)) {
			sub1 <- d.all[,col.X[15]] >= x15[i1-1] & d.all[,col.X[15]] < x15[i1]
			sub2 <- d.all[,col.X[18]] >= x18[i2-1] & d.all[,col.X[18]] < x18[i2]
			if (sum(sub1&sub2) > 0) tmat[i1,i2] <- mean(d.all[ sub1&sub2, col.Y ])
		}
	}
}

	# b,l,t,r
	mi <- 0.132
	mb <- c(2*mi/(7.2/3),0.8,0.528/2.25,1-1.5*mi/2.25)
	ms <- c(0.82,0.85,0.528/2.25,1-1.5*mi/2.25)

	library(scales)
	pdf("pdf/ex/most.pdf", height=2.25, width=7.2)
		par(mfrow=c(1,3))

#par(fig=c(0.00,0.30,0,1)) #, new=TRUE )
		#par(mar = c(4,4.1,0.1,0.1)) # b,l,t,r
#print( par("mai")/par("mar") )
		#plot(d.all[,col.X[15]], d.all[,col.Y], pch=16, cex=0.5, col=alpha("black", 0.5), xlab="", ylab="")
		#mtext(expression(X[15]), side=1, line=2.5, cex=0.75)
		#mtext(expression(X[18]), side=2, line=2, cex=0.75)
		par(mar = c(4,2.0,1.5,3.1)) # b,l,t,r
		image.plot(x=x15, y=y, z=hmat1, legend.width=0.5, legend.mar=0, legend.shrink=0.25,
			main=expression(paste("Heatmap of ",x[15]," vs ",y)), xlab="", ylab="", #xlab=expression(X[15]), ylab=expression(Y),
			bigplot=mb, smallplot=ms)

if (TRUE) {
#par(fig=c(0.30,0.60,0,1), new=TRUE )
		#par(mar = c(4,3.0,0.1,1)) # b,l,t,r
		#plot(d.all[,col.X[18]], d.all[,col.Y], pch=16, cex=0.5, col=alpha("black", 0.5), xlab="", ylab="")
		#mtext(expression(X[15]), side=1, line=2.5, cex=0.75)
		#mtext(expression(X[18]), side=2, line=2, cex=0.75)
		par(mar = c(4,2.0,1.5,3.1)) # b,l,t,r
		image.plot(x=x18, y=y, z=hmat2, legend.width=0.5, legend.mar=0, legend.shrink=0.25,
			main=expression(paste("Heatmap of ",x[18]," vs ",y)), xlab="", ylab="", #xlab=expression(X[15]), ylab=expression(Y),
			bigplot=mb, smallplot=ms)
}

#par(fig=c(0.60,1.00,0,1), new=TRUE )
		par(mar = c(4,2.0,1.5,3.1)) # b,l,t,r
		image.plot(x=x15, y=x18, z=tmat, legend.width=0.5, legend.mar=0, legend.shrink=0.25,
			main=expression(paste("Mean of ",y," for ",x[15]," vs ",x[18])), xlab="", ylab="", #xlab=expression(X[15]), ylab=expression(Y),
			bigplot=mb, smallplot=ms)
		#image(x=x15, y=x18, z=tmat, xlab="", ylab="")
		#title("")
		#mtext(expression(X[15]), side=1, line=2.5, cex=0.75)
		#mtext(expression(X[18]), side=2, line=2, cex=0.75)

#par(fig=c(0.90,1.00,0,1), new=TRUE )
#	image.plot(x=x15, y=x18, z=tmat, legend.only=TRUE, horizontal=FALSE, legend.width=5, smallplot=c(0.20,0.30,0.20,0.80))

	graphics.off()
}

if (FALSE) { # estimate prediction error with CV

"cv_error" <- function(k) {
	# subset data
	subsets <- split(sample(1:n.train),ceiling( (1:n.train)/(n.train/k) ))

	res <- sapply(1:k, function(fold) {
		# setup data for this fold
		d.other <- d.train[-subsets[[fold]],]
		n.other <- nrow(d.other)
		d.fold  <- d.train[subsets[[fold]],]
		n.fold  <- nrow(d.fold)

		# fit model
# IND C / 100
		# assign to blocks
		NperB <- 100
		nb    <- ceiling(n.fold/NperB)
		cat("# per block = ",NperB," => # of blocks = ",nb,"\n",sep="")

		# cluster on locations
		km <- kmeans(d.other[,col.X], nb, iter.max=100)
		B  <- km$cluster

		# fit the model
		t1   <- proc.time()
		fit  <- hd.estimate(data=list(X=d.other[,col.X], Yobs=d.other[,col.Y]), nb, B, cbind(1:nb,1:nb), nX, rep(log(1),nX), rep(FALSE,nX), ce_cov, ce_partial, verbose=TRUE, parallel=TRUE, maxIter=5)
		time <- proc.time()-t1
print(time)

		preds <- compute_preds(fit, d.other, d.fold)
print(preds)

		list(fit=fit, time=time, preds=preds)
	})


}

set.seed(1983)
cv_error(5)

}
