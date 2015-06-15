source("R/cov.R")
source("R/estimate.R")
library(fields)

options(cores=3); options(mc.cores=3)
set.seed(1983)

# read in sarcos data
load("data/sarcos.RData")

d.train <- sarcos$sarcos.inv
#d.train <- d.train[sample(nrow(d.train), 10000),] # sample the data
d.test  <- sarcos$sarcos.inv.test

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

# scale Y to have mean 0 and SD 1
Y.mean <- mean(d.train[,col.Y])
Y.sd   <- sd(d.train[,col.Y])
d.train[,col.Y] <- (d.train[,col.Y] - Y.mean)/Y.sd
d.test[,col.Y]  <- (d.test[,col.Y] - Y.mean)/Y.sd

n.train <- nrow(d.train)
n.test <- nrow(d.test)

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
}

if (TRUE) { # fit models

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


}
