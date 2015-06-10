# run simulation study on number of blocks
library(fields)
library(lhs)
library(MASS)
library(multicore)

source("R/sim_study.R")

# fixed design parameters
sim.design <- list(
	# number of replications
	Nreps=100,
	# number used in sensitivity analysis
	Nsens=5000
)

sim.factors <- expand.grid(
	pred=TRUE,  # prediction study?
	# generate data for this number of inputs
	#p=c(10,20,30),
	#p=c(1), b=c(1),
	#p=c(2), b=c(1),
	#p=c(3),  b=c(1),
	#p=c(4),  b=c(1),
	#p=c(5),  b=c(1),
	p=c(10), b=c(1.5,4),
	#p=c(20), b=c(7.55),
	#p=c(30), b=c(2),
	# tau: difficulty of problem
	#tau=c(3,10,20),
	tau=c(10),
	# number of observations per input
	Nper_p=c(1000,2500,5000,10000),
	#Nper_p=c(10000),
	#Nper_p=c(2000),
	# number of locations to predict at
	Npred=1000,
	# number of zeros
	Fzero=0,
	# observation variance?
	sigma2=1
)

options(cores=3)

# run the experiment for each combination of factors
#res <- lapply(1:nrow(sim.factors), function(i) {
#res <- lapply(1:1, function(i) {
res <- lapply(which_exp, function(i) {
#	try({ rm("predXobs", "predXpred", "predSigma", "predInvSigma", pos = ".GlobalEnv") })
	print(sim.factors[i,])
	exp_res <- sim_exp_pred(sim.design, sim.factors[i,], i, which_part)
	save(exp_res, file=paste0("output/nblocks_exp_pred_",i,"_",which_part,".RData"))

print(head(exp_res))

  exp_res
})

mres <- do.call(rbind, res)
