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
	pred=FALSE,  # prediction study?
	# generate data for this number of inputs
	# p: number of inputs, b: sparsity
	p=c(10), b=c(1.5,4),
	#p=c(10), b=c(1,3,9),
	# tau: difficulty of problem
	tau=c(10),
	# number of observations per input
	Nper_p=c(50,100,250),
	#Nper_p=c(50),
	# number of locations to predict at
	Npred=1000,
	# % of observations per block
	Nblock_obs_ind=100, # this is a placeholder modified during the sim study
	# which inputs are not important?
	Fzero=0,
	# observation variance?
	sigma2=1
)

if (FALSE) {
	for (i in 1:nrow(sim.factors)) {
		dat <- generate_data(sim.design, sim.factors[i,])
		print(round(dat$theta,2))
	}
	done
}

options(cores=3)

# run the experiment for each combination of factors
#res <- lapply(1:nrow(sim.factors), function(i) {
#res <- lapply(1:1, function(i) {
res <- lapply(which_exp, function(i) {
	print(sim.factors[i,])
	exp_res <- sim_exp_est(sim.design, sim.factors[i,], i, which_part)
	save(exp_res, file=paste0("output/nblocks_exp_",i,"_",which_part,".RData"))

print(head(exp_res))

  exp_res
})

mres <- do.call(rbind, res)
