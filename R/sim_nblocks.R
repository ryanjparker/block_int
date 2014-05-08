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
	Nsens=1000
)

sim.factors <- expand.grid(
	# generate data for this number of inputs
	#p=c(5,10,15),
	p=c(10),
	#p=c(1),
	# b: sparsity
	b=c(1,3,9),
	# tau: difficulty of problem
	tau=c(3),
	# number of observations per input
	Nper_p=c(50,100,250),
	# number of locations to predict at
	Npred=100,
	# % of observations per block
	Nblock_obs_ind=100,
	# which inputs are not important?
	Fzero=1,
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

options(cores=4)

# run the experiment for each combination of factors
#res <- lapply(1:nrow(sim.factors), function(i) {
#res <- lapply(1:1, function(i) {
res <- lapply(which_exp, function(i) {
  print(sim.factors[i,])
  exp_res <- sim_exp(sim.design, sim.factors[i,], i)
	save(exp_res, file=paste0("output/nblocks_exp_",i,".RData"))

print(head(exp_res))

  exp_res
})

mres <- do.call(rbind, res)
