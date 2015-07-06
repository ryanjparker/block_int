print(.libPaths())
.libPaths("/mnt/home/rjparker/Rlib")
print(.libPaths())
setwd("/mnt/home/rjparker/git/block_int")

which_exp <- 3
which_part <- 2
source("R/sim_nblocks_pred.R")
