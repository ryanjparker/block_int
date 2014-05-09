# analyze output

"analyze_sims" <- function(exp) {

load(paste0("output/nblocks_exp_",exp,".RData"))
with(exp_res, {

cat("Timing\n")
cat("Full:",mean(full.time),"\n")
cat("250:",mean(ir250.time),mean(full.time/ir250.time),mean(ic250.time),mean(full.time/(ir250.time+ic250.time)),"\n")
cat("100:",mean(ir100.time),mean(full.time/ir100.time),mean(ic100.time),mean(full.time/(ir100.time+ic100.time)),"\n")
cat("50:",mean(ir50.time),mean(full.time/ir50.time),mean(ic50.time),mean(full.time/(ir50.time+ic50.time)),"\n")
cat("25:",mean(ir25.time),mean(full.time/ir25.time),mean(ic25.time),mean(full.time/(ir25.time+ic25.time)),"\n")

cat("RMSE theta\n")
cat("Full:",mean(full.rmse_t),"\n")
cat("250:",mean(ir250.rmse_t,na.rm=TRUE),mean(ir250.rmse_t/full.rmse_t,na.rm=TRUE),mean(ic250.rmse_t,na.rm=TRUE),mean(ic250.rmse_t/full.rmse_t,na.rm=TRUE),"\n")
cat("100:",mean(ir100.rmse_t,na.rm=TRUE),mean(ir100.rmse_t/full.rmse_t,na.rm=TRUE),mean(ic100.rmse_t,na.rm=TRUE),mean(ic100.rmse_t/full.rmse_t,na.rm=TRUE),"\n")
cat("50:",mean(ir50.rmse_t,na.rm=TRUE),mean(ir50.rmse_t/full.rmse_t,na.rm=TRUE),mean(ic50.rmse_t,na.rm=TRUE),mean(ic50.rmse_t/full.rmse_t,na.rm=TRUE),"\n")
cat("25:",mean(ir25.rmse_t,na.rm=TRUE),mean(ir25.rmse_t/full.rmse_t,na.rm=TRUE),mean(ic25.rmse_t,na.rm=TRUE),mean(ic25.rmse_t/full.rmse_t,na.rm=TRUE),"\n")

cat("RMSE pred\n")
cat("Oracle:",mean(orac.rmse_p),"\n")
cat("Full:",mean(full.rmse_p),"\n")
cat("250:",mean(ir250.rmse_full_p,na.rm=TRUE),mean(ir250.rmse_full_p/full.rmse_p,na.rm=TRUE),mean(ic250.rmse_full_p,na.rm=TRUE),mean(ic250.rmse_full_p/full.rmse_p,na.rm=TRUE),"\n")
cat("100:",mean(ir100.rmse_full_p,na.rm=TRUE),mean(ir100.rmse_full_p/full.rmse_p,na.rm=TRUE),mean(ic100.rmse_full_p,na.rm=TRUE),mean(ic100.rmse_full_p/full.rmse_p,na.rm=TRUE),"\n")
cat("50:",mean(ir50.rmse_full_p,na.rm=TRUE),mean(ir50.rmse_full_p/full.rmse_p,na.rm=TRUE),mean(ic50.rmse_full_p,na.rm=TRUE),mean(ic50.rmse_full_p/full.rmse_p,na.rm=TRUE),"\n")
cat("25:",mean(ir25.rmse_full_p,na.rm=TRUE),mean(ir25.rmse_full_p/full.rmse_p,na.rm=TRUE),mean(ic25.rmse_full_p,na.rm=TRUE),mean(ic25.rmse_full_p/full.rmse_p,na.rm=TRUE),"\n")

cat("RMSE sensitivity\n")
cat("Oracle:",mean(orac.rmse_s),"\n")
cat("Full:",mean(full.rmse_s),"\n")
cat("250:",mean(ir250.rmse_s,na.rm=TRUE),mean(ir250.rmse_s/full.rmse_s,na.rm=TRUE),mean(ic250.rmse_s,na.rm=TRUE),mean(ic250.rmse_s/full.rmse_s,na.rm=TRUE),"\n")
cat("100:",mean(ir100.rmse_s,na.rm=TRUE),mean(ir100.rmse_s/full.rmse_s,na.rm=TRUE),mean(ic100.rmse_s,na.rm=TRUE),mean(ic100.rmse_s/full.rmse_s,na.rm=TRUE),"\n")
cat("50:",mean(ir50.rmse_s,na.rm=TRUE),mean(ir50.rmse_s/full.rmse_s,na.rm=TRUE),mean(ic50.rmse_s,na.rm=TRUE),mean(ic50.rmse_s/full.rmse_s,na.rm=TRUE),"\n")
cat("25:",mean(ir25.rmse_s,na.rm=TRUE),mean(ir25.rmse_s/full.rmse_s,na.rm=TRUE),mean(ic25.rmse_s,na.rm=TRUE),mean(ic25.rmse_s/full.rmse_s,na.rm=TRUE),"\n")

#11-apply(do.call("rbind",res[[1]]$orac.Si),1,rank)

})

}

for (i in c(1)) {
	analyze_sims(i)
}
