# analyze output

"analyze_sims" <- function(exp) {

load(paste0("output/nblocks_exp_",exp,".RData"))
with(exp_res, {

	se.mult.factor <- 1
	s.mult.factor <- 1

	good.orac  <- orac.status
	good.full  <- full.status
	good.ir250 <- ir250.status
	good.ic250 <- ic250.status
	good.ir100 <- ir100.status
	good.ic100 <- ic100.status
	good.ir50 <- ir50.status
	good.ic50 <- ic50.status
	good.ir25 <- ir25.status
	good.ic25 <- ic25.status

	print(c(
mean(good.orac),
mean(good.full),
mean(good.ir250),
mean(good.ic250),
mean(good.ir100),
mean(good.ic100),
mean(good.ir50),
mean(good.ic50),
mean(good.ir25),
mean(good.ic25)
))

return(NA)

	# oracle
	orac <- list()
	orac$mse_t_mu <- 0
	orac$mse_t_se <- 0
	orac$rmse_p_mu <- mean(orac.rmse_p[good.orac])
	orac$rmse_p_se <- sd(orac.rmse_p)/sqrt(length(good.orac))*se.mult.factor

	# full
	sub <- which(good.orac&good.full)
	full <- list()
	full$rel_rmse_p_mu <- mean(full.rmse_p[sub]/orac.rmse_p[sub])
	full$rel_rmse_p_se <- sd(full.rmse_p[sub]/orac.rmse_p[sub])/sqrt(length(sub))*se.mult.factor
	full$rmse_s_mu <- mean(full.rmse_s_T[sub])*s.mult.factor
	full$rmse_s_se <- sd(full.rmse_s_T[sub])/sqrt(length(sub))*se.mult.factor
	full$top3 <- mean( unlist( lapply(sub, function(i) {
		a <- sort(sort(orac.Si_T[[i]],dec=T,index.return=TRUE)$ix[1:3])
		b <- sort(sort(full.Si_T[[i]],dec=T,index.return=TRUE)$ix[1:3])
		sum(a-b)==0
	}) ) )

	# ir 250
	sub <- which(good.orac&good.full&good.ir250)
	ir250 <- list()
	ir250$rel_rmse_t_mu <- mean(ir250.rmse_t[sub]/full.rmse_p[sub])
	ir250$rel_rmse_t_se <- sd(ir250.rmse_t[sub]/full.rmse_p[sub])/sqrt(length(sub))*se.mult.factor
	ir250$rel_rmse_p_mu <- mean(ir250.rmse_full_p[sub]/orac.rmse_p[sub])
	ir250$rel_rmse_p_se <- sd(ir250.rmse_full_p[sub]/orac.rmse_p[sub])/sqrt(length(sub))*se.mult.factor
	ir250$rmse_s_mu <- mean(ir250.rmse_s_T[sub])*s.mult.factor
	ir250$rmse_s_se <- sd(ir250.rmse_s_T[sub])/sqrt(length(sub))*se.mult.factor
	ir250$top3 <- mean( unlist( lapply(sub, function(i) {
		a <- sort(sort(orac.Si_T[[i]],dec=T,index.return=TRUE)$ix[1:3])
		b <- sort(sort(ir250.Si_T[[i]],dec=T,index.return=TRUE)$ix[1:3])
		sum(a-b)==0
	}) ) )
	ir250$speedup_mu <- mean(full.time[sub]/ir250.time[sub])
	ir250$speedup_se <- sd(full.time[sub]/ir250.time[sub])/sqrt(length(sub))*se.mult.factor

	# ic 250
	sub <- which(good.orac&good.full&good.ir250&good.ic250)
	ic250 <- list()
	ic250$rel_rmse_t_mu <- mean(ic250.rmse_t[sub]/full.rmse_p[sub])
	ic250$rel_rmse_t_se <- sd(ic250.rmse_t[sub]/full.rmse_p[sub])/sqrt(length(sub))*se.mult.factor
	ic250$rel_rmse_p_mu <- mean(ic250.rmse_full_p[sub]/orac.rmse_p[sub])
	ic250$rel_rmse_p_se <- sd(ic250.rmse_full_p[sub]/orac.rmse_p[sub])/sqrt(length(sub))*se.mult.factor
	ic250$rmse_s_mu <- mean(ic250.rmse_s_T[sub])*s.mult.factor
	ic250$rmse_s_se <- sd(ic250.rmse_s_T[sub])/sqrt(length(sub))*se.mult.factor
	ic250$top3 <- mean( unlist( lapply(sub, function(i) {
		a <- sort(sort(orac.Si_T[[i]],dec=T,index.return=TRUE)$ix[1:3])
		b <- sort(sort(ic250.Si_T[[i]],dec=T,index.return=TRUE)$ix[1:3])
		sum(a-b)==0
	}) ) )
	ic250$speedup_mu <- mean(full.time[sub]/(ir250.time[sub]+ic250.time[sub]))
	ic250$speedup_se <- sd(full.time[sub]/(ir250.time[sub]+ic250.time[sub]))/sqrt(length(sub))*se.mult.factor

	# ir 100
	sub <- which(good.orac&good.full&good.ir100)
	ir100 <- list()
	ir100$rel_rmse_t_mu <- mean(ir100.rmse_t[sub]/full.rmse_p[sub])
	ir100$rel_rmse_t_se <- sd(ir100.rmse_t[sub]/full.rmse_p[sub])/sqrt(length(sub))*se.mult.factor
	ir100$rel_rmse_p_mu <- mean(ir100.rmse_full_p[sub]/orac.rmse_p[sub])
	ir100$rel_rmse_p_se <- sd(ir100.rmse_full_p[sub]/orac.rmse_p[sub])/sqrt(length(sub))*se.mult.factor
	ir100$rmse_s_mu <- mean(ir100.rmse_s_T[sub])*s.mult.factor
	ir100$rmse_s_se <- sd(ir100.rmse_s_T[sub])/sqrt(length(sub))*se.mult.factor
	ir100$top3 <- mean( unlist( lapply(sub, function(i) {
		a <- sort(sort(orac.Si_T[[i]],dec=T,index.return=TRUE)$ix[1:3])
		b <- sort(sort(ir100.Si_T[[i]],dec=T,index.return=TRUE)$ix[1:3])
		sum(a-b)==0
	}) ) )
	ir100$speedup_mu <- mean(full.time[sub]/ir100.time[sub])
	ir100$speedup_se <- sd(full.time[sub]/ir100.time[sub])/sqrt(length(sub))*se.mult.factor

	# ic 100
	sub <- which(good.orac&good.full&good.ir100&good.ic100)
	ic100 <- list()
	ic100$rel_rmse_t_mu <- mean(ic100.rmse_t[sub]/full.rmse_p[sub])
	ic100$rel_rmse_t_se <- sd(ic100.rmse_t[sub]/full.rmse_p[sub])/sqrt(length(sub))*se.mult.factor
	ic100$rel_rmse_p_mu <- mean(ic100.rmse_full_p[sub]/orac.rmse_p[sub])
	ic100$rel_rmse_p_se <- sd(ic100.rmse_full_p[sub]/orac.rmse_p[sub])/sqrt(length(sub))*se.mult.factor
	ic100$rmse_s_mu <- mean(ic100.rmse_s_T[sub])*s.mult.factor
	ic100$rmse_s_se <- sd(ic100.rmse_s_T[sub])/sqrt(length(sub))*se.mult.factor
	ic100$top3 <- mean( unlist( lapply(sub, function(i) {
		a <- sort(sort(orac.Si_T[[i]],dec=T,index.return=TRUE)$ix[1:3])
		b <- sort(sort(ic100.Si_T[[i]],dec=T,index.return=TRUE)$ix[1:3])
		sum(a-b)==0
	}) ) )
	ic100$speedup_mu <- mean(full.time[sub]/(ir100.time[sub]+ic100.time[sub]))
	ic100$speedup_se <- sd(full.time[sub]/(ir100.time[sub]+ic100.time[sub]))/sqrt(length(sub))*se.mult.factor

	# ir 50
	sub <- which(good.orac&good.full&good.ir50)
	ir50 <- list()
	ir50$rel_rmse_t_mu <- mean(ir50.rmse_t[sub]/full.rmse_p[sub])
	ir50$rel_rmse_t_se <- sd(ir50.rmse_t[sub]/full.rmse_p[sub])/sqrt(length(sub))*se.mult.factor
	ir50$rel_rmse_p_mu <- mean(ir50.rmse_full_p[sub]/orac.rmse_p[sub])
	ir50$rel_rmse_p_se <- sd(ir50.rmse_full_p[sub]/orac.rmse_p[sub])/sqrt(length(sub))*se.mult.factor
	ir50$rmse_s_mu <- mean(ir50.rmse_s_T[sub])*s.mult.factor
	ir50$rmse_s_se <- sd(ir50.rmse_s_T[sub])/sqrt(length(sub))*se.mult.factor
	ir50$top3 <- mean( unlist( lapply(sub, function(i) {
		a <- sort(sort(orac.Si_T[[i]],dec=T,index.return=TRUE)$ix[1:3])
		b <- sort(sort(ir50.Si_T[[i]],dec=T,index.return=TRUE)$ix[1:3])
		sum(a-b)==0
	}) ) )
	ir50$speedup_mu <- mean(full.time[sub]/ir50.time[sub])
	ir50$speedup_se <- sd(full.time[sub]/ir50.time[sub])/sqrt(length(sub))*se.mult.factor

	# ic 50
	sub <- which(good.orac&good.full&good.ir50&good.ic50)
	ic50 <- list()
	ic50$rel_rmse_t_mu <- mean(ic50.rmse_t[sub]/full.rmse_p[sub])
	ic50$rel_rmse_t_se <- sd(ic50.rmse_t[sub]/full.rmse_p[sub])/sqrt(length(sub))*se.mult.factor
	ic50$rel_rmse_p_mu <- mean(ic50.rmse_full_p[sub]/orac.rmse_p[sub])
	ic50$rel_rmse_p_se <- sd(ic50.rmse_full_p[sub]/orac.rmse_p[sub])/sqrt(length(sub))*se.mult.factor
	ic50$rmse_s_mu <- mean(ic50.rmse_s_T[sub])*s.mult.factor
	ic50$rmse_s_se <- sd(ic50.rmse_s_T[sub])/sqrt(length(sub))*se.mult.factor
	ic50$top3 <- mean( unlist( lapply(sub, function(i) {
		a <- sort(sort(orac.Si_T[[i]],dec=T,index.return=TRUE)$ix[1:3])
		b <- sort(sort(ic50.Si_T[[i]],dec=T,index.return=TRUE)$ix[1:3])
		sum(a-b)==0
	}) ) )
	ic50$speedup_mu <- mean(full.time[sub]/(ir50.time[sub]+ic50.time[sub]))
	ic50$speedup_se <- sd(full.time[sub]/(ir50.time[sub]+ic50.time[sub]))/sqrt(length(sub))*se.mult.factor

	# ir 25
	sub <- which(good.orac&good.full&good.ir25)
	ir25 <- list()
	ir25$rel_rmse_t_mu <- mean(ir25.rmse_t[sub]/full.rmse_p[sub])
	ir25$rel_rmse_t_se <- sd(ir25.rmse_t[sub]/full.rmse_p[sub])/sqrt(length(sub))*se.mult.factor
	ir25$rel_rmse_p_mu <- mean(ir25.rmse_full_p[sub]/orac.rmse_p[sub])
	ir25$rel_rmse_p_se <- sd(ir25.rmse_full_p[sub]/orac.rmse_p[sub])/sqrt(length(sub))*se.mult.factor
	ir25$rmse_s_mu <- mean(ir25.rmse_s_T[sub])*s.mult.factor
	ir25$rmse_s_se <- sd(ir25.rmse_s_T[sub])/sqrt(length(sub))*se.mult.factor
	ir25$top3 <- mean( unlist( lapply(sub, function(i) {
		a <- sort(sort(orac.Si_T[[i]],dec=T,index.return=TRUE)$ix[1:3])
		b <- sort(sort(ir25.Si_T[[i]],dec=T,index.return=TRUE)$ix[1:3])
		sum(a-b)==0
	}) ) )
	ir25$speedup_mu <- mean(full.time[sub]/ir25.time[sub])
	ir25$speedup_se <- sd(full.time[sub]/ir25.time[sub])/sqrt(length(sub))*se.mult.factor

	# ic 25
	sub <- which(good.orac&good.full&good.ir25&good.ic25)
	ic25 <- list()
	ic25$rel_rmse_t_mu <- mean(ic25.rmse_t[sub]/full.rmse_p[sub])
	ic25$rel_rmse_t_se <- sd(ic25.rmse_t[sub]/full.rmse_p[sub])/sqrt(length(sub))*se.mult.factor
	ic25$rel_rmse_p_mu <- mean(ic25.rmse_full_p[sub]/orac.rmse_p[sub])
	ic25$rel_rmse_p_se <- sd(ic25.rmse_full_p[sub]/orac.rmse_p[sub])/sqrt(length(sub))*se.mult.factor
	ic25$rmse_s_mu <- mean(ic25.rmse_s_T[sub])*s.mult.factor
	ic25$rmse_s_se <- sd(ic25.rmse_s_T[sub])/sqrt(length(sub))*se.mult.factor
	ic25$top3 <- mean( unlist( lapply(sub, function(i) {
		a <- sort(sort(orac.Si_T[[i]],dec=T,index.return=TRUE)$ix[1:3])
		b <- sort(sort(ic25.Si_T[[i]],dec=T,index.return=TRUE)$ix[1:3])
		sum(a-b)==0
	}) ) )
	ic25$speedup_mu <- mean(full.time[sub]/(ir25.time[sub]+ic25.time[sub]))
	ic25$speedup_se <- sd(full.time[sub]/(ir25.time[sub]+ic25.time[sub]))/sqrt(length(sub))*se.mult.factor

cat("Full:\n");print(as.data.frame(full))
cat("IR 250:\n");print(as.data.frame(ir250))
cat("IC 250:\n");print(as.data.frame(ic250))
cat("IR 100:\n");print(as.data.frame(ir100))
cat("IC 100:\n");print(as.data.frame(ic100))
cat("IR 50:\n");print(as.data.frame(ir50))
cat("IC 50:\n");print(as.data.frame(ic50))
cat("IR 25:\n");print(as.data.frame(ir25))
cat("IC 25:\n");print(as.data.frame(ic25))
done

if (TRUE) {
cat("RMSE theta\n")
cat("Full:",mean(full.rmse_t),"\n")
cat("250:",mean(ir250.rmse_t,na.rm=TRUE),mean(ir250.rmse_t/full.rmse_t,na.rm=TRUE),mean(ic250.rmse_t,na.rm=TRUE),mean(ic250.rmse_t/full.rmse_t,na.rm=TRUE),"\n")
cat("100:",mean(ir100.rmse_t,na.rm=TRUE),mean(ir100.rmse_t/full.rmse_t,na.rm=TRUE),mean(ic100.rmse_t,na.rm=TRUE),mean(ic100.rmse_t/full.rmse_t,na.rm=TRUE),"\n")
cat("50:",mean(ir50.rmse_t,na.rm=TRUE),mean(ir50.rmse_t/full.rmse_t,na.rm=TRUE),mean(ic50.rmse_t,na.rm=TRUE),mean(ic50.rmse_t/full.rmse_t,na.rm=TRUE),"\n")
cat("25:",mean(ir25.rmse_t,na.rm=TRUE),mean(ir25.rmse_t/full.rmse_t,na.rm=TRUE),mean(ic25.rmse_t,na.rm=TRUE),mean(ic25.rmse_t/full.rmse_t,na.rm=TRUE),"\n")
}

if (TRUE) {
cat("RMSE pred\n")
cat("Oracle:",mean(orac.rmse_p),"\n")
cat("Full:",mean(full.rmse_p),"\n")
cat("250:",mean(ir250.rmse_full_p,na.rm=TRUE),mean(ir250.rmse_full_p/full.rmse_p,na.rm=TRUE),mean(ic250.rmse_full_p,na.rm=TRUE),mean(ic250.rmse_full_p/full.rmse_p,na.rm=TRUE),"\n")
cat("100:",mean(ir100.rmse_full_p,na.rm=TRUE),mean(ir100.rmse_full_p/full.rmse_p,na.rm=TRUE),mean(ic100.rmse_full_p,na.rm=TRUE),mean(ic100.rmse_full_p/full.rmse_p,na.rm=TRUE),"\n")
cat("50:",mean(ir50.rmse_full_p,na.rm=TRUE),mean(ir50.rmse_full_p/full.rmse_p,na.rm=TRUE),mean(ic50.rmse_full_p,na.rm=TRUE),mean(ic50.rmse_full_p/full.rmse_p,na.rm=TRUE),"\n")
cat("25:",mean(ir25.rmse_full_p,na.rm=TRUE),mean(ir25.rmse_full_p/full.rmse_p,na.rm=TRUE),mean(ic25.rmse_full_p,na.rm=TRUE),mean(ic25.rmse_full_p/full.rmse_p,na.rm=TRUE),"\n")
}

if (FALSE) {
cat("RMSE sensitivity (main)\n")
cat("Full:",mean(full.rmse_s_S),"\n")
cat("250:",mean(ir250.rmse_s_S,na.rm=TRUE),mean(ir250.rmse_s_S/full.rmse_s_S,na.rm=TRUE),mean(ic250.rmse_s_S,na.rm=TRUE),mean(ic250.rmse_s_S/full.rmse_s_S,na.rm=TRUE),"\n")
cat("100:",mean(ir100.rmse_s_S,na.rm=TRUE),mean(ir100.rmse_s_S/full.rmse_s_S,na.rm=TRUE),mean(ic100.rmse_s_S,na.rm=TRUE),mean(ic100.rmse_s_S/full.rmse_s_S,na.rm=TRUE),"\n")
cat("50:",mean(ir50.rmse_s_S,na.rm=TRUE),mean(ir50.rmse_s_S/full.rmse_s_S,na.rm=TRUE),mean(ic50.rmse_s_S,na.rm=TRUE),mean(ic50.rmse_s_S/full.rmse_s_S,na.rm=TRUE),"\n")
cat("25:",mean(ir25.rmse_s_S,na.rm=TRUE),mean(ir25.rmse_s_S/full.rmse_s_S,na.rm=TRUE),mean(ic25.rmse_s_S,na.rm=TRUE),mean(ic25.rmse_s_S/full.rmse_s_S,na.rm=TRUE),"\n")
}

if (TRUE) {
cat("RMSE sensitivity (total)\n")
cat("Full:",mean(full.rmse_s_T),"\n")
cat("250:",mean(ir250.rmse_s_T,na.rm=TRUE),mean(ir250.rmse_s_T/full.rmse_s_T,na.rm=TRUE),mean(ic250.rmse_s_T,na.rm=TRUE),mean(ic250.rmse_s_T/full.rmse_s_T,na.rm=TRUE),"\n")
cat("100:",mean(ir100.rmse_s_T,na.rm=TRUE),mean(ir100.rmse_s_T/full.rmse_s_T,na.rm=TRUE),mean(ic100.rmse_s_T,na.rm=TRUE),mean(ic100.rmse_s_T/full.rmse_s_T,na.rm=TRUE),"\n")
cat("50:",mean(ir50.rmse_s_T,na.rm=TRUE),mean(ir50.rmse_s_T/full.rmse_s_T,na.rm=TRUE),mean(ic50.rmse_s_T,na.rm=TRUE),mean(ic50.rmse_s_T/full.rmse_s_T,na.rm=TRUE),"\n")
cat("25:",mean(ir25.rmse_s_T,na.rm=TRUE),mean(ir25.rmse_s_T/full.rmse_s_T,na.rm=TRUE),mean(ic25.rmse_s_T,na.rm=TRUE),mean(ic25.rmse_s_T/full.rmse_s_T,na.rm=TRUE),"\n")
}

if (TRUE) {
cat("Timing\n")
cat("Full:",round(mean(full.time),2),"\n")
cat("250:",round(mean(ir250.time),2),round(mean(full.time/ir250.time),2),round(mean(ic250.time),2),round(mean(full.time/(ir250.time+ic250.time)),2),"\n")
cat("100:",round(mean(ir100.time),2),round(mean(full.time/ir100.time),2),round(mean(ic100.time),2),round(mean(full.time/(ir100.time+ic100.time)),2),"\n")
cat("50:",round(mean(ir50.time),2),round(mean(full.time/ir50.time),2),round(mean(ic50.time),2),round(mean(full.time/(ir50.time+ic50.time)),2),"\n")
cat("25:",round(mean(ir25.time),2),round(mean(full.time/ir25.time),2),round(mean(ic25.time),2),round(mean(full.time/(ir25.time+ic25.time)),2),"\n")
}


#11-apply(do.call("rbind",res[[1]]$orac.Si),1,rank)

})

}

for (i in c(1,2,3,4,5,6)) {
	analyze_sims(i)
}
