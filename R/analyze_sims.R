# analyze output

"round_df" <- function(df) {
	for (i in 1:ncol(df)) {
		if (is.numeric(df[,i])) df[,i] <- round(df[,i], 2)
	}

	df
}

"analyze_sims" <- function(exp) {

	Nparts <- 20
	if (Nparts > 1) {
		cer <- data.frame()
		for (part in 1:Nparts) {
		#for (part in (1:Nparts)[-c(9,20)]) {
			load(paste0("output/nblocks_exp_",exp,"_",part,".RData"))
			cer <- rbind(cer, exp_res)
		}
		exp_res <- cer
	} else {
		load(paste0("output/nblocks_exp_",exp,".RData"))
	}

with(exp_res, {

	se.mult.factor <- 1
	s.mult.factor <- 1

	good.orac  <- orac.status
	good.full  <- full.status
	# ind
	good.ir250 <- ir250.status
	good.ic250 <- ic250.status
	good.ir100 <- ir100.status
	good.ic100 <- ic100.status
	good.ir50 <- ir50.status
	good.ic50 <- ic50.status
	good.ir25 <- ir25.status
	good.ic25 <- ic25.status
	# dep
	good.dr50 <- dr50.status
	good.dc50 <- dc50.status
	good.dr25 <- dr25.status
	good.dc25 <- dc25.status

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
mean(good.ic25),
mean(good.dr50),
mean(good.dc50),
mean(good.dr25),
mean(good.dc25)
))

	# oracle
	orac <- list()
	orac$mse_t_mu <- 0
	orac$mse_t_se <- 0
	orac$rmse_p_mu <- mean(orac.rmse_p[good.orac])
	orac$rmse_p_se <- sd(orac.rmse_p)/sqrt(length(good.orac))*se.mult.factor

	# full
	sub <- which(good.orac&good.full)
print(length(sub))
	full <- list()
	full$rmse_t_mu <- mean(full.rmse_t[sub])
	full$rmse_p_mu <- mean(full.rmse_p[sub])
	full$rmse_s_mu <- mean(full.rmse_s_T[sub])*s.mult.factor
	full$rel_rmse_t_mu <- 1
	full$rel_rmse_p_mu <- 1 # mean(full.rmse_p[sub]/orac.rmse_p[sub])
	full$rel_rmse_s_mu <- 1 #mean(full.rmse_s_T[sub])*s.mult.factor
	#full$rel_rmse_p_se <- sd(full.rmse_p[sub]/orac.rmse_p[sub])/sqrt(length(sub))*se.mult.factor
	#full$rel_rmse_s_se <- sd(full.rmse_s_T[sub])/sqrt(length(sub))*se.mult.factor
	full$top3 <- mean( unlist( lapply(sub, function(i) {
		a <- sort(sort(orac.Si_T[[i]],dec=T,index.return=TRUE)$ix[1:3])
		b <- sort(sort(full.Si_T[[i]],dec=T,index.return=TRUE)$ix[1:3])
		sum(a-b)==0
	}) ) )
	full$speedup_mu <- 1
	full$time <- mean(full.time[sub])/60

	# ir 250
	sub <- which(good.orac&good.full&good.ir250)
	ir250 <- list()
	ir250$rmse_t_mu <- mean(ir250.rmse_t[sub])
	ir250$rmse_p_mu <- mean(ir250.rmse_full_p[sub])
	ir250$rmse_s_mu <- mean(ir250.rmse_s_T[sub])*s.mult.factor
	ir250$rel_rmse_t_mu <- mean(ir250.rmse_t[sub]/full.rmse_t[sub])
	ir250$rel_rmse_p_mu <- mean(ir250.rmse_full_p[sub]/full.rmse_p[sub])
	ir250$rel_rmse_s_mu <- mean(ir250.rmse_s_T[sub]/full.rmse_s_T[sub])*s.mult.factor
	ir250$top3 <- mean( unlist( lapply(sub, function(i) {
		a <- sort(sort(orac.Si_T[[i]],dec=T,index.return=TRUE)$ix[1:3])
		b <- sort(sort(ir250.Si_T[[i]],dec=T,index.return=TRUE)$ix[1:3])
		sum(a-b)==0
	}) ) )
	ir250$speedup_mu <- mean(full.time[sub]/ir250.time[sub])
	ir250$time <- mean(ir250.time[sub])/60
	#ir250$rel_rmse_t_se <- sd(ir250.rmse_t[sub]/full.rmse_p[sub])/sqrt(length(sub))*se.mult.factor
	#ir250$rel_rmse_p_se <- sd(ir250.rmse_full_p[sub]/orac.rmse_p[sub])/sqrt(length(sub))*se.mult.factor
	#ir250$rel_rmse_s_se <- sd(ir250.rmse_s_T[sub]/full.rmse_s_T[sub])/sqrt(length(sub))*se.mult.factor
	#ir250$speedup_se <- sd(full.time[sub]/ir250.time[sub])/sqrt(length(sub))*se.mult.factor

	# ic 250
	sub <- which(good.orac&good.full&good.ir250&good.ic250)
	ic250 <- list()
	ic250$rmse_t_mu <- mean(ic250.rmse_t[sub])
	ic250$rmse_p_mu <- mean(ic250.rmse_full_p[sub])
	ic250$rmse_s_mu <- mean(ic250.rmse_s_T[sub])*s.mult.factor
	ic250$rel_rmse_t_mu <- mean(ic250.rmse_t[sub]/full.rmse_t[sub])
	ic250$rel_rmse_p_mu <- mean(ic250.rmse_full_p[sub]/full.rmse_p[sub])
	ic250$rel_rmse_s_mu <- mean(ic250.rmse_s_T[sub]/full.rmse_s_T[sub])*s.mult.factor
	ic250$top3 <- mean( unlist( lapply(sub, function(i) {
		a <- sort(sort(orac.Si_T[[i]],dec=T,index.return=TRUE)$ix[1:3])
		b <- sort(sort(ic250.Si_T[[i]],dec=T,index.return=TRUE)$ix[1:3])
		sum(a-b)==0
	}) ) )
	#ic250$speedup_mu <- mean(full.time[sub]/(ir250.time[sub]+ic250.time[sub]))
	ic250$speedup_mu <- mean(full.time[sub]/(ic250.time[sub]))
	ic250$time <- mean(ic250.time[sub])/60

	# ir 100
	sub <- which(good.orac&good.full&good.ir100)
	ir100 <- list()
	ir100$rmse_t_mu <- mean(ir100.rmse_t[sub])
	ir100$rmse_p_mu <- mean(ir100.rmse_full_p[sub])
	ir100$rmse_s_mu <- mean(ir100.rmse_s_T[sub])*s.mult.factor
	ir100$rel_rmse_t_mu <- mean(ir100.rmse_t[sub]/full.rmse_t[sub])
	ir100$rel_rmse_p_mu <- mean(ir100.rmse_full_p[sub]/full.rmse_p[sub])
	ir100$rel_rmse_s_mu <- mean(ir100.rmse_s_T[sub]/full.rmse_s_T[sub])*s.mult.factor
	ir100$top3 <- mean( unlist( lapply(sub, function(i) {
		a <- sort(sort(orac.Si_T[[i]],dec=T,index.return=TRUE)$ix[1:3])
		b <- sort(sort(ir100.Si_T[[i]],dec=T,index.return=TRUE)$ix[1:3])
		sum(a-b)==0
	}) ) )
	ir100$speedup_mu <- mean(full.time[sub]/ir100.time[sub])
	ir100$time <- mean(ir100.time[sub])/60

	# ic 100
	sub <- which(good.orac&good.full&good.ir100&good.ic100)
	ic100 <- list()
	ic100$rmse_t_mu <- mean(ic100.rmse_t[sub])
	ic100$rmse_p_mu <- mean(ic100.rmse_full_p[sub])
	ic100$rmse_s_mu <- mean(ic100.rmse_s_T[sub])*s.mult.factor
	ic100$rel_rmse_t_mu <- mean(ic100.rmse_t[sub]/full.rmse_t[sub])
	ic100$rel_rmse_p_mu <- mean(ic100.rmse_full_p[sub]/full.rmse_p[sub])
	ic100$rel_rmse_s_mu <- mean(ic100.rmse_s_T[sub]/full.rmse_s_T[sub])*s.mult.factor
	ic100$top3 <- mean( unlist( lapply(sub, function(i) {
		a <- sort(sort(orac.Si_T[[i]],dec=T,index.return=TRUE)$ix[1:3])
		b <- sort(sort(ic100.Si_T[[i]],dec=T,index.return=TRUE)$ix[1:3])
		sum(a-b)==0
	}) ) )
	#ic100$speedup_mu <- mean(full.time[sub]/(ir100.time[sub]+ic100.time[sub]))
	ic100$speedup_mu <- mean(full.time[sub]/(ic100.time[sub]))
	ic100$time <- mean(ic100.time[sub])/60

	# ir 50
	sub <- which(good.orac&good.full&good.ir50)
	ir50 <- list()
	ir50$rmse_t_mu <- mean(ir50.rmse_t[sub])
	ir50$rmse_p_mu <- mean(ir50.rmse_full_p[sub])
	ir50$rmse_s_mu <- mean(ir50.rmse_s_T[sub])*s.mult.factor
	ir50$rel_rmse_t_mu <- mean(ir50.rmse_t[sub]/full.rmse_t[sub])
	ir50$rel_rmse_p_mu <- mean(ir50.rmse_full_p[sub]/full.rmse_p[sub])
	ir50$rel_rmse_s_mu <- mean(ir50.rmse_s_T[sub]/full.rmse_s_T[sub])*s.mult.factor
	ir50$top3 <- mean( unlist( lapply(sub, function(i) {
		a <- sort(sort(orac.Si_T[[i]],dec=T,index.return=TRUE)$ix[1:3])
		b <- sort(sort(ir50.Si_T[[i]],dec=T,index.return=TRUE)$ix[1:3])
		sum(a-b)==0
	}) ) )
	ir50$speedup_mu <- mean(full.time[sub]/ir50.time[sub])
	ir50$time <- mean(ir50.time[sub])/60

	# ic 50
	sub <- which(good.orac&good.full&good.ir50&good.ic50)
	ic50 <- list()
	ic50$rmse_t_mu <- mean(ic50.rmse_t[sub])
	ic50$rmse_p_mu <- mean(ic50.rmse_full_p[sub])
	ic50$rmse_s_mu <- mean(ic50.rmse_s_T[sub])*s.mult.factor
	ic50$rel_rmse_t_mu <- mean(ic50.rmse_t[sub]/full.rmse_t[sub])
	ic50$rel_rmse_p_mu <- mean(ic50.rmse_full_p[sub]/full.rmse_p[sub])
	ic50$rel_rmse_s_mu <- mean(ic50.rmse_s_T[sub]/full.rmse_s_T[sub])*s.mult.factor
	ic50$top3 <- mean( unlist( lapply(sub, function(i) {
		a <- sort(sort(orac.Si_T[[i]],dec=T,index.return=TRUE)$ix[1:3])
		b <- sort(sort(ic50.Si_T[[i]],dec=T,index.return=TRUE)$ix[1:3])
		sum(a-b)==0
	}) ) )
	#ic50$speedup_mu <- mean(full.time[sub]/(ir50.time[sub]+ic50.time[sub]))
	ic50$speedup_mu <- mean(full.time[sub]/(ic50.time[sub]))
	ic50$time <- mean(ic50.time[sub])/60

	# ir 25
	sub <- which(good.orac&good.full&good.ir25)
	ir25 <- list()
	ir25$rmse_t_mu <- mean(ir25.rmse_t[sub])
	ir25$rmse_p_mu <- mean(ir25.rmse_full_p[sub])
	ir25$rmse_s_mu <- mean(ir25.rmse_s_T[sub])*s.mult.factor
	ir25$rmse_t_mu <- mean(ir25.rmse_t[sub])
	ir25$rmse_p_mu <- mean(ir25.rmse_full_p[sub])
	ir25$rmse_s_mu <- mean(ir25.rmse_s_T[sub])*s.mult.factor
	ir25$rel_rmse_t_mu <- mean(ir25.rmse_t[sub]/full.rmse_t[sub])
	ir25$rel_rmse_p_mu <- mean(ir25.rmse_full_p[sub]/full.rmse_p[sub])
	ir25$rel_rmse_s_mu <- mean(ir25.rmse_s_T[sub]/full.rmse_s_T[sub])*s.mult.factor
	ir25$top3 <- mean( unlist( lapply(sub, function(i) {
		a <- sort(sort(orac.Si_T[[i]],dec=T,index.return=TRUE)$ix[1:3])
		b <- sort(sort(ir25.Si_T[[i]],dec=T,index.return=TRUE)$ix[1:3])
		sum(a-b)==0
	}) ) )
	ir25$speedup_mu <- mean(full.time[sub]/ir25.time[sub])
	ir25$time <- mean(ir25.time[sub])/60

	# ic 25
	sub <- which(good.orac&good.full&good.ir25&good.ic25)
	ic25 <- list()
	ic25$rmse_t_mu <- mean(ic25.rmse_t[sub])
	ic25$rmse_p_mu <- mean(ic25.rmse_full_p[sub])
	ic25$rmse_s_mu <- mean(ic25.rmse_s_T[sub])*s.mult.factor
	ic25$rel_rmse_t_mu <- mean(ic25.rmse_t[sub]/full.rmse_t[sub])
	ic25$rel_rmse_p_mu <- mean(ic25.rmse_full_p[sub]/full.rmse_p[sub])
	ic25$rel_rmse_s_mu <- mean(ic25.rmse_s_T[sub]/full.rmse_s_T[sub])*s.mult.factor
	ic25$top3 <- mean( unlist( lapply(sub, function(i) {
		a <- sort(sort(orac.Si_T[[i]],dec=T,index.return=TRUE)$ix[1:3])
		b <- sort(sort(ic25.Si_T[[i]],dec=T,index.return=TRUE)$ix[1:3])
		sum(a-b)==0
	}) ) )
	#ic25$speedup_mu <- mean(full.time[sub]/(ir25.time[sub]+ic25.time[sub]))
	ic25$speedup_mu <- mean(full.time[sub]/(ic25.time[sub]))
	ic25$time <- mean(ic25.time[sub])/60

	# dr 50
	sub <- which(good.orac&good.full&good.dr50)
	dr50 <- list()
	dr50$rmse_t_mu <- mean(dr50.rmse_t[sub])
	dr50$rmse_p_mu <- mean(dr50.rmse_full_p[sub])
	dr50$rmse_s_mu <- mean(dr50.rmse_s_T[sub])*s.mult.factor
	dr50$rel_rmse_t_mu <- mean(dr50.rmse_t[sub]/full.rmse_t[sub])
	dr50$rel_rmse_p_mu <- mean(dr50.rmse_full_p[sub]/full.rmse_p[sub])
	dr50$rel_rmse_s_mu <- mean(dr50.rmse_s_T[sub]/full.rmse_s_T[sub])*s.mult.factor
	dr50$top3 <- mean( unlist( lapply(sub, function(i) {
		a <- sort(sort(orac.Si_T[[i]],dec=T,index.return=TRUE)$ix[1:3])
		b <- sort(sort(dr50.Si_T[[i]],dec=T,index.return=TRUE)$ix[1:3])
		sum(a-b)==0
	}) ) )
	dr50$speedup_mu <- mean(full.time[sub]/dr50.time[sub])
	dr50$time <- mean(dr50.time[sub])/60

	# dc 50
	sub <- which(good.orac&good.full&good.dr50&good.dc50)
	dc50 <- list()
	dc50$rmse_t_mu <- mean(dc50.rmse_t[sub])
	dc50$rmse_p_mu <- mean(dc50.rmse_full_p[sub])
	dc50$rmse_s_mu <- mean(dc50.rmse_s_T[sub])*s.mult.factor
	dc50$rel_rmse_t_mu <- mean(dc50.rmse_t[sub]/full.rmse_t[sub])
	dc50$rel_rmse_p_mu <- mean(dc50.rmse_full_p[sub]/full.rmse_p[sub])
	dc50$rel_rmse_s_mu <- mean(dc50.rmse_s_T[sub]/full.rmse_s_T[sub])*s.mult.factor
	dc50$top3 <- mean( unlist( lapply(sub, function(i) {
		a <- sort(sort(orac.Si_T[[i]],dec=T,index.return=TRUE)$ix[1:3])
		b <- sort(sort(dc50.Si_T[[i]],dec=T,index.return=TRUE)$ix[1:3])
		sum(a-b)==0
	}) ) )
	#dc50$speedup_mu <- mean(full.time[sub]/(dr50.time[sub]+dc50.time[sub]))
	dc50$speedup_mu <- mean(full.time[sub]/(dc50.time[sub]))
	dc50$time <- mean(dc50.time[sub])/60

	# dr 25
	sub <- which(good.orac&good.full&good.dr25)
	dr25 <- list()
	dr25$rmse_t_mu <- mean(dr25.rmse_t[sub])
	dr25$rmse_p_mu <- mean(dr25.rmse_full_p[sub])
	dr25$rmse_s_mu <- mean(dr25.rmse_s_T[sub])*s.mult.factor
	dr25$rel_rmse_t_mu <- mean(dr25.rmse_t[sub]/full.rmse_t[sub])
	dr25$rel_rmse_p_mu <- mean(dr25.rmse_full_p[sub]/full.rmse_p[sub])
	dr25$rel_rmse_s_mu <- mean(dr25.rmse_s_T[sub]/full.rmse_s_T[sub])*s.mult.factor
	dr25$top3 <- mean( unlist( lapply(sub, function(i) {
		a <- sort(sort(orac.Si_T[[i]],dec=T,index.return=TRUE)$ix[1:3])
		b <- sort(sort(dr25.Si_T[[i]],dec=T,index.return=TRUE)$ix[1:3])
		sum(a-b)==0
	}) ) )
	dr25$speedup_mu <- mean(full.time[sub]/dr25.time[sub])
	dr25$time <- mean(dr25.time[sub])/60

	# dc 25
	sub <- which(good.orac&good.full&good.dr25&good.dc25)
	dc25 <- list()
	dc25$rmse_t_mu <- mean(dc25.rmse_t[sub])
	dc25$rmse_p_mu <- mean(dc25.rmse_full_p[sub])
	dc25$rmse_s_mu <- mean(dc25.rmse_s_T[sub])*s.mult.factor
	dc25$rel_rmse_t_mu <- mean(dc25.rmse_t[sub]/full.rmse_t[sub])
	dc25$rel_rmse_p_mu <- mean(dc25.rmse_full_p[sub]/full.rmse_p[sub])
	dc25$rel_rmse_s_mu <- mean(dc25.rmse_s_T[sub]/full.rmse_s_T[sub])*s.mult.factor
	dc25$top3 <- mean( unlist( lapply(sub, function(i) {
		a <- sort(sort(orac.Si_T[[i]],dec=T,index.return=TRUE)$ix[1:3])
		b <- sort(sort(dc25.Si_T[[i]],dec=T,index.return=TRUE)$ix[1:3])
		sum(a-b)==0
	}) ) )
	#dc25$speedup_mu <- mean(full.time[sub]/(dr25.time[sub]+dc25.time[sub]))
	dc25$speedup_mu <- mean(full.time[sub]/(dc25.time[sub]))
	dc25$time <- mean(dc25.time[sub])/60

	# round numeric fields
	full <- round_df(as.data.frame(full))
	ir250 <- round_df(as.data.frame(ir250))
	ic250 <- round_df(as.data.frame(ic250))
	ir100 <- round_df(as.data.frame(ir100))
	ic100 <- round_df(as.data.frame(ic100))
	ir50  <- round_df(as.data.frame(ir50))
	ic50  <- round_df(as.data.frame(ic50))
	ir25  <- round_df(as.data.frame(ir25))
	ic25  <- round_df(as.data.frame(ic25))
	dr50  <- round_df(as.data.frame(dr50))
	dc50  <- round_df(as.data.frame(dc50))
	dr25  <- round_df(as.data.frame(dr25))
	dc25  <- round_df(as.data.frame(dc25))

"tex_show" <- function(t,r) {
	str <- paste("& &",t)
	for (i in 1:length(r)) {
		if (length(grep("rel_", colnames(r[i]))) == 1) next
		if (colnames(r[i]) == "time") next
		f <- 1
		if (colnames(r[i]) == "top3") f <- 100

		str <- paste(str, "& ",r[i]*f)

		# does this have a rel?
		has.rel <- grep(paste0("rel_",colnames(r[i])), colnames(r))

		if (length(has.rel) > 0) {
			str <- paste0(str, " (", r[has.rel], ")")
		}

	}
	cat(str,"\\\\ \n")
	str
}

	print(full)
	tex_show("Full",full)
	tex_show("Ind R 250",ir250)
	tex_show("Ind C 250",ic250)
	tex_show("Ind R 100",ir100)
	tex_show("Ind C 100",ic100)
	tex_show("Ind R 50",ir50)
	tex_show("Ind C 50",ic50)
	tex_show("Ind R 25",ir25)
	tex_show("Ind C 25",ic25)
	tex_show("Dep R 50",dr50)
	tex_show("Dep C 50",dc50)
	tex_show("Dep R 25",dr25)
	tex_show("Dep C 25",dc25)
return(NA)

done
	print(full)
	print(rbind(
		ir250,ic250,ir100,ic100,ir50,ic50,
		ir25,ic25,dr50,dc50,dr25,dc25
	))

return(NA)

cat("Full:\n");print(full)
cat("IR 250:\n");print(ir250)
cat("IC 250:\n");print(ic250)
cat("IR 100:\n");print(ir100)
cat("IC 100:\n");print(ic100)
cat("IR 50:\n");print(ir50)
cat("IC 50:\n");print(ic50)
cat("IR 25:\n");print(ir25)
cat("IC 25:\n");print(ic25)
cat("DR 50:\n");print(dr50)
cat("DC 50:\n");print(dc50)
cat("DR 25:\n");print(dr25)
cat("DC 25:\n");print(dc25)

return(NA)

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

for (i in c(6)) {#,8,9)) {
	analyze_sims(i)
}
