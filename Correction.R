################################################################################
## PINK POPULATION ANALYSIS
################################################################################
## The following code explores an error in the reported methods of Peacock et al.
## (2013) regarding the subsetting of populations with < 20 spawner-recruit pairs.
## Specifically, the paper reports results using all SR data, rather than the
## subsetted dataset as reported.
## Here, we explore the implications of this error for the results.
## Date: June 28, 2022
## Contact: stephanie.j.peacock at gmail.com
## -------------------------------------------------------------------------------
## This code is provided to make all treatment of pink salmon population data and 
## subsequent analyses transparent. The authors give full permission to use and 
## modify this code, but ask that appropriate acknowledgement be given.
################################################################################
################################################################################

library(parallel)
library(lme4)
library(gplots)

#-------------------------------------------------------------------------------
# 1) Read in data
# ** Note that you must run "1_DataCompilation.R" 
# to generate the datafile used in this analysis**
#-------------------------------------------------------------------------------

# Spawner-recruit data (from 1_DataCompilation.R)
Z<-read.csv("StockRecruitData_PINK.csv", header=TRUE)
if(names(Z)[1]=="X") Z<-Z[,2:dim(Z)[2]]

Z$Yr.numeric<-Z$Yr
Z$Yr<-as.factor(Z$Yr)
Z$Area<-as.factor(Z$Area)
Z$S<-as.numeric(Z$S)
Z$Popn<-as.factor(Z$Popn)

#-------------------------------------------------------------------------------
# 2) Setting W -> NA for return years 1991-2001
#-------------------------------------------------------------------------------
Z$WildLice2<-Z$WildLice
Z$WildLice2[is.element(Z$Yr.numeric, c(1991:2001))&Z$Area==12]<-NA
Z1<-subset(Z, is.na(WildLice2)==FALSE)

#-------------------------------------------------------------------------------
# 3) Model fitting
#-------------------------------------------------------------------------------
#install.packages("lme4")
library(lme4)

# Results as a list with
#  1) seen in paper, using all populations
#  2) using subsetted populations with >= 20 SR pairs

mod.null <- list(
	lmer(Survival~S:Population+(1|Yr/Area), data=Z1, REML=TRUE),
	lmer(Survival~S:Population+(1|Yr/Area), data=Z1[which(Z1$LongEnough == 1), ], REML=TRUE))

mod.alt <- list(
	lmer(Survival~S:Population+WildLice+(1|Yr/Area), data=Z1, REML=TRUE),
	lmer(Survival~S:Population+WildLice+(1|Yr/Area), data=Z1[which(Z1$LongEnough == 1), ], REML=TRUE))

anova(mod.null[[1]], mod.alt[[1]])
anova(mod.null[[2]], mod.alt[[2]])
# Effect of lice still significant

cat("Point estimates on parameters in paper: \n r = ", fixef(mod.alt[[1]])[1], "\n c = ", fixef(mod.alt[[1]])[2])

cat("Point estimates on parameters using subsetted data: \n r = ", fixef(mod.alt[[2]])[1], "\n c = ", fixef(mod.alt[[2]])[2])

# Minimal change in growth rate and c parameters

#-------------------------------------------------------------------------------
# 4) Bootstrap confidence intervals on parameter estimates
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# a) Function to bootstrap:
sim <- function(x){
	require(lme4)
	a <- params[[1]]
	b.i <- params[[2]]
	c <- params[[3]]
	sigma.ya<-params[[4]]
	sigma.y<-params[[5]]
	sigma.e<-params[[6]]
	
	set.seed(job.seeds[x,2]) #Set random seed for the job
	
	R<-numeric(length(Z$S))	#Simulated Survival
	theta.y<-rnorm(length(levels(Z$Yr)), 0, sigma.y)
	theta.ya<-rnorm(max(Z$ya), 0, sigma.ya)
	epsilon<-rnorm(length(Z$S), 0, sigma.e)
	
	R<-Z$S*exp(a+b.i[Z$r]*Z$S+c*Z$WildLice2+theta.ya[Z$ya]+theta.y[as.numeric(Z$Yr)]+epsilon)
	
	SS<-log(R/Z$S)
	
	Z2<-data.frame(Area=as.factor(Z$Area), Population=as.factor(Z$Population), Yr=as.factor(Z$Yr), S=Z$S, Survival=SS, WildLice2=Z$WildLice2)
	
	lice.mod <- lmer(Survival~S:Population+WildLice2+(1|Yr/Area), data=Z2, REML=TRUE)
	
	out.params <- fixef(lice.mod)[1:2]
	
	return(out.params)
	
} #end function

#------------------------------------------------------------------------------
# Apply to both full and subsetted datasets
CI <- list(); length(CI) <- 2
for(s in 1:2){
	
	if(s == 1) Z <- Z1 # full
	if(s == 2) Z <- Z1[Z1$LongEnough == 1, ] # subsetted
	
	Z$ya <- integer(dim(Z)[1])
	for(i in 1:dim(ranef(mod.alt[[s]])$`Area:Yr`)[1]){
		j <- strsplit(dimnames(ranef(mod.alt[[s]])$`Area:Yr`)[[1]], ":")[[i]]
		Z$ya[Z$Yr == j[2] & Z$Area == j[1]] <- i
	}
	
	Z$r<-integer(dim(Z)[1])
	for(i in 1:length(unique(Z$Population))){
		j<-strsplit(names(fixef(mod.alt[[s]]))[3:length(fixef(mod.alt[[s]]))], "S:Population")[[i]]
		Z$r[Z$Population==j[2]]<-i
	}
	
	
	# b) Set up estimated variances for bootstrap algorithm
	b <-as.numeric(fixef(mod.alt[[s]])[3:length(fixef(mod.alt[[s]]))])
	a <- as.numeric(fixef(mod.alt[[s]])[1])
	c <- as.numeric(fixef(mod.alt[[s]])[2])
	
	# Standard deviations!! Not variances!!
	sigma.ya <- sqrt(VarCorr(mod.alt[[s]])$'Area:Yr'[1])
	sigma.y <- sqrt(VarCorr(mod.alt[[s]])$Yr[1]) 
	sigma.e <- attr(VarCorr(mod.alt[[s]]), "sc")
	
	params <- list(a, b, c, sigma.ya, sigma.y, sigma.e)	
	
	#---------------------------------------------------
	# c) Parallel computation (uncomment to run)
	#    **************************WARNING**************************
	#	 	This takes approx. 30 mins with n.cores=3 on my MacBook Pro, 
	#		so it may tie up your computer for a while. 
	#	************************************************************ 
	#---------------------------------------------------
	
	detectCores() #Number of cores available
	n.cores <- 6
	if(n.cores > detectCores()) cat("Attempting to parallelize on more cores than you have available!  Set n.cores < detectCores()")
	
	#Figure out random number sequences for different chains
	n.jobs <- 1000 # How many iterations to bootstrap
	RNGkind("L'Ecuyer-CMRG")
	set.seed(1234)
	job.seeds <- matrix(nrow = n.jobs, ncol = 7)
	job.seeds[1,] <- .Random.seed
	for(i in 2:n.jobs) job.seeds[i,]<-nextRNGStream(job.seeds[i-1,])
	
	t0 <- proc.time()
	cl <- makeCluster(n.cores)
	clusterExport(cl, varlist=list("job.seeds", "Z", "params"))
	X <- clusterApply(cl, x = c(1:1000), fun = sim)
	stopCluster(cl)
	cat("Process time (minutes) = ", (proc.time()-t0)[3]/60)
	
	# Unlist results
	p.all2 <- matrix(nrow=1000, ncol=2)
	for(i in 1:1000) p.all2[i,] <- as.numeric(X[[i]])
	
	CI.s <- apply(p.all2, 2, quantile, c(0.025, 0.975))
	CI.s <- rbind(CI.s[1,], fixef(mod.alt[[s]])[1:2], CI.s[2,])
	colnames(CI.s)<-c("r", "c")
	rownames(CI.s)<-c("2.5%", "MLE", "97.5%")
	# print(CI)
	
	CI[[s]] <- CI.s
	
} # end s

# saveRDS(CI, file = "paramEstimates.rds")
#---------------------------------------------------
# d) Computation of the percent mortality due to sea lice
#---------------------------------------------------

# # Parameter confidenece intervals 
# CI <- readRDS("paramEstimates.rds")

# Wild lice estimates
W <- data.frame(return.year=c(2001:2009)+1, WildLice=c(12.1739692,6.228714679,0.692532577,6.226036908,2.657006199,0.906920939,0.869532124,0.390081339,0.202908529))

mort <- list(
	cbind(CI[[1]][2,2]*W$WildLice, CI[[1]][1,2]*W$WildLice, CI[[1]][3,2]*W$WildLice),
	cbind(CI[[2]][2,2]*W$WildLice, CI[[2]][1,2]*W$WildLice, CI[[2]][3,2]*W$WildLice))
	
p.mort <- list(
	100*(1-exp(mort[[1]])),
	100*(1-exp(mort[[2]])))
	
for(s in 1:2){
	colnames(p.mort[[s]])<-c("MLE", "97.5%", "2.5%")
	rownames(p.mort[[s]])<-c(2002:2010)
}

#---------------------------------------------------
# e) Plotting figure 6 
#---------------------------------------------------

quartz(width = 8, height = 6)
par(mfcol=c(2,2), mar=c(5,5,2,1), oma=c(0,0,2,0), cex=0.9)

for(s in 1:2){
	
	if(s == 1) Z <- Z1 # full
	if(s == 2) Z <- Z1[Z1$LongEnough == 1, ] # subsetted
	
	# a) Pink salmon survival
	plot(Z$S[Z$Area!=12]*10^(-6), Z$Survival[Z$Area!=12], col=c(grey(0.7)), ylim=c(-7,8), bty="l", las=1, ylab="Pink salmon survival", xlab=expression(paste("Spawner abundance (x",10^6,")")), pch=8)
	
	points(Z$S[Z$Area==12&Z$Yr.numeric<2002]*10^(-6), Z$Survival[Z$Area==12&Z$Yr.numeric<2002], col=grey(0.7), pch=8)
	points(Z$S[Z$Area==12&Z$Yr.numeric>=2002]*10^(-6), Z$Survival[Z$Area==12&Z$Yr.numeric>=2002], pch=21, bg="white")
	points(Z$S[Z$Area==12&Z$Yr.numeric==2004]*10^(-6), Z$Survival[Z$Area==12&Z$Yr.numeric==2004], pch=19)
	
	mtext(side=3, line=1, adj=0, c("(a)", "(c)")[s])
	mtext(side = 3, line = 2.5, c("Peacock et al. (2013)", "Revised subsetted data")[s])
	
	text(2.5, 7, paste("Rivers: ", length(unique(Z$River)), "\nPop'ns: ", length(unique(Z$Population)), "\nSR pairs: ", nrow(Z), sep = ""), adj = 0, cex = 0.9, xpd = NA)
	
	# b) Estimated percent mortality
	plotCI(2002:2010, p.mort[[s]][1:9,1], ui=p.mort[[s]][1:9,2], li=p.mort[[s]][1:9,3], bty="l", xlab="Return year", ylab="Estimated percent mortality", gap=0.3, las=1, yaxs = "i", ylim = c(0, 100))
	abline(h = seq(0, 100, 10), lty = 3, col = grey(0.8))
	points(2004, p.mort[[s]][3,1], pch=19)
	mtext(side=3, line=1, adj=0, c("(b)", "(d)")[s])
	
} # end s

################################################################################
# Compare datasets
################################################################################
lu <- function(x){length(unique(x))}

length(unique(Z1$River))
length(unique(Z1$River[Z1$LongEnough == 1]))

d <- data.frame(
	Area = c(as.character(unique(Z1$Area)), "all"),
	nRivers = c(
		paste0(lu(Z1$River[Z1$LongEnough == 1 & Z1$Area == 7]), " (", lu(Z1$River[Z1$Area == 7]), ")"),
		paste0(lu(Z1$River[Z1$LongEnough == 1 & Z1$Area == 8]), " (", lu(Z1$River[Z1$Area == 8]), ")"),
		paste0(lu(Z1$River[Z1$LongEnough == 1 & Z1$Area == 9]), " (", lu(Z1$River[Z1$Area == 9]), ")"),
		paste0(lu(Z1$River[Z1$LongEnough == 1 & Z1$Area == 10]), " (", lu(Z1$River[Z1$Area == 10]), ")"),
		paste0(lu(Z1$River[Z1$LongEnough == 1 & Z1$Area == 12]), " (", lu(Z1$River[Z1$Area == 12]), ")"),
		paste0(lu(Z1$River[Z1$LongEnough == 1]), " (", lu(Z1$River), ")")),
	
	nPopns = c(
		paste0(lu(Z1$Population[Z1$LongEnough == 1 & Z1$Area == 7]), " (", lu(Z1$Population[Z1$Area == 7]), ")"),
		paste0(lu(Z1$Population[Z1$LongEnough == 1 & Z1$Area == 8]), " (", lu(Z1$Population[Z1$Area == 8]), ")"),
		paste0(lu(Z1$Population[Z1$LongEnough == 1 & Z1$Area == 9]), " (", lu(Z1$Population[Z1$Area == 9]), ")"),
		paste0(lu(Z1$Population[Z1$LongEnough == 1 & Z1$Area == 10]), " (", lu(Z1$Population[Z1$Area == 10]), ")"),
		paste0(lu(Z1$Population[Z1$LongEnough == 1 & Z1$Area == 12]), " (", lu(Z1$Population[Z1$Area == 12]), ")"),
		paste0(lu(Z1$Population[Z1$LongEnough == 1]), " (", lu(Z1$Population), ")")),
	
	nData = c(
		paste0(length(Z1$Population[Z1$LongEnough == 1 & Z1$Area == 7]), " (", length(Z1$Population[Z1$Area == 7]), ")"),
		paste0(length(Z1$Population[Z1$LongEnough == 1 & Z1$Area == 8]), " (", length(Z1$Population[Z1$Area == 8]), ")"),
		paste0(length(Z1$Population[Z1$LongEnough == 1 & Z1$Area == 9]), " (", length(Z1$Population[Z1$Area == 9]), ")"),
		paste0(length(Z1$Population[Z1$LongEnough == 1 & Z1$Area == 10]), " (", length(Z1$Population[Z1$Area == 10]), ")"),
		paste0(length(Z1$Population[Z1$LongEnough == 1 & Z1$Area == 12]), " (", length(Z1$Population[Z1$Area == 12]), ")"),
		paste0(length(Z1$Population[Z1$LongEnough == 1]), " (", length(Z1$Population), ")"))
)

write.csv(d, file = "Summary_of_dataset_changes.csv")
unique(Z1$River[Z1$LongEnough == 1 & Z1$Area == 12])
unique(Z1$River[Z1$Area == 12])

# Viner Sound Creek removed
quartz(width = 3.2, height = 2.5, pointsize = 10)
par(mar = c(3,4,1,1))
plot(as.numeric(as.character(Z1$Yr[Z1$River == "VINER SOUND CREEK"])), Z1$Survival[Z1$River == "VINER SOUND CREEK"], "n", xlim = c(1960, 2010), xlab = "", ylab = "Survival (log R/S)", ylim = range(Z1$Survival), yaxt = "n")
axis(side = 2, at = c(-5, 0, 5), las = 1)
abline(h = 0, col = grey(0.8))
points(as.numeric(as.character(Z1$Yr[Z1$River == "VINER SOUND CREEK"])), Z1$Survival[Z1$River == "VINER SOUND CREEK"], pch = as.numeric(as.factor(Z1$EO[Z1$River == "VINER SOUND CREEK"])))
