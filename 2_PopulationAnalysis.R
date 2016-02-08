################################################################################
## PINK POPULATION ANALYSIS
################################################################################
## The following code is the second of two code supplements for the article:
## Cessation of a salmon decline with control of parasites
## Stephanie Peacock*, Martin Krkosek, Stan Proboscz, Craig Orr & Mark Lewis
## *stephanie.peacock@ualberta.ca
## -------------------------------------------------------------------------------
## This code is provided to make all treatment of pink salmon population data and 
## subsequent analyses transparent. The authors give full permission to use and 
## modify this code, but ask that appropriate aknowledgement be given.
################################################################################
################################################################################
rm(list=ls()) #Clear workspace

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
require(lme4)

mod.null<-lmer(Survival~S:Population+(1|Yr/Area), data=Z1, REML=TRUE)
mod.alt<-lmer(Survival~S:Population+WildLice+(1|Yr/Area), data=Z1, REML=TRUE)
anova(mod.null, mod.alt)

cat("Point estimates on parameters: \n r = ", fixef(mod.alt)[1], "\n c = ", fixef(mod.alt)[2])

#-------------------------------------------------------------------------------
# 4) Bootstrap confidence intervals on parameter estimates
#-------------------------------------------------------------------------------

Z1$ya<-integer(dim(Z1)[1])
for(i in 1:dim(ranef(mod.alt)$`Area:Yr`)[1]){
	j<-strsplit(dimnames(ranef(mod.alt)$`Area:Yr`)[[1]], ":")[[i]]
	Z1$ya[Z1$Yr==j[2]&Z1$Area==j[1]]<-i
	}
	
Z1$r<-integer(dim(Z1)[1])
for(i in 1:length(unique(Z1$Population))){
	j<-strsplit(names(fixef(mod.alt))[3:length(fixef(mod.alt))], "S:Population")[[i]]
	Z1$r[Z1$Population==j[2]]<-i
	}

#-------------------------------------------------------------------------------
# a) Function to bootstrap:
sim<-function(x){
	require(lme4)
	a<-params[[1]]
	b.i<-params[[2]]
	c<-params[[3]]
	sigma.ya<-params[[4]]
	sigma.y<-params[[5]]
	sigma.e<-params[[6]]
	
	set.seed(job.seeds[x,2]) #Set random seed for the job
	
	R<-numeric(length(Z1$S))	#Simulated Survival
	theta.y<-rnorm(length(levels(Z1$Yr)), 0, sigma.y)
	theta.ya<-rnorm(max(Z1$ya), 0, sigma.ya)
	epsilon<-rnorm(length(Z1$S), 0, sigma.e)
	
	R<-Z1$S*exp(a+b.i[Z1$r]*Z1$S+c*Z1$WildLice2+theta.ya[Z1$ya]+theta.y[as.numeric(Z1$Yr)]+epsilon)

	SS<-log(R/Z1$S)
	
	Z2<-data.frame(Area=as.factor(Z1$Area), Population=as.factor(Z1$Population), Yr=as.factor(Z1$Yr), S=Z1$S, Survival=SS, WildLice2=Z1$WildLice2)
	
	lice.mod<-lmer(Survival~S:Population+WildLice2+(1|Yr/Area), data=Z2, REML=TRUE)
		
	out.params<-fixef(lice.mod)[1:2]
	
	return(out.params)
	
} #end function
#-------------------------------------------------------------------------------

# b) Set up estimated variances for bootstrap algorithm
b <-as.numeric(fixef(mod.alt)[3:length(fixef(mod.alt))])
a<-as.numeric(fixef(mod.alt)[1])
c<-as.numeric(fixef(mod.alt)[2])

# Standard deviations!! Not variances!!
sigma.ya<-sqrt(VarCorr(mod.alt)$'Area:Yr'[1])
sigma.y<-sqrt(VarCorr(mod.alt)$Yr[1]) 
sigma.e<-attr(VarCorr(mod.alt), "sc")

params<-list(a, b, c, sigma.ya, sigma.y, sigma.e)	

#---------------------------------------------------
# c) Parallel computation (uncomment to run)
#    **************************WARNING**************************
#	 	This takes approx. 30 mins with n.cores=3 on my MacBook Pro, 
#		so it may tie up your computer for a while. 
#	************************************************************ 
#---------------------------------------------------

# require(parallel)
# detectCores() #Number of cores available
# n.cores<-3
# if(n.cores>detectCores()) cat("Attempting to parallelize on more cores than you have available!  Set n.cores < detectCores()")

# #Figure out random number sequences for different chains
# n.jobs<-1000 # How many iterations to bootstrap
# RNGkind("L'Ecuyer-CMRG")
# set.seed(1234)
# job.seeds<-matrix(nrow=n.jobs, ncol=7)
# job.seeds[1,]<-.Random.seed
# for(i in 2:n.jobs) job.seeds[i,]<-nextRNGStream(job.seeds[i-1,])
	
# t0<-proc.time()
# cl<-makeCluster(n.cores)
# clusterExport(cl, varlist=list("job.seeds", "Z1", "params"))
# X<-clusterApply(cl, x=c(1:1000), fun=sim)
# stopCluster(cl)
# cat("Process time (minutes) = ", (proc.time()-t0)[3]/60)
	
# # Unlist results
# p.all2<-matrix(nrow=1000, ncol=2)
# for(i in 1:1000) p.all2[i,]<-as.numeric(X[[i]])

# CI<-apply(p.all2, 2, quantile, c(0.025, 0.975))
# CI<-rbind(CI[1,], fixef(mod.alt)[1:2], CI[2,])
# colnames(CI)<-c("r", "c")
# rownames(CI)<-c("2.5%", "MLE", "97.5%")
# print(CI)

#---------------------------------------------------
# d) Computation of the percent mortality due to sea lice
#---------------------------------------------------
# Parameter confidenece intervals 
CI<-matrix(c(0.873361,1.088313,1.301514,-0.29934107,-0.19028109,-0.08669928), nrow=3, ncol=2)
colnames(CI)<-c("r", "c")
rownames(CI)<-c("2.5%", "MLE", "97.5%")

# Wild lice estimates
W<-data.frame(return.year=c(2001:2009)+1, WildLice=c(12.1739692,6.228714679,0.692532577,6.226036908,2.657006199,0.906920939,0.869532124,0.390081339,0.202908529))

mort<-cbind(CI[2,2]*W$WildLice, CI[1,2]*W$WildLice, CI[3,2]*W$WildLice)
p.mort<-100*(1-exp(mort))
colnames(p.mort)<-c("MLE", "97.5%", "2.5%")
rownames(p.mort)<-c(2002:2010)

#---------------------------------------------------
# e) Plotting figure 6 
#---------------------------------------------------
#install.packages("gplots")
require(gplots)
par(mfrow=c(2,1), mar=c(5,5,2,2), oma=c(0,0,0,0), cex=0.9)

# a) Pink salmon survival
plot(Z1$S[Z1$Area!=12]*10^(-6), Z1$Survival[Z1$Area!=12], col=c(grey(0.7)), ylim=c(-7,8), bty="l", las=1, ylab="Pink salmon survival", xlab=expression(paste("Spawner abundance (x",10^6,")")), pch=8)

points(Z1$S[Z1$Area==12&Z1$Yr.numeric<2002]*10^(-6), Z1$Survival[Z1$Area==12&Z1$Yr.numeric<2002], col=grey(0.7), pch=8)
points(Z1$S[Z1$Area==12&Z1$Yr.numeric>=2002]*10^(-6), Z1$Survival[Z1$Area==12&Z1$Yr.numeric>=2002], pch=21, bg="white")
points(Z1$S[Z1$Area==12&Z1$Yr.numeric==2004]*10^(-6), Z1$Survival[Z1$Area==12&Z1$Yr.numeric==2004], pch=19)

mtext(side=3, line=1, adj=0, "(a)")

# b) Estimated percent mortality
plotCI(2002:2010, p.mort[1:9,1], ui=p.mort[1:9,2], li=p.mort[1:9,3], bty="n", xlab="Return year", ylab="Estimated percent mortality", gap=0.3, las=1)
points(2004, p.mort[3,1], pch=19)
mtext(side=3, line=1, adj=0, "(b)")
