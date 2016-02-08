################################################################################
## DATA COMPILATION
################################################################################
## The following code is the first of two code supplements for the article:
## Cessation of a salmon decline with control of parasites
## Stephanie Peacock*, Martin Krkosek, Stan Proboscz, Craig Orr & Mark Lewis
## *stephanie.peacock@ualberta.ca
## -------------------------------------------------------------------------------
## This code is provided to make all treatment of pink salmon population data and 
## subsequent analyses transparent. The authors give full permission to use and 
## modify this code, but ask that appropriate aknowledgement be given.
################################################################################
################################################################################

################################################################################
## A. Calculation of exploitation rates using Pmax technique
################################################################################

#-------------------------------------------------------------------------------
# 1) Read in data
#-------------------------------------------------------------------------------
z<-read.csv("nuSEDs_PINK.csv")

Yr<-rep(1950:2010, length(unique(z$River)))
River<-rep(unique(z$River), each=length(1950:2010))
Escapement<-numeric(length(Yr)); Area<-numeric(length(Yr)); L<-numeric(length(Yr))
for(i in 1:length(Yr)){
	L[i]<-length(z[z$Yr==Yr[i]&z$River==River[i],1])
	# If there is no escapement value, fill in NA
	if(L[i]==0) Escapement[i]<-NA 
	# If there is a single, unique value, use it!
	if(L[i]==1) Escapement[i]<-z$Escapement[z$Yr==Yr[i]&z$River==River[i]]
	# If there are more than one value, then take the sum (usually one is NA)
	if(L[i]==2) Escapement[i]<-sum(z$Escapement[z$Yr==Yr[i]&z$River==River[i]], na.rm=TRUE)
	Area[i]<-z$Area[z$River==River[i]][1]
}

z<-data.frame(Area, River=as.numeric(River), RiverName=River, Yr, Escapement)

cat("Total number of rivers in nuSEDS data: ", length(unique(z$River)))

#-------------------------------------------------------------------------------
# 2) Calculate adjusted escapement using Pmax (Appendix A)
#-------------------------------------------------------------------------------

zz<-data.frame(AdjustedEsc=numeric(length(unique(Area))*length(1950:2010)))
zz$Area<-rep(unique(z$Area),each=length(1950:2010))
zz$Yr<-rep(c(1950:2010), length(unique(z$Area)))

for(a in 1:length(unique(z$Area))){ # For each area
	z1<-subset(z, Area==unique(z$Area)[a]) #subset data
	
	#choose index streams as those that have at least 20 years data
	ind<-as.numeric(tapply(is.na(z1$Escapement), z1$River, sum))
	index<-unique(z1$River)[which(ind<20)]
	
	#The observed sum is the sum of escapement for index streams each year
	obs.sum<-as.numeric(tapply(z1$Escapement[is.element(z1$River, index)], z1$Yr[is.element(z1$River, index)], sum, na.rm=TRUE))
	
	#The decade mean escapement for each stream 
	decade.mean.esc<-matrix(nrow=length(unique(z1$River)), ncol=7)
	decade.mean.esc[,1]<-as.numeric(tapply(z1$Escapement[z1$Yr<1960], z1$River[z1$Yr<1960], mean, na.rm=TRUE))
	decade.mean.esc[,2]<-as.numeric(tapply(z1$Escapement[z1$Yr>=1960&z1$Yr<1970], z1$River[z1$Yr>=1960&z1$Yr<1970], mean, na.rm=TRUE))
	decade.mean.esc[,3]<-as.numeric(tapply(z1$Escapement[z1$Yr>=1970&z1$Yr<1980], z1$River[z1$Yr>=1970&z1$Yr<1980], mean, na.rm=TRUE))
	decade.mean.esc[,4]<-as.numeric(tapply(z1$Escapement[z1$Yr>=1980&z1$Yr<1990], z1$River[z1$Yr>=1980&z1$Yr<1990], mean, na.rm=TRUE))
	decade.mean.esc[,5]<-as.numeric(tapply(z1$Escapement[z1$Yr>=1990&z1$Yr<2000], z1$River[z1$Yr>=1990&z1$Yr<2000], mean, na.rm=TRUE))
	decade.mean.esc[,6]<-as.numeric(tapply(z1$Escapement[z1$Yr>=2000], z1$River[z1$Yr>=2000], mean, na.rm=TRUE))
	decade.mean.esc[,7]<-as.numeric(tapply(z1$Escapement, z1$River, mean, na.rm=TRUE))
	
	dimnames(decade.mean.esc)<-list(unique(z1$River), c("1950s", "1960s", "1970s", "1980s", "1990s", "2000s", "all"))
	
	#Calculate the sum of escapement for each decade, from index streams only
	if(length(index)==1) decade.sum<-decade.mean.esc[which(is.element(unique(z1$River),index)),]
	if(length(index)>1) decade.sum<-apply(decade.mean.esc[which(is.element(unique(z1$River),index)),], 2, sum, na.rm=TRUE)
	
	#Porportion of decade.sum that each stream contributed, for each decade
	prop.obs.total<-matrix(nrow=length(index), ncol=6)
	for(i in 1:6){prop.obs.total[,i]<-decade.mean.esc[which(is.element(unique(z1$River),index)),i]/decade.sum[i]}
	apply(prop.obs.total, 2, sum) 
	
	#Calculate annual escapement contribution from each index stream
	esc.contribution<-matrix(nrow=length(index), ncol=length(1950:2010))
	for(i in 1:length(index)){
		for(y in 1:length(1950:2010)){
			if(sum(z1$Escapement[which(z1$River==index[i]&z1$Yr==c(1950:2010)[y])], na.rm=TRUE)==0){esc.contribution[i,y]<-0}else{
				if(y<=10) esc.contribution[i,y]<-prop.obs.total[i,1]
				if(y>10&y<=20) esc.contribution[i,y]<-prop.obs.total[i,2]
				if(y>20&y<=30) esc.contribution[i,y]<-prop.obs.total[i,3]
				if(y>30&y<=40) esc.contribution[i,y]<-prop.obs.total[i,4]
				if(y>40&y<=50) esc.contribution[i,y]<-prop.obs.total[i,5]
				if(y>50) esc.contribution[i,y]<-prop.obs.total[i,6]
				}
			}
		}
	
	#Determine expansion factor: how much of annual escapement for index streams is NA (i.e., infilling of index streams)
	expansion<-apply(esc.contribution, 2, sum)
	
	#Adjusted total escapement of the area for the year
	adj.sum<-obs.sum/expansion
	
	#plot(1950:2010, obs.sum, "l", ylab="Escapement", xlab="Year", bty="n", main=paste("Area ", unique(z$Area)[a]))
	#lines(1950:2010, adj.sum, col=2)
	
	#Analysis of 1980s baseline period for index streams contribution
	if(length(index)>1) proportion<-as.numeric(apply(decade.mean.esc[which(is.element(unique(z1$River),index)),], 2, sum, na.rm=TRUE)/apply(decade.mean.esc, 2, sum, na.rm=TRUE)) else proportion<-decade.mean.esc[which(is.element(unique(z1$River),index)),]/apply(decade.mean.esc, 2, sum, na.rm=TRUE)
	
	proportion<-c(rep(proportion[1:6], each=10), proportion[6])
	
	#Observed sum of all streams (not just index) in area (i.e., infilling of non-index streams)
	obs.sum.all<-adj.sum/proportion
	#lines(1950:2010, obs.sum.all, col=3)
	
	#Brian S. multiplies the observed sum by a 1.5 expansion factor, presumably for the streams that are not monitored at all.  It would be goo to see how sensitive the results are to this.
	adj.sum.all<-obs.sum.all*1.5
	#lines(1950:2010, adj.sum.all, col=4)

	
	zz$AdjustedEsc[zz$Area==unique(zz$Area)[a]]<-adj.sum.all
	
	} #end area

#-------------------------------------------------------------------------------
# 3) Calculate exploitation rates
#-------------------------------------------------------------------------------
x<-read.csv("Catch_PINK.csv")

zz$Catch<-numeric(dim(zz)[1])
for(a in unique(zz$Area)){
	for(j in 1950:2010){
		if(length(x$Catch[x$Area==a&x$Year==j])>0){
			zz$Catch[zz$Area==a&zz$Yr==j]<-x$Catch[x$Area==a&x$Year==j] 
			}else{
			zz$Catch[zz$Area==a&zz$Yr==j]<-NA
			}
		}
	}
zz$Exploitation<-zz$Catch/(zz$AdjustedEsc+zz$Catch)

#-------------------------------------------------------------------------------
# 4) Substitute DFO estimates for Area 12 (includes estimates of #fish caught 
#    that were returning to more southern rivers (e.g., the Fraser River))
#-------------------------------------------------------------------------------
Area12Exploitation<-data.frame(Exploitation=c(NA, NA, NA, NA,0.382085422, NA, 0.577450409,NA, 0.608614425, 0.609438028,0.581290831,0.624488262,0.511452327,0.53910599, 0.559106731, 0.407737058,0.709759091,0.68330814,0.702556339,0.355568067, 0.661001531, 0.417602208,0.495093597,0.466323627,0.557995733,0.307493388, 0.69227994, 0.642881521,0.535438691,0.597624993,0.445960806,0.445568795, 0.239768737, 0.439214304,0.363600223,0.521287118,0.267575097,0.503395014, 0.350754263, 0.724103169,0.664221614,0.481413686,0.315936987,0.518137595, 0.13366282, 0.337189177,0.018862085,0.316506018,0.031687856,0.001688128, 0.233618801, 0.150603921,0.108160486,0.074265328,0.06162374,0.071483675, 0.074, 0.011066875, 0.0262,0.049,0.129), ReturnYear=c(1950:2010))

for(j in 1950:2010){
	zz$Exploitation[zz$Yr==j&zz$Area==12]<-Area12Exploitation$Exploitation[Area12Exploitation$ReturnYear==j]
}
################################################################################
## B. Selection of "best" rivers and connection to louse data
################################################################################

#-------------------------------------------------------------------------------
# 1) Read in data
#-------------------------------------------------------------------------------
z<-read.csv("nuSEDS_PINK.csv")
x<-zz

#-------------------------------------------------------------------------------
# 2) Set up the overall S-R database 1960-2010, all Rivers
#-------------------------------------------------------------------------------

# a) Enter a row of NA for each missing year
L<-levels(z$River);
yr1<-c(); area1<-c(); river1<-c(); 

for(i in 1:length(L)){
	river<-L[i]; zz<-subset(z,River==river); area<-zz$Area[1];
	for(t in 1960:2010){yr1<-c(yr1,t); area1<-c(area1,area); river1<-c(river1,river)};
	}
	
# b) Create dataframe
Z<-data.frame(Yr=yr1, Area=area1, River=river1); Z$Escapement<-NA; Z$Expl<-NA;

# c) Add exploitation rates to dataframe
for(a in 7:12){
	for(y in 1960:2010){
		Z$Expl[which(Z$Yr==y & Z$Area==a)]<-x$Exploitation[which(x$Yr==y & x$Area==a)]
	}
}

# d) Add the nuSEDs escapement values (takes a minute)
for(i in 1:length(Z$Yr)){
	E<-z$Escapement[which(z$Yr==Z$Yr[i] & z$River==Z$River[i])]; 
	if(length(E)==1) Z$Escapement[i]<-E;
	if(length(E)>1) E<-E[which(E>0)]; if(length(E)!=0) Z$Escapement[i]<-E;
	}

# e) Replace 0s with NAs
Z$Escapement[which(Z$Escapement==0)]<-NA;	

# f) Recruitment estimates R = N/(1-u)
Z$R<-Z$Escapement/(1-Z$Expl)

# g) Spawner estimates to pair with recruitment (S(t-2) corresponds to R(t))
Z$S<-NA; Z$S[3:length(Z$S)]<-Z$Escapement[1:(length(Z$S)-2)];
Z$S[which(Z$Yr==1960 | Z$Yr==1961)]<-NA;

# h) Calculate survival as log(R/S)
Z$Survival<-log(Z$R/Z$S);

#-------------------------------------------------------------------------------
# 3) Remove ambiguous or confounding farm exposure, enhancement, etc
#-------------------------------------------------------------------------------

# Area 12
Z12<-subset(Z,Area==12)
N12<-unique(Z12$River)
N12
N.BA<-c("AHNUHATI RIVER", "AHTA RIVER", "GLENDALE CREEK", "KAKWEIKEN RIVER", "KINGCOME RIVER", "LULL CREEK", "VINER SOUND CREEK", "WAKEMAN RIVER");

BA<-c(); for(i in 1:length(N.BA)) BA<-rbind(BA,subset(Z,River==N.BA[i]))

# Area 7
Z7<-subset(Z,Area==7)
N7<-unique(Z7$River)
N7
N.7<-c("PINE RIVER", "NEEKAS CREEK", "TANKEEAH RIVER", "KWAKUSDIS RIVER", "BULLOCK CHANNEL CREEKS", "QUARTCHA CREEK", "LEE CREEK", "ROSCOE CREEK", "CLATSE CREEK", "WALKER LAKE CREEK", "GOAT BUSHU CREEK", "DEER PASS LAGOON CREEKS", "KUNSOOT RIVER", "KADJUSDIS RIVER", "MCLOUGHLIN CREEK", "COOPER INLET CREEKS");

A7<-c(); for(i in 1:length(N.7)) A7<-rbind(A7,subset(Z,River==N.7[i]))

ZZ<-subset(Z,Area!=12 & Area!=7); ZZ<-rbind(A7,ZZ,BA);

#-------------------------------------------------------------------------------
# 4) Specify odd and even year populations within the same river
#-------------------------------------------------------------------------------

ZZ$EO<-NA; ZZ$EO[which(ZZ$Yr %% 2 ==0)]<-"E"; ZZ$EO[which(ZZ$Yr %% 2 ==1)]<-"O";

# Populations
ZZ$Popn<-NA; N<-unique(ZZ$River); j<-1;
for(i in 1:length(N)){
	ZZ$Popn[which(ZZ$River==N[i] & ZZ$EO=="O")]<-j; j<-j+1;
	ZZ$Popn[which(ZZ$River==N[i] & ZZ$EO=="E")]<-j; j<-j+1;
	}
ZZ$Popn<-factor(ZZ$Popn);

#-------------------------------------------------------------------------------
# 4) Only keep populations with a minimum of 20 spawner-recruit pairs
#-------------------------------------------------------------------------------

Min.N<-20; Z.lengths<-c(); L<-unique(ZZ$Popn); 
for(i in 1:length(L)){
	z1<-subset(ZZ,Popn==L[i]);	
	z2<-which(is.na(z1$Survival)); 
	Z.lengths[i]<-length(z1$S)-length(z2);
	}

R.LongEnough<-which(Z.lengths>=20)
R.River<-c(); R.Area<-c();
for(i in 1:length(R.LongEnough)){
	z1<-subset(ZZ,Popn==R.LongEnough[i])
	R.Area<-c(R.Area,z1$Area[1])
	R.River<-c(R.River,z1$Popn[1])
	}

#Separate out even and odd year populations 
ZZ$Population<-as.factor(paste(ZZ$River, ZZ$EO, sep=" "))

#Only keep those entries with estimates of survival
Z1<-subset(ZZ, is.na(ZZ$Survival)==FALSE)

cat("Final dataset: \n Total number of populations (even/odd): ", length(unique(Z1$Population)), "\n Total number of S-R pairs: ", dim(Z1)[1], "\n Total number of rivers: ", length(unique(Z1$River)))

################################################################################
## C. Inclusion of louse covariate data
################################################################################

# Wild lice estimates
W<-data.frame(return.year=c(2001:2010)+1, mean=c(12.1739692,6.228714679,0.692532577,6.226036908,2.657006199,0.906920939,0.869532124,0.390081339,0.202908529,0.627256331))

Z1$WildLice<-rep(0,dim(Z1)[1])
for(i in 1:dim(W)[1]){
	Z1$WildLice[Z1$Area==12&Z1$Yr==W$return.year[i]]<-rep(W$mean[i], length(which(Z1$Area==12&Z1$Yr==W$return.year[i])))
}

################################################################################
## D. Write final database to file
################################################################################
write.csv(Z1, "StockRecruitData_PINK.csv")
