rm(list=ls())

# require(MCMCglmm)
# require(pedantics)

source("~/github/pedigree_simulations/R/extract_pedigree.R")
source("~/github/pedigree_simulations/R/simulate_pedigree.R")

#ped<-read.csv("~/Dropbox/0_blue_tits/skew/Data/Intermediate/adult_pedigree.csv")
load("~/Dropbox/0_blue_tits/skew/Data/Intermediate/chick_data.Rdata")

head(ped)
head(THBW_noRep)

ped<-THBW_noRep[,c("bird_id","dam_P","sire_P","year")]
names(ped) <- c("animal","dam","sire","year")
ped$cohort <- as.numeric(as.character(ped$year))
head(ped)

extract_pedigree(ped)

# js <- sum(ped2[!(is.na(ped2[,2])&is.na(ped2[,3])),1] %in% unique(c(ped2[,2],ped2[,3]))) / length(ped2[!(is.na(ped2[,2])&is.na(ped2[,3])),1])
# # individuals caught as adults, are in pedigree, and so inflate chick survival
# # use just those with at least one parent


out<- simulate_pedigree(
	years = 5,
	n_females = 100,
	p_breed = 0.5, 
	fecundity = 4,
	p_sire = 0.8, 
	p_retain = 0.5,
	juv_surv = 0.2,
	adult_surv = 0.9,
	immigration = 0,
	constant_pop = TRUE ,
	known_age_structure = FALSE)

ped <- out$ped

extract_pedigree(ped)




# rm(list=ls())

# install.packages("compoisson")

# out <- simulate_pedigree(years = 5, n_females = 200, constant_pop = TRUE )
# ped <- out$ped
# js <- 

# wo_last_year <- subset(ped, cohort!=max(cohort, na.rm=TRUE))










######
## ----- male skew
######

k
# https://stats.stackexchange.com/questions/4659/relationship-between-binomial-and-beta-distributions

# a=k
# b= n-k+1
# n= popsize
# k= mean sirings


N=100
par(mfcol=c(2,5))
for(k in c(1,2,5,10)){
	x<-rbinom(10000,N,k/N)/(N)
	 hist(x,breaks=seq(0,1,0.01))
	print(sd(x)/mean(x))
	y<-rbeta(10000,k,N-k +1)
	hist(y,breaks=seq(0,1,0.01))
	print(sd(y)/mean(y))
}

N=100
par(mfcol=c(2,5))
for(k in c(1,2,5,10)){
	x<-rbinom(10000,N,k/N)
	hist(x)
	print(sd(x)/mean(x))
	y<-N*rbeta(10000,k,N-k +1)
	hist(y)
	print(sd(y)/mean(y))
}


cohort_split<-split(ped,~cohort)
sire_skew<-lapply(cohort_split,function(x){
	table(x$sire)/length(unique(x$sire))
})
### need to include any males that were known to be alive but didn't get any paternity


hist(c(sire_skew, recursive=TRUE), breaks=50, freq=FALSE)

beta_param<-EnvStats::ebeta(c(sire_skew, recursive=TRUE))$parameters
x<-seq(0,1,length.out=1000)
y<-dbeta(x, shape1=beta_param["shape1"], shape2=beta_param["shape2"])
lines(y~x)


## make a dataframe with all years alive for each individual 
## then match with some form of sire_skew to get offspring per year 


## make a dataframe with all years alive for each individual 
## then match with some form of sire_skew to get offspring per year 


