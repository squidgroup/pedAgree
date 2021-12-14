rm(list=ls())

require(MCMCglmm)
require(pedantics)

#ped<-read.csv("~/Dropbox/0_blue_tits/skew/Data/Intermediate/adult_pedigree.csv")
load("~/Dropbox/0_blue_tits/skew/Data/Intermediate/chick_data.Rdata")

head(ped)
head(THBW_noRep)

ped<-THBW_noRep[,c("bird_id","dam_P","sire_P","year")]
names(ped) <- c("animal","dam","sire","year")
ped$cohort <- as.numeric(as.character(ped$year))
head(ped)

# js <- sum(ped2[!(is.na(ped2[,2])&is.na(ped2[,3])),1] %in% unique(c(ped2[,2],ped2[,3]))) / length(ped2[!(is.na(ped2[,2])&is.na(ped2[,3])),1])
# # individuals caught as adults, are in pedigree, and so inflate chick survival
# # use just those with at least one parent





# rm(list=ls())

# install.packages("compoisson")

# out <- simulate_pedigree(years = 5, n_females = 200, constant_pop = TRUE )
# ped <- out$ped
# js <- 

# wo_last_year <- subset(ped, cohort!=max(cohort, na.rm=TRUE))


######
## ----- Juvenile survival
######

## in order not to downward bias survival too much, maybe first work out modal age at first repro, and then exclude that many years off the end?

chicks_to_recruit <- subset(ped,!(is.na(dam)&is.na(sire)) & cohort!=max(cohort, na.rm=TRUE))[,1]
all_adults <- unique(c(ped[,"dam"],ped[,"sire"]))

js <- sum(chicks_to_recruit %in% all_adults) / length(chicks_to_recruit)


######
## ----- Adult survival
######

cohorts<-sort(unique(ped$cohort))

## overall adult survival rate
as_year<- sapply(1:(length(cohorts)-1), function(i){
	ped_t0 <- subset(ped,cohort==cohorts[i])
	ped_t1 <- subset(ped,cohort %in% cohorts[cohorts>cohorts[i]])

	dams_t0 <- unique(ped_t0$dam)
	dams_t1 <- unique(ped_t1$dam)
	sires_t0 <- unique(ped_t0$sire)
	sires_t1 <- unique(ped_t1$sire)

	adults_t0 <- c(dams_t0,sires_t0)
	adults_t1 <- c(dams_t1, sires_t1)

	c(survived=sum(adults_t0 %in% adults_t1) ,total=length(adults_t0))
})
as_year
as_total <- rowSums(as_year)
as <- as_total[1]/as_total[2]
as_year[1,]/as_year[2,]



######
## ----- Probability of mating, Fecundity, EPP and fidelity
######


juv_df <- subset(ped,!(is.na(dam)&is.na(sire)))[,c("animal","cohort")]
juv_df$age <- 0


male_rs <- aggregate(cbind(animal,dam)~sire+cohort,ped,function(x)length(unique(x)))
names(male_rs)[3:4] <- c("fecundity","n_dams")
male_cohort <- lapply(split(male_rs, ~sire), function(x){
	min_year <- if(x$sire[1] %in% juv_df$animal){
		juv_df[match(x$sire[1],juv_df$animal),"cohort"] + 1 	
		##assumes that adults breed at age 1
	}else{
		min(x$cohort)
	}
	years <- min_year:max(x$cohort)
	missing_years <-years[!years %in% x$cohort]
	if(length(missing_years)>0){
		x<- rbind(x,data.frame(sire=unique(x$sire),cohort=missing_years,fecundity = 0, n_dams=0))
	}
	x[order(x$cohort),]

})


female_cohort<-lapply(split(subset(ped,!is.na(dam))[,c("dam","sire","cohort")], ~dam), function(x){
	z <- do.call(rbind,lapply(split(x,~cohort), function(y){
			males<-table(y$sire)
			data.frame(	dam=unique(y$dam),
				cohort=unique(y$cohort),
				social_male=sample(names(which(males==max(males))),1),
				## what if its split equally?? at the moment just sample random male with equal paternity
				n_males = length(males),
				prop_paternity = males[which(males==max(males))[1]]/nrow(y),
				fecundity = nrow(y)

			)
		}))
	min_year <- if(z$dam[1] %in% juv_df$animal){
		juv_df[match(z$dam[1],juv_df$animal),"cohort"] + 1 	
		##assumes that adults breed at age 1
	}else{
		min(z$cohort)
	}

	years <- min_year:max(z$cohort)
	missing_years <-years[!years %in% z$cohort]
	if(length(missing_years)>0){
		z<-rbind(z,data.frame(dam=unique(z$dam),cohort=missing_years,social_male=NA,n_males = NA,prop_paternity = NA,fecundity = 0))
	}
	z<- z[order(z$cohort),]

	z$last_male_alive <- NA
	z$last_male_kept <- NA

	if(nrow(z)>1){
		for(i in 2:nrow(z)){
			z$last_male_kept[i] <- z$social_male[i-1]==z$social_male[i]
			z$last_male_alive[i] <- z$cohort[i] %in% male_cohort[[z$social_male[i-1]]]$cohort
		}}
	z$divorce <- ifelse(!z$last_male_alive,NA, 
					ifelse(z$last_male_alive & z$last_male_kept,0,1)
					)
	z

})

f_df<-do.call(rbind,female_cohort)

f <- sum(f_df$fecundity)/sum(table(f_df$fecundity)[-1])

p_breed <- table(f_df$fecundity)[1]/sum(table(f_df$fecundity))

p_sire <- mean(f_df$prop_paternity, na.rm=TRUE)

p_retain <- table(f_df$divorce)[1]/sum(table(f_df$divorce))


# plot(table(f_df$n_males))
# hist(f_df$prop_paternity, breaks=50)


######
## ----- Immigration
######

1 - ((js*f)/2 + as)



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




