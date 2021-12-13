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
## ----- Fecundity
######

f <- mean(aggregate(animal~dam+cohort,ped, function(x)length(x))[,3])

######
## ----- Immigration
######

1 - ((js*f)/2 + as)

#var(aggregate(animal~dam+cohort,ped, function(x)length(x))[,3])


######
## ----- Probability of mating
######

juv_df <- subset(ped,!(is.na(dam)&is.na(sire)))[,c("animal","cohort")]
juv_df$age <- 0

# dam_df <- aggregate(animal~dam+cohort,ped,length)[1]

female_cohort<-lapply(split(subset(ped,!is.na(dam))[,c("dam","sire","cohort")], ~dam), function(x){
	ee1 <- do.call(rbind,lapply(split(x,~cohort), function(y){
			males<-table(y$sire)
			data.frame(	dam=unique(y$dam),
				cohort=unique(y$cohort),
				social_male=names(which(males==max(males)))[1],
				## what if its split equally?? at the moment just take first male in list
				n_males = length(males),
				prop_paternity = males[which(males==max(males))[1]]/nrow(y),
				fecundity = nrow(y)

			)
		}))

})

female_cohort[[100]] %in% juv_df$animal



## split for each dam, then split for each cohort, from this return fecundity, number of males 
## mate switch y/n, previous male alive y/n, 


######
## ----- EPP and fidelity
######


ba<-split(ped,~dam+cohort)
social_male <- lapply(ba, function(x){
	if(nrow(x)==0){NULL}else{
	males<-table(x$sire)
	c(dam=unique(x$dam),cohort=unique(x$cohort),social_male=names(which(males==max(males))))	
	}
	
	## what if its split equally??
	## add in total chicks and proportion from social male
})

table(ba[[3]]$sire)
ba[[20]]



######
## ----- male skew
######


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




