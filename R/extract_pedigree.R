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
## ----- Chick survival
######

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

ba<-split(ped,~dam+cohort)
social_male <- sapply(ba, function(x){
	males<-table(x$sire)
	c(dam=unique(x$dam),cohort=unique(x$cohort),social_male=names(which(males==max(males))))
	## what if its split equally??
	## add in total chicks and proportion from social male
})



table(ba[[3]]$sire)
ba[[20]]






