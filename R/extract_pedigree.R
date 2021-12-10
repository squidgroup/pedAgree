rm(list=ls())

require(MCMCglmm)
require(pedantics)

#ped<-read.csv("~/Dropbox/0_blue_tits/skew/Data/Intermediate/adult_pedigree.csv")
load("~/Dropbox/0_blue_tits/skew/Data/Intermediate/chick_data.Rdata")

head(ped)
head(THBW_noRep)

ped2<-THBW_noRep[,c("bird_id","dam_P","sire_P","year")]

js <- sum(ped2[!(is.na(ped2[,2])&is.na(ped2[,3])),1] %in% unique(c(ped2[,2],ped2[,3]))) / length(ped2[!(is.na(ped2[,2])&is.na(ped2[,3])),1])
# individuals caught as adults, are in pedigree, and so inflate chick survival
# use just those with at least one parent





rm(list=ls())
out <- simulate_pedigree(years = 5, n_females = 200, constant_pop = TRUE )
ped <- out$ped
# js <- 

# wo_last_year <- subset(ped, cohort!=max(cohort, na.rm=TRUE))


chicks_to_recruit <- subset(ped,!(is.na(dam)&is.na(sire)) & cohort!=max(cohort, na.rm=TRUE))[,1]
all_adults <- unique(c(ped[,"dam"],ped[,"sire"]))

js <- sum(chicks_to_recruit %in% all_adults) / length(chicks_to_recruit)

unique(ped$cohort)

## overall adult survival rate
as_year<- sapply(1:4, function(i){
	adults_t0 <- unique(c(subset(ped,cohort==i)$dam,subset(ped,cohort==i)$sire))
	adults_t1 <- unique(c(subset(ped,cohort==i+1)$dam,subset(ped,cohort==i+1)$sire))
	c(survived=sum(adults_t0 %in% adults_t1),total=length(adults_t0))
})
as_total <- rowSums(as_year)
as <- as_total[1]/as_total[2]

mean(aggregate(animal~dam+cohort,ped, function(x)length(x))[,3])

var(aggregate(ped[,1]~ped[,2]+ped[,5],ped, function(x)length(x))[,3])




