
fixedPed_rou <- read.csv("Data/Pedigree_rouviere_SQuID.txt", sep="", stringsAsFactors=TRUE)
fixedPed_city <- read.csv("Data/Pedigree_city_SQuID.txt", sep="", stringsAsFactors=TRUE)

ped=fixedPed_city
names(ped)[1]="animal"
names(ped)[4]="cohort"

#juvenile survival
chicks_to_recruit <- subset(ped,!is.na(cohort) & cohort!=max(cohort, na.rm=TRUE))[,1]
all_adults <- unique(c(ped[,"dam"],ped[,"sire"]))
js <- sum(chicks_to_recruit %in% all_adults) / length(chicks_to_recruit)
js

fem_rec=length(which(chicks_to_recruit %in% ped$dam))/length(which(chicks_to_recruit %in% all_adults))
hist(replicate(1000,mean(sample(c(0,1),length(which(chicks_to_recruit %in% all_adults)),replace=TRUE))),xlim=c(0,1),xlab="prob female",main=" ")#can this happen by chance?
arrows(fem_rec,0,fem_rec,1000,col='red')


#female probability of breeding (mean)
dam_years=unique(ped[!is.na(ped$dam),c("dam","cohort")])
juv_df <- subset(ped,!(is.na(cohort)))[,c("animal","cohort")]
dams=unique(dam_years$dam)

female_breeding <- lapply(dams, function(x){
  min_year <- if(x %in% juv_df$animal){
    juv_df[match(x,juv_df$animal),"cohort"] + 1 	
  }else{
    min(dam_years[dam_years$dam==x,"cohort"])
  }
  years <- min_year:max(dam_years[dam_years$dam==x,"cohort"])
  years_breed=length(dam_years[dam_years$dam==x,"cohort"])
  missing_years <-length(years[!years %in% dam_years[dam_years$dam==x,"cohort"]])
  data.frame("dam"=x,"prob_breed"=years_breed/(years_breed+missing_years))
})

female_breeding <-do.call("rbind",female_breeding)
mean(female_breeding$prob_breed)

#female fecundity (mean and variance)

female_rs <- aggregate(animal~dam+cohort,ped,function(x)length(unique(x)))
mean(female_rs[,3])
var(female_rs[,3])

#female juvenile survival (mean):Juv sex not available

#female age at first reproduction (mode?)

juv_df <- subset(ped,!(is.na(cohort)))[,c("animal","dam","cohort")]
fem_rec<- subset(juv_df,juv_df$animal %in% ped$dam)[,c("animal","cohort")]#year of birth
fem_rec=fem_rec[order(fem_rec$animal,fem_rec$cohort),]
dam_rec<-unique(subset(ped,ped$dam %in% juv_df$animal)[,c("dam","cohort")])#years breeding
dam_rec<-dam_rec[order(dam_rec$dam,dam_rec$cohort),]
dam_rec<-aggregate(cohort~dam,dam_rec,function(x) x[1])
female_breeding= data.frame("female"=fem_rec$animal,"age_1stbreed"=dam_rec$cohort-fem_rec$cohort)

as.numeric(dimnames(table(female_breeding$age_1stbreed))[[1]][which(tabulate(female_breeding$age_1stbreed)==max(tabulate(female_breeding$age_1stbreed)))])#mode

#female adult survival (mean)
cohorts<-sort(unique(ped$cohort))
as_year<- sapply(1:(length(cohorts)-1), function(i){
  ped_t0 <- subset(ped,cohort==cohorts[i])
  ped_t1 <- subset(ped,cohort %in% cohorts[cohorts>cohorts[i]])
  dams_t0 <- unique(ped_t0$dam)
  dams_t1 <- unique(ped_t1$dam)
  c(survived=sum(dams_t0 %in% dams_t1) ,total=length(dams_t0))
})
as_year
as_total <- rowSums(as_year)
as <- as_total[1]/as_total[2]
as
as_year[1,]/as_year[2,]

#female immigration (mean)
no_cohort=ped[which(is.na(ped$cohort)),]
cohorts<-sort(unique(ped$cohort))

im_year<- sapply(2:(length(cohorts)), function(i){
  dams<-unique(ped[which(ped$cohort==cohorts[i]& !is.na(ped$dam)),"dam"])
  n_fem=length(dams)
  dams=dams[which(dams %in% no_cohort$animal)]#remove recruits
  old_fems= unique(ped[which(ped$cohort %in% cohorts[cohorts<cohorts[i]]),"dam"])
  dams=dams[-which(dams %in% old_fems)]#remove older immigrants
  n_im_fem=length(dams)
  c(immgrant=n_im_fem ,total=n_fem)
})

im_total <- rowSums(im_year)
im <- im_total[1]/im_total[2]
im
im_year[1,]/im_year[2,]
