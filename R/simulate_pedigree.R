
## what do we output - how do we incorporate alive individuals that don't breed or multiple measurements
## output pedigree and data structure?

	# years = 5
	# n_females = 50
	# p_breed = 1
	# fecundity = 4
	# p_sire = 1
	# p_retain = 0.8
	# polgyny_rate = 0
	# juv_surv = 0.25
	# adult_surv = 0.5
	# immigration = 0
	# constant_pop = TRUE 
	# known_age_structure = FALSE

simulate_pedigree <- function(
	years = 5,
	n_females = 50,
	
	p_breed = 1, #probability that a female breeds
	fecundity = 4, #number of zygotes a female produces each year
	p_sire = 1, # probability that 'social' male sires all offspring 
	p_retain = 0, #probability that social partnership is retained
	# polgyny_rate = 0,
	juv_surv = 0.25,
	adult_surv = 0.5,
	immigration = 0,
	constant_pop = TRUE ,
	known_age_structure = FALSE){

options(stringsAsFactors=FALSE)

	# det_growth_rate <- (juv_surv * fecundity)/2 + adult_surv + immigration 

	# v_as <- adult_surv * (1-adult_surv)
	# v_js <- juv_surv * (1-juv_surv)
	# v_fec <- fecundity
	# v_rec <- fecundity^2*v_js + juv_surv^2*v_fec + v_js*v_fec

	# stoch_growth_rate <- det_growth_rate - (v_as + v_rec)/(2*n_females)

	# immigration <- 1- stoch_growth_rate


	# make pedigree for base population
	pedigree <- data.frame(
		animal = paste0("0_",1:(n_females*2)),
		dam = NA,
		sire = NA,
		sex = rep(c("F","M"),each=n_females),
		cohort=NA
	)

	# make list that stores who is alive in each year
	dat <- list()
	dat[[1]] <-  data.frame(animal = pedigree$animal, sex = pedigree$sex, age=NA)
	## think about using rgeom here = should it be different for year 1
	## do we start with age unknown as that would be realistic to the sampling?
	# plot(table(rgeom(1000,0.5)+1))
	# dgeom(0:10,0.5)

	pairs <- list()

	# year=2
	for(year in 1:years){
		## import individuals that are around as adults pre-breeding (dat[[year]])

		# probability of breeding just on females? assume that male breeding is dependent on females?

		# get vectors of females and males available to breed, accounting for the probability of females breeding
		females <- subset(dat[[year]],sex=="F")$animal
		breeding_females <- females[as.logical(rbinom(length(females),1,p_breed))]

		#work out number of pairs that can be formed 
		n_pair <- length(breeding_females)#min(length(breeding_females),length(males))
		# print(length(males))
		# print(length(breeding_females))


		# number of offspring per female
		# n_juv <- rpois(n_pair,fecundity)
		n_juv <- rep(fecundity,n_pair)


### for each female
### start with 'social' male from previous year
### how many eggs does he sire of total, with probability p_polyandry
### randomly choose a different male, with probability p_polyandry
### how many of the remaining eggs does he sire
### continue until all eggs are sired		

		
		males <- subset(dat[[year]],sex=="M")$animal

		
		### assign 'social' male
		if(year==1){
			## should this be with replacement?
			social_male<-sample(males,n_pair, replace=TRUE)	
		}else{
			social_male<-sapply(breeding_females,function(bf){
				if(bf %in% pairs[[year-1]]$female 
				& pairs[[year-1]]$male[match(bf,pairs[[year-1]]$female)] %in% males
				& rbinom(1,1,p_retain)==1 ){
				  pairs[[year-1]]$male[match(bf,pairs[[year-1]]$female)]
				}else{
					sample(males,1, replace=TRUE)
				}
			})
		}
		
		
# if female existed last year, who was male. else sample new male
		# if male alive, then rbinom(1,1,p_retain), 
		  # if retain then male, 
		  # else sample new male
		pairs[[year]] <- data.frame(
			female = breeding_females,
			male = social_male
			)		

		## make ped incorporating EPP and fecundity
		ped <- data.frame(
			animal=paste0(year,"_",1:sum(n_juv)),
			# fecundity
			dam=rep(breeding_females,n_juv),
			# EPP
			sire=c(lapply(1:n_pair,function(i){
					n_sired <- rbinom(1,n_juv[i],p_sire)
					c(rep(social_male[i], n_sired), sample(males,n_juv[i]-n_sired,replace=TRUE))
				}), recursive=TRUE),
			#equal sex ratio
			sex=sample(c("M","F"),sum(n_juv),replace=TRUE),
			cohort=year)
		pedigree <- rbind(pedigree,ped)

		# print(table(aggregate(animal~dam+cohort,pedigree, function(x)length(x))[,3]))
		# print(table(aggregate(animal~dam+cohort,pedigree, function(x)length(x))[,2]))

		## create individuals present in the next year
		next_year_ind <- if(constant_pop){
			rbind(
				### need to ensure equal sex ratio of recruits, otherwise population size fluctuations
				ped[sample(which(ped[,"sex"]=="F"), juv_surv*fecundity*n_females/2, replace=FALSE),c(1,4)],
				ped[sample(which(ped[,"sex"]=="M"), juv_surv*fecundity*n_females/2, replace=FALSE),c(1,4)],
				if(adult_surv!=0) cbind(animal=sample(females, adult_surv*n_females, replace=FALSE), sex="F"),
				if(adult_surv!=0) cbind(animal=sample(males, adult_surv*n_females, replace=FALSE), sex="M")
			)
		}else{
			rbind(
				ped[as.logical(rbinom(nrow(ped),1,juv_surv)),c(1,4)],
				cbind(animal=females[as.logical(rbinom(length(females),1,adult_surv))], sex="F"),
				cbind(animal=males[as.logical(rbinom(length(males),1,adult_surv))], sex="M")
			)
		}
	  next_year_ind$age<-(year+1) - pedigree[match(next_year_ind$animal,pedigree$animal),"cohort"]
	  dat[[year+1]] <- next_year_ind
	}

	return(list(pedigree=pedigree,data_str=do.call(rbind,dat)))

}
