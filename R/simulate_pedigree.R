

#' @title simulate_pedigree
#' @description Individual based simulation based on specified demographic parameters.
#' @param years number of time steps
#' @param n_females starting number of breeding females
#' @param afr age at first reproduction
#' @param p_breed probability that a female breeds
#' @param fecundity number of juveniles a female produces each year
#' @param fixed_fecundity logical. is fecundity fixed or drawn from a Poisson distribution
#' @param juv_surv survival of juveniles until local recruitment, where recruitment is defined as having genetic offspring
#' @param adult_surv survival of adults across years
#' @param immigration yearly immigration, as a proportion of starting number of females (n_females)
#' @param p_polyandry probability that a female has any polyandry
#' @param p_sire probability that 'social' male sires all offspring 
#' @param p_retain probability that social partnership is retained
#' @param	constant_pop Logical. Should there be stochastic variation in population size
#' @param known_age_structure Currently not in use
#' @param verbose logical - print simulation progress? useful for debugging


#' @details ...
#' 
#' @author Joel Pick - joel.l.pick@gmail.com
#' @return A list with two elements: population data and a pedigree
#' @examples
#' \dontrun{
#' }
#' @export


## what do we output - how do we incorporate alive individuals that don't breed or multiple measurements
## output pedigree and data structure?

	# years = 5
	# n_females = 10

	# p_breed = 1
	# juv_surv = 0.25
	# adult_surv = 0.5
	# immigration = 0
	# afr=1

	# fecundity = 4
	# p_sire = 1
	# p_retain = 0.8
	# polgyny_rate = 0

	# constant_pop = TRUE 
	# known_age_structure = FALSE

##  todo
# sex specific rates
# juvenile survival

simulate_pedigree <- function(
	years = 5,
	n_females = 50,
	afr=1,
	p_breed = 1,
	fecundity = 4,
	fixed_fecundity = TRUE,
	p_polyandry = 0,
	p_sire = 1, 
	p_retain = 0, #
	juv_surv = 0.25,
	adult_surv = 0.5,
	immigration = 0,
	# polgyny_rate = 0,
	constant_pop = TRUE ,
	known_age_structure = FALSE,
	verbose =FALSE){

  options(stringsAsFactors=FALSE) # as long as later version of R - dont need


	## get environment, for transform det function
	Renv <- environment()
	
	# check and generate all sex specific variables

	general_check("afr", env=Renv, rate=FALSE, sex_specific=TRUE)

	lapply(c("p_breed", "adult_surv","juv_surv","immigration"),general_check, env=Renv, rate=TRUE, sex_specific=TRUE)

	lapply(c("years", "n_females", "fecundity"), general_check, env=Renv, rate=FALSE, sex_specific=FALSE)

	lapply(c("p_retain", "p_sire"), general_check, env=Renv, rate=TRUE, sex_specific=FALSE)

	det_growth_rate_f <- (p_breed_f * juv_surv_f * fecundity)/2 + adult_surv_f + immigration_f

	# male growth rate probably doesn't matter as long as all eggs are fertilised? although at some point would run out of males
	# det_growth_rate_m <- (p_breed_m *juv_surv_m * fecundity)/2 + adult_surv_m + immigration_m

	if(det_growth_rate_f>1) warning("growth rate is more than 1") 
  if(det_growth_rate_f<1) warning("growth rate is less than 1") 

	# v_as <- adult_surv * (1-adult_surv)
	# v_imm <- immigration * (1-immigration)
	# v_js <- juv_surv * (1-juv_surv)
	# v_fec <- fecundity
	# v_rec <- fecundity^2*v_js + juv_surv^2*v_fec + v_js*v_fec

	# stoch_growth_rate <- det_growth_rate - (v_as + v_rec + v_imm)/(2*n_females)

	# immigration <- 1- stoch_growth_rate


###
# MAKE STARTING POPULATION
###

	# make pedigree for base population
	# pedigree <- data.frame(
	# 	animal = paste0("0_",1:(n_females*2)),
	# 	dam = NA,
	# 	sire = NA,
	# 	sex = rep(c("F","M"),each=n_females),
	# 	## starting age structure for female and males, based on constant survival rate
	# 	## we could delete age/cohort of these individuals in the output as that would be realistic to a real pedigree?
	# 	cohort= -1 * c(
	# 		if(adult_surv_f==0){
	# 			rep(0,n_females)
	# 	  }else{
	# 		  rgeom(n_females,adult_surv_f)
	# 	  },
	# 	  if(adult_surv_m==0){
	# 	  	rep(0,n_females)
	# 	  }else{
	# 	  	rgeom(n_females,adult_surv_m)
	# 	  })
	# )

## this needs to be sex specific with JS
 	yearly_recruits <- n_females*p_breed*fecundity*juv_surv/2
	starting_n <- c(n_females,rep(yearly_recruits,afr-1))

	pedigree <- data.frame(
		animal = paste0("0_",1:sum(starting_n*2)),# or could code them with their cohort. using 0 means all founder are the same
		dam = NA,
		sire = NA,
		# equal sex ratio to start
		sex = rep(c("F","M"),sum(starting_n)),
		## starting age structure for female and males, based on constant survival rate
		## we could delete age/cohort of these individuals in the output as that would be realistic to a real pedigree?
		cohort= rep(-1*(afr-1):0,c(n_females,rep(yearly_recruits,afr-1))*2)
		## need to start sims with some age structure, otherwise with afr>1 the population size will drop after first year, as new offspring wont have recruited yet. So this way has founder that are recruited for the first afr-1 years. the new number of recruits per year is n_females*p_breed*fecundity*juv_surv for each sex
	)
## could make some stochasticity in the starting number if constant_pop=FALSE


	# make list that stores who is available to breed(or alive??) in each year
	dat <- list()
	# dat[[1]] <-  data.frame(animal = pedigree$animal, sex = pedigree$sex, age=NA, year=1)
	dat[[1]] <-  data.frame(
		animal = pedigree$animal, 
		sex = pedigree$sex, 
		age= 1-pedigree$cohort, 
		year= 1
		)

	if(verbose) cat("starting pop created \n")


	## stores male-female pairings across years
	pairs <- list()

	# year=1
	for(year in 1:years){
		## import individuals that are around as adults pre-breeding (dat[[year]])

	if(verbose) cat("year ",year,": ")

	####
	# BREEDING FEMALE AND MALE, AND PAIRING
	####

		# get vectors of females and males available to breed, accounting for the probability of females breeding
#		females <- subset(dat[[year]],sex=="F")$animal
		females <- subset(dat[[year]],sex=="F" & age>=afr_f)$animal
		## maybe just dont include these individuals in dat?

		breeding_females <- females[as.logical(stats::rbinom(length(females),1,p_breed_f))]

		# dat[[year]][(match(females,dat[[year]]$animal)),"age"]

		#work out number of pairs that can be formed 
		n_pair <- length(breeding_females)#min(length(breeding_females),length(males))
		# print(length(males))
		# print(length(breeding_females))

		
		males <- subset(dat[[year]],sex=="M"& age>=afr_m)$animal
		breeding_males <- males[as.logical(stats::rbinom(length(males),1,p_breed_m))]

		if(verbose) cat("Breeding Individuals, ")
### maybe need something that says that if no breeding female skip offspring creation?
		
		### assign 'social' male
		if(p_retain==0 || adult_surv==0 || year==1){ 
			## when there is no pairing to be done, just take random males, :
			social_male <- sample_males(breeding_females,breeding_males)
		}else{
			social_male<-sapply(breeding_females,function(bf){
				## MAKE VECTOR OF THOSE PAIRING UP AGAIN, then assign remaining males 
				if(bf %in% pairs[[year-1]]$female
				# female bred in the last year 
				& pairs[[year-1]]$male[match(bf,pairs[[year-1]]$female)] %in% breeding_males
				# and paired male is breeding this year
				& stats::rbinom(1,1,p_retain)==1
				# and they retain each other
				 ){
				  pairs[[year-1]]$male[match(bf,pairs[[year-1]]$female)]
				}else{
					NA
				}
			})
			social_male[is.na(social_male)] <- sample_males(
				females = breeding_females[is.na(social_male)], 
				males = breeding_males, 
				unpaired_males = breeding_males[!breeding_males%in%social_male]
				)
			# for the ones that havent been assigned, sample the males. 
			#there is a potential issue here if the number of females is larger than the males, then the males being chosen to sire the unpaired females will get lots of pairings, but the ones already paired wont get any extra.
		}
		## what do do with probability of breeding and mate retention?! If mate retention is 1, and one member of the pair doesnt breed once, then they will end up swapping, so there will be some level of divorce
		
		
# if female existed last year, who was male. else sample new male
		# if male alive, then rbinom(1,1,p_retain), 
		  # if retain then male, 
		  # else sample new male
		pairs[[year]] <- data.frame(
			female = breeding_females,
			male = social_male
			)		

		if(verbose) cat("Breeding pairs created, ")

	####
	# FEMALE FECUNDITY
	####

		# number of offspring per female
		n_juv <- if(fixed_fecundity) {
			rep(fecundity,n_pair)
		}else{
			stats::rpois(n_pair,fecundity)
		}
		
		if(verbose) cat("Eggs created, ")


	####
	# PATERNITY AND PEDIGREE CREATION
	# CREATION OF NEW OFFSPRING
	####

### for each female
### start with 'social' male 
### how many eggs does he sire of total, with probability p_polyandry
### randomly choose a different male, with probability p_polyandry
### how many of the remaining eggs does he sire
### continue until all eggs are sired		

# print(year)
		## make ped incorporating EPP and fecundity
		ped <- if(sum(n_juv)>0){
			data.frame(
				animal=paste0(year,"_",1:sum(n_juv)),
				# fecundity
				dam=rep(breeding_females,n_juv),
				# EPP
				### can make this more efficient - if p_polyandry=0 then rep(social_male,n_juv)
				#if p_polyandry==1 & p_sire==0 then sample(males,n_juv*length(breeding_females))
				sire=c(
					lapply(1:n_pair,function(i){
					## probability of any EPP
					polyandry <- stats::rbinom(1,1,p_polyandry)
					if(polyandry){
						## if there is EPP, how much
						## this is calculated by sampling how many of the offspring the paired male sired, and then giving the same probability to subsequent males. This means that extra pair males will be few, and have several offspring if p_sire if high - think this is more realistic
						if(p_sire==0) {
							sample(breeding_males,n_juv[i])
						}else{
							within_sires <- fill_sires(n_juv[i],p_sire)
							if(length(within_sires)>0){
								c(social_male[i], sample(breeding_males,max(within_sires)-1,replace=TRUE))[within_sires]
								## have put replace=TRUE as if there are few males this might not work - some males might get chosen twice, but this will likely only happen when N is low, and so isnt unrealistic anyway
							}else{ NULL }
						}
						# n_sired <- rbinom(1,n_juv[i],p_sire)
						# c(rep(social_male[i], n_sired), sample(males,n_juv[i]-n_sired,replace=TRUE))
					}else{
						rep(social_male[i],n_juv[i])
					}
				})
					, recursive=TRUE),
				#equal sex ratio
				sex=sample(c("M","F"),sum(n_juv),replace=TRUE),
				cohort=year
			)
		}else{
			NULL
		}
		# print(table(aggregate(animal~dam+cohort,pedigree, function(x)length(x))[,3]))
		# print(table(aggregate(animal~dam+cohort,pedigree, function(x)length(x))[,2]))
		

	if(verbose) cat("Offspring created, ")



		## create individuals present in the next year

	##----------
	## Juvenile Survival
	##----------
	
		if(is.null(ped)){
			next_year_juvF <- NULL
			next_year_juvM <- NULL
		}else{
			pedM <- ped[ped[,"sex"]=="M",]
			pedF <- ped[ped[,"sex"]=="F",]

			next_year_juvF <- if(nrow(pedF)==0){
				NULL	
			}else if(constant_pop){
				sample(pedF[,"animal"], round(juv_surv_f*fecundity*n_females/2), replace=FALSE)
			}else{
				pedF[as.logical(stats::rbinom(nrow(pedF),1,juv_surv_f)),"animal"]
			}

			
			next_year_juvM <- if(nrow(pedM)==0){
				NULL
			}else if(constant_pop){
				### need to ensure equal sex ratio of recruits, otherwise population size fluctuations
				sample(pedM[,"animal"], round(juv_surv_m*fecundity*n_females/2), replace=FALSE)

			}else{
				pedM[as.logical(stats::rbinom(nrow(pedM),1,juv_surv_m)),"animal"]

			}
		}
		
		if(verbose) cat("Juvenile survival, ")


		## or could take the new recruits from the whole pedigree, subset by cohort = year-afr - maybe this is better, because otherwise dat always has a load of pre afr individuals in it

	##----------
	## Adult Survival
	##----------

		next_year_AF <- if(adult_surv_f==0){
			NULL
		}else if(constant_pop){
			sample(females, round(adult_surv_f*n_females), replace=FALSE)
		}else{
			females[as.logical(stats::rbinom(length(females),1,adult_surv_f))]
		}
		
		next_year_AM <- if(adult_surv_m==0){
	    NULL
		}else if(constant_pop){
			sample(males, round(adult_surv_m*n_females), replace=FALSE)
		}else{
			males[as.logical(stats::rbinom(length(males),1,adult_surv_m))]
		}
		
	if(verbose) cat("Adult survival, ")

	##----------
	## IMMIGRATION
	##----------

		## need to sort out with constant pop
		## might be worth giving the immigrants age=afr, so cohort=year-afr

		n_imm <- if(constant_pop){
			c(immigration_f,immigration_m)*n_females
		}else{
			stats::rbinom(2,n_females,c(immigration_f,immigration_m))
		}
		
		imm_females <- if(n_imm[1]>0){
			data.frame(
				animal=paste(year+1,"IF",seq_len(n_imm[1]),sep="_"),
				dam=NA,
				sire=NA,
				sex="F",#rep(c("F","M"),c(n_imm)),
				cohort=year-afr +1
				)
			# paste(year+1,"IF",seq_len(n_imm[1]),sep="_")	
		}else{NULL}
		
		imm_males <- if(n_imm[2]>0){
			data.frame(
				animal=paste(year+1,"IM",seq_len(n_imm[2]),sep="_"),
				dam=NA,
				sire=NA,
				sex="M",#rep(c("F","M"),c(n_imm)),
				cohort=year-afr +1
				)
			# paste(year+1,"IM",seq_len(n_imm[2]),sep="_")
		}else{NULL}


		immigrants <- if(year!=years){
			rbind(imm_females,imm_males)
			# data.frame(
			# 	animal=c(imm_females,imm_males),
			# 	dam=NA,
			# 	sire=NA,
			# 	sex=rep(c("F","M"),c(n_imm)),
			# 	cohort=year-afr +1
			# 	)
		}else{
			NULL
		}

	if(verbose) cat("immigration, ")


	##----------
	## OUTPUTS
	##----------
		
		#update pedigree
		pedigree <- rbind(pedigree,ped,immigrants)

		# put together all individuals surviving to the next time step
		next_year_ind <- data.frame(
			animal = c(next_year_juvF,
			 	next_year_juvM,
			 	next_year_AF,
			 	next_year_AM,
			 	immigrants[,"animal"]),
			sex= c(rep(c("F","M","F","M"),c(length(next_year_juvF),length(next_year_juvM),length(next_year_AF),length(next_year_AM))),immigrants[,"sex"])
		)

	  if(nrow(next_year_ind)==0){
	  	message("Population went extinct in year ", year+1)
	  	break
	  }

		## maybe we don't need this - redundant with cohort?
	  next_year_ind$age <- (year+1) - pedigree[match(next_year_ind$animal,pedigree$animal),"cohort"]
	  next_year_ind$year <- year+1
	  dat[[year+1]] <- next_year_ind

	  
		
	  if(length(unique(next_year_ind$sex))==1){
	  	if(unique(next_year_ind$sex)=="F"){
	  		message("Males went extinct in year ", year+1)	
	  	}
	  	if(unique(next_year_ind$sex)=="M"){
	  		message("Females went extinct in year ", year+1)	
	  	}
	  	break
	  }
	  if(length(unique(next_year_ind$sex))==0){
	  	message("Population went extinct in year ", year+1)
	  	break
	  }
	if(verbose) cat("\n")

	}

	return(list(pedigree=pedigree,data_str=do.call(rbind,dat)))

}

