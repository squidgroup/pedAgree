
is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol

between <- function(x,min,max) x>=min & x<=max

fill_sires <- function(n,p){
	n_remaining<-n
	sire <- 1
	sires <- rep(NA,length=n)
	while(n_remaining>0){
		n_sired <- (rbinom(1,n_remaining,p))
		if(n_sired>0){	 
			sires[(n-sum(is.na(sires)) +1):(n-sum(is.na(sires)) + n_sired)] <- sire
			n_remaining <- sum(is.na(sires))
				 sire <- sire+1}
	}
	sires
}

general_check <- function(name,env, rate=TRUE, sex_specific=TRUE){
	x <- get(name,pos=env)

	if(!is.vector(x)) stop(paste(name,"is not a vector"))
	
	if(rate){
		if(!all(between(x,0,1))) stop(paste(name,"must be between 0 and 1"))	
	}else{
		if(!all(is.wholenumber(x) & x>0)) stop(paste(name,"must be a whole number and >0"))	
	}
	
	
	if(sex_specific){
		if(!length(x) %in% 1:2) stop(paste(name,"should be length 1 or 2"))
	
		assign(paste0(name,"_f"),x[1],pos=env)
		if(length(x)==2){	
			assign(paste0(name,"_m"),x[2],pos=env)
		}else{
			assign(paste0(name,"_m"),x[1],pos=env)
		}
	}else{
		if(!length(x) == 1) stop(paste(name,"should be length 1"))
	}
}


#' @title simulated_variance
#' @description Calculate simulated mean and variance in response variable(s)
#' @param years number of time steps
#' @param n_females starting number of breeding females
#' @param afr age at first reproduction
#' @param p_breed probability that a female breeds
#' @param fecundity number of juveniles a female produces each year
#' 
#' @param juv_surv survival of juveniles until local recruitment, where recruitment is defined as having genetic offspring
#' @param adult_surv survival of adults across years
#' @param immigration yearly immigration, as a proportion of starting number of females (n_females)
#' @param p_polyandry probability that a female has any polyandry
#' @param p_sire probability that 'social' male sires all offspring 
#' @param p_retain probability that social partnership is retained

#' @details 
#' 
#' @author Joel Pick - joel.l.pick@gmail.com
#' @return A list with means and variance at each level
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
	juv_surv = 0.25,
	adult_surv = 0.5,
	immigration = 0,
	p_polyandry = 0,
	p_sire = 1, 
	p_retain = 0, #
	# polgyny_rate = 0,
	constant_pop = TRUE ,
	known_age_structure = FALSE){

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
	# det_growth_rate_m <- (juv_surv_m * fecundity)/2 + adult_surv_m + immigration_m

	if(det_growth_rate_f>1) warning("growth rate is more than 1") 
  if(det_growth_rate_f<1) warning("growth rate is less than 1") 

	# v_as <- adult_surv * (1-adult_surv)
	# v_imm <- immigration * (1-immigration)
	# v_js <- juv_surv * (1-juv_surv)
	# v_fec <- fecundity
	# v_rec <- fecundity^2*v_js + juv_surv^2*v_fec + v_js*v_fec

	# stoch_growth_rate <- det_growth_rate - (v_as + v_rec + v_imm)/(2*n_females)

	# immigration <- 1- stoch_growth_rate




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

	## stores male-female pairings across years
	pairs <- list()

	# year=1
	for(year in 1:years){
		## import individuals that are around as adults pre-breeding (dat[[year]])

		# probability of breeding just on females? assume that male breeding is dependent on females?

		# get vectors of females and males available to breed, accounting for the probability of females breeding
#		females <- subset(dat[[year]],sex=="F")$animal
		females <- subset(dat[[year]],sex=="F" & age>=afr_f)$animal
		## maybe just dont include these individuals in dat?

		breeding_females <- females[as.logical(rbinom(length(females),1,p_breed_f))]

		# dat[[year]][(match(females,dat[[year]]$animal)),"age"]

		#work out number of pairs that can be formed 
		n_pair <- length(breeding_females)#min(length(breeding_females),length(males))
		# print(length(males))
		# print(length(breeding_females))

		# number of offspring per female
		n_juv <- if(fixed_fecundity) {
			rep(fecundity,n_pair)
		}else{
			rpois(n_pair,fecundity)
		}
		
		males <- subset(dat[[year]],sex=="M"& age>=afr_m)$animal
		breeding_males <- males[as.logical(rbinom(length(males),1,p_breed_m))]

### for each female
### start with 'social' male from previous year
### how many eggs does he sire of total, with probability p_polyandry
### randomly choose a different male, with probability p_polyandry
### how many of the remaining eggs does he sire
### continue until all eggs are sired		

		
		### assign 'social' male
		if(year==1 || adult_surv==0){ # or adult_surv=0?
			## should this be with replacement?
			social_male<-sample(breeding_males, n_pair, replace=FALSE)	#TURE?
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
		## what do do with probability of breeding and mate retention?!
		
		
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
			### can make this more efficient - if p_polyandry=0 then rep(social_male,n_juv)
			#if p_polyandry==1 & p_sire==0 then sample(males,n_juv*length(breeding_females))

			sire=c(lapply(1:n_pair,function(i){
				## probability of any EPP
				polyandry <- rbinom(1,1,p_polyandry)
				if(polyandry){
					## if there is EPP, how much
					## this is calculated by sampling how many of the offspring the paired male sired, and then giving the same probability to subsequent males. This means that extra pair males will be few, and have several offspring if p_sire if high - thin this is more realistic
					if(p_sire==0) {
						sample(males,n_juv[i])
					}else{
						within_sires <- fill_sires(n_juv[i],p_sire)
						c(social_male[i], sample(males,max(within_sires)-1,replace=FALSE))[within_sires]
					}
					# n_sired <- rbinom(1,n_juv[i],p_sire)
					# c(rep(social_male[i], n_sired), sample(males,n_juv[i]-n_sired,replace=TRUE))
				}else{
					rep(social_male[i],n_juv[i])
				}
			}), recursive=TRUE),
			#equal sex ratio
			sex=sample(c("M","F"),sum(n_juv),replace=TRUE),
			cohort=year)

		# print(table(aggregate(animal~dam+cohort,pedigree, function(x)length(x))[,3]))
		# print(table(aggregate(animal~dam+cohort,pedigree, function(x)length(x))[,2]))

		## immigrants
		## need to sort out with constant pop
		## might be worth giving the immigrants age=afr, so cohort=year-afr

		n_imm <- if(constant_pop){
			c(immigration_f,immigration_m)*n_females
		}else{
			rbinom(2,n_females,c(immigration_f,immigration_m))
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




		## create individuals present in the next year


## or could take the new recruits from the whole pedigree, subset by cohort = year-afr - maybe this is better, because otherwise dat always has a load of pre afr individuals in it

		next_year_ind <- if(constant_pop){
			rbind(
				### need to ensure equal sex ratio of recruits, otherwise population size fluctuations
				ped[sample(which(ped[,"sex"]=="F"), juv_surv_f*fecundity*n_females/2, replace=FALSE),c(1,4)],
				ped[sample(which(ped[,"sex"]=="M"), juv_surv_m*fecundity*n_females/2, replace=FALSE),c(1,4)],
				if(adult_surv>0){	
					cbind(animal=sample(females, adult_surv_f*n_females, replace=FALSE), sex="F")},
				if(adult_surv>0){	
					cbind(animal=sample(males, adult_surv_m*n_females, replace=FALSE), sex="M")},
				immigrants[,c(1,4)]
			)
		}else{
			rbind(
				## need to make sex-specific
				ped[as.logical(rbinom(nrow(ped),1,juv_surv)),c(1,4)],
				cbind(animal=females[as.logical(rbinom(length(females),1,adult_surv_f))], sex="F"),
				cbind(animal=males[as.logical(rbinom(length(males),1,adult_surv_m))], sex="M"),
				immigrants[,c(1,4)]
			)
		}

		pedigree <- rbind(pedigree,ped,immigrants)

		## maybe we don't need this - redundant with cohort?
	  next_year_ind$age<-(year+1) - pedigree[match(next_year_ind$animal,pedigree$animal),"cohort"]
	  next_year_ind$year <- year+1
	  dat[[year+1]] <- next_year_ind
	}

	return(list(pedigree=pedigree,data_str=do.call(rbind,dat)))

}
