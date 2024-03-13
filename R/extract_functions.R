
	######
	## ----- Function to extract 'age at first reproduction' from a pedigree
	######

AFR <- function(ped,sex_specific=TRUE){
  dams <- as.character(na.omit(unique(ped[,"dam"])))
	sires <- as.character(na.omit(unique(ped[,"sire"])))

	## AFR for females
	  ## for each dam, work out her cohort and her first year of repro - difference between them in AFR
  # i=dams[900]
  AFR_d <- t(sapply(dams,function(i){
	  cohort <- ped[ped[,1]==i, "cohort"]
	  parent_cohort <- na.omit(subset(ped, dam==i)$cohort)
		if(length(parent_cohort)>0){
	  	YFR <- min(parent_cohort)
	  	YFR-cohort	
	  }else{
	  	NA
	  }	
  }))
	
	## AFR for males
	AFR_s <- sapply(sires,function(i){
	  cohort <- ped[ped[,1]==i, "cohort"]
	  parent_cohort <- na.omit(subset(ped, sire==i)$cohort)
	  if(length(parent_cohort)>0){
	  	YFR <- min(parent_cohort)
	  	YFR-cohort	
	  }else{
	  	NA
	  }	  
	})

	## return modal AFR
	if(sex_specific){
		counts_d <- table(AFR_d)
  	counts_s <- table(AFR_s)

  	c(
  		f=as.numeric(names(counts_d[which(counts_d==max(counts_d))])),
  		m=as.numeric(names(counts_s[which(counts_s==max(counts_s))]))
  	)
	}else{
		counts <- table(c(AFR_d,AFR_s))
  	as.numeric(names(counts[which(counts==max(counts))]))
	}
}



	######
	## ----- Function to extract 'juvenile survival' from a pedigree
	######

juvenile_survival <- function(ped,sex_specific=TRUE){

	## chicks to include in the calculation of J
	  ## exclude those without parents (immigrants and founder) - these will inflate J
	  ## miss off last cohort (at least), and these cant recruit, so will bias J downwards
		## in order not to downward bias survival too much, maybe first work out modal age at first repro, and then exclude that many years off the end? Or correct for it?  
	chicks_to_recruit <- subset(ped,!(is.na(dam)&is.na(sire)) & cohort!=max(cohort, na.rm=TRUE))[,1]

	## all individuals observed as parents
  dams <- na.omit(unique(ped[,"dam"]))
	sires <- na.omit(unique(ped[,"sire"]))

	if(sex_specific){
		c(
			f = 2* sum(chicks_to_recruit %in% dams) / length(chicks_to_recruit),
			m = 2* sum(chicks_to_recruit %in% sires) / length(chicks_to_recruit)
		)
	}else{
		sum(chicks_to_recruit %in% c(dams,sires)) / length(chicks_to_recruit)
	}
}







## essentially these give paramters of a zero inflated process 

p_breed <- function(ped){
	breed_attempts <- aggregate(cohort~dam,ped,function(x)length(unique(x)))$cohort	
	last_year <- aggregate(cohort~dam,ped,max)$cohort
	first_year <- aggregate(cohort~dam,ped,min)$cohort

	mean(breed_attempts/(last_year - first_year +1))
}


fecundity <- function(ped){
	
	female_rs <- aggregate(animal~dam+cohort,ped,length)
	mean(female_rs[,3])
	# var(female_rs[,3])

}



adult_survival <- function(ped, sex_specific=TRUE){

	cohorts<-sort(unique(ped$cohort))

	## overall adult survival rate
	as_year_d <- sapply(1:(length(cohorts)-1), function(i){
		ped_t0 <- subset(ped,cohort==cohorts[i])
		ped_t1 <- subset(ped,cohort %in% cohorts[cohorts>cohorts[i]])

		dams_t0 <- unique(ped_t0$dam)
		dams_t1 <- unique(ped_t1$dam)

 		c(survived=sum(dams_t0 %in% dams_t1) ,total=length(dams_t0))
	})

	as_year_s <- sapply(1:(length(cohorts)-1), function(i){
		ped_t0 <- subset(ped,cohort==cohorts[i])
		ped_t1 <- subset(ped,cohort %in% cohorts[cohorts>cohorts[i]])

		sires_t0 <- unique(ped_t0$sire)
		sires_t1 <- unique(ped_t1$sire)

 		c(survived=sum(sires_t0 %in% sires_t1) ,total=length(sires_t0))
	})


	as_total_d <- rowSums(as_year_d)
	as_total_s <- rowSums(as_year_s)
	as_total <- as_total_s + as_total_d
	# as_d <- as_total_d[1]/as_total_d[2]
	# as_s <- as_total_s[1]/as_total_s[2]

	if(sex_specific){

	  c(f=as_total_d[1]/as_total_d[2], m=as_total_s[1]/as_total_s[2])
	}else{
	  as_total[1]/as_total[2]
	}


}


######
#	immigration
######
immigration <- function(ped, sex_specific=TRUE){
	founders <- ped$animal[is.na(ped$dam) & is.na(ped$sire)]
	year_d <- aggregate(cohort~dam,subset(ped,dam%in%founders),min)$cohort
	year_s <- aggregate(cohort~sire,subset(ped,sire%in%founders),min)$cohort
	all_d <- table(year_d)
	all_s <- table(year_s)
	c(
		f=mean(all_d[-1]/all_d[1]),
		m=mean(all_s[-1]/all_s[1])
		)
}




#   last_cohort <- max((unique(ped$cohort)), na.rm=TRUE)

# 	last_year_f <- aggregate(cohort~dam,ped,max)$cohort
# 	first_year_f <- aggregate(cohort~dam,ped,min)$cohort

# 	lifespan_f <- last_year_f-first_year_f +1

# # mean((last_year_f-first_year_f)/(last_year_f-first_year_f+1))

# 1-1/mean(lifespan_f)



#### immigration 



