extract_pedigree <- function(ped){

	######
	## ----- Juvenile survival
	######

	## in order not to downward bias survival too much, maybe first work out modal age at first repro, and then exclude that many years off the end? Or correct of it?  

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

	f_table <- table(f_df$fecundity)
	f <- sum(f_df$fecundity)/sum(	f_table[names(f_table)!=0])

	p_breed <- 1-f_table[1]/sum(f_table)

	p_sire <- mean(f_df$prop_paternity, na.rm=TRUE)

	p_retain <- table(f_df$divorce)[1]/sum(table(f_df$divorce))


	# plot(table(f_df$n_males))
	# hist(f_df$prop_paternity, breaks=50)


	######
	## ----- Immigration
	######

	immigration <- 1 - ((js*f)/2 + as)

	out <- list(juv_surv=js,
		adult_surv=as,
		fecundity=f,
		p_breed = p_breed,
		p_sire = p_sire,
		p_retain=p_retain,
		immigration=immigration
		)
		lapply(out, function(x) {
			names(x)<-NULL
			x
		})

}