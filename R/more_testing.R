x <- matrix(0,ncol=6,nrow=6)

rownames(x) <- colnames(x) <- c("p_breed","fecundity","p_sire","p_retain","juv_surv","adult_surv")


rm(list=ls())

det<-NULL

G1 <- diag(c(3,4))
rownames(G1) <- colnames(G1) <- c("fecundity","juv_surv")
det<-list(G=G1)


transform_det <- function(det, env){
	# print(environment())
	# variance components
	components <- c("G", "PE", "MG", "ME", "C", "Y", "E")
	# vital rates that can vary 
	rates <- c("p_breed","fecundity","p_sire","p_retain","juv_surv","adult_surv")
	empty_mat <- matrix(0,ncol=length(rates),nrow=length(rates),dimnames=list(rates,rates))
# G <- PE <- MG <- ME <- C <- Y <- E <- empty_mat


	for(component in components){
		x <- empty_mat
		x[rownames(det[[component]]),colnames(det[[component]])] <- det[[component]]
		assign(component,x,pos=env)
	}
}
transform_det(det)

G
PE

test <- function(det){
	Renv <- environment()
	transform_det(det, env=Renv)
	G
}

test(det)


G <- diag(c(3,4))
rownames(G) <- colnames(G) <- c("fecundity","juv_surv")

x[rownames(G),colnames(G)] <- G

x

bv <- MASS::mvrnorm(n=10, mu=rep(0,ncol(x)), x)

bv <- MASS::mvrnorm(n=10, mu=matrix(rnorm(60),10,6), x)

a=1:10
b=1:10
c=1:10
d=1:10
e= FALSE
cbind(a,b,if(e){c},d)
cbind(
	a,
	b,
	if(e)c,
	if(e)d
	)

det <-  NULL
det[["G"]]
fecundity = 4
fecundity_l <- log(fecundity) - (0.2)/2
hist(exp(rnorm(100,fecundity_l,0.2 )))

hist()
rm(list=ls())

source("/Users/joelpick/github/pedigree_simulations/R/functions.R")
source("/Users/joelpick/github/pedigree_simulations/R/simulate_pedigree_var.R")

pop<-simulate_pedigree(
	years = 50,
	n_females = 100,
	fecundity = 4,
	juv_surv = 0.5,
	adult_surv = 0,
	det = list(rates=c("fecundity"), G=0.01),
	constant_pop=FALSE,
	verbose=FALSE
	)

fec<-aggregate(animal~cohort,aggregate(animal~cohort+dam,pop$ped,length),mean)
pop_size <- aggregate(animal~cohort,pop$ped,length)

plot(animal~cohort,fec, ylab="mean annual fecundity", xlab="year")
plot(animal~cohort,pop_size, ylab="population size", xlab="year")
