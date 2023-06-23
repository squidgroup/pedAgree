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

source("/Users/joelpick/github/squidPed/R/functions.R")
source("/Users/joelpick/github/squidPed/R/simulate_pedigree_var.R")

pop<-simulate_pedigree(
	years = 10,
	n_females = 1000,
	fecundity = 4,
	fixed_fecundity=FALSE,
	juv_surv = 0.25,
	p_sire = 0,
	adult_surv = 0.5,
	det = list(rates=c("adult_surv"), G=0.1),
	constant_pop=FALSE,
	verbose=FALSE
	)

fec<-aggregate(animal~cohort,aggregate(animal~cohort+dam,pop$ped,length),mean)
pop_size <- aggregate(animal~cohort,pop$ped,length)

# plot(animal~cohort,fec, ylab="mean annual fecundity", xlab="year")
plot(animal~cohort,pop_size, ylab="population size", xlab="year")
diff(pop_size$animal) / pop_size$animal[1:9]

ped <- pop$ped

dams <- unique(ped$dam)
recruits <- c(dams,unique(ped$sire))
ped$recruited <- ped$animal%in%recruits
diff(table(ped$recruited,ped$cohort)[2,]/colSums(table(ped$recruited,ped$cohort)))

pnorm(0,qnorm(0.5,0,sqrt(1.05)),sqrt(1.05))^2*0.05


dat<-aggregate(animal~dam,ped,	function(x) sum(x%in%recruits))
names(dat) <- c("animal","recruits")
dat$cohort <- ped[match(dat$animal,ped$animal),"cohort"]
dat<-subset(dat,cohort!=9)

hist(dat$recruits)
mean(dat$recruits);var(dat$recruits)

ped2<-MCMCglmm::prunePed(ped[,1:3],keep=dat$animal,make.base=TRUE)
Ainv <- MCMCglmm::inverseA(ped2)$Ainv

prior <- list(
	R = list(V = 1, nu=0.0001),
	G=list(g1=list(V = 1, nu = 1, alpha.mu=0, alpha.V=100))
)

system.time(
	mod <- MCMCglmm::MCMCglmm(recruits~1,random=~animal,ginverse=list(animal=Ainv),family="poisson",data=dat, prior=prior)
)
summary(mod)
plot(mod)
 verbose = FALSE, nitt=50000, thin=40, burnin=10000

# sapply(dams,function(x){
# 	y=subset(ped,dam==x)
# 	c(zygotes=nrow(y),recruits=sum(y$animal%in%recruits))
# })



out2<-replicate(50,{
	ped<-simulate_pedigree(
		years = 10,
		n_females = 500,
		fecundity = 4,
		fixed_fecundity=TRUE,
		juv_surv = 0.5,
		adult_surv = 0,
		det = list(rates=c("juv_surv"), G=0.023),
		constant_pop=FALSE,
		verbose=FALSE
		)$ped

	# fec<-aggregate(animal~cohort,aggregate(animal~cohort+dam,ped,length),mean)
	pop_size <- aggregate(animal~cohort,ped,length)$animal


})

out2
plot(NA,xlim=c(1,10),ylim=range(out2))
apply(out2,2,lines)



G=0.02
MG=0.01

out<-do.call(rbind,parallel::mclapply(1:200,function(x){
	ped<-simulate_pedigree(
		years = 4,
		n_females = 10000,
		fecundity = 4,
		fixed_fecundity=FALSE,
		juv_surv = 0.5,
		adult_surv = 0,
		det = list(rates=c("juv_surv"), G=G,MG=MG),
		constant_pop=FALSE,
		verbose=FALSE
		)$ped

	# fec<-aggregate(animal~cohort,aggregate(animal~cohort+dam, ped,length),mean)
	# pop_size <- aggregate(animal~cohort, ped,length)$animal

	dams <- unique(ped$dam)
	recruits <- c(dams,unique(ped$sire))
	ped$recruited <- ped$animal%in%recruits

	diff(table(ped$recruited,ped$cohort)[2,]/colSums(table(ped$recruited,ped$cohort)))
},mc.cores=8))
colMeans(out)
apply(out,2,sd)/sqrt(1000)


trans_func <- pnorm(0,qnorm(0.5,0,sqrt(1+G+MG)),sqrt(1+G+MG))^2
trans_func*(G+0.5*MG)


pnorm(0,qnorm(0.5,0,sqrt(1.023)),sqrt(1.023))^2*0.023
(pnorm(0,qnorm(0.5,0,sqrt(1.023)),sqrt(1.023))^2*0.023)/0.25

