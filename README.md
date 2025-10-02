# pedAgree (version 0.0.1.0)
R package for demographically realistic pedigree generation using individual based simulations

## Installation

```
devtools::install_github("squidgroup/pedAgree")
library(pedAgree)
```

## Use

Just one function exported at the moment `simulate_pedigree()`. This takes a bunch of parameters (many can be sex specific), to generate stochastic pedigrees via individual based simulations.

### Pedigree Size
- **years**: Number of time steps/synchronous reproductive events
- **n_females**: Starting number of breeding females
-	**constant_pop**: Should there be stochastic variation in population size? (Logical)

### Demography
- **afr**: Age at first reproduction
- **p_breed**: Probability that a female breeds. Can be sex-specific if entered as a vector of length 2 (female,male) 
- **fecundity**: Mean number of juveniles a female produces each year.
- **fixed_fecundity**: Is fecundity fixed or drawn from a Poisson distribution? (Logical)
- **juv_surv**: Probability of juveniles surviving until local recruitment, where recruitment is defined as having genetic offspring. Can be sex-specific if entered as a vector of length 2 (female,male) 
- **adult_surv**: Probability of survival of adults across years. Can be sex-specific if entered as a vector of length 2 (female,male) 
- **immigration**:: Yearly immigration rate, as a proportion of starting number of females (n_females). Can be sex-specific if entered as a vector of length 2 (female,male) 

### Mating System
- **p_polyandry**: Probability that a female mates with multiple males.
- **p_sire**: Probability that 'social' male sires all offspring. Conditional on p_polyandry being >0
- **p_retain**: Probability that social partnership is retained across years. 

