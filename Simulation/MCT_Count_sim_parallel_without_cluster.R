## Alternative code for parallelization without using a cluster computer
### using future_lapply instead of slurm_apply()
#### Example for simulation for Type I Error for 4 groups (function mct_count_4.

library(future.apply)
library(parallel)
library(devtools)
#devtools::install_github("cran/cotram")
library(writexl)
library(multcomp)
library(sandwich)
library(cotram)
library(tram)
library(nparcomp)
#devtools::install_github("thomas-fung/mpcmp")
library(mpcmp)
#devtools::install_github("cran/emmeans")
library(emmeans)
library(dplyr)

source("./MCT_Count_functions.R")

###------------------------------ Settings ----------------------------------### 
nsim=1e+4 
nboot=5e+03
a=c(4) 
n1=c(seq(6,30,2))
n2=c(seq(10,34,2))
n3=c(seq(16,40,2))
m=c(seq(0,24,2))
lambda1=c(1,6,10) 
delta1=delta2=delta3=c(0)
Contrast=c("Dunnett","Tukey","GrandMean")

############## POI ##############
distribution=c("POI")
size1=c(rep(NA,13)) 

#Balanced
POI1 <- data.frame(n1=n1, n2=n1, n3=n1, n4=n1, m, size1=size1, size2=size1, size3=size1, size4=size1, no=1)
#Unbalanced
POI1 <- rbind(POI1, data.frame(n1=n1, n2=n2, n3=n2, n4=n3, m, size1=size1, size2=size1, size3=size1, size4=size1, no=2))

expand.grid.df <- function(...) Reduce(function(...) merge(..., by=NULL), list(...))
scenarios_POI <- expand.grid.df(data.frame(nsim),data.frame(nboot),
                                data.frame(a),
                                data.frame(lambda1),
                                data.frame(POI1),
                                data.frame(delta1),
                                data.frame(delta2),
                                data.frame(delta3),
                                data.frame(distribution),
                                data.frame(Contrast))

############## CMP ###############
distribution=c("CMP") 
size1=c(rep(0.2,length(n1)))
size2=c(rep(0.5,length(n1)))
size3=c(rep(2,length(n1)))

#Balanced homoscedatic
CMP1 <- data.frame(n1=n1, n2=n1, n3=n1, n4=n1, m, size1=size2, size2=size2, size3=size2, size4=size2, no=3)
CMP1 <- rbind(CMP1, data.frame(n1=n1, n2=n1, n3=n1, n4=n1, m, size1=size3, size2=size3, size3=size3, size4=size3, no=4))
#Balanced heteroscedastic
CMP1 <- rbind(CMP1, data.frame(n1=n1, n2=n1, n3=n1, n4=n1, m, size1=size1, size2=size1, size3=size2, size4=size2, no=5))
CMP1 <- rbind(CMP1, data.frame(n1=n1, n2=n1, n3=n1, n4=n1, m, size1=size1, size2=size2, size3=size3, size4=size3, no=6))
#Unbalanced positive pairing
CMP1 <- rbind(CMP1, data.frame(n1=n1, n2=n2, n3=n2, n4=n3, m, size1=size2, size2=size2, size3=size1, size4=size1, no=7))
CMP1 <- rbind(CMP1, data.frame(n1=n1, n2=n2, n3=n2, n4=n3, m, size1=size3, size2=size3, size3=size2, size4=size1, no=8))
#Unbalanced negative pairing
CMP1 <- rbind(CMP1, data.frame(n1=n1, n2=n2, n3=n2, n4=n3, m, size1=size1, size2=size1, size3=size2, size4=size2, no=9))
CMP1 <- rbind(CMP1, data.frame(n1=n1, n2=n2, n3=n2, n4=n3, m, size1=size1, size2=size2, size3=size3, size4=size3, no=10))

scenarios_CMP <- expand.grid.df(data.frame(nsim),data.frame(nboot),
                                data.frame(a),
                                data.frame(lambda1),
                                data.frame(CMP1),
                                data.frame(delta1),
                                data.frame(delta2),
                                data.frame(delta3),
                                data.frame(distribution),
                                data.frame(Contrast))


############## NB ##############
distribution=c("NB") 
size1=c(rep(0.75,length(n1)))
size2=c(rep(2,length(n1)))
size3=c(rep(3,length(n1)))
size4=c(rep(5,length(n1)))


#Balanced homoscedastic
NB1 <- data.frame(n1=n1, n2=n1, n3=n1, n4=n1, m, size1=size3, size2=size3, size3=size3, size4=size3, no=11)
#Unbalanced positive pairing
NB1 <- rbind(NB1, data.frame(n1=n1, n2=n2, n3=n2, m, n4=n3, size1=size4, size2=size3, size3=size3, size4=size2, no=12))
#Unbalanced negative pairing
NB1 <- rbind(NB1, data.frame(n1=n1, n2=n2, n3=n2, n4=n3, m, size1=size1, size2=size2, size3=size2, size4=size3, no=13))
NB1 <- rbind(NB1, data.frame(n1=n1, n2=n2, n3=n2, n4=n3, m, size1=size2, size2=size3, size3=size3, size4=size4, no=14))

scenarios_NB <- expand.grid.df(data.frame(nsim),data.frame(nboot),
                               data.frame(a),
                               data.frame(lambda1),
                               data.frame(NB1),
                               data.frame(delta1),
                               data.frame(delta2),
                               data.frame(delta3),
                               data.frame(distribution),
                               data.frame(Contrast))

scenarios <- rbind(scenarios_POI, scenarios_CMP, scenarios_NB)

###############################################################################
### --------------- parallelization with future_apply --------------------- ###
# function to run one simulation scenario
run_scenario <- function(scn){
  res_scn <- mct_count_4(nsim=scn$nsim,
                         nboot=scn$nboot,
                         a=scn$a,
                         lambda1=scn$lambda1,
                         n1=scn$n1,n2=scn$n2,n3=scn$n3,n4=scn$n4,
                         m=scn$m,
                         size1=scn$size1,size2=scn$size2,
                         size3=scn$size3,size4=scn$size4,
                         no=scn$no,
                         delta1=scn$delta1,
                         delta2=scn$delta2,
                         delta3=scn$delta3,
                         distribution=scn$distribution,
                         Contrast=scn$Contrast)
  
  return(res_scn)
}

# parallelize scenarios
plan(multisession, workers=detectCores()*.75)
sim_out <- future_lapply(1:nrow(scenarios), function(i) run_scenario(scn=scenarios[i,]), future.seed = T)
plan(sequential)

# Merge results across all simulation scenarios
sim_results <- do.call(rbind, sim_out)
# save - insert the suitable Name for Excel file
write_xlsx(sim_results, "MCT_Count_typeI_4groups.xlsx")
