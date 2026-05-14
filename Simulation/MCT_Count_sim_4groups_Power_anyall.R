library(devtools)
#devtools::install_github("cran/cotram")
library(writexl)
library(multcomp)
library(sandwich)
library(cotram)
library(tram)
library(nparcomp)
library(rslurm)
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
n1=c(12)
n2=c(16)
n3=c(22)
m=c(0)
lambda1=c(6) 
delta1=c(0,3)
delta2=c(1.5)
delta3=c(2.5)
Contrast=c("Dunnett","Tukey")

############## POI ##############
distribution=c("POI")
size1=c(NA) 

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
size1=c(0.2)
size2=c(0.5)
size3=c(2)

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
size1=c(0.75)
size2=c(2)
size3=c(3)
size4=c(5)

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


###------------------------- Parallel on Cluster ----------------------------###
sjob <- slurm_apply(mct_count_4, scenarios, jobname = "Stat_Complete",
                    nodes = ceiling(nrow(scenarios)/4), cpus_per_node = 4, submit = TRUE,
                    slurm_options = list(time = "3-00", partition = "medium", "mem-per-cpu"="4G", 
                                         output="%x-%j.log"),
                    global_objects = c("mct_count_sim"))
results <- get_slurm_out(sjob, outtype = 'table', wait = TRUE)

write_xlsx(results, "MCT_Count_power_anyall_4groups.xlsx")

version
