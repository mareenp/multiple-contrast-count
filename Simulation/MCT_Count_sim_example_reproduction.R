### examplary simulation of one setting to check reproducibility of type-I error
### and one setting to check reproducibility of power

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
library(readxl)

source("./MCT_Count_functions.R")

### ------------------------- Type I error --------------------------------- ###
## saved results from cluster
results_type1 <- read_xlsx("../2_results/Results_Simulation_Excel/MCT_Count_typeI_4groups.xlsx", guess_max = 2000)

## one setting for NB distributed data 
example_type1 <- mct_count_4(nsim=10000,
            nboot=5000,
            a=4,
            lambda1=1,
            n1=26, n2=30, n3=30, n4=36,
            m=20,
            size1=2, size2=3,
            size3=3, size4=5,
            no=14,
            delta1=0,
            delta2=0,
            delta3=0,
            distribution="NB",
            Contrast="Dunnett")

#### Test equality of results
results_type1[10537:10544,]
example_type1 


### ------------------------------ Power ----------------------------------- ###
## saved results from cluster
results_power <- read_xlsx("../2_results/Results_Simulation_Excel/MCT_Count_power_4groups.xlsx", guess_max = 2200)

## one setting for POI distributed data 
example_power <- mct_count_4(nsim=10000,
                       nboot=5000,
                       a=4,
                       lambda1=6,
                       n1=12, n2=12, n3=12, n4=12,
                       m=0,
                       size1=0, size2=0,
                       size3=0, size4=0,
                       no=14,
                       delta1=2.4,
                       delta2=0,
                       delta3=0,
                       distribution="POI",
                       Contrast="Dunnett")

#### Test equality of results
results_power[537:544,]
example_power 


