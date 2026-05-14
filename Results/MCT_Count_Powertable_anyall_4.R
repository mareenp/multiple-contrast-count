setwd("../Results/Results_Simulation_Excel/")

pacman::p_load(ggplot2, writexl, readxl, dplyr, tidyr, purrr, janitor, 
               forcats, gridExtra, patchwork, cowplot, xtable)

Paa4 <- read_excel("MCT_Count_power_anyall_4groups.xlsx")

## change directory to the folder for tables
setwd("../../2_results/Results_Tables/")

############################# Power types ####################################
Paa4_tab <- Paa4 %>%  
  pivot_longer(cols = global:all,
               names_to = "Powertype",
               values_to = "Power") %>%
  pivot_wider(names_from = "Method",
              values_from = "Power") %>%
  subset(select=c(Powertype,boot,hom,het,nb,q_poi,lambda,delta1,no,Dist,Contrast))


##------------------------- Tables S1 - S4 --------------------------------##
##### Dunnett
### Table S1
S1 <- Paa4_tab %>% filter(lambda==6, delta1==0, Contrast=="Dunnett") %>% 
  subset(select=-c(Contrast, delta1)) %>% 
  rename(Setting=no) %>%
  mutate(across(2:6, ~ round(.x, 2)))
## export Table S1 as Excel
write_xlsx(S1, "TableS1.xlsx")
## Latex Code for Table S1
S1_latex <- xtable(S1,
                   caption = "Power-types for contrast Dunnett with $\\tilde{\\lambda}_1=\\tilde{\\lambda}_2=\\lambda, \\tilde{\\lambda}_3=\\lambda + 1.5, \\tilde{\\lambda}_4=\\lambda + 2.5$ and $\\lambda=6$, resulting in one true and two false null hypotheses.")

print(S1_latex, include.rownames = FALSE)

### Table S2
S2 <- Paa4_tab %>% filter(lambda==6, delta1==3, Contrast=="Dunnett") %>% 
  subset(select=-c(Contrast, delta1)) %>% 
  rename(Setting=no) %>%
  mutate(across(2:6, ~ round(.x, 2)))
## export Table S2 as Excel
write_xlsx(S2, "TableS2.xlsx")
## Latex Code for Table S2
S2_latex <- xtable(S2,
                   caption = "Power-types for contrast Dunnett, $\\tilde{\\lambda}_1=\\lambda, \\tilde{\\lambda}_2=\\lambda + 1.5, \\tilde{\\lambda}_3=\\lambda + 2.5, \\tilde{\\lambda}_4=\\lambda + 3$ and $\\lambda=6$, all null hypotheses are false.")

print(S2_latex, include.rownames = FALSE)

##### Tukey
### Table S3
S3 <- Paa4_tab %>% filter(lambda==6, delta1==0, Contrast=="Tukey") %>% 
  subset(select=-c(Contrast, delta1)) %>% 
  rename(Setting=no) %>%
  mutate(across(2:6, ~ round(.x, 2)))
## export Table S3 as Excel
write_xlsx(S3, "TableS3.xlsx")
## Latex Code for Table S3
S3_latex <- xtable(S3,
                   caption = "Power-types for contrast Tukey, $\\tilde{\\lambda}_1=\\tilde{\\lambda}_2=\\lambda, \\tilde{\\lambda}_3=\\lambda + 1.5, \\tilde{\\lambda}_4=\\lambda + 2.5$ and $\\lambda=6$, resulting in one true and five false null hypotheses.")

print(S3_latex, include.rownames = FALSE)

### Table S4
S4 <- Paa4_tab %>% filter(lambda==6, delta1==3, Contrast=="Tukey") %>% 
  subset(select=-c(Contrast, delta1)) %>% 
  rename(Setting=no) %>%
  mutate(across(2:6, ~ round(.x, 2)))
## export Table S4 as Excel
write_xlsx(S4, "TableS4.xlsx")
## Latex Code for Table S4
S4_latex <- xtable(S4,
                   caption = "Power-types for contrast Tukey, $\\tilde{\\lambda}_1=\\lambda, \\tilde{\\lambda}_2=\\lambda + 1.5, \\tilde{\\lambda}_3=\\lambda + 2.5, \\tilde{\\lambda}_4=\\lambda + 3$ and $\\lambda=6$, all null hypotheses are false.")

print(S4_latex, include.rownames = FALSE)
