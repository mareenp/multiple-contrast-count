setwd("../Results/Results_Tables/")

devtools::install_github("thomas-fung/mpcmp")
pacman::p_load(dplyr, stringr, xtable, mpcmp)

###------------------------------- Table 1 ----------------------------------###
table1 <- read.delim(text="Setting & Distrib. & Sample sizes & Dispersion & Dispersion parameter & Interpretation\
1 & POI & n=n_1 + m & equidispersed & - & Balanced homosc. \
2 & POI & n=n_2 + m & equidispersed & - & Unbalanced homosc. \
3 & CMP & n=n_1 + m & overdispersed & nu=(0.5,0.5,0.5,0.5)' & Balanced homosc. \
4 & CMP & n=n_1 + m & underdispersed & nu=(2,2,2,2)' & Balanced homosc. \
5 & CMP & n=n_1 + m & overdispersed & nu=(0.2,0.2,0.5,0.5)' & Balanced heterosc. \
6 & CMP & n=n_1 + m &  over- and underd. & nu=(0.2,0.5,2,2)' & Balanced heterosc. \
7 & CMP & n=n_2 + m & overdispersed & nu=(0.5,0.5,0.2,0.2)' & Positive pairing \
8 & CMP & n=n_2 + m & over- and underd. & nu=(2,2,0.5,0.2)' & Positive pairing \
9 & CMP & n=n_2 + m & overdispersed & nu=(0.2,0.2,0.5,0.5)' & Negative pairing \
10 & CMP & n=n_2 + m & over- and underd. & nu=(0.2,0.5,2,2)' & Negative pairing \
11 & NB & n=n_1 + m & overdispersed & size=(3,3,3,3)' & Balanced homosc. \
12 & NB & n=n_2 + m & overdispersed & size=(5,3,3,2)' & Positive pairing \
13 & NB & n=n_2 + m & overdispersed & size=(0.75,2,2,3)' & Negative pairing \
14 & NB & n=n_2 + m & overdispersed & size=(2,3,3,5)' & Negative pairing\ "
,sep="&")

### save Table 1 as Excel
write_xlsx(table1, "Table1.xlsx")

###------------------------------- Table 2 ----------------------------------###
table2 <- read.delim(text="Characteristic_Sensitive_to & boot & hom & het & het_log & het_sqrt & poi & nb  & q_poi \
# groups                                       & No  & No & No & No & No & No & No & No \
Type of contrast                                  & No  & Yes & No  & Yes & Yes & Yes & Yes & Yes \
Data Distribution                                 & No  & Yes & No  & Yes & Yes & Yes & Yes & Yes \
Overdispersion                                    & No  & Yes & No  & Yes & Yes & Yes & No  & No \
Underdispersion                                   & No  & Yes & No  & Yes & Yes & Yes & Yes & No \ 
Small Sample Sizes combined with small lambda     & Yes & Yes & Yes & Yes & Yes & Yes & Yes & Yes \
Small Sample Sizes combined with higher lambda    & Yes  & Yes & Yes & Yes & Yes & Yes & Yes & Yes", 
sep="&")

### save Table 2 as Excel
write_xlsx(table2, "Table2.xlsx")

###------------------------------- Table A1 ---------------------------------###
tableA1 <- tibble(Distribution = c(paste("POI", 1:3), paste("CMP", 1:9),
                                   paste("NB", 1:12)),
                  Lambda = c(1,6,10,rep(1,3),rep(6,3),rep(10,3),
                             rep(1,4),rep(6,4),rep(10,4)),
                  Dispersion = c(NA,NA,NA,rep(c(0.2,0.5,2),3),rep(c(0.75,2,3,5),3)),
                  Variance = rep(NA, 24),
                  PropZero = rep(NA, 24),
                  VarByMean = rep(NA, 24))

tableA1 <- tableA1 %>% mutate(Variance=round(if_else(str_detect(Distribution, "POI"),Lambda,
                                               if_else(str_detect(Distribution, "NB"),Lambda+Lambda^2/Dispersion,NA)),1))

set.seed(123456)
x1=lapply(1:3,function(arg){matrix(rpois(n=100000,lambda=tableA1$Lambda[1:3][arg]),ncol=1)})
x2=lapply(1:9,function(arg){matrix(rcomp(n=100000,mu=tableA1$Lambda[4:12][arg],nu=tableA1$Dispersion[4:12][arg]),ncol=1)})
x3=lapply(1:12,function(arg){matrix(rnbinom(n=100000,mu=tableA1$Lambda[13:24][arg],size=tableA1$Dispersion[13:24][arg]),ncol=1)})

  
tableA1$Variance[4:12] <- round(unlist(lapply(x2, function(mat) apply(mat, 2, var))),1)
tableA1$PropZero[1:3] <- round(unlist(lapply(x1, function(arg) {sum(arg==0) / 100000})),2)
tableA1$PropZero[4:12] <- round(unlist(lapply(x2, function(arg) {sum(arg==0) / 100000})),2)
tableA1$PropZero[13:24] <- round(unlist(lapply(x3, function(arg) {sum(arg==0) / 100000})),2)

tableA1 <- tableA1 %>% mutate(VarByMean=round(Variance/Lambda,1)) 

### save Table A1 as Excel
write_xlsx(tableA1, "TableA1.xlsx")

### print for Latex
A1_latex <- xtable(tableA1)
print(A1_latex)

###------------------------------- Table S5 ----------------------------------###
tableS5 <- read.delim(text="Setting & Distrib. & Sample sizes & Dispersion & Dispersion parameter & Interpretation\
1 & POI & n=n_1 + m & equidispersed & - & Balanced homosc. \
2 & POI & n=n_2 + m & equidispersed & - & Unbalanced homosc. \
3 & CMP & n=n_1 + m & overdispersed & nu=(0.5,0.5,0.5)' & Balanced homosc. \
4 & CMP & n=n_1 + m & underdispersed & nu=(2,2,2)' & Balanced homosc. \
5 & CMP & n=n_1 + m & overdispersed & nu=(0.2,0.5,0.5)' & Balanced heterosc. \
6 & CMP & n=n_1 + m &  over- and underd. & nu=(0.2,0.5,2)' & Balanced heterosc. \
7 & CMP & n=n_2 + m & overdispersed & nu=(0.5,0.5,0.2)' & Positive pairing \
8 & CMP & n=n_2 + m & over- and underd. & nu=(2,0.5,0.2)' & Positive pairing \
9 & CMP & n=n_2 + m & overdispersed & nu=(0.2,0.5,0.5)' & Negative pairing \
10 & CMP & n=n_2 + m & over- and underd. & nu=(0.2,0.5,2)' & Negative pairing \
11 & NB & n=n_1 + m & overdispersed & size=(3,3,3)' & Balanced homosc. \
12 & NB & n=n_2 + m & overdispersed & size=(5,3,2)' & Positive pairing \
13 & NB & n=n_2 + m & overdispersed & size=(0.75,2,3)' & Negative pairing \
14 & NB & n=n_2 + m & overdispersed & size=(2,3,5)' & Negative pairing\ "
                     ,sep="&")

### save Table S5 as Excel
write_xlsx(tableS5, "TableS5.xlsx")
