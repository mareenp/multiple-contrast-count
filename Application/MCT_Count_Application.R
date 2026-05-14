setwd("../Application")

pacman::p_load(nparcomp, dplyr, emmeans, COUNT, janitor, ggplot2, xtable,
               multcomp, sandwich, tram, cotram)

###------------------------ Application function ----------------------------###
mct_count <- function(dat, nboot, count, group, Contrast){
  #-----------------------------definitions------------------------------------#
  set.seed(839264)
  x_mat <- as_tibble(dat)
  x_mat <- x_mat %>% rename(count = !!count, group = !!group)
  levels_group <- levels(x_mat$group)
  a = length(c(levels(x_mat$group)))
  n = c(as.data.frame(table(x_mat$group))[,2])
  N = nrow(x_mat)
  nmat = matrix(rep(n,1),ncol=1)
  nc = nrow(contrMat(n=n,type=Contrast))
  erglog_boot = CIlow_boot = CIup_boot = c(rep(NA,nc))
  x = lapply(1:a,function(arg){as.matrix(unname(x_mat[x_mat$group==levels_group[arg],]["count"]))})
  result=data.frame(Method=c(rep("boot",nc), rep("hom",nc), rep("het",nc), 
                             rep("het_log",nc), rep("het_sqrt",nc), 
                             rep("poi",nc), rep("nb",nc), rep("q_poi",nc)),
                    Contrast=rep(row.names(contrMat(n=n,type=Contrast)),8),
                    Estimate=rep(NA,nc*8),
                    SE=rep(NA,nc*8),
                    Statistic=rep(NA,nc*8),
                    p_value=rep(NA,nc*8),
                    CI_95_low=rep(NA,nc*8),
                    CI_95_up=rep(NA,nc*8))
  
  #----------------------Contrast Matrices-------------------------------------#
  CC = contrMat(n=n,type=Contrast)
  CC2 = CC^2
  nc = nrow(CC)
  #-------------------------Resampling matrix B--------------------------------#
  ## List of a matrices with dimension n*nboot
  B = lapply(1:a,function(arg){(matrix(sample(matrix(1:n[arg],n=n[arg],ncol=nboot), replace = T),n=n[arg],ncol=nboot))})
  
  #--------------------------Estimate the rates--------------------------------#
  lambda = lapply(x,colMeans)
  lambdas = matrix(unlist(lambda),ncol=1,byrow=TRUE)
  loglambdas = log(lambdas)
  Clambda = CC%*%lambdas
  Cloglambdas = CC%*%loglambdas
  
  #---------------------------Estimate the Variances---------------------------#
  x2 = lapply(1:a,function(arg){x[[arg]]^2})
  sigma = lapply(1:a,function(arg){(colSums(x2[[arg]])-n[arg]*lambdas[arg,]^2)/(n[arg]-1)})
  sigman = lapply(1:a,function(arg){(colSums(x2[[arg]])-n[arg]*lambdas[arg,]^2)/(n[arg]*(n[arg]-1))}) #sigma/n
  sigmas = matrix(unlist(sigma),ncol=1,byrow=TRUE)
  sigmasn = matrix(unlist(sigman),ncol=1,byrow=TRUE)
  sigmaslambda2 = sigmas/lambdas^2
  sigmasnlambda2 = sigmasn/lambdas^2
  Csigmas = CC2%*%sigmas
  Csigmasn = CC2%*%sigmasn
  Csigmaslambda2 = CC2%*%sigmaslambda2
  Csigmasnlambda2 = CC2%*%sigmasnlambda2
  Cninv=CC%*%diag(1/n)%*%t(CC)
  vpool=colSums((nmat-1)*sigmas)/(N-a)
  vC=sapply(vpool,function(arg){(diag(Cninv*arg))})
  
  #-----------------------------Compute the Test Statistics--------------------#
  Tlog = Cloglambdas/sqrt(Csigmasnlambda2)
  
  i <- 1 #only 1 data set instead of simulation study 
  Tlog0 = max(abs(Tlog[,i]))
  
  if(!is.na(Tlog0) & Tlog0!=Inf & Tlog0!=-Inf){
    
    #---------------------------- Bootstrap -----------------------------------#
    x_i  = lapply(1:a,function(arg){matrix(x[[arg]][,i],ncol=1)}) #list of matrices for a groups
    #list of matrices of resampled X-vectors (X*), column = bootstrap sample
    x_star = lapply(1:a,function(arg){matrix(x_i[[arg]][B[[arg]]],ncol=nboot)}) 
    
    #-------------------- Estimate the rates - lambda_star --------------------#
    lambda_star = lapply(x_star,colMeans)
    lambdas_star = matrix(unlist(lambda_star),ncol=nboot,byrow=TRUE)
    loglambdas_star = log(lambdas_star)
    Clambda_star = CC%*%lambdas_star
    Cloglambdas_star = CC%*%loglambdas_star
    
    #------------------ Estimate the Variances - sigma_star -------------------#
    x2_star = lapply(1:a,function(arg){x_star[[arg]]^2})
    sigma_star = lapply(1:a,function(arg){(colSums(x2_star[[arg]])-n[arg]*lambdas_star[arg,]^2)/(n[arg]-1)})
    sigman_star = lapply(1:a,function(arg){(colSums(x2_star[[arg]])-n[arg]*lambdas_star[arg,]^2)/(n[arg]*(n[arg]-1))}) #sigma/n
    sigmas_star = matrix(unlist(sigma_star),ncol=nboot,byrow=TRUE)
    sigmasn_star = matrix(unlist(sigman_star),ncol=nboot,byrow=TRUE)
    sigmaslambda2_star = sigmas_star/lambdas_star^2
    sigmasnlambda2_star = sigmasn_star/lambdas_star^2
    Csigmas_star = CC2%*%sigmas_star
    Csigmasn_star = CC2%*%sigmasn_star
    Csigmaslambda2_star = CC2%*%sigmaslambda2_star
    Csigmasnlambda2_star = CC2%*%sigmasnlambda2_star
    
    #----------------------- Compute the Test Statistics ----------------------#
    Tlog_star = (Cloglambdas_star-Cloglambdas[,i])/sqrt(Csigmasnlambda2_star)
    Tlog0_star = lapply(1:nboot, function(arg){max(abs(Tlog_star[,arg]))})
    #-----------critical values of bootstrap distribution----------------------# 
    critlog_boot = quantile(unlist(Tlog0_star), probs = c(0.95), na.rm = T)
    
    erglog_boot <- sapply(1:nc, function(h) mean(unlist(lapply(1:nboot, function(arg){max(abs(na.omit(Tlog_star[,arg]))) > abs(Tlog[h])}))))
    CIlow_boot = Cloglambdas - critlog_boot * sqrt(Csigmasnlambda2)
    CIup_boot = Cloglambdas + critlog_boot * sqrt(Csigmasnlambda2)
    
  }
  
  #----------------------------- competitors LM -------------------------------# 
  ### homoscedastic variances
  mod1 <- lm(count~group, data=x_mat)
  hom <- glht(mod1, linfct = mcp(group = Contrast))
  ### heteroscedastic variances using sandwich variance estimator 
  het <- glht(mod1, linfct = mcp(group = Contrast),vcov=vcovHC)
  ###  Data transformation log(y+1), heteroscedastic variances 
  x_mat$lcount<-log(x_mat$count+1)
  mod2 <- lm(lcount~group, data=x_mat)
  het_log <- glht(mod2, linfct = mcp(group = Contrast), vcov=vcovHC)
  ### Data transformation sqrt(x+0.5), heteroscedastic variances 
  x_mat$scount<-sqrt(x_mat$count+0.5)
  mod3 <- lm(scount~group, data=x_mat)
  het_sqrt <- glht(mod3, linfct = mcp(group = Contrast), vcov=vcovHC)
  
  #----------------------- competitors GLM ------------------------------------#
  emm1 <- glm(count ~ group, data = x_mat, family = poisson)
  emm2 <- glm.nb(count ~ group, data = x_mat)
  emm3 <- glm(count ~ group, data = x_mat, family = quasipoisson(link="log"))
  
  ### emmeans uses different spelling for Contrasts
  contr <- if_else(Contrast=="Dunnett", "dunnett",
                   if_else(Contrast=="Tukey", "tukey", 
                           if_else(Contrast=="GrandMean", "eff", NA)))
  
  emm11 <- emmeans(emm1, specs="group", contr=contr, adjust="mvt")
  emm21 <- emmeans(emm2, specs="group", contr=contr, adjust="mvt")
  emm31 <- emmeans(emm3, specs="group", contr=contr, adjust="mvt")
  
  #####
  # Methods
  # 1 - boot; 2 - hom; 3 - het; 4 - het_log; 5 - het_sqrt
  # 6 - poi; 7 - nb; 8 - q_poi
  # Estimate
  result$Estimate[result$Method=="boot"] <- Cloglambdas
  result$Estimate[result$Method=="hom"] <- summary(hom)$test$coefficients
  result$Estimate[result$Method=="het"] <- summary(het)$test$coefficients
  result$Estimate[result$Method=="het_log"] <- summary(het_log)$test$coefficients
  result$Estimate[result$Method=="het_sqrt"] <- summary(het_sqrt)$test$coefficients
  result$Estimate[result$Method=="poi"] <- summary(emm11$contrasts)$estimate
  result$Estimate[result$Method=="nb"] <- summary(emm21$contrasts)$estimate
  result$Estimate[result$Method=="q_poi"] <- summary(emm31$contrasts)$estimate
  
  # SE
  result$SE[result$Method=="boot"] <-  sqrt(Csigmasnlambda2)
  result$SE[result$Method=="hom"] <- summary(hom)$test$sigma
  result$SE[result$Method=="het"] <- summary(het)$test$sigma
  result$SE[result$Method=="het_log"] <- summary(het_log)$test$sigma
  result$SE[result$Method=="het_sqrt"] <- summary(het_sqrt)$test$sigma
  result$SE[result$Method=="poi"] <- summary(emm11$contrasts)$SE
  result$SE[result$Method=="nb"] <- summary(emm21$contrasts)$SE
  result$SE[result$Method=="q_poi"] <- summary(emm31$contrasts)$SE
  
  # Statistic
  result$Statistic[result$Method=="boot"] <-  abs(Tlog)
  result$Statistic[result$Method=="hom"] <- abs(summary(hom)$test$tstat)
  result$Statistic[result$Method=="het"] <- abs(summary(het)$test$tstat)
  result$Statistic[result$Method=="het_log"] <- abs(summary(het_log)$test$tstat)
  result$Statistic[result$Method=="het_sqrt"] <- abs(summary(het_sqrt)$test$tstat)
  result$Statistic[result$Method=="poi"] <- abs(summary(emm11$contrasts)$z.ratio)
  result$Statistic[result$Method=="nb"] <- abs(summary(emm21$contrasts)$z.ratio)
  result$Statistic[result$Method=="q_poi"] <- abs(summary(emm31$contrasts)$z.ratio)
  
  # p-value
  result$p_value[result$Method=="boot"] <- erglog_boot
  result$p_value[result$Method=="hom"] <- summary(hom)$test$pvalues
  result$p_value[result$Method=="het"] <- summary(het)$test$pvalues
  result$p_value[result$Method=="het_log"] <- summary(het_log)$test$pvalues
  result$p_value[result$Method=="het_sqrt"] <- summary(het_sqrt)$test$pvalues
  result$p_value[result$Method=="poi"] <- summary(emm11$contrasts)$p.value
  result$p_value[result$Method=="nb"] <- summary(emm21$contrasts)$p.value
  result$p_value[result$Method=="q_poi"] <- summary(emm31$contrasts)$p.value
  
  # CI lower
  result$CI_95_low[result$Method=="boot"] <- CIlow_boot
  result$CI_95_low[result$Method=="hom"] <- confint(hom)$confint[,2]
  result$CI_95_low[result$Method=="het"] <- confint(het)$confint[,2]
  result$CI_95_low[result$Method=="het_log"] <- confint(het_log)$confint[,2]
  result$CI_95_low[result$Method=="het_sqrt"] <- confint(het_sqrt)$confint[,2]
  result$CI_95_low[result$Method=="poi"] <- confint(emm11$contrasts)$asymp.LCL
  result$CI_95_low[result$Method=="nb"] <- confint(emm21$contrasts)$asymp.LCL
  result$CI_95_low[result$Method=="q_poi"] <- confint(emm31$contrasts)$asymp.LCL
  
  # CI upper
  result$CI_95_up[result$Method=="boot"] <- CIup_boot
  result$CI_95_up[result$Method=="hom"] <- confint(hom)$confint[,3]
  result$CI_95_up[result$Method=="het"] <- confint(het)$confint[,3]
  result$CI_95_up[result$Method=="het_log"] <- confint(het_log)$confint[,3]
  result$CI_95_up[result$Method=="het_sqrt"] <- confint(het_sqrt)$confint[,3]
  result$CI_95_up[result$Method=="poi"] <- confint(emm11$contrasts)$asymp.UCL
  result$CI_95_up[result$Method=="nb"] <- confint(emm21$contrasts)$asymp.UCL
  result$CI_95_up[result$Method=="q_poi"] <- confint(emm31$contrasts)$asymp.UCL
  
  return(as.data.frame(c(result[1:2],round(result[,3:8],3))))
}


#####----------------------------- data example -------------------------------#####
data(badhealth)
#?badhealth
# group with "bad health" n=112
group_badh <- badhealth[badhealth$badh==1,]

# Agegroups
group_badh$age_cat <- cut(group_badh$age,
                    breaks = c(-Inf, 29, 39, 49, Inf),
                    labels = c("20-29", "30-39", "40-49", "50-60"))
tabyl(group_badh$age_cat)

### Table 3
tab3 <- group_badh %>% group_by(age_cat) %>% summarize(sum(numvisit>=0), mean(numvisit), var(numvisit), sum(numvisit==0)/sum(numvisit>=0) ,var(numvisit)/ mean(numvisit))
tab3
#print(xtable(tab3)) # Table as Latex Code
## export Table 3 as Excel
write_xlsx(tab3, "Table3.xlsx")

### use mct_count function - creates Table 4
out <- mct_count(dat=group_badh, nboot=5000, count="numvisit", group="age_cat", Contrast="Dunnett")
out 
#print(xtable(out)) # Table as Latex Code
## export Table 4 as Excel
write_xlsx(out, "Table4.xlsx")

#### Boxplot - creates Figure 4
p <- ggplot(group_badh, aes(y=numvisit, color=age_cat))+
  geom_boxplot()+
  theme_minimal(base_size = 12) +
  labs(y="Number of visits",
       color = "Age")+
  scale_x_continuous(breaks = c())+
  theme(legend.position = "bottom") +
  guides(color = guide_legend(nrow = 1))
p

png(filename = paste("Application_Boxplot.png"), 
    width = 2000, height = 2000, res = 400)
p
dev.off()
