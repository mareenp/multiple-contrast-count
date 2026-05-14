# dat - dataset
# nboot - number of bootstraps
# group - group-variable
# Contrast - type of contrast ("Dunnett" / "Tukey" / "GrandMean")
# offset - offset-variable, if used. Default is "none" indicating that none is used. Weighted version is used.

mct_count <- function(dat, nboot, count, group, Contrast, offset="none"){
  #-----------------------------definitions------------------------------------#
  set.seed(839264)
  x_mat <- as_tibble(dat)
  if(offset=="none") x_mat <- x_mat %>% rename(count = !!count, group = !!group)
  if(offset!="none") x_mat <- x_mat %>% rename(count = !!count, group = !!group, offset = !!offset)
  if(offset!="none") x_mat <- x_mat %>% mutate(x_t = count / offset)
  levels_group <- levels(x_mat$group)
  a = length(c(levels(x_mat$group)))
  n = as.vector(table(x_mat$group))
  N = nrow(x_mat)
  nmat = matrix(rep(n,1),ncol=1)
  nc = nrow(contrMat(n=n,type=Contrast))
  erglog_boot = CIlow_boot = CIup_boot = c(rep(NA,nc))
  x = lapply(1:a,function(arg){as.matrix(unname(x_mat[x_mat$group==levels_group[arg],]["count"]))})
  if(offset!="none"){
    x_t_off = lapply(1:a,function(arg){as.matrix(unname(x_mat[x_mat$group==levels_group[arg],]["x_t"]))})
    t_offset = lapply(1:a,function(arg){as.matrix(unname(x_mat[x_mat$group==levels_group[arg],]["offset"]))})  
  }
  
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
  
  ### without offset
  if(offset=="none"){
    #--------------------------Estimate the rates--------------------------------#
    lambda = lapply(x,colMeans)
    lambdas = matrix(unlist(lambda),ncol=1,byrow=TRUE)
    loglambdas = log(lambdas)
    Cloglambdas = CC%*%loglambdas
    
    #---------------------------Estimate the Variances---------------------------#
    x2 = lapply(1:a,function(arg){x[[arg]]^2})
    sigman = lapply(1:a,function(arg){(colSums(x2[[arg]])-n[arg]*lambdas[arg,]^2)/(n[arg]*(n[arg]-1))}) #sigma/n
    sigmasn = matrix(unlist(sigman),ncol=1,byrow=TRUE)
    sigmasnlambda2 = sigmasn/lambdas^2
    Csigmasnlambda2 = CC2%*%sigmasnlambda2
    
    #-----------------------------Compute the Test Statistics--------------------#
    Tlog = Cloglambdas/sqrt(Csigmasnlambda2)
    
    i <- 1 #only 1 data set instead of simulation study 
    Tlog0 = max(abs(Tlog[,i]))
    
    #---------------------------- Bootstrap -----------------------------------#
    x_i  = lapply(1:a,function(arg){matrix(x[[arg]][,i],ncol=1)}) #list of matrices for a groups
    #list of matrices of resampled X-vectors (X*), column = bootstrap sample
    x_star = lapply(1:a,function(arg){matrix(x_i[[arg]][B[[arg]]],ncol=nboot)}) 
    
    #-------------------- Estimate the rates - lambda_star --------------------#
    lambda_star = lapply(x_star,colMeans)
    lambdas_star = matrix(unlist(lambda_star),ncol=nboot,byrow=TRUE)
    loglambdas_star = log(lambdas_star)
    Cloglambdas_star = CC%*%loglambdas_star
    
    #------------------ Estimate the Variances - sigma_star -------------------#
    x2_star = lapply(1:a,function(arg){x_star[[arg]]^2})
    sigman_star = lapply(1:a,function(arg){(colSums(x2_star[[arg]])-n[arg]*lambdas_star[arg,]^2)/(n[arg]*(n[arg]-1))}) #sigma/n #####
    sigmasn_star = matrix(unlist(sigman_star),ncol=nboot,byrow=TRUE)
    sigmasnlambda2_star = sigmasn_star/lambdas_star^2
    Csigmasnlambda2_star = CC2%*%sigmasnlambda2_star
  }
  
  ### with offset - weighted
  if(offset!="none"){
    T_i = as.matrix(unname(tapply(x_mat$offset, x_mat$group, sum)), nsol=1)
    #--------------------------Estimate the rates--------------------------------#
    x_colsums = lapply(x,colSums)
    x_colsums = matrix(unlist(x_colsums),ncol=1,byrow=TRUE)
    lambdas = x_colsums / T_i
    loglambdas = log(lambdas)
    Cloglambdas = CC%*%loglambdas
    
    #---------------------------Estimate the Variances---------------------------#
    x_bar_w <- mapply(function(xi, ti) {sum(xi) / sum(ti)}, x, t_offset, SIMPLIFY = FALSE)
    Z_ik <- mapply(function(xi, ti, x_q) {xi - ti * x_q}, x, t_offset, x_bar_w, SIMPLIFY = FALSE)
    sigman_i <- mapply(function(ti, zi) {
      T_i <- sum(ti)
      K_i <- sum( (ti^2) / ((T_i - 2 * ti) * T_i) )
      sigma_sum <- sum( (T_i * zi^2) / (T_i - 2 * ti) )
      sigman_i <- (1 / ((1 + K_i) * T_i^2)) * sigma_sum
      return(c(sigman_i))
    }, t_offset, Z_ik)
    sigmasn = matrix(sigman_i,ncol=1,byrow=TRUE) 
    sigmasnlambda2 = sigmasn/lambdas^2
    Csigmasnlambda2 = CC2%*%sigmasnlambda2
    
    #-----------------------------Compute the Test Statistics--------------------#
    Tlog = Cloglambdas/sqrt(Csigmasnlambda2)
    
    i <- 1 #only 1 data set instead of simulation study 
    Tlog0 = max(abs(Tlog[,i]))
    
    #---------------------------- Bootstrap -----------------------------------#
    x_i  = lapply(1:a,function(arg){matrix(x[[arg]][,i],ncol=1)}) #list of matrices for a groups #######
    #list of matrices of resampled X-vectors (X*), column = bootstrap sample
    x_star = lapply(1:a,function(arg){matrix(x_i[[arg]][B[[arg]]],ncol=nboot)})
    T_offset_star = lapply(1:a,function(arg){matrix(t_offset[[arg]][B[[arg]]],ncol=nboot)}) 
    T_i_star = lapply(T_offset_star,colSums) 
    T_i_star = matrix(unlist(T_i_star),ncol=nboot,byrow=TRUE)
    
    #-------------------- Estimate the rates - lambda_star --------------------#
    x_colsums_star = lapply(x_star,colSums)
    x_colsums_star = matrix(unlist(x_colsums_star),ncol=nboot,byrow=TRUE)
    lambdas_star = x_colsums_star / T_i_star
    loglambdas_star = log(lambdas_star)
    Cloglambdas_star = CC%*%loglambdas_star
    
    #------------------ Estimate the Variances - sigma_star -------------------#
    x_bar_w_star = mapply(function(xi, ti) {colSums(xi) / colSums(ti)}, x_star, T_offset_star, SIMPLIFY = FALSE)
    Z_ik_star = mapply(function(xi, ti, x_q) {xi - sweep(ti, 2, x_q, "*")}, x_star, T_offset_star, x_bar_w_star, SIMPLIFY = FALSE)
    sigman_i_star = mapply(function(ti, zi) {
      T_i = colSums(ti)
      T_i_mat = matrix(rep(T_i, each = nrow(ti)), nrow = nrow(ti))
      denominator = T_i_mat - 2 * ti
      K_i <- colSums( (ti^2) / (denominator * T_i_mat) )
      sigma_sum <- colSums( (T_i_mat * zi^2) / denominator )
      sigman_i <- (1 / ((1 + K_i) * T_i^2)) * sigma_sum
      return(c(sigman_i))
    }, T_offset_star, Z_ik_star)
    sigmasn_star = matrix(sigman_i_star,ncol=nboot,byrow=TRUE) 
    sigmasnlambda2_star = sigmasn_star/lambdas_star^2
    Csigmasnlambda2_star = CC2%*%sigmasnlambda2_star
  }
  
  
  #----------------------- Compute the Test Statistics ----------------------#
  Tlog_star = (Cloglambdas_star-Cloglambdas[,i])/sqrt(Csigmasnlambda2_star)
  
    #----------------------- Compute the Test Statistics ----------------------#
    Tlog0_star = lapply(1:nboot, function(arg){max(abs(Tlog_star[,arg]))})
    #----------------------- Compute p-values ----------------------#
    erglog_boot <- sapply(1:nc, function(h) mean(unlist(lapply(1:nboot, function(arg){max(abs(na.omit(Tlog_star[,arg]))) > abs(Tlog[h])}))))
    #-----------critical values of bootstrap distribution & confidence intervals ----------------------# 
    critlog_boot = quantile(unlist(Tlog0_star), probs = c(0.95), na.rm = T)
    CIlow_boot = Cloglambdas - critlog_boot * sqrt(Csigmasnlambda2)
    CIup_boot = Cloglambdas + critlog_boot * sqrt(Csigmasnlambda2)

  if(offset!="none"){
    #----------------------------- competitors LM -------------------------------# 
    ### homoscedastic variances
    mod1 <- lm(count~group, data=x_mat, offset = offset)
    hom <- glht(mod1, linfct = mcp(group = Contrast))
    ### heteroscedastic variances using sandwich variance estimator 
    het <- glht(mod1, linfct = mcp(group = Contrast),vcov=vcovHC)
    ###  Data transformation log(y+1), heteroscedastic variances 
    x_mat$lcount<-log(x_mat$count+1)
    mod2 <- lm(lcount~group, data=x_mat, offset = offset)
    het_log <- glht(mod2, linfct = mcp(group = Contrast), vcov=vcovHC)
    ### Data transformation sqrt(x+0.5), heteroscedastic variances 
    x_mat$scount<-sqrt(x_mat$count+0.5)
    mod3 <- lm(scount~group, data=x_mat, offset = offset)
    het_sqrt <- glht(mod3, linfct = mcp(group = Contrast), vcov=vcovHC)
    
    #----------------------- competitors GLM ------------------------------------#
    emm1 <- glm(count ~ group, data = x_mat, family = poisson, offset = offset)
    emm2 <- glm.nb(count ~ group, data = x_mat, offset(offset))
    emm3 <- glm(count ~ group, data = x_mat, family = quasipoisson(link="log"), offset = offset)
    
    ### emmeans uses different spelling for Contrasts
    contr <- if_else(Contrast=="Dunnett", "dunnett",
                     if_else(Contrast=="Tukey", "tukey", 
                             if_else(Contrast=="GrandMean", "eff", NA)))
    
    emm11 <- emmeans(emm1, specs="group", contr=contr, adjust="mvt")
    emm21 <- emmeans(emm2, specs="group", contr=contr, adjust="mvt")
    emm31 <- emmeans(emm3, specs="group", contr=contr, adjust="mvt")
  }
    
  if(offset=="none"){
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
  }
  
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



