
###------------------------------ Main Function ------------------------------------###
mct_count_sim <-function(nsim,nboot,a,param,n,m,no,delta1,delta2,delta3,distribution,Contrast){
  #--------------------------save results into-----------------------------------#
  set.seed(12345)
  
  #-----------------------------sample size--------------------------------------#
  N = sum(n)
  nmat = matrix(rep(n,nsim),ncol=nsim)
  #----------------------Contrast Matrices---------------------------------------#
  CC = contrMat(n=n,type=Contrast)
  CC2 = CC^2
  nc = nrow(CC)
  ## for computation of any-pairs and all-pairs power check which H0 are true
  if(Contrast=="Dunnett" & is.na(param[[1]][4])==F){
    H0=rep(c((param[[1]][2]==param[[1]][1]), (param[[1]][3]==param[[1]][1]), (param[[1]][4]==param[[1]][1])),8*nsim)}
  if(Contrast=="Dunnett" & is.na(param[[1]][4])==T){
    H0=rep(c((param[[1]][2]==param[[1]][1]), (param[[1]][3]==param[[1]][1])),8*nsim)}
  if(Contrast=="Tukey" & is.na(param[[1]][4])==F){
    H0=rep(c((param[[1]][2]==param[[1]][1]), (param[[1]][3]==param[[1]][1]), (param[[1]][4]==param[[1]][1]), 
             (param[[1]][3]==param[[1]][2]), (param[[1]][4]==param[[1]][2]), (param[[1]][4]==param[[1]][3])),8*nsim)}
  if(Contrast=="Tukey" & is.na(param[[1]][4])==T){
    H0=rep(c((param[[1]][2]==param[[1]][1]), (param[[1]][3]==param[[1]][1]), 
             (param[[1]][3]==param[[1]][2])),8*nsim)}
  if(Contrast=="GrandMean"){ #any-pairs and all-pairs power is not computable for GrandMean
    # to only calculate type-I or global power, all H0 are set to "TRUE" (as construct to not calculate any- and all-pairs power)
    H0=rep("TRUE",8*nsim*nc)
  }
  result_sim=data.frame(Method=c(rep("boot",nc), rep("hom",nc), 
                                 rep("het",nc), rep("het_log",nc), 
                                 rep("het_sqrt",nc), rep("poi",nc), 
                                 rep("nb",nc), rep("q_poi",nc)),
                        Contrast=rep(row.names(contrMat(n=n,type=Contrast)),8*nsim),
                        p_value=rep(NA,nc*8*nsim),
                        i=rep(1:nsim, each = 8*nc),
                        H0=H0) 
  
  #-------------------------Resampling matrix B----------------------------------#
  ## List of a matrices with dimension n*nboot
  B = lapply(1:a,function(arg){(matrix(sample(matrix(1:n[arg],n=n[arg],ncol=nboot), replace = T),n=n[arg],ncol=nboot))})
  
  #----------------------------Data Generation-----------------------------------#
  if(distribution=="POI"){
    x=lapply(1:a,function(arg){matrix(rpois(n=n[arg]*nsim,lambda=param[[1]][arg]),ncol=nsim)})}
  if(distribution=="NB"){  #x is list of a matrices
    x=lapply(1:a,function(arg){matrix(rnbinom(n=n[arg]*nsim,mu=param[[1]][arg],size=param[[2]][arg]),ncol=nsim)})}
  if(distribution=="CMP"){  #x is list of a matrices
    x=lapply(1:a,function(arg){matrix(rcomp(n=n[arg]*nsim,mu=param[[1]][arg],nu=param[[2]][arg]),ncol=nsim)})}
  
  #----------------------------Estimate the rates--------------------------------#
  lambda = lapply(x,colMeans)
  lambdas = matrix(unlist(lambda),ncol=nsim,byrow=TRUE)
  loglambdas = log(lambdas)
  Clambda = CC%*%lambdas
  Cloglambdas = CC%*%loglambdas
  
  #-----------------------------Estimate the Variances---------------------------#
  x2 = lapply(1:a,function(arg){x[[arg]]^2})
  sigma = lapply(1:a,function(arg){(colSums(x2[[arg]])-n[arg]*lambdas[arg,]^2)/(n[arg]-1)})
  sigman = lapply(1:a,function(arg){(colSums(x2[[arg]])-n[arg]*lambdas[arg,]^2)/(n[arg]*(n[arg]-1))}) #sigma/n
  sigmas = matrix(unlist(sigma),ncol=nsim,byrow=TRUE)
  sigmasn = matrix(unlist(sigman),ncol=nsim,byrow=TRUE)
  sigmaslambda2 = sigmas/lambdas^2
  sigmasnlambda2 = sigmasn/lambdas^2
  Csigmas = CC2%*%sigmas
  Csigmasn = CC2%*%sigmasn
  Csigmaslambda2 = CC2%*%sigmaslambda2
  Csigmasnlambda2 = CC2%*%sigmasnlambda2
  Cninv=CC%*%diag(1/n)%*%t(CC)
  vpool=colSums((nmat-1)*sigmas)/(N-a)
  vC=sapply(vpool,function(arg){(diag(Cninv*arg))})
  
  #------------------------------- Compute the Test Statistic --------------------#
  Tlog = Cloglambdas/sqrt(Csigmasnlambda2)
  
  #------------------------------ Start of Simulation -----------------------------#
  for(i in 1:nsim){
    Tlog0 = max(abs(Tlog[,i]))
    
    if(!is.na(Tlog0) & Tlog0!=Inf & Tlog0!=-Inf){
      #---------------------- competitors using LM (multcomp) --------------------# 
      ### dataframe of i. column of x & column for group 
      x_mat <- as.data.frame(matrix(c(unlist(lapply(1:a,function(arg){matrix(x[[arg]][,i],ncol=1)})), c(rep(c(seq(1,a,1)),n))), ncol = 2))
      colnames(x_mat) <- c("count", "group")
      x_mat$group <- as.factor(x_mat$group)
      ### homoscedastic variances
      mod1 <- lm(count~group, data=x_mat)
      hom <- summary(glht(mod1, linfct = mcp(group = Contrast)))
      ### Approx. normal model using sandwich variance estimator 
      het <- summary(glht(mod1, linfct = mcp(group = Contrast),vcov=vcovHC))
      ### Data transformation log(y+1) 
      x_mat$lcount<-log(x_mat$count+1)
      mod2 <- lm(lcount~group, data=x_mat)
      het_log <- summary(glht(mod2, linfct = mcp(group = Contrast), vcov=vcovHC))
      ### Data transformation sqrt(x+0.5)
      x_mat$scount<-sqrt(x_mat$count+0.5)
      mod3 <- lm(scount~group, data=x_mat)
      het_sqrt <- summary(glht(mod3, linfct = mcp(group = Contrast), vcov=vcovHC))
      
      #------------------------ p-values LM ---------------------------------#
      result_sim$p_value[result_sim$Method=="hom" & result_sim$i==i] <- summary(hom)$test$pvalues
      result_sim$p_value[result_sim$Method=="het" & result_sim$i==i] <- summary(het)$test$pvalues
      result_sim$p_value[result_sim$Method=="het_log" & result_sim$i==i] <- summary(het_log)$test$pvalues
      result_sim$p_value[result_sim$Method=="het_sqrt" & result_sim$i==i] <- summary(het_sqrt)$test$pvalues
      
      #----------------- competitors using GLM (emmeans) ----------------------#
      poi <- glm(count ~ group, data = x_mat, family = poisson)
      nb <- glm.nb(count ~ group, data = x_mat)
      q_poi <- glm(count ~ group, data = x_mat, family = quasipoisson(link="log"))
      
      # emmeans uses different spelling for Contrasts
      contr <- dplyr::if_else(Contrast=="Dunnett", "dunnett",
                              dplyr::if_else(Contrast=="Tukey", "tukey", 
                                             dplyr::if_else(Contrast=="GrandMean", "eff", NA)))
      
      poi1 <- emmeans(poi, specs="group", contr=contr, adjust="mvt")
      nb1 <- emmeans(nb, specs="group", contr=contr, adjust="mvt")
      q_poi1 <- emmeans(q_poi, specs="group", contr=contr, adjust="mvt")
      
      #------------------------ p-values GLM ---------------------------------#
      result_sim$p_value[result_sim$Method=="poi" & result_sim$i==i] <- summary(poi1$contrasts)$p.value
      result_sim$p_value[result_sim$Method=="nb" & result_sim$i==i] <- summary(nb1$contrasts)$p.value
      result_sim$p_value[result_sim$Method=="q_poi" & result_sim$i==i] <- summary(q_poi1$contrasts)$p.value
      
      #------------------------------ Bootstrap -------------------------------------#
      x_i  = lapply(1:a,function(arg){matrix(x[[arg]][,i],ncol=1)}) 
      #list of matrices of resampled X-vectors (X*), column = bootstrap sample
      x_star = lapply(1:a,function(arg){matrix(x_i[[arg]][B[[arg]]],ncol=nboot)}) 
      
      #------------------------- Estimate the rates -----------------------------#
      lambda_star = lapply(x_star,colMeans)
      lambdas_star = matrix(unlist(lambda_star),ncol=nboot,byrow=TRUE)
      loglambdas_star = log(lambdas_star)
      Clambda_star = CC%*%lambdas_star
      Cloglambdas_star = CC%*%loglambdas_star
      
      #------------------------ Estimate the Variances ---------------------------#
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
      
      #--------------------------- Compute the Test Statistics --------------------#
      Tlog_star = (Cloglambdas_star-Cloglambdas[,i])/sqrt(Csigmasnlambda2_star)
      Tlog0_star = lapply(1:nboot, function(arg){max(abs(Tlog_star[,arg]))})
      
      #------------------------ p-values boot ---------------------------------#
      result_sim$p_value[result_sim$Method=="boot" & result_sim$i==i] <- sapply(1:nc, function(h) mean(unlist(lapply(1:nboot, function(arg){max(abs(na.omit(Tlog_star[,arg]))) > abs(Tlog[h,i])}))))
    }
    
    else{
      result_sim$p_value[result_sim$i==i] = NA 
    }
    
    #----------------------------End of Simulation-------------------------------#
  }
  
  #------------------------- results (type-I / power) ---------------------------#
  result_global <- result_sim %>% group_by(Method,i) %>% summarise(p_global=(min(p_value)<0.05))
  
  if(nrow(result_sim %>% filter(H0=="FALSE"))>0){
    result_aa <- result_sim %>% filter(H0=="FALSE") %>% group_by(Method,i) %>% 
      summarize(p_any=if_else(any((p_value<0.05)),TRUE,FALSE),
                p_all=if_else(all((p_value<0.05)),TRUE,FALSE))
    result <- merge(result_global, result_aa)}
  if(nrow(result_sim %>% filter(H0=="FALSE"))==0){
    result <- result_global %>% mutate(p_any=NA,
                                       p_all=NA)
  }
  
  results <- result %>% group_by(Method) %>% summarise(global=mean(na.omit(p_global)),
                                                       any=mean(na.omit(p_any)),
                                                       all=mean(na.omit(p_all)),
                                                       percent_NA=mean(is.na(p_global)))
  results <- results %>% mutate(nsim=nsim,
                                nboot=nboot,
                                Dist=distribution,
                                Contrast=Contrast,
                                delta1=delta1,
                                delta2=delta2,
                                delta3=delta3,
                                m=m,
                                no=no)
  return(results)
}


###------------------------- Simulation with 4 groups -------------------------------###
mct_count_4 <- function(nsim,nboot,a,lambda1,n1,n2,n3,n4,m,size1,size2,size3,size4,
                        no,delta1,delta2,delta3,distribution,Contrast){
  n=c(n1,n2,n3,n4)
  param=list(c(lambda1,lambda1 + delta1, lambda1 + delta2, lambda1 + delta3),c(size1,size2,size3,size4))
  
  ret <- mct_count_sim(nsim,nboot,a,param,n,m,no,delta1,delta2,delta3,distribution,Contrast)
  return(cbind(ret,data.frame(lambda=lambda1,n1=n1,n2=n2,n3=n3,n4=n4,
                              size1=size1,size2=size2,size3=size3,size4=size4))
  )
}


###-------------------------- Simulation with 3 groups ------------------------------###
mct_count_3 <- function(nsim,nboot,a,lambda1,n1,n2,n3,m,size1,size2,size3,
                        no,delta1,delta2,delta3,distribution,Contrast){
  n=c(n1,n2,n3)
  param=list(c(lambda1,lambda1 + delta1, lambda1 + delta2),c(size1,size2,size3))
  
  ret <- mct_count_sim(nsim,nboot,a,param,n,m,no,delta1,delta2,delta3,distribution,Contrast)
  return(cbind(ret,data.frame(lambda=lambda1,n1=n1,n2=n2,n3=n3,
                              size1=size1,size2=size2,size3=size3))
  )
}

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
