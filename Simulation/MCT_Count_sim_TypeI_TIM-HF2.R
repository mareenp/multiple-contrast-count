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

###------------------------------ Main Function ------------------------------------###
mct_count_sim <-function(nsim,nboot,a,param,n,m,no,zeros){
  #--------------------------save results into-----------------------------------#
  set.seed(12345)
  
  #-----------------------------sample size--------------------------------------#
  N = sum(n)
  nmat = matrix(rep(n,nsim),ncol=nsim)
  #----------------------Contrast Matrices---------------------------------------#
  ## individual Contrast matrix 
  CC <- matrix(c(
    1,  0,  -1,  0,  # 1 vs 3
    0,  1,  0,  -1  # 2 vs 4
    
  ), nrow = 2, byrow = TRUE)
  Contrast = CC
  
  CC2 = CC^2
  nc = nrow(CC)
  
  result_sim=data.frame(Method=c(rep("boot",nc), rep("hom",nc), 
                                 rep("het",nc), rep("het_log",nc), 
                                 rep("het_sqrt",nc), rep("poi",nc), 
                                 rep("nb",nc), rep("q_poi",nc)),
                        p_value=rep(NA,nc*8*nsim),
                        i=rep(1:nsim, each = 8*nc),
                        H0=rep("TRUE",8*nsim*nc)) 
  
  #-------------------------Resampling matrix B----------------------------------#
  ## List of a matrices with dimension n*nboot
  B = lapply(1:a,function(arg){(matrix(sample(matrix(1:n[arg],n=n[arg],ncol=nboot), replace = T),n=n[arg],ncol=nboot))})
  
  #----------------------------Data Generation-----------------------------------#
  #x is list of a matrices
  is_zero = lapply(1:a,function(arg){matrix(rbinom(n=n[arg]*nsim,size=1, prob = zeros),ncol=nsim)})
  x_count=lapply(1:a,function(arg){matrix(rnbinom(n=n[arg]*nsim,mu=param[[1]][arg],size=param[[2]][arg]),ncol=nsim)})
  x=mapply(function(binary, count){ifelse(binary == 1, 0, count)}, is_zero, x_count, SIMPLIFY = FALSE)
  
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
      
      # emmeans uses list for Contrasts
      contr <-  list(
        "1 - 3" = CC[1,],
        "2 - 4" = CC[2,]
      )
      
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
                                                       percent_NA=mean(is.na(p_global)))
  results <- results %>% mutate(nsim=nsim,
                                nboot=nboot,
                                m=m,
                                no=no)
  return(results)
}


###------------------------- Simulation with 4 groups -------------------------------###
mct_count_4 <- function(nsim,nboot,a,lambda1,n1,n2,n3,n4,m,size1,size2,size3,size4,
                        no,zeros){
  n=c(n1,n2,n3,n4)
  param=list(c(lambda1,lambda1, lambda1, lambda1),c(size1,size2,size3,size4))
  
  ret <- mct_count_sim(nsim,nboot,a,param,n,m,no,zeros)
  return(cbind(ret,data.frame(lambda=lambda1,n1=n1,n2=n2,n3=n3,n4=n4,
                              size1=size1,size2=size2,size3=size3,size4=size4))
  )
}

###------------------------------ Settings ----------------------------------### 
nsim=1e+4 
nboot=5e+03
a=c(4) 
n1=c(seq(60,400,20))
m=c(seq(0,340,20))
lambda1=c(20)
zeros=c(0.65)

############## NBZ ##############
size1=c(rep(0.4,length(n1)))
size2=c(rep(1.35,length(n1)))
size3=c(rep(0.9,length(n1)))
size4=c(rep(5.2,length(n1)))


scenarios <- data.frame(nsim, nboot, a, lambda1, 
                   n1=n1, n2=n1, n3=n1, n4=n1, m, size1, size2, 
                   size3, size4, no=20, zeros)


###------------------------- Parallel on Cluster ----------------------------###
sjob <- slurm_apply(mct_count_4, scenarios, jobname = "Stat_Complete",
                    nodes = ceiling(nrow(scenarios)/4), cpus_per_node = 4, submit = TRUE,
                    slurm_options = list(time = "3-00", partition = "medium", "mem-per-cpu"="4G", 
                                         output="%x-%j.log"),
                    global_objects = c("mct_count_sim"))
results <- get_slurm_out(sjob, outtype = 'table', wait = TRUE)

write_xlsx(results, "MCT_Count_typeI_TIM-HF2.xlsx")

version
