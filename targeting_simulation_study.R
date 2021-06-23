rm(list = ls())
library(devtools)
install_github("LendieFollett/CBART")
library(CBART)
library(ggplot2)
library(dplyr)
library(coda)
library(xtable)
library(reshape2)
library(randomForest)
library(gridExtra)
library(caret)


nd=5000        ## number of draws to keep
burn=5000       ## number of draws to discard


#parameters varying
# baseline parameterization: P=10, sigma_e=1, lambda=2, n=2000, rho=0.7 u_dist = "Exponential".
# Experiment 1 (beta0): Use beta0 values -10,0,1,2 
# Experiment 2 (beta1): Use beta1 values  -2,0,0.5,1
# Experiment 3 (rho): Use rho values 0.1, 0.3, 0.5, 0.7, 0.9
# Experiment 4 (sigma_e): Use sigma_e values 1, 2, 3, 4, 5
# Experiment 5 (P): Use P values 5, 10, 20, 50, 100
# Experiment 6 (n): Use n values 500, 1000, 2000, 5000, 10000


#SET DISTRIBUTION OF ASYMMETRIC TERM
#"Normal" or "Exponential"
true_udist <- "Normal"
nrep = 10
base <- data.frame(true_P = 10,
                n = 2000,
                true_rho = 0.7,
                true_beta0 = 1,
                true_beta1 = 0,
                true_sigma_e = 1,
                true_udist = true_udist)
ex1 <- data.frame(rbind(base, base, base, base))%>%mutate(
                  true_beta0 = c(-10,0,1,2 ),
                  experiment = 1)
ex2 <- data.frame(rbind(base, base, base, base))%>%mutate(
                  true_beta1 = c(-2,0,0.5,1),
                  experiment = 2)
ex3 <- data.frame(rbind(base, base, base, base, base))%>%mutate(
                  true_rho = c(0.1, 0.3, 0.5, 0.7, 0.9),
                  experiment = 3)
ex4 <- data.frame(rbind(base, base, base, base, base))%>%mutate(
                  true_sigma_e = c(1, 2, 3, 4, 5),
                  experiment = 4)
ex5 <- data.frame(rbind(base, base, base, base, base))%>%mutate(
                  true_P = c(5, 10, 20, 50, 100),
                  experiment = 5)
ex6 <- data.frame(rbind(base, base, base, base, base))%>%mutate(
                  n = c(500, 1000, 2000, 5000, 10000),
                  experiment = 6)

all <- rbind(ex1, ex2, ex3, ex4, ex5, ex6)
G <- all[rep(1:nrow(all), each = nrep),]
rownames(G) <- 1:nrow(G)
G$idx <- 1:nrow(G)


all_preds <- rep(list(NA),nrow(G))

hold <- data.frame(OLR_MSE = rep(NA, nrow(G)),
                   RF_MSE = rep(NA, nrow(G)),
                   BART_MSE = rep(NA, nrow(G)),
                   FMT_MSE = rep(NA, nrow(G)),
                   
                   OLR_fcor = rep(NA, nrow(G)),
                   RF_fcor = rep(NA, nrow(G)),
                   BART_fcor = rep(NA, nrow(G)),
                   
                   BART_time = rep(NA, nrow(G)),
                   experiment = G$experiment,
                   idx = G$idx)

#START SIMULATIONS
for (idx in 1:nrow(G)){
  #set varying levels
  beta0 <- G$true_beta0[idx] #1
  beta1 <- G$true_beta1[idx] #2
  udist <- G$true_udist[idx] #3
  rho <- G$true_rho[idx] #4
  sig <- G$true_sigma_e[idx] #5
  p <- G$true_P[idx] #6
  n <- G$n[idx] #7
  rep <- G$rep[idx]
  
  train_idx <- 1:(n)
  test_idx <- (n+1):(n + 500)
  
  set.seed(352833 + idx)
  x <- matrix(rbeta((n + 500)*p, 1, 1), ncol = p, byrow = FALSE)
  #true_lambda <- exp(true_beta0)
  true_f <- rep(NA, n)
  true_f2 <- rep(NA, n)
  true_u <- rep(NA, n)
  if (udist == "Exponential"){
    for (i in 1:(n + 500)){
      lam <- exp(beta0 + beta1*x[i,1])
      set.seed(238520 + i + idx*500)
      true_u[i] <- rexp(1, rate = 1/lam)
      set.seed(2352 + i + idx*500)
      epsilon <- rnorm(1, 0, sig)
      true_f[i] <- 10*sin(pi*x[i,1]*x[i,2]) + 20*(x[i,3] - 0.5)^2 + 10*x[i,4] + 5*x[i,5]  
      true_f2[i] <- true_f[i] + epsilon
    }
  }else{ #i.e., if it's "Normal"
    for (i in 1:(n + 500)){
      lam <- exp(beta0 + beta1*x[i,1])
      set.seed(238520 + i + idx*500)
      true_u[i] <- abs(rnorm(1,mean =0, sd = sqrt(pi/2)*lam))
      set.seed(2352 + i + idx*500)
      epsilon <- rnorm(1, 0, sig)
      true_f[i] <- 10*sin(pi*x[i,1]*x[i,2]) + 20*(x[i,3] - 0.5)^2 + 10*x[i,4] + 5*x[i,5]  
      true_f2[i] <- true_f[i] + epsilon
    } 
  }
  set.seed(23589  + idx*500)
  y <- true_f2  - true_u*rbinom((n + 500), size = 1, prob = rho)
  
  x_train <- x[train_idx,]
  x_test <- x[test_idx,]
  nt <- nrow(x_test)
  
  p1 <- ggplot() +
    geom_point(aes(x=1:length(y), y[order(true_f)]),alpha = I(.5)) +
    geom_line(aes(x = 1:length(y), y = true_f[order(true_f)])) +
    ggtitle("Test & Training") +
    labs(x = "ID (ordered by true potential)", y = "y")
  
  #to estimate sigma, fit a random forest
  temp <- data.frame(y = y[train_idx], x_train)
  rf <- randomForest(y ~ ., data = temp)
  rf_pred <- predict(rf, data.frame(x_test))
  #olr for baseline
  olr <- lm(y ~ ., data = temp)
  olr_pred <- predict(olr, data.frame(x_test))

  temp <- data.frame(y = log(pmax(predict(rf)-y[train_idx],0)+.2),x =x_train )#model.matrix(pmin(predict(olr, data.frame(x_train))-y_train,0) ~ cbind(pc_train))
  lambda_start <- (lm(y ~., data= temp) %>%coef())[c(1,2)]
  lambda_sd <- ((lm(y ~., data= temp) %>%vcov())%>%diag()%>%sqrt())[c(1,2)]*2 
  
  
  set.seed(36726231 + idx)  ## MCMC posterior sampling: set seed for reproducibility
  t1 <- Sys.time()
  bf = cbart(x_train,y[train_idx],nskip=burn,
                   sigdf=3, #default is 3
                   ntree = 100,
                   sigquant = .99,#default is .9
                   sigest = sqrt(mean(rf$mse))/2,
                   ndpost=nd,
                   x.test = x_test,
                   printevery=1000L,
                   sparse = TRUE,
                   ux.train = cbind(1, x_train[,1]), #training data for asymmetric betas
                   ux.test =cbind(1, x_test[,1]),
                   lambda_prop_sd = lambda_sd,
                   lambda_start = lambda_start,
                   gamma_prop_sd = 1) #testing data for asymmetric betas
  t2 <- Sys.time()
  bart_pred <- apply(bf$potential.test, 2, mean)
  hold$BART_time[idx] <- difftime(t2, t1, units='secs')
  hold[idx,c("OLR_MSE",   "RF_MSE",    "BART_MSE","FMT_MSE",
             "OLR_fcor",   "RF_fcor",   "BART_fcor")] <- 
    c(mean((true_f[test_idx]-olr_pred)^2),
      mean((true_f[test_idx]-rf_pred)^2),
      mean((true_f[test_idx]-bart_pred)^2),
      mean((true_f[test_idx]-y[-train_idx])^2),
      
      cor( true_f[test_idx], y=olr_pred),
      cor( true_f[test_idx], y=rf_pred),
      cor( true_f[test_idx], y=bart_pred)
    )
  all_preds[[idx]] <- data.frame(true_potential = true_f[test_idx],
                                 BART = bart_pred,
                                 OLR = olr_pred,
                                 RF = rf_pred,
                                 y = y[test_idx])
  
  p2 <- ggplot() + geom_point(aes(x = true_f[test_idx], y=bart_pred)) +
    geom_abline(aes(slope = 1, intercept = 0)) +
    geom_text(aes(x = mean(true_f), y = mean(true_f) + 2*sd(true_f), 
                  label = round(cor( true_f[test_idx], y=bart_pred),3))) +
    labs(x = "True Potential", y = "Posterior Mean") +
    ggtitle("BART")
  
  p3 <- ggplot() + geom_point(aes(x = true_f[test_idx], y=rf_pred)) +
    geom_abline(aes(slope = 1, intercept = 0)) +
    geom_text(aes(x = mean(true_f), y = mean(true_f) + 2*sd(true_f), 
                  label = round(cor( true_f[test_idx], y=rf_pred),3))) +
    labs(x = "True Potential", y = "Predicted") +
    ggtitle("RF")
  
  p4 <- ggplot() + geom_point(aes(x = true_f[test_idx], y=olr_pred)) +
    geom_abline(aes(slope = 1, intercept = 0)) +
    geom_text(aes(x = mean(true_f), y = mean(true_f) + 2*sd(true_f), 
                  label = round(cor( true_f[test_idx], y=olr_pred),3))) +
    labs(x = "True Potential", y = "Predicted") +
    ggtitle("OLR")
  
  
  p5 <- grid.arrange(p1,p2,p3,p4, nrow = 1)
  if(idx %% 5 == 0){
  ggsave(paste0("Simulation_Plots/sim",true_udist, idx, ".pdf"),plot = p5, width = 14, height = 6)
  }
  print(hold)
}

#merge on experiment parameters
hold2 <- merge(G, hold, by = "idx")
#save results
saveRDS(hold2, paste0(true_udist, "_results.RDS"))
saveRDS(all_preds, paste0(true_udist,"_all_preds.RDS"))



