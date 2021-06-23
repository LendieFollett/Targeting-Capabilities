rm(list = ls())
library(devtools)
install_github("LendieFollett/CBART")
library(CBART)
library(randomForest)
library(ggplot2)
library(dplyr)
library(reshape2)
library(caret)
library(xtable)
source("plotting_functions.R")

idx = 10
beta <-c(-2, 4)
rho <- .7
sig <-1
p <-1
n <- 1000
ntest = 500

train_idx <- 1:(n)
test_idx <- (n+1):(n + ntest)

set.seed(352833 + idx)
x <- matrix(rbeta((n + ntest)*p, 1,1), ncol = p, byrow = FALSE)
#true_lambda <- exp(true_beta0)
true_f <- rep(NA, n)
true_f2 <- rep(NA, n)
true_u <- rep(NA, n)

  for (i in 1:(n + ntest)){
    set.seed(2388520 + i + idx*500)
    true_u[i] <- rexp(1, rate = 1/exp(beta[1] + beta[2]*x[i,1]))
    set.seed(2352 + i + idx*500)
    epsilon <- rnorm(1, 0, sig)
    true_f[i] <- sin(pi*(x[i,1]-0.5)) 
    true_f2[i] <- true_f[i] + epsilon
  }

set.seed(23589  + idx*500)
y <- true_f2  - true_u*rbinom((n + ntest), size = 1, prob = rho)

x_train <- x[train_idx]%>%as.matrix()
x_test <- x[-train_idx]%>%as.matrix()
nt <- nrow(x_test)

p1 <- ggplot() +
  geom_point(aes(x=x, y, colour = (1:length(y)%in%train_idx)),alpha = I(.5)) +
  geom_line(aes(x = x, y = true_f)) +
  ggtitle("Test & Training") +
  labs(x = "ID (ordered by true potential)", y = "y")
p1
#exp(beta[1] + beta[2]*x[,1])
#to estimate sigma, fit a random forest
temp <- data.frame(y = y[train_idx], x=x_train,x2 = x_train^2,x3 = x_train^3)
rf <- randomForest(y ~ x, data = temp, ntree = 1000)
plot(rf)
rf_pred <- predict(rf, data.frame(x=x_test))
mean((rf_pred- y[-train_idx])^2)

#quantile rf
#qrf <- quantregForest(x_train, y[train_idx],keep.inbag=TRUE)
#rf_pred1 <- predict(qrf, quantiles = c(.45,.50, .55))

#olr for baseline
olr <- lm(y ~ x, data = temp)
olr_pred <- predict(olr, data.frame(x=x_test,x2 = x_test^2,x3 = x_test^3))
mean((olr_pred- y[-train_idx])^2)
temp <- data.frame(y = log(pmax(predict(rf, data.frame(x=x_train))-y[train_idx],0)+.2),x =x_train )#model.matrix(pmin(predict(olr, data.frame(x_train))-y_train,0) ~ cbind(pc_train))
#lambda_start <- solve(t(temp)%*%temp)%*%t(temp)%*%log(pmax(predict(rf, data.frame(x_train))-y_train,0)+.2)
lambda_start <- lm(y ~., data= temp) %>%coef()
lambda_sd <- ((lm(y ~., data= temp) %>%vcov())%>%diag()%>%sqrt())*2 

set.seed(367231 + idx)  ## MCMC posterior sampling: set seed for reproducibility

bf = cbart(x_train,y[train_idx],nskip=50000,
                 sigdf=3, #default is 3
                 ntree = 50,
                 sigquant = .9,#default is .9
                 sigest = sqrt(mean(rf$mse))/2,
                 ndpost=50000,
                 x.test = x_test,
                 printevery=1000L,
                 sparse = TRUE,
                 ux.train = cbind(1, x_train), #training data for asymmetric betas
                 ux.test =cbind(1, x_test),
                 lambda_prop_sd = lambda_sd,
                 #lambda_prior_sd = 2.5,
                 lambda_start = lambda_start,
                 gamma_prop_sd = 1) #testing data for asymmetric betas


bart_pred <- apply(bf$potential.test, 2, mean)


temp2 <- data.frame(x = x_test, y = y[test_idx], BART = bart_pred, RF = rf_pred, OLR = olr_pred,True = true_f[test_idx]) %>%
  melt(id.vars = c("x", "y", "True")) %>%
  mutate(variable = factor(variable, levels = c("BART","RF", "OLR"),
                           labels = c("CBART","RF", "OLR"))) %>%
  group_by(variable) %>%
  mutate(poor =  value < quantile(value, .3)) %>% ungroup

temp3 <- rbind(temp2, 
               data.frame(x = x_test, y = y[test_idx], True = true_f[test_idx], variable = "FMT", value = y[test_idx], poor = (y[test_idx] < quantile(y[test_idx], .2))))

temp3$poor <- factor(temp3$poor, levels = c("TRUE", "FALSE"),
                                  labels = c("Yes", "No"))

#FIGURE 4: The functionings-capabilities distinction with a single explanatory variable
get_onex_simulation_study_plot("CBART") + theme(legend.position = c(.11,.2), legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid'), 
                          legend.text = element_text(size = 18), legend.title = element_text(size = 18),
                          text = element_text(size = 18))
ggsave(paste0("Simulation_Plots/one_x_BART.pdf")) 
get_onex_simulation_study_plot("OLR")+ theme(text = element_text(size = 18))
ggsave(paste0("Simulation_Plots/one_x_OLR.pdf"))
get_onex_simulation_study_plot("RF")+ theme(text = element_text(size = 18))
ggsave(paste0("Simulation_Plots/one_x_RF.pdf"))
get_onex_simulation_study_plot("FMT")+ theme(text = element_text(size = 18))
ggsave(paste0("Simulation_Plots/one_x_FMT.pdf"))

caps <- data.frame(
bart_pred_capl = apply(bf$cap.test, 2, quantile, .025),
bart_pred_capu = apply(bf$cap.test, 2, quantile, .975),
potential = bart_pred,
y = y[-train_idx],
x = x_test)
p3 <- ggplot(data = caps) + 
  geom_errorbar(aes(x = x, ymin = bart_pred_capl, ymax = bart_pred_capu),colour = "grey")+
  geom_point(aes(x = x, y=y),alpha = I(.3)) +
  geom_line(aes(x = x, y = potential), size = .75)+
  labs(x = "x", y = "y") +theme_bw() +
  scale_linetype("") #+scale_y_continuous(limits = c(-15, 15))
p3
ggsave(paste0("Simulation_Plots/one_x_cap.pdf"),plot = p3)


apply(bf$lambda, 2, mean)
plot(bf$lambda[,1])
plot(bf$lambda[,2])
plot(bf$rho)
mean(bf$rho[-c(1:1000)])
apply(bf$lambda, 2,function(x){length(unique(x))/length(x)})


temp2 <- data.frame(x = x_test, y = y[test_idx], BART = bart_pred, RF = rf_pred, OLR = olr_pred,True = true_f[test_idx])

all_summary<-rep(list(NA),4)
power = 1 #power in Foster-Greer-Thorbecke
j = 0
for(percent in c(.1,.2,.3,.4)){
  j=j+1
  true_cutoff <- quantile(temp2$True, percent)
  true_class <- ifelse(temp2$True < true_cutoff, 1, 0)%>%as.factor()
  
  BART_cutoff <- quantile(temp2$BART, percent)
  BART_class <- ifelse(temp2$BART < BART_cutoff, 1, 0)%>%as.factor()
  
  OLR_cutoff <- quantile(temp2$OLR, percent)
  OLR_class <- ifelse(temp2$OLR < OLR_cutoff, 1, 0)%>%as.factor()
  
  RF_cutoff <- quantile(temp2$RF, percent)
  RF_class <- ifelse(temp2$RF < RF_cutoff, 1, 0)%>%as.factor()
  
  FMT_cutoff <- quantile(temp2$y, percent)
  FMT_class <- ifelse(temp2$y < FMT_cutoff, 1, 0)%>%as.factor()

  
  all_summary[[j]] <- data.frame(model = c("BART", "OLR",  "RF", "FMT"),
                                 percent = percent,
                                 rbind(confusionMatrix(BART_class, true_class,positive = "1")$byClass,
                                       confusionMatrix(OLR_class, true_class,positive = "1")$byClass,
                                       confusionMatrix(RF_class, true_class,positive = "1")$byClass,
                                       confusionMatrix(FMT_class, true_class,positive = "1")$byClass)
  )
  
}


all_summaryd <- do.call("rbind",all_summary)

all_summaryd$ER <- 1-all_summaryd$Sensitivity #Y
all_summaryd$IE <- 1-all_summaryd$Precision#SAME AS EE
all_summaryd$TD <- all_summaryd$Sensitivity - (1-all_summaryd$Specificity)

print(xtable(all_summaryd[,c("model","percent","ER", "TD")]), include.rownames = FALSE)




