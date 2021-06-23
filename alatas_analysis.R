rm(list = ls())
library(Rcpp)
library(devtools)
install_github("LendieFollett/CBART")
library(CBART)
library(MASS)
library(ggplot2)
library(tibble)
library(GGally)
library(hrbrthemes)
library(viridis)
library(randomForest)
library(dplyr)
library(gridExtra)
library(reshape2)
library(tidytext)
library(coda)
library(BART)
library(caret)
library(xtable)
library(tidyr)
source("plotting_functions.R")

full_data <- read.csv("Updated Data/alatas.csv")
full_data$nonfood <- 100*full_data$nonfood
full_data$nonstaple <- 100*full_data$nonstaple

#groups of x variables
m1 <- c("hhsize","hhsize_ae","hhage","hhmale","hhmarried","hhage2", "hhsize2", "hhmalemarr",
        "hheduc2","hheduc3","hheduc4",
        "age04","higheduc2","higheduc3","higheduc4","depratio")
m2.1 <- c("pcfloor", "tfloor","twall", "toilet","water","lighting", "troof",
          "fcook","house", "ac","computer","radio","tv", "dvd","satellite", 
          "gas", "refrigerator", "bicycle", "motorcycle", "auto", "hp", 
          "jewelry","chicken","cow")
m2 <- c(m1, m2.1)
m3.1 <- c("credit","hhsector1", "hhsector2","hhsector3",
          "formal","informal", "eschild","jschild","sschild")
m3 <- c(m2,m3.1) #full collection


nskip <- 90000 #how many iterations for burn-in
ndpost <- 10000 #how many iterations to keep after burn-in

#run get_results() for each providence individually
#province:  12   33   73 
r12_1 <- get_results(12, m3)

r33_1 <- get_results(33, m3)

r73_1 <- get_results(73, m3)


x_vars <- full_data[,m3]
adata_preds_all <- rbind(r12_1$adata_preds,
                         r33_1$adata_preds,
                         r73_1$adata_preds)

#Add FMT test based on consumption: original and adult equivalency adjusted
adata_preds_all <- adata_preds_all %>%
  group_by(village, province, district, subdistrict) %>%
  mutate(FMT_ranking_consumption =rank(consumption)/length(village),
         FMT_ranking_consumption_ae =rank(consumption_ae)/length(village),
        FMT_rankingnon_food =rank(nonfood)/length(village),
        FMT_rankingnon_staple =rank(nonfood)/length(village),
        FMT_ranking_diversity =rank(nonfood)/length(village),
        prop_treated = mean(treated) ) %>%ungroup()


#Diagnostic checking
grid.arrange(qplot(1:ndpost,r33_1[[1]]$lambda[,1], geom = "line"),qplot(1:ndpost,r33_1[[1]]$lambda[,2], geom = "line"))
grid.arrange(qplot(1:ndpost,r33_1[[2]]$lambda[,1], geom = "line"),qplot(1:ndpost,r33_1[[2]]$lambda[,2], geom = "line"))
grid.arrange(qplot(1:ndpost,r33_1[[3]]$lambda[,1], geom = "line"),qplot(1:ndpost,r33_1[[3]]$lambda[,2], geom = "line"))
grid.arrange(qplot(1:ndpost,r33_1[[4]]$lambda[,1], geom = "line"),qplot(1:ndpost,r33_1[[4]]$lambda[,2], geom = "line"))
grid.arrange(qplot(1:ndpost,r33_1[[5]]$lambda[,1], geom = "line"),qplot(1:ndpost,r33_1[[5]]$lambda[,2], geom = "line"))

summary(apply(r12_1[[1]]$u.train, 2, mean))
summary(apply(r33_1[[5]]$u.train, 2, mean))
summary(apply(r73_1[[5]]$u.train, 2, mean))

#how do predictions compare
p1 <- qplot(OLR_consumption, BART_consumption,alpha=I(.5), data = adata_preds_all) +
  geom_abline(aes(slope = 1, intercept = 0)) +labs(y = "CBART consumption")
p2 <- qplot(OLR_food, BART_food,alpha=I(.5), data = adata_preds_all) +
  geom_abline(aes(slope = 1, intercept = 0))+labs(y = "CBART food")
p3 <- qplot(OLR_diversity, BART_diversity,alpha=I(.5), data = adata_preds_all) +
  geom_abline(aes(slope = 1, intercept = 0))+labs(y = "CBART diversity")
p4 <- qplot(OLR_nonfood, BART_nonfood,alpha=I(.5), data = adata_preds_all) +
  geom_abline(aes(slope = 1, intercept = 0))+labs(y = "CBART nonfood")
p5 <- qplot(OLR_nonstaple, BART_nonstaple,alpha=I(.5), data = adata_preds_all) +
  geom_abline(aes(slope = 1, intercept = 0))+labs(y = "CBART nonstaple")
p6 <- grid.arrange(p1,p2,p3,p4,p5)
p6
ggsave(paste0("OLRBART_pred_compare", "m3", ".pdf"),plot = p6, width = 8, height = 8)


#------inclusion/exclusion analysis---
#FIGURES 8, 9
binary_comparison_plot(list = "HH Characteristic", 
                       method_list = c("FMT", "CBT", "OLR", "RF"))
ggsave("sd_characteristics.pdf", height = 8, width = 14)

binary_comparison_plot(list = "HH Characteristic", 
                       method_list = c("FMT", "FMT (AE)", "OLR", "OLR (AE)"))
ggsave("sd_characteristics_ae.pdf", height = 8, width = 14)
variables <- read.csv("variables.csv")
binary_comparison_plot(list = unique(variables$category),
                       method_list = c("FMT", "CBT", "OLR", "RF"))
ggsave("sd_characteristics_all.pdf", height = 14, width = 20)


get_binary_table()
