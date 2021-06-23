source("plotting_functions.R")

true_udist <- "Normal"#"Exponential" #
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
hold2 <- readRDS(paste0(true_udist, "_results.RDS"))
all_preds <- readRDS(paste0(true_udist,"_all_preds.RDS"))

hold2_long <- hold2 %>% 
  dplyr::select(-c(idx, experiment.y, BART_time, OLR_fcor, RF_fcor, BART_fcor)) %>%
  melt(id.var = c("true_beta0", "true_beta1" ,"true_rho", "true_sigma_e", "true_P", "n", "true_udist", "experiment.x")) %>%
  melt(id.var = c("variable", "value", "experiment.x")) 
colnames(hold2_long) <- c("Model", "MSE","experiment", "parameter", "parameter_value")
MSE_data <- hold2_long%>%
  subset((experiment == 1 & parameter == "true_beta0") |
           (experiment == 2 & parameter == "true_beta1") |
           (experiment == 3 & parameter == "true_rho") | 
           (experiment == 4 & parameter == "true_sigma_e")|
           (experiment == 5 & parameter == "true_P")|
           (experiment == 6 & parameter == "n")) %>%
  mutate(parameter_value = as.numeric(parameter_value),
         Model = factor(Model, levels = c("BART_MSE","FMT_MSE", "OLR_MSE", "RF_MSE"),
                        labels = c("CBART","FMT", "OLR", "RF")))


m1 <- get_simulation_study_plot(1, expression(phi[0]), 1)+ labs(y = "MSE")+
  theme(legend.position = c(.35,.7), 
        legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid'), 
         legend.text = element_text(size = 18), legend.title = element_text(size = 18),
          text = element_text(size = 18)) 
m2 <- get_simulation_study_plot(2, expression(phi[1]), 2)+ labs(y = " ")+
  theme(axis.title.y = element_blank(), text = element_text(size = 18),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) 
m3 <- get_simulation_study_plot(3, expression(rho), 3) +  labs(y = " ")+
  theme(axis.title.y = element_blank(), text = element_text(size = 18),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) 
m4 <- get_simulation_study_plot(4, expression(sigma), 4)+  labs(y = " ")+
  theme(axis.title.y = element_blank(), text = element_text(size = 18),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) 
m5 <- get_simulation_study_plot(5, "P", 5)+  labs(y = " ")+
  theme(axis.title.y = element_blank(), text = element_text(size = 18),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) 
m6 <- get_simulation_study_plot(6, "n", 6)+  labs(y = " ")+
  theme(axis.title.y = element_blank(), text = element_text(size = 18),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) 


legend <- (qplot(parameter_value, MSE,geom = "line", linetype = Model, data = MSE_data) + theme_bw() ) %>%get_legend()
#FIGURE 5: AVERAGE MSE ACROSS PARAMETERIZATIONS FOR MULTIVARIATE EXPERIMENTS
MSE_plot <- grid.arrange(m1,m2,m3,m4,m5,m6, nrow = 1,
                         widths = c(1.1,1,1,1,1,1))
ggsave(paste0("Simulation_Plots/MSE", true_udist, ".pdf"),MSE_plot, width = 20, height = 5)


#------- classification --------

all_summary<-rep(list(NA),length(all_preds)*4)
power = 2 #power in Foster-Greer-Thorbecke
j  = 0
for (idx in 1:length(all_preds)){
  for(percent in c(.1,.2,.3,.4)){
    j = j + 1
    
    true_cutoff <- quantile(all_preds[[idx]]$true_potential, percent)
    true_class <- ifelse(all_preds[[idx]]$true_potential < true_cutoff, 1, 0)%>%as.factor()
    
    BART_cutoff <- quantile(all_preds[[idx]]$BART, percent)
    BART_class <- ifelse(all_preds[[idx]]$BART < BART_cutoff, 1, 0)%>%as.factor()
    
    OLR_cutoff <- quantile(all_preds[[idx]]$OLR, percent)
    OLR_class <- ifelse(all_preds[[idx]]$OLR < OLR_cutoff, 1, 0)%>%as.factor()
    
    RF_cutoff <- quantile(all_preds[[idx]]$RF, percent)
    RF_class <- ifelse(all_preds[[idx]]$RF < RF_cutoff, 1, 0)%>%as.factor()
    
    FMT_cutoff <- quantile(all_preds[[idx]]$y, percent)
    FMT_class <- ifelse(all_preds[[idx]]$y < FMT_cutoff, 1, 0)%>%as.factor()
    all_summary[[j]] <- data.frame(model = c("BART", "OLR",  "RF", "FMT"),
                                   percent = percent,
                                   sim = idx,
                                   rbind(confusionMatrix(BART_class, true_class,positive = "1")$byClass,
                                         confusionMatrix(OLR_class, true_class,positive = "1")$byClass,
                                         confusionMatrix(RF_class, true_class,positive = "1")$byClass,
                                         confusionMatrix(FMT_class, true_class,positive = "1")$byClass)
    )
    
  }
}


all_summaryd <- do.call("rbind",all_summary)

all_summaryd <- merge(all_summaryd, G, all.x = TRUE, all.y = FALSE,by.x = "sim", by.y = "idx")
#Sensitivity= P(beneficiary | true poor)
#Specificity= P(non-beneficiary | true non-poor)
#Precision = P(true poor|beneficiary)
#1-Precision = P(true non-poor|beneficiary)
#EE = E2/P = C/(A+C) = 1-sensitivity
#IE = E1/B = B/(A+B) = 1-precision
#TD = C1/P - E1/NP = A/(A+C)-B/(B+D)
#                  =sensitivity-(1-specificity)

all_summaryd$ER <- 1-all_summaryd$Sensitivity #Y
all_summaryd$IE <- 1-all_summaryd$Precision#SAME AS EE
all_summaryd$TD <- all_summaryd$Sensitivity - (1-all_summaryd$Specificity)

all_summaryd2 <- all_summaryd %>%
  dplyr::select(c(model, ER,IE,TD,true_P,    n, true_rho, true_beta0, true_beta1, true_sigma_e,  true_udist, experiment,percent)) %>%
  melt(id.var = c("model", "ER","IE",  "TD", "experiment","percent"))%>%
  subset((experiment == 1 & variable == "true_beta0") |
           (experiment == 2 & variable == "true_beta1") |
           (experiment == 3 & variable == "true_rho") | 
           (experiment == 4 & variable == "true_sigma_e")|
           (experiment == 5 & variable == "true_P")|
           (experiment == 6 & variable == "n")) %>%
  melt(id.var = c("model", "experiment", "variable", "value", "percent")) 
colnames(all_summaryd2) <- c("model", "experiment", "parameter", "parameter_value","percent", "class_measure", "value")

try_percent <- .4

p1 <- all_summaryd2 %>% 
  mutate(experiment = factor(experiment, levels = c(1,2,3,4,5,6),
                             labels = c("1", "2", "3","4", "5", "6")),
         model = factor(model, levels = c("BART","FMT", "OLR", "RF"),
                labels = c("CBART","FMT", "OLR", "RF"))) %>%
  group_by(class_measure, experiment, parameter_value, model, percent) %>%
  summarise(value = mean(value)) %>%
  ungroup %>%
  subset(percent == try_percent & class_measure %in% c("ER") )%>%
  ggplot(aes(x = as.numeric(parameter_value), y = value, linetype = model, group =interaction(percent, model))) + 
  facet_grid(.~experiment, scales = "free")+
  geom_line() +
  geom_point() +
  scale_colour_grey()+
  scale_linetype_manual("Model",values=c("solid","longdash", "dotted", "dotdash"))+
  #ggtitle(paste0("Cutoff = ", try_percent*100, "%")) +
  #labs(x = param_label, y = "") +
  theme_bw() +
  theme(legend.position = c(.07,.7), legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid'), 
        legend.text = element_text(size = 18), legend.title = element_text(size = 18),
        text = element_text(size = 18)) +
  labs(x = "", y= "ER")


p2 <-  all_summaryd2 %>% 
  mutate(experiment = factor(experiment, levels = c(1,2,3,4,5,6),
                             labels = c("1", "2", "3","4", "5", "6")),
         model = factor(model, levels = c("BART","FMT", "OLR", "RF"),
                        labels = c("CBART","FMT", "OLR", "RF"))) %>%
  group_by(class_measure, experiment, parameter_value, model, percent) %>%
  summarise(value = mean(value)) %>%
  ungroup %>%
  subset(percent == try_percent & class_measure %in% c("TD") )%>%
  ggplot(aes(x = as.numeric(parameter_value), y = value, linetype = model, group =interaction(percent, model))) + 
  facet_grid(.~experiment, scales = "free")+
  scale_colour_grey()+
  scale_linetype_manual(values=c("solid","longdash", "dotted", "dotdash"))+
  geom_line() +
  geom_point() +
  #ggtitle(paste0("Cutoff = ", try_percent*100, "%")) +
  #labs(x = param_label, y = "") +
  theme_bw() +
  theme(legend.position = "none",text = element_text(size = 18)) +
  labs(y = "TD",
       x = expression(paste(phi[0],  "                                       ",
                            phi[1],  "                                       ",
                            rho,     "                                       ", 
                            sigma,   "                                       ",
                            "P",     "                                       ",
                            "n"," ")),parse = TRUE) 
#FIGURE 6: ERROR RATES AND TARGETING DIFFERENTIALS FOR MULTIVARIATE EXPERIMENTS
p3 <- grid.arrange(p1,p2, nrow = 2)  
  
ggsave(paste0("Simulation_Plots/class",true_udist, try_percent*100,".pdf"),plot=p3, width = 18, height = 10)



