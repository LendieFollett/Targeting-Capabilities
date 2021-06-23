
#EMPIRICAL STUDY FUNCTIONS

get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

#rank1 = first rank (involves CBART), rank2 = second rank
#e.g., rank1 = "BART_ranking_food", rank2 = "OLR_ranking_consumption"
#r = poverty rate
#e.g., r = 0.3
get_binary_comparison <- function(adata_preds, rank1, rank2, r){
  
  #within each subvillage,
  #is the rank less than or equal to the prop treated
  adata_preds <- adata_preds %>%
    group_by(village, province, district, subdistrict) %>%
    mutate(m1_cutoff = (!! sym(rank1) <= prop_treated),
           m2_cutoff = (!! sym(rank2) <= prop_treated)) %>%
    ungroup()
  
  
  binary_xs <- setdiff(m3[m3 != "hhsize_ae"], c("depratio", "sschild", "jschild", 
                "eschild", "age04", "hhage", "hhsize", "pcfloor"))
  numeric_xs <- c("depratio", "sschild", "jschild", 
                  "eschild", "age04", "hhage", "hhsize", "pcfloor")
  
  adata_preds <- adata_preds %>%
    mutate_at(binary_xs, function(x){(x - mean(x))/(1*sd(x))}) %>%
    mutate_at(numeric_xs, function(x){x*2}) %>%
    mutate(category=paste0(m1_cutoff,  m2_cutoff))%>%
    mutate(category = factor(category, levels = c("FALSEFALSE",  "FALSETRUE",  "TRUEFALSE",   "TRUETRUE"),
                             labels =           c("Neither",     "Not CBART",   "CBART Only",   "Both")))
          

  
  
  adata_preds2 <- adata_preds %>%
    select(c(m3[m3!="hhsize_ae"], m1_cutoff, m2_cutoff, category))%>%
  #adata_preds[,!colnames(adata_preds)%in%c("hhsize2","hhage2","hybrid", "community", "pmt")] %>%
    group_by(category, m1_cutoff, m2_cutoff)%>%
    summarise_at(vars(m3[!m3%in%c("hhsize2","hhage2","hhsize_ae")]), mean)%>%ungroup()
  
  order <- adata_preds2 %>%
    subset(category %in% c("CBART Only", "Not CBART")) %>% 
    select(m3[!m3%in%c("hhsize2","hhage2", "hhsize_ae")]) %>%
    apply(2, diff)%>%melt() %>% select(value) %>%
    order
  
  p2 <- adata_preds2%>%
     mutate(subset = paste0("BART=",m1_cutoff,"other=",m2_cutoff),
           agree = category %in% c("Both", "Neither"))%>%
    subset(agree == FALSE) %>%
    mutate(Error = factor(subset, levels = c( paste0("BART=","FALSE","other=","TRUE"), paste0("BART=","TRUE","other=","FALSE")),
                          labels = c("IE", "EE"))) %>%
    dplyr::select(-c(agree, subset)) %>%
    melt(id.vars = c("category", "m1_cutoff", "m2_cutoff", "Error")) %>%
    select(c(Error, variable, value)) %>%
    tidyr::spread(Error, value) %>%
    mutate(diff = EE-IE) #CBART ONLY GROUP MEAN - OTHER ONLY GROUP MEAN
    
  
  return(p2)
}

#list = variable classification from variables.csv to be shown
#method_list = models to be shown
binary_comparison_plot <- function(list, method_list){
ie_results <- list()
j=0
for (method in c( "OLR", "RF", "CBT", "FMT", "FMT_AE", "OLR_AE")){
  for (response in c("nonfood", "nonstaple", "diversity")){
    j=j+1
    if(method == "FMT"){
      ie_results[[j]]<- get_binary_comparison(adata_preds=adata_preds_all, 
                                              rank1=paste0("BART_ranking_",response),
                                              rank2 ="FMT_ranking_consumption",
                                              r = 0.3)
    }else if(method == "FMT_AE"){
      ie_results[[j]]<- get_binary_comparison(adata_preds=adata_preds_all, 
                                              rank1=paste0("BART_ranking_",response),
                                              rank2 ="FMT_ranking_consumption_ae",
                                              r = 0.3) 
    }else if(method == "OLR_AE"){
      ie_results[[j]]<- get_binary_comparison(adata_preds=adata_preds_all, 
                                              rank1=paste0("BART_ranking_",response),
                                              rank2 = paste0(method,"_ranking_consumption_ae"),
                                              r = 0.3) 
    #}else if(method == "FMT_same"){
    #  ie_results[[j]]<- get_binary_comparison(adata_preds=adata_preds_all, 
    #                                           rank1=paste0("BART_ranking_",response),
    #                                           rank2 = paste0("FMT_ranking_",response),
    #                                           r = 0.3) 
    # }else if(method == "OLR_same"){
    #   ie_results[[j]]<- get_binary_comparison(adata_preds=adata_preds_all, 
    #                                           rank1=paste0("BART_ranking_",response),
    #                                           rank2 = paste0("OLR_ranking_",response),
    #                                           r = 0.3) 
      }else{
      ie_results[[j]]<- get_binary_comparison(adata_preds=adata_preds_all, 
                                              rank1=paste0("BART_ranking_",response),
                                              rank2 = paste0(method,"_ranking_consumption"),
                                              r = 0.3) 
    }

    ie_results[[j]]$method <- method
    ie_results[[j]]$response <- response
    
  }
}

ie_results <- do.call(rbind, ie_results)

variables <- read.csv("variables.csv")

names_order <- ie_results %>% 
  merge(variables, by.x = "variable",by.y = "Variable.Name", all.x = TRUE, all.y = FALSE) %>%
  subset(method %in% c('OLR',"RF", "FMT" ) & variable != "hhmalemarr") %>% 
  group_by(Variable.Definition2)%>%
  summarise(diff = mean(diff))%>%ungroup()

order <- names_order%>% 
dplyr::select(diff) %>%order
  
names <-names_order$Variable.Definition2[!names_order$Variable.Definition2%in%c("HH size squared","Age of head of HH squared", "Head of HH is male and married")]



ie_results <- ie_results %>%
  merge(variables, by.x = "variable",by.y = "Variable.Name", all.x = TRUE, all.y = FALSE) %>%
  subset(variable != "hhmalemarr") %>%
  mutate(variable = factor(variable, levels = names[order]),
         response = factor(response, levels = c("nonfood", "nonstaple", "diversity"),
                           labels = c("Non-food", "Non-staple", "Dietary diversity"))) %>% 

  mutate(method = factor(method, levels = c("CBT", "RF", "FMT", "FMT_AE", "OLR", "OLR_AE", "FMT_same", "OLR_same"),
                         labels = c("CBT", "RF", "FMT", "FMT (AE)", "OLR", "OLR (AE)", "FMT (same)", "OLR (same)")))

  p <- ie_results %>% subset(category %in% c(list) &
                             method %in% c(method_list))%>%
    mutate(Variable.Definition2 = factor(Variable.Definition2, levels = names[order]))%>%
    ggplot()  +
  geom_hline(aes(yintercept = 0), colour = "grey")+
   geom_point(aes(x = Variable.Definition2, y = diff, colour = method, shape = method))+
   geom_line(aes(x = Variable.Definition2, y = diff, colour = method, group = method))+
   facet_grid(~response)+#, scales = "free_y", space ="free_y")+
   theme_bw()+
   labs(x = "Characteristic", 
        y = "Difference in Averages Between EE and IE Classes")+
   scale_colour_grey("Method")  +
  scale_shape("Method")+
  coord_flip()+theme(legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid'), 
                    legend.text = element_text(size = 18), legend.title = element_text(size = 18),
                    text = element_text(size = 18)) 
#p
#write.csv(ie_results, "sd_characteristics.csv")
return(p)
}

#rank1 is CBART based on food security
#rank2 is FMT or PMT method
#rank2_name is name of FMT or PMT method for labelling
#fsi is food security indicator column name used to assess source of disagreement
ie_decomp <- function(adata_preds, fsi, rank1, rank2, rank2_name, r){
  adata_preds <- adata_preds %>%
    group_by(village, province, district, subdistrict) %>%
    mutate(m1_cutoff = (!! sym(rank1) <= prop_treated),
           m2_cutoff = (!! sym(rank2) <= prop_treated),
           fsi_cutoff = as.logical(!! sym(fsi) <= quantile(!! sym(fsi),prop_treated))) %>%
    ungroup()
  

  adata_preds <- adata_preds %>%
    mutate(category=paste0(m1_cutoff,  m2_cutoff))%>%
    mutate(category = factor(category, levels = c("FALSEFALSE",  "FALSETRUE",  "TRUEFALSE",   "TRUETRUE"),
                             labels =           c("Neither",     "Not CBART",   "CBART Only",   "Both")))# %>%


#for denominators
  PR <- nrow(adata_preds%>%subset(m1_cutoff == TRUE)) #number of truly poor based on CBART
  BE <- nrow(adata_preds%>%subset(m2_cutoff == TRUE)) #number beneficiaries based on FMT
  
  adata_preds2 <- rbind(adata_preds %>% 
                      subset(category %in% c("CBART Only")) %>%
                        group_by(category) %>% 
                      summarise(fsi_poor = sum(fsi_cutoff)/PR,#proportion classified as poor via FSI
                                fsi_nonpoor = sum(!fsi_cutoff)/PR,
                                er = length(fsi_cutoff)/PR) %>% #exclusion error rate EER
                        ungroup(),
                      adata_preds %>% 
                        subset(category %in% c("Not CBART")) %>%
                        group_by(category) %>% 
                        summarise(fsi_poor = sum(fsi_cutoff)/BE,#proportion classified as poor via FSI
                                  fsi_nonpoor = sum(!fsi_cutoff)/BE,
                                  er = length(fsi_cutoff)/BE)%>% #inclusion error rate IER
                        ungroup())
  
  adata_preds3 <- adata_preds %>% 
    group_by(m1_cutoff) %>%
    summarise(Bene = sum(m2_cutoff),
              NonBene = sum(!m2_cutoff)) 
  
  CI <- adata_preds3 %>% subset(m1_cutoff == TRUE) %>%select(Bene)
  IE <- adata_preds3 %>% subset(m1_cutoff == FALSE) %>%select(Bene)
  
  PR <- adata_preds3 %>% subset(m1_cutoff == TRUE) %>%sum
  NP <- adata_preds3 %>% subset(m1_cutoff == FALSE) %>%sum
  
  TD <- CI/PR - IE/NP
  
  final <- data.frame(IE_diversity = adata_preds2 %>% subset(category == "Not CBART") %>%select(fsi_nonpoor) %>% as.numeric,
                      IE_agency    = adata_preds2 %>% subset(category == "Not CBART") %>%select(fsi_poor)%>% as.numeric,
                      ER1          = adata_preds2 %>% subset(category == "Not CBART") %>%select(er)%>% as.numeric,
                      EE_diversity = adata_preds2 %>% subset(category == "CBART Only") %>%select(fsi_poor)%>% as.numeric,
                      EE_agency    = adata_preds2 %>% subset(category == "CBART Only") %>%select(fsi_nonpoor)%>% as.numeric,
                      ER2          = adata_preds2 %>% subset(category == "CBART Only") %>%select(er)%>% as.numeric,
             TD = TD %>% as.numeric)
  
  return(final)
}

get_binary_table <- function(){

ie_results <- list()
j=0
for (method in c( "OLR", "RF", "CBT", "FMT")){
  for (response in c("nonfood", "nonstaple", "diversity")){
    j=j+1
    ie_results[[j]]<- ie_decomp(adata_preds=adata_preds_all, fsi=response, 
              rank1=paste0("BART_ranking_", response), 
              rank2=paste0(method,"_ranking_consumption"), 
              rank2_name=method, 
              r=0.3) 
    ie_results[[j]]$method <- method
    ie_results[[j]]$response <- response
    
  }
}

ie_results <- do.call(rbind, ie_results)
ie_results_mean <- ie_results%>%select(-c(response))%>% group_by(method)  %>% summarise(across(everything(), mean))

rbind(ie_results, data.frame(ie_results_mean, response = "Average")) %>% 
  arrange(method, response) %>%
  mutate(response = factor(response,
                           levels = c("nonfood", "nonstaple", "diversity", "Average"),
                           labels = c("Non-food consumption", "Non-staple consumption", "Dietary diversity", "Average"))) %>%
  mutate(row = paste0("(", 1:16, ")")) %>%
  relocate(where(is.character),where(is.factor), .before = where(is.numeric)) %>%
  xtable() %>% print(include.rownames=FALSE)

}

#SIMULATION STUDY FUNCTIONS
get_onex_simulation_study_plot<- function(which){
  p2 <- ggplot() + 
    geom_point(aes(x = x, y=y,fill = poor),colour = "grey30", shape=21, stroke=1.1, size = 2, alpha = I(.6),data = subset(temp3, variable == which)) +
    geom_line(aes(x = x, y = value), linetype = 1, size = .75,data = subset(temp3, variable !="FMT" & variable == which), colour = "black")+
    geom_line(aes(x = x, y = True), linetype = 2, size = .75, data = subset(temp3,variable == which))+
    labs(x = "x", y = "y") +theme_bw() +
    scale_y_continuous(limits = c(-20, 5))+
    scale_fill_manual("Inclusion\nstatus",values = c("grey30", "white"))+
    geom_vline(aes(xintercept = .30), data = temp3) +
    theme(legend.position = "none")
  print(p2)
}

get_simulation_study_plot <- function(ex, parameter, ex_label){
  MSE_data %>% subset(experiment == ex)%>%
    mutate(experiment = ex_label)%>%
    ggplot(aes(x = (parameter_value), y = MSE, linetype = Model, group =Model)) + 
    facet_grid(.~experiment, scales = "free")+
    stat_summary(geom = "line", fun = mean) +
    stat_summary(geom = "point", fun = mean) +
    theme_bw() +
    labs(x = parameter) +
    theme(legend.position = "none") +
    scale_y_continuous(limits = c(0, max(MSE_data$MSE))) +
    scale_linetype_manual("Model",values=c("solid","longdash", "dotted", "dotdash"))
  #scale_x_discrete(breaks = unique(MSE_data$parameter_value[MSE_data$experiment == ex]),
  #                  labels = labels)
}