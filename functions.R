
#makes predictions for olr, olr_ae, rf, bart, cbart
#uses x_train, x_test, y_train (for 1 particular response)
#returns list of (1) data frame of characteristics, all model predictions, and (2) the cbart object
get_preds <- function(y_name,const){
  
  temp <- data.frame(y = y_train, x_train)
  olr <- lm(y ~ ., data = temp)
  olr_pred <- predict(olr, data.frame(x_test))
  
  rf <- randomForest(y ~., data = temp)
  rf_pred <- predict(rf, data.frame(x_test))
  
  temp <- data.frame(y = y_train, x_train_ae)
  olr_ae <- lm(y ~ ., data = temp)
  olr_pred_ae <- predict(olr_ae, data.frame(x_test_ae))
  
  temp <- data.frame(y = log(pmax(predict(rf, data.frame(x_train))-y_train,0)+.2),x =pc_train )
  lambda_start <- lm(y ~., data= temp) %>%coef()
  lambda_sd <- ((lm(y ~., data= temp) %>%vcov())%>%diag()%>%sqrt())*const
  
  bf = cbart(x_train[,bart_idx],y_train,nskip=nskip,
             sigdf=3, #default is 3
             ntree = 50,
             sigquant = .99,#default is .9
             sigest = sqrt(mean(rf$mse))/2,
             ndpost=ndpost,
             x.test = x_test[,bart_idx],
             printevery=1000L,
             sparse = TRUE,
             ux.train = cbind(1, pc_train), #array(1, dim = c(nrow(x_train),2)),##training data for asymmetric betas
             ux.test = cbind(1, pc_test),#array(1, dim = c(nrow(x_test),2)),##testing data for asymmetric betas
             lambda_prop_sd=lambda_sd, #random walk jump size
             lambda_start = lambda_start,
             gamma_prop_sd = sqrt(mean(rf$mse))/5
  )
  
  
  
  bf2 = wbart(x_train,y_train,nskip=nskip,
              sigdf=3, #default is 3
              ntree = 50,
              sigquant = .9,#default is .9
              sigest = sqrt(mean(rf$mse)),
              ndpost=ndpost,
              printevery=1000L,
              x.test = x_test)
  
  bart_pred <- apply(bf$potential.test, 2, mean)
  bart2_pred <- bf2$yhat.test.mean
  
  temp <- data.frame(addl_test,
                     BART = bart_pred,
                     BART2 = bart2_pred,
                     OLR = olr_pred, 
                     RF = rf_pred,
                     OLR_AE = olr_pred_ae) %>%
    group_by(village, province, district, subdistrict) %>% #within village
    mutate(BART_ranking = rank(BART)/length(BART),       #CBART
           BART2_ranking = rank(BART2)/length(BART),     #BART
           OLR_ranking = rank(OLR)/length(BART),         #OLR
           OLR_AE_ranking = rank(OLR_AE)/length(OLR_AE), #OLR adult equivalency
           RF_ranking = rank(RF)/length(BART),           #RF
           CBT_ranking = rank(rank)/length(rank)        #CBT
    ) %>%ungroup()
  #modify names of newly added predictions, rankings to include y variable
  colnames(temp)[-c(1:ncol(addl_test))]<-paste(colnames(temp)[-c(1:ncol(addl_test))],y_name , sep = "_")
  return(list(temp = temp,bf=bf))
}

#sets up training/teseting data sets
#runs get_preds() for all y variables
get_results <- function(prov, x_names){
  data <- full_data%>%subset(province == prov)
  #training rows
  train_idx <- which(data$community != 1 & data$hybrid !=1)
  
  x <- data[,x_names]%>%mutate_at(intersect(x_names, c("depratio", "sschild", "jschild", 
                                                       "eschild", "age04", "hhage", "hhsize", "pcfloor")),
                                  function(x){(x - mean(x))/(2*sd(x))})%>%as.matrix()
  #adult equivalence adjustment
  x_ae <- data[,x_names]%>%mutate_at(intersect(x_names, c("depratio", "sschild", "jschild", 
                                                          "eschild", "age04", "hhage", "hhsize_ae", "pcfloor")),
                                     function(x){(x - mean(x))/(2*sd(x))})%>%as.matrix()
  
  x_train <- x[train_idx,];x_train_ae <- x_ae[train_idx,]
  x_test <- x[-train_idx,];x_test_ae <- x_ae[-train_idx,]
  
  #to remove squared terms, interactions for CBART model
  bart_idx <- which(!colnames(x_train) %in% c("hhage2", "hhsize2", "hhmalemarr"))
  addl_test <- data[-train_idx,(1:ncol(data))[!colnames(data) %in% x_names]]
  
  #PCA to summarise x matrix for asymmetric mean linear predictor
  pc_train1 <- prcomp(x_train[,bart_idx])$x[,c(1)]%>%as.data.frame() 
  m <- apply(pc_train1,2,mean); s <-  apply(pc_train1,2,sd)
  pc_train <- pc_train1%>%mutate_all(function(x){(x - mean(x))/(2*sd(x))}) %>%as.matrix()
  pc_test <- predict( prcomp(x_train[,bart_idx]),x_test)[,c(1)]%>%as.data.frame()%>%apply(1,function(x){(x - m)/(2*s)})
  
  
  y_train <- log(data[train_idx,"consumption"])
  y_test <- log(data[-train_idx,"consumption"])
  cons_out <- get_preds("consumption", 10)
  
  y_train <- log(data[train_idx,"food"])
  y_test <- log(data[-train_idx,"food"])
  food_out <- get_preds("food",10)
  
  y_train <- (data[train_idx,"nonfood"])
  y_test <- (data[-train_idx,"nonfood"])
  nonfood_out <- get_preds("nonfood",2)
  
  y_train <- (data[train_idx,"nonstaple"])
  y_test <- (data[-train_idx,"nonstaple"])
  nonstaple_out <- get_preds("nonstaple", 2)
  
  y_train <- (data[train_idx,"diversity"])
  y_test <- (data[-train_idx,"diversity"])
  diversity_out <- get_preds("diversity", 2)
  
  y_train <- log(data[train_idx,"consumption_ae"])
  y_test <- log(data[-train_idx,"consumption_ae"])
  cons_ae_out <- get_preds("consumption_ae", 2)
  
  adata_preds <- cbind(x_test,
                       cons_out$temp,
                       cons_ae_out$temp[,-c(1:ncol(addl_test))],
                       food_out$temp[,-c(1:ncol(addl_test))],
                       nonfood_out$temp[,-c(1:ncol(addl_test))],
                       nonstaple_out$temp[,-c(1:ncol(addl_test))],
                       diversity_out$temp[,-c(1:ncol(addl_test))])
  
  
  
  return(list(bfcons=cons_out[[2]],
              bffood = food_out[[2]],
              bfnonfood = nonfood_out[[2]],
              bfnonstaple = nonstaple_out[[2]],
              bfdiversity = diversity_out[[2]],
              adata_preds=adata_preds))
}

