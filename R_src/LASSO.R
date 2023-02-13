LASSO = function(fulldata, export = F, fit_svm = F, fit_rf = F, foldfreq = 0.6) {
  retenv = new.env()
  #############################################################

  lasso_models = list()
  
  auc_RF = c()
  auc_RF_permute = c()
  auc_svm = c()
  auc_svm_permute = c()
  
  
  #This part is process for the lasso regression
  auc_data = c()
  auc_data_permute = c()
  
  #Create the variable to hold values:
  
  # this is the default for createFolds anyway
  k = 10
  numRuns = 10
  
  for(runCount in 1:numRuns) {
    
    selFeat1 <- c()
    compSelFeat1 <- c()
    chosenFeat1 <- c()
    total_oriY <- c()
    total_oriYhat <- c()
    
    selFeat2 <- c()
    compSelFeat2 <- c()
    chosenFeat2 <- c()
    total_oriY <- c()
    total_permYhat <- c()
    
    append.Y = c()
    append.Y_permute = c()
    
    append_RF = c()
    append_RF_permute = c()
    append_svm = c()
    append_svm_permute = c()
    
    trueY = c()
    permuteY = c()
    selectedVars = c()
    total_permY = c()
    
    # scramble response vector
    permutedY = sample(fulldata$Group)
    
    #create kfold associate to a data frame.
    kfolds <- createFolds(y = fulldata$Group, k=k, list = FALSE, returnTrain = FALSE)
    myData <- cbind(fulldata, kfolds)
    
    ##############################################
    
    for (select in 1:k) {
      #Create different subdata by different fold.
      chooseFold <- which(kfolds == select)
      
      # copy original data
      newData = myData
      
      for (n_run in 1:2) {
        
        # select training data
        
        #take out the random row in newData and put into train, 
        #and the test data only contain the one which the train data does not have
        train <- newData[-chooseFold, ]
        test <- newData[chooseFold, ]
        
        #remove the kfold column which is the last column out of the data to calculate
        train <- train[, -ncol(train)]
        test <- test[, -ncol(test)]
        
        #Since X data will be constant
        X_train <- as.matrix(train[, -1])
        
        # not used until later
        X_test <- as.matrix(test[, -1])
        
        y = train$Group
        
        if (n_run == 2) {
          newData$Y = permutedY
        }
        
        all_signif_vars = c()
        
        for (num_FS in 1:numRuns) {
          
          # run cv lasso to find best lambda
          lasso_run <- cv.glmnet(x=X_train, y=y, alpha=1, nfolds = k,
                                 standardize = TRUE)
          
          #run lasso on training data with best lambda
          lasso_run_lmin <- glmnet(x=X_train, y=y, alpha=1, 
                                   lambda = lasso_run$lambda.min)
          
          c = coef(lasso_run_lmin)
          inds = which(c!=0)
          
          #remove intercept from the list
          tmp_variables <- tail(row.names(c)[inds], -1)
          #print(tmp_variables)
          
          all_signif_vars<-c(all_signif_vars,tmp_variables)
          # save features
          
          # lasso_predict1 = predict(lasso_run_lmin, newx = X_test)
        }
        
        # select features based on frequency
        freq = sort(table(all_signif_vars),decreasing=TRUE)/numRuns
        
        # only pick features which show up above frequency threshold
        variables = names(which(freq > foldfreq))
        
        # variables = all_signif_vars
        #subset train and test data correspond to the chosen features
        selTrain <- train[, (names(train) %in% variables)]
        newTrain <- cbind.data.frame(train$Group, selTrain)
        
        selTest <- test[, (names(test) %in% variables)]
        newTest <- cbind.data.frame(test$Group, selTest)
        
        colnames(newTrain)[1] <- "Y"
        colnames(newTest)[1] <- "Y"
        
        # for true Y and permuteY, put in separate lists but both named Y
        if (n_run == 1) {
          trueY = append(trueY, newTest$Y)
          Y = trueY
        } else {
          permuteY = append(permutedY, newTest$Y)
          Y = permutedY
        }
        
        if (fit_svm){
        # svm
          svmfit = svm(y= newTrain$Y, x=as.matrix(newTrain[,-1]), kernel="linear", cost=10, scale=FALSE)
          yhat.SVM = predict(svmfit, newdata = newTest[, -1], type = "response")
        }
        
        if (fit_rf) {
        # random forest
          RFfit <- randomForest(y= newTrain$Y, x=as.matrix(newTrain[,-1]), importance=TRUE, ntree = 100)
          yhat.RF = predict(RFfit, newdata = newTest[, -1], type = "response")
        }
        
        # measuring yhat
        lasso_run2 <- glmnet(y= newTrain$Y, x=as.matrix(newTrain[,-1]), alpha=1,
                             lambda = lasso_run$lambda.min)
        
        yhat <-predict(lasso_run2, newx=as.matrix(newTest[,-1], type = "response"))
        
        # store yhat in appropriate location depending on mode
        if (n_run == 1) {
          append.Y = append(append.Y, yhat)
          selFeat1 <- c(selFeat1, variables)
          # compSelFeat1 <- c(compSelFeat1,selFeat1[select == k])
          compSelFeat1 <- c(compSelFeat1,selFeat1)
          
          
          lasso_models[[runCount]] = list(fit = lasso_run2, variables = variables)
          #svm and rf
          if (fit_svm) {
            append_svm = c(append_svm, yhat.SVM)
            
          }
          if (fit_rf) {
            append_RF = c(append_RF, yhat.RF)
          }
          
        } else {
          append.Y_permute = append(append.Y_permute, yhat)
          selFeat2 <- c(selFeat2, variables)
          # compSelFeat2 <- c(compSelFeat2,selFeat2[select == k])
          compSelFeat2 <- c(compSelFeat2,selFeat2)
          
          #svm and rf
          if (fit_svm) {
            append_svm_permute = c(append_svm_permute, yhat.SVM)
            
          }
          if (fit_rf) {
            append_RF_permute = c(append_RF_permute, yhat.RF)
          }
        }
      }
      chosenFeat1 <- c(chosenFeat1,compSelFeat1)
      chosenFeat2 <- c(chosenFeat2,compSelFeat2)
      
    }
    if(length(unique(permutedY))>2){
      #measure the classification accuracy
      total_oriY <- c(total_oriY, append.Y)
      total_oriYhat <- c( total_oriYhat, trueY)
      
      total_permY <- c(total_permY, append.Y_permute)
      total_permYhat <- c(total_permYhat, permutedY)
    }
    
    # save AUC scores
    auc_data[runCount] = auc(trueY, append.Y)
    auc_data_permute[runCount] = auc(permutedY, append.Y_permute)
    
    if (fit_rf) {
      auc_RF[runCount] = auc(trueY, append_RF)
      auc_RF_permute[runCount] = auc(permutedY, append_RF_permute)
    }
    
    if (fit_svm) {
      auc_svm[runCount] = auc(trueY, append_svm)
      auc_svm_permute[runCount] = auc(permutedY, append_svm_permute)
    }
  }
  
  classifiedAccuracy <- confusionMatrix(as.factor(trueY), 
                                        as.factor(round(append.Y)))
  
  classifiedAccuracyPerm <- confusionMatrix(as.factor(permutedY), 
                                            as.factor(round(append.Y_permute)))
  
  # to get freq features as needed
  cf_ori = sort(table(compSelFeat1), decreasing = T)
  # cf_ori = names(which(cf_ori > foldfreq * ))
  
  chosenFeat_ori <- table(chosenFeat1) %>% as.data.frame() %>% arrange(desc(Freq))
  chosenFeat_perm <- table(chosenFeat2) %>% as.data.frame() %>% arrange(desc(Freq))
    if (export) {
      
    write.table(chosenFeat_ori, 
                file="output/chosenFeat_ori", 
                append = FALSE, sep = " ",
                row.names = FALSE, col.names = TRUE) 
    write.table(chosenFeat_perm, 
                file="output/chosenFeat_perm2", 
                append = FALSE, sep = " ",
                row.names = FALSE, col.names = TRUE) 
    
    # make dataframes for classification info
    class_accur = data.frame(cbind(t(classifiedAccuracy$byClass)))
    
    class_accurperm = data.frame(cbind(t(classifiedAccuracyPerm$byClass)))
    
    write.table(class_accur, 
                file="output/classifiedAccuracy", 
                append = FALSE, sep = " ",
                row.names = FALSE, col.names = TRUE) 
    write.table(class_accurperm, 
                file="output/classifiedAccuracyPerm", 
                append = FALSE, sep = " ",
                row.names = FALSE, col.names = TRUE) 
    
  }
  
  retenv$feats = list(chosenFeat = chosenFeat_ori, compSelFeat = cf_ori, neg = chosenFeat_perm)
  retenv$auc = list(lasso_auc = list(pred = auc_data, perm = auc_data_permute), 
                    svm_auc = list(pred = auc_svm, perm = auc_svm_permute), 
                    rf_auc = list(pred = auc_RF, perm = auc_RF_permute))
  
  retenv$models = lasso_models
  return(retenv)
}