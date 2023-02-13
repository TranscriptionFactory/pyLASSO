LASSO_SVM = function(fulldata, export = F, fit_rf = F, foldfreq = 0.6, alphaTune = 1.0) {
  retenv = new.env()
  #############################################################
  
  # save svm models
  svm_models = list()
  
  # auc values for random forest and svm
  auc_RF = c()
  auc_RF_permute = c()
  auc_svm = c()
  auc_svm_permute = c()
  
  # Confusion matrices
  cfm = list()
  cfm_perm = list()

  #Create the variable to hold values:
  
  # number of folds to create
  k = 10
  
  # number of repeats 
  numRuns = 10
  
  for(runCount in 1:numRuns) {
    #create kfold associate to a data frame.
    kfolds <- createFolds(y = fulldata$Group, k=k, list = FALSE, returnTrain = FALSE)
    myData <- cbind(fulldata, kfolds)
    
    # scramble response vector
    permutedY = sample(fulldata$Group)
    
    # variables chosen for true and permuted data
    chosenFeat <- c()
    chosenFeat_perm <- c()
    
    # temp variables for saving predicted groups
    append_RF = c()
    append_RF_permute = c()
    append_svm = c()
    append_svm_permute = c()
    
    # temp variables for true groups and permuted groups
    trueY = c()
    permuteY = c()
    
    ##############################################
    
    for (fold in 1:k) {
      #Create different subdata by different fold.
      currentFold <- which(kfolds == fold)
      
      # copy original data
      newData = myData
      
      for (n_run in 1:2) {
        
        # select training data
        
        #take out the random row in newData and put into train, 
        #and the test data only contain the one which the train data does not have
        train <- newData[-currentFold, ]
        test <- newData[currentFold, ]
        
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
        
        # run cv lasso to find best lambda
        lasso_run <- cv.glmnet(x=X_train, y=y, alpha=1, nfolds = k,
                               standardize = TRUE)
        
        #run lasso on training data with best lambda
        lasso_run_lmin <- glmnet(x=X_train, y=y, alpha=1, 
                                 lambda = lasso_run$lambda.min * alphaTune)
        
        c = coef(lasso_run_lmin)
        inds = which(c!=0)
        
        #remove intercept from the list
        variables <- tail(row.names(c)[inds], -1)
        
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
        
          # svm
        svmfit = svm(y= newTrain$Y, x=as.matrix(newTrain[,-1]), kernel="linear", cost=10, scale=FALSE)
        yhat.SVM = predict(svmfit, newdata = newTest[, -1], type = "response")
        
        if (fit_rf) {
          # random forest
          RFfit <- randomForest(y= newTrain$Y, x=as.matrix(newTrain[,-1]), importance=TRUE, ntree = 100)
          yhat.RF = predict(RFfit, newdata = newTest[, -1], type = "response")
        }
        
        # store yhat in appropriate location depending on mode
        if (n_run == 1) {
          chosenFeat <- c(chosenFeat, variables)
          
          #svm and rf
          append_svm = c(append_svm, yhat.SVM)
          svm_models[[runCount]] = list(fit = svmfit, variables = variables)
          
          if (fit_rf) {
            append_RF = c(append_RF, yhat.RF)
          }
          
        } else {
          chosenFeat_perm <- c(chosenFeat_perm, variables)
          
          #svm and rf
          append_svm_permute = c(append_svm_permute, yhat.SVM)
            
          if (fit_rf) {
            append_RF_permute = c(append_RF_permute, yhat.RF)
          }
        }
      }
    }

    # save AUC scores
    if (fit_rf) {
      auc_RF[runCount] = auc(trueY, append_RF)
      auc_RF_permute[runCount] = auc(permutedY, append_RF_permute)
    }
    
    auc_svm[runCount] = auc(trueY, append_svm)
    auc_svm_permute[runCount] = auc(permutedY, append_svm_permute)
  
  cfm[[runCount]] <- confusionMatrix(as.factor(trueY), 
                                        as.factor(round(append_svm)))
  
  cfm_perm[[runCount]] <- confusionMatrix(as.factor(permutedY), 
                                            as.factor(round(append_svm_permute)))
  }
  # to get freq features as needed
  cf_ori = sort(table(chosenFeat), decreasing = T)

  chosenFeat <- table(chosenFeat) %>% as.data.frame() %>% arrange(desc(Freq))
  chosenFeat_perm <- table(chosenFeat_perm) %>% as.data.frame() %>% arrange(desc(Freq))
  if (export) {
    
    write.table(chosenFeat, 
                file="output/chosenFeat", 
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
  
  retenv$feats = list(chosenFeat = chosenFeat, neg = chosenFeat_perm)
  retenv$auc = list(svm_auc = list(pred = auc_svm, perm = auc_svm_permute), 
                    rf_auc = list(pred = auc_RF, perm = auc_RF_permute))
  retenv$cf = list(cfm = cfm, cfm_perm = cfm_perm)
  
  retenv$models = svm_models
  return(retenv)
}