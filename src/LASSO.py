import pandas as pd
import numpy as np
from sklearn import metrics
from sklearn.metrics import roc_curve, auc
import math
from sklearn.metrics import confusion_matrix
import types
import random
from sklearn.model_selection import KFold
from sklearn.linear_model import Lasso
from sklearn import svm

def calculateAUC(trueY, permutedY,append_model, append_model_permute, unique_levels = 2):

   if len(pd.unique(trueY)) == 2:
    
    auc_model = auc(trueY, append_model)
  
    # auc_model_permute = auc(permutedY, append_model_permute)
    auc_model_permute = metrics.auc(permutedY, append_model_permute)
  
   else:
      auc_model = auc(roc_curve(trueY, append_model))
    
    # auc_model_permute = auc(multiclass.roc(permutedY, append_model_permute))
      auc_model_permute = auc(roc_curve(permutedY, append_model_permute))

   def normbetween(x, unique_levels) :
    # the number of levels is the range we want
    
      if (unique_levels <= 2) :
        return( 1 / (1 + math.exp(-x)))
    
      else :
        norm_factor = 1/(unique_levels - 1)
      
        return( (1 / (norm_factor + math.exp(-x))) + 1 )

   cfm = confusion_matrix(pd.Categorical(round(append_model)),
                       pd.Categorical(trueY))

   cfm_permute = confusion_matrix(pd.Categorical(round(append_model_permute)),
                               pd.Categorical(permutedY))

   return(list(auc_model = auc_model, auc_model_permute = auc_model_permute,
              cfm = cfm, cfm_permute = cfm_permute))

def LASSO_Grid(fulldata, export = F, foldfreq = 0.6, alphaValues = [1.0]) :
  retenv = types.SimpleNamespace()
  #############################################################

  random.seed(100)
  alphaValues = float(alphaValues)
  
  # number of folds to create
  k = 10
  
  # number of repeats 
  numRuns = 10
  
  # initialize results list
  gridResults = list(numRuns)
  ##################################################################################################################################
  for runCount in 1:numRuns:
    
    gridValues = list()
    for val in alphaValues:
      gridValues[[len(gridValues) + 1]] = list(alpha = val)
    
    
#create kfold associate to a data frame.
  z = fulldata['Group']

# specify the number of folds
  k = 5

# create KFold object
  kf = KFold(n_splits=k, shuffle=True, random_state=42)

# create an empty list to store the folds
  folds = []

# loop through the folds
  for train_index, test_index in kf.split(z):
    # append the test indices to the folds list
    folds.append(test_index)

# print the folds
  print(folds)

  mydata = pd.concat([fulldata, folds], axis=1)    
    # scramble response vector
  permutedY = random.sample(z)
    
  unique_levels = len(pd.unique(z))
    
    ##############################################
    
  for fold in range(1,k) :
      #Create different subdata by different fold.
      currentFold = np.where(folds == fold)
      
      # copy original data
      newData = myData
      variables = list()
  ################################################################################################################################### 
      for entry in range(1,len(gridValues)) :

        for n_run in range(1,2) :
          #take out the random row in newData and put into train, 
          #and the test data only contain the one which the train data does not have
          ########################################################################################################################
          train = newData[-currentFold, ]
          test = newData[currentFold, ]
          ########################################################################################################################
          #remove the kfold column which is the last column out of the data to calculate
          train = train[:, :-1]
          test = test[:, :-1]
          ########################################################################################################################
          #Since X data will be constant
          X_train = train.iloc[:, 1:].values

          # not used until later
          X_test = test.iloc[:, 1:].values
          
          y = train['Group']
          ######################################################################################################################
          # run cv lasso to find best lambda
        
          lasso_run = LassoCV(cv=k, normalize=True, fit_intercept=False)
          lasso_run.fit(X_train, y)

          alpha_value = lasso_run.alpha_

          lasso_run_lmin = Lasso(alpha=alpha_value*gridValues[entry]['alpha'], fit_intercept=False, normalize=False)
          lasso_run_lmin.fit(X_train, y)

          c = lasso_run_lmin.coef_
          inds = np.where(c != 0)[0]
          #######################################################################################################################
           

            #remove intercept from the list
          variables[entry] = list(c.index[inds][1:])

          selTrain = train.loc[:, train.columns.isin(variables[entry])]
          newTrain = pd.concat([train['Group'], selTrain], axis=1)
            
          selTest = test.loc[:, test.columns.isin(variables[entry])]
          newTest = pd.concat([test['Group'], selTest], axis=1)
            
          newTrain = newTrain.rename(columns={newTrain.columns[0]: "Y"})
          newTest = newTest.rename(columns={newTest.columns[0]: "Y"})
          
          # for true Y and permuteY, put in separate lists but both named Y
          if (n_run == 1) :
              gridValues[entry]['trueY'] = gridValues[entry]['trueY'] + newTest['Y'].tolist()
              
          
          else : 
              gridValues[entry]['permuteY'].append(newTest['Y'])
              
            
              newTrain = newTrain.droplevels()
              newTest = newTest.droplevels()

          svmfit = svm.SVC(kernel='linear', C=10, probability=True)
          svmfit.fit(np.array(newTrain.drop('Y', axis=1)), np.array(newTrain['Y']))

          # make predictions on new data
          yhat_SVM = svmfit.predict_proba(np.array(newTest.drop('Y', axis=1)))[:, 1]
          # random forest
          RFfit = RandomForestClassifier(n_estimators=100, random_state=0)
          RFfit.fit(X=newTrain.iloc[:, 1:], y=newTrain['Y'])

          # Making predictions on the test set
          yhat_RF = RFfit.predict(newTest.iloc[:, 1:])

            # store yhat in appropriate location depending on mode
          if (n_run == 1) {
              gridValues[entry]['chosenFeats'] = variables[entry]
              
              #svm and rf
              gridValues[entry]['append_svm'].extend(yhat_SVM)
              gridValues[entry]['svm_models'] = {'fit': svmfit, 'variables': variables[entry]}

              # Appending the predictions of the random forest and the fitted model to the corresponding lists in the dictionary
              gridValues[entry]['append_RF'].extend(yhat_RF)
              gridValues[entry]['rf_models'] = {'fit': RFfit, 'variables': variables[entry]}

          else :
             # Storing the selected variables in the dictionary
             gridValues[entry]['chosenFeat_perm'] = variables[entry]

             # Appending the predictions of the SVM and the random forest to the corresponding lists in the dictionary
             gridValues[entry]['append_svm_permute'].extend(yhat_SVM)
             gridValues[entry]['append_RF_permute'].extend(yhat_RF)  

          for (entry in range(1,len(gridValues)):
      
# Calculating the AUC for the SVM and the random forest
              gridValues[entry]['svm'] = calculateAUC(gridValues[entry]['trueY'],
                                         permutedY,
                                         gridValues[entry]['append_svm'],
                                         gridValues[entry]['append_svm_permute'],
                                         unique_levels)

              gridValues[entry]['rf'] = calculateAUC(gridValues[entry]['trueY'],
                                        permutedY,
                                        gridValues[entry]['append_RF'],
                                        gridValues[entry]['append_RF_permute'],
                                        unique_levels)

    

    gridResults[[runCount]] = gridValues
  

  retenv['gridResults'] = gridResults

  return(retenv)
   
