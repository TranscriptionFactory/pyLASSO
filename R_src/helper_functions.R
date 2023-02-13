#######################END###############################

#There are 5 functions to process the QC
#      update_GeneName()
#      QC_process()
#      transpose_data()
#      rm_stableVar()
#      znorm_data()
##########This part is filteriing process for the X_axis data################


#Updating GeneName function
update_GeneName <- function(data){
  #update the GeneName
  data[c('A','B','C')] <- str_split_fixed(toupper(data$GeneName), ";", 3)
  for (i in 1: nrow(data)){
    #if (data$A[i] == ""){
    if (is.na(data$A[i])){
      data$GeneName[i] <- paste0("NoName")
    }
    #else if (str_detect(data$A[i], 'HCG')){
    else if (str_count(data$A[i], "HCG") > 0){
      data$GeneName[i] <- data$B[i]
    }
    else {
      data$GeneName[i] <- data$A[i]
    }
  }
  new_data <- data %>% filter(GeneName!="NoName") %>% dplyr::select(-c(A,B,C))
  return(new_data)
}



#filtering function
QC_process <- function(Xdata) {
  
  Xdata$Unique_peptide <- as.numeric(Xdata$Unique_peptide)
  
  # filter potential contaminants
  Xdata = Xdata %>% filter(is.na(Potential_contaminant) == TRUE) %>%
    select(-Potential_contaminant)
  
  # not filtering potential contaminants, just low peptides
  Xdata <- Xdata %>% filter(Unique_peptide >= 1) %>%
    dplyr::select(-Unique_peptide)
  

  Xdata[,2:ncol(Xdata)] <- apply(Xdata[,2:ncol(Xdata)], 2, as.numeric)
  
  impute.mean <- function(x){
    bad_cells = which( (is.na(x) | is.nan(x) | is.infinite(x) | !is.numeric(x)))

    return(replace(x, bad_cells, 0))
  }

  safe_scaling = function(x) {
    if(sum(x) != 0){
      return(scale(x, T, T))
    }
    else {
      return(x)
    }
  }

  # Xdata_scaled = Xdata
  # center data
  # Xdata_scaled[, 2:ncol(Xdata)] = lapply(Xdata[, 2:ncol(Xdata)], safe_scaling)

  # Xdata_imputed = impute.mean(Xdata_scaled)
  # Xdata = Xdata_scaled
  # Xdata[,2:ncol(Xdata)] <- apply(Xdata[,2:ncol(Xdata)], 2, impute.mean)
  
  # NOT removing proteins with missing data
  #Remove protein which has more than 50% missing data in all samples:
  # Xdata$Misscount <- rowSums(Xdata==0, na.rm=FALSE)*100/length(Xdata)
  # Xdata <- Xdata %>% filter(Misscount <= 75) %>% dplyr::select(!Misscount)
  return(Xdata)
}
# 
# #transpose the dataframe => column is GeneName and row is sampleID, the data is the intensity data. 
# transpose_data <- function(Xdata) {
#   gene_names = Xdata$GeneName
#   t_Xdata <- as.data.frame(as.numeric(t(Xdata[,-1])),
#                            col.names = Xdata$GeneName)
#   # names(t_Xdata) <- gene_names
#   # t_Xdata$genes = factor(row.names(t_Xdata))
#   return(t_Xdata)
# }
# 
# 
# 
# ##Unsupervised, filter all the features cross data (1 feature cross all sample, if there is no variants, remove), 
# # take standard variation per column (per protein/gene across sample), 
# #Filter the feature if the feature values are closed to the mean. This is called variance base filter.
# #This is easy to understand because if there is no feature variant cross sample, it mean that feature has to affect
# #on the sample, we should filter it because the feature give an equal affect to the sample of disease. It also mean
# #this protein or gene does not affect!!!. Drop anything below 10%, plot boxplot and remove the bottom quartile. 
# #Variance = sqrt(standard deviation)
# #Low Variance Filter method
# #t_dataX, transpose form of data and the dataX is the regular form of the data.
# #this is function remove stable variant cross all sampleID
# 
# rm_stableVar <- function(t_dataX, dataX) {
#   #create the empty var to whole value 
#   temp_var <- data.frame(matrix(ncol = length(t_dataX), nrow = 1))
#   colnames(temp_var) <-names(t_dataX)
#   
#   for (i in 1:length(t_dataX)) {
#     temp_var[,i] <- var(t_dataX[,i])
#   }
#   
#   #vector trans_temp is holding the variance value for the t_Xdata
#   trans_temp <- t(temp_var)
#   colnames(trans_temp) <- "variance"
#   add_dataX <- cbind(dataX,trans_temp)
#   #create the quantile variable measurement, then remove by filtering the variance value 
#   #which less than .25
#   
#   #NOT removing based on variance
#   # res1<-quantile(trans_temp, .25)
#   add2_dataX <- add_dataX %>% dplyr::select(!variance)
#   new_dataX <- as.data.frame(t(add2_dataX[,-1]))
#   colnames(new_dataX) <- add2_dataX$GeneName
#   return(new_dataX)
# }
# 
# 
# 
# #Normalzation process by Z-score by function znorm_data
# znorm_data <- function(x){
#   newset <- data.frame(matrix(ncol = length(x), nrow = nrow(x)))
#   colnames(newset) <- names(x)
#   newset[is.na(newset)] <- 0
#   for (i in 2:length(x)) {
#     newset[,i] <- scale(x[,i])
#   }
#   return(newset)
# }