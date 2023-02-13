# #######################END###############################
# 
# #There are 5 functions to process the QC
# #      update_GeneName()
# #      QC_process()
# #      transpose_data()
# #      rm_stableVar()
# #      znorm_data()
# ##########This part is filteriing process for the X_axis data################
# 
# 
# #Updating GeneName function
# update_GeneName <- function(data){
#   #update the GeneName
#   data[c('A','B','C')] <- str_split_fixed(toupper(data$GeneName), ";", 3)
#   for (i in 1: nrow(data)){
#     #if (data$A[i] == ""){
#     if (is.na(data$A[i])){
#       data$GeneName[i] <- paste0("NoName")
#     }
#     #else if (str_detect(data$A[i], 'HCG')){
#     else if (str_count(data$A[i], "HCG") > 0){
#       data$GeneName[i] <- data$B[i]
#     }
#     else {
#       data$GeneName[i] <- data$A[i]
#     }
#   }
#   new_data <- data %>% filter(GeneName!="NoName") %>% dplyr::select(-c(A,B,C))
#   return(new_data)
# }
# 
# 
# 
# #filtering function
# QC_process <- function(Xdata) {
#   
#   Xdata = Xdata %>% dplyr::select(-Potential_contaminant)
#   Xdata$Unique_peptide <- as.numeric(Xdata$Unique_peptide)
#   # Xdata <- as.numeric(Xdata)
#   
#   
#   # previous
#   # Xdata <- Xdata %>% 
#   #   #filtering the positive comtamninant and Unique peptide < 1
#   #   filter(is.na(Potential_contaminant)==TRUE) %>%
#   #   filter(is.na(Unique_peptide)==FALSE & Unique_peptide >= 1) %>%
#   #   dplyr::select(!c(Potential_contaminant,Unique_peptide))
#   # 
#   # 
#   # not filtering potential contaminants, just low peptides
#   # Xdata <- Xdata %>% filter(Unique_peptide >= 1) %>%
#   #   dplyr::select(!c(Potential_contaminant,Unique_peptide))
#   Xdata <- Xdata %>% filter(Unique_peptide >= 1) %>%
#     dplyr::select(-Unique_peptide)
#   
#   
#   Xdata[,2:ncol(Xdata)] <- apply(Xdata[,2:ncol(Xdata)], 2, as.numeric)
#   
#   # Xdata = as.data.frame(as.numeric(Xdata))
#   #cast entire df as numeric
#   # Xdata[,1:ncol(Xdata)] <- sapply(Xdata[,1:ncol(Xdata)],as.numeric)
#   # Xdata <- sapply(Xdata,as.numeric)
#   
#   
#   impute.mean <- function(x){
#     
#     bad_cells = which( (is.na(x) | is.nan(x) | is.infinite(x) | !is.numeric(x)))
#     return(replace(x, bad_cells, 0))
#   }
#   # } replace(x, is.na(x) | is.nan(x) | is.infinite(x) | !is.numeric(x), 0)
#   
#   Xdata[,2:ncol(Xdata)] <- apply(Xdata[,2:ncol(Xdata)], 2, impute.mean)
#   
#   # Xdata[,1:ncol(Xdata)] <- apply(Xdata[,1:ncol(Xdata)],as.numeric)
#   
#   
#   #Xdata$CountNA <- rowSums(is.na(Xdata), na.rm=T)
#   #Xdata <- Xdata %>% filter(CountNA == 0) %>% dplyr::select(!CountNA)
#   # print(paste0("after filtering the protein has NA record cross all the samples is ",dim(Xdata)))
#   
#   
#   # Xdata = Xdata %>% mutate(across(2:ncol(Xdata), ~replace_na(.x, 0.001)))
#   
#   # impute.mean <- function(x) replace(x, is.na(x) | is.nan(x) | is.infinite(x), base::mean(x))
#   # Xdata <- sapply(Xdata, impute.mean) %>% data.frame(stringsAsFactors = F)
#   
#   
#   #After remove NA, how has to update the data-frame as numeric
#   # Xdata[,2:ncol(Xdata)] <- sapply(Xdata[,2:ncol(Xdata)],as.numeric)
#   
#   # impute.mean <- function(x) replace(x, is.na(x) | is.nan(x) | is.infinite(x), base::mean(x[!is.na(x) & !is.nan(x) & !is.infinite(x)]))
#   # Xdata <- sapply(Xdata, impute.mean) %>% data.frame(stringsAsFactors = F)
#   # 
#   
#   #Remove protein which has more than 50% missing data in all samples:
#   # Xdata$Misscount <- rowSums(Xdata==0, na.rm=FALSE)*100/length(Xdata)
#   # Xdata <- Xdata %>% filter(Misscount <= 75) %>% dplyr::select(!Misscount)
#   return(Xdata)
# }
# 
# 
# 
# #transpose the dataframe => column is GeneName and row is sampleID, the data is the intensity data. 
# transpose_data <- function(Xdata) {
#   t_Xdata <- as.data.frame(as.numeric(t(Xdata[,-1])))
#   colnames(t_Xdata) <- t_Xdata$GeneName
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
#     #var() comput the variance of data
#     temp_var[,i] <- var(t_dataX[,i])
#   }
#   
#   #vector trans_temp is holding the variance value for the t_Xdata
#   trans_temp <- t(temp_var)
#   colnames(trans_temp) <- "variance"
#   add_dataX <- cbind(dataX,trans_temp)
#   #create the quantile variable measurement, then remove by filtering the variance value 
#   #which less than .25
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