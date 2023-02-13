processFile = function(filepath, df, fullModules = F) {
  results = list()
  con = file(filepath, "r")
  gcount = 0
  allModules = list()
  while ( TRUE ) {
    line = readLines(con, n = 1)
    if ( length(line) == 0 ) {
      break
    }
    
    genes = str_split(line, "\t|\n")[[1]]
    allModules[length(allModules) + 1] = list(l = genes)
    gcount = gcount + 1
    genelist = c()
    for (g in genes) {
      if (g %in% names(df[, -1])) {
        genelist = c(genelist, g)
      }
    }
    if ( length(genelist) >= 3) {
      if (fullModules) {
        # get all genes in the modules
        results[length(results) + 1] = list(l = genes)
      } else {
        # get only genes in data from module
        results[length(results) + 1] = list(l = genelist)
      }
    }
  }
  
  close(con)
  final_results = list(results = results, allModules = allModules)
  return(final_results)
}

getModules = function(dfs, cluster_filepath) {

  modules = processFile(cluster_filepath, dfs)
  
  modulesInData = modules$results
  tempModules <- modulesInData
  data = dfs
  
  M <- matrix(0, ncol=nrow(data), nrow = length(modulesInData))
  i <- 1
  Y = data$Group
  moduleNames = list()
  for (i in 1:length(tempModules)) {
    
    module <- modulesInData[[i]]
    # if(unique(module[module != ""]))
    module <- unique(module[module != ""])
    
    # module <- module[2:length(module)]
    subDF <- data %>% dplyr::select(all_of(module))
    # moduleNames = c(moduleNames, module)
    moduleNames[[length(moduleNames) + 1]] = module
    # # removing features with 0's
    # limit <- 0.1 * nrow(subDF)
    # subDF <- subDF[, which(apply(subDF, 2, function(col) !any(table(col) > limit)))]
    subDF <- cbind.data.frame(Y, subDF)
    
    # Zr = base::colMeans(subDF[,-1])
    
    grp_means = matrix(0, nrow = length(unique(Y)), ncol = ncol(subDF))
    
    for (val in 1:length(unique(Y))) {
      grp_means[val,] = colSums(subset(subDF, Y == unique(Y)[val]))
    }
    
    grp_means = grp_means[, -1]
    # Zr <- subset(subDF, Y==1)
    # Znr <- subset(subDF, Y==0)
    # muZr <- colMeans(Zr)
    # muZnr <- colMeans(Znr)
    # muZr <- muZr[c(2:length(muZr))]
    # muZnr <- muZnr[c(2:length(muZnr))]
    
    tempSubDF <- subDF[,-1]
    m <- c()
    N <- c()
    SN <- c()
    tempSign <- 1
    for (k in 1:nrow(as.matrix(tempSubDF))) {
      for (j in 1:ncol(as.matrix(tempSubDF))) {
        
        
        # group means for col j
        j_means = sum(grp_means[, j]^2)
        
        # calculate sample vs group scores for each group
        # std 
        # get value of Y for sample k
        y_k = Y[k]
        j_score = (tempSubDF[k, j] - grp_means[, j]) / sum(1 + grp_means[, j]^2)

        # get abs values

        m[j] = sum(j_score)
        
        
        # m[j] = log2( (tempSubDF[k, j]^2) / (Zr[j]^2) )
        
        # m[j] <- 2^( (tempSubDF[k, j]^2)/( (muZr[j]^2 + muZnr[j]^2) ) )
        
        SN[j] <- sign(m[j])
        
        
      }
      
      m[!is.finite(m)] <- 0.0001    
      
      posSN <- 0
      negSN <- 0
      posSN <- sum(SN > 0, na.rm = T)
      negSN <- sum(SN < 0, na.rm = T)
      
      m  <- abs(m)
      
      if (posSN > negSN){
        N[k] <- sum(m)*1
      } else{
        N[k] <- sum(m)*(-1)
      }
    }
    
    M[i, ] <- N
  }
  
  modulePA <- t(M)  #%>% data.frame(stringsAsFactors = F)
  
  # module names correspond to the genes in each module
  # names(modulePA) = moduleNames
  names(modulePA) = 1:length(moduleNames)
  # modulePA = modulePA[, !duplicated(names(modulePA))]
  Group = data$Group
  modulePA <- cbind.data.frame(Group, modulePA)
  
  
  retenv = new.env()
  retenv$modulePA = modulePA
  retenv$modules = modules
  return(retenv)

}
