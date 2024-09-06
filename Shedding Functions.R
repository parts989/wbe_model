create_infectioncounter <- function(n){
  infectioncounter <- data.frame(1:n, vector(mode = 'numeric', length = n),vector(mode = 'numeric', length = n),vector(mode = 'numeric', length = n))
  colnames(infectioncounter) <- c('ID', 'Day_of_infection', 'Trajectory', 'traj_counter')
  infectioncounter
}


day_infection_counts <- function(I, R, shed_fits, n, d, infectioncounter, posprob){
  
  oldinfects <- nrow(subset(infectioncounter, infectioncounter[,2] != 0))
  totalinfects <- I+R
  
  #if(max(oldinfects +1) <= n){
  #newinfects <- seq(min(oldinfects+1,n),totalinfects)
  #} else {newinfects <- c()}
  if(oldinfects < totalinfects){
    newinfects <- totalinfects - oldinfects
    for (i in (oldinfects+1):(oldinfects+newinfects)){
      infectioncounter[i,3] <- rbinom(1,1,posprob)
      if (infectioncounter[i,3] == 1){
        infectioncounter[i,3] <- sample(1:length(names(shed_fits)), 1)
      }
      infectioncounter[i,2] <- d
    }
  }
  
  dayshed <- c()
  
  for (j in 1:totalinfects){
    infectioncounter[j,4] <- infectioncounter[j,4] + 1
    
    if (infectioncounter[j,3] > 0){
      shed <- shed_fits[[infectioncounter[j,3]]][infectioncounter[j,4],3] #indexing the shedding value for each day
    } else if (infectioncounter[j,3] == 0) {
      shed <- 0
    }
    
    dayshed <- append(dayshed, shed)
  }
  zeros <- seq(0,0,length.out = n-length(dayshed))
  dayshed <- append(dayshed,zeros)
  
  
  dayshed[which(is.na(dayshed) == T)] <- 0
  
  infectioncounter <- infectioncounter %>%
    mutate(dayshed)
  
  names(infectioncounter)[names(infectioncounter) == 'dayshed'] <- d
  
  infectioncounter
}

#build function that divides data into 10 quantiles

day_quant <- function(vec, n = 10, Pop, Mf = 28000, Uw = 1E5, LOD = 0){
  
  c <- data.frame('a' = ntile(-vec[which(vec > LOD)], n = n), 'b' = vec[which(vec > LOD)])
  
  c <- c %>%
    split(~.$a)
  
  quant_split <- data.frame('Quant' =c(), 'ww_cont' = c())
  for(i in names(c)){
    len <- nrow(c[[i]])
    quant_split <- bind_rows(quant_split,data.frame('Quant' = i, 'ww_cont' = sum(total_dil(Cf = c[[i]]$b, Pshed = 1, Ptotal = Pop))))
  }
  
  quant_split
}