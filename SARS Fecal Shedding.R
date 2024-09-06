cohort1 <- read.table('Data/092722_allvirusdata_PSO.csv', header = T, sep = ',')
IDs <- sort(unique(cohort1$ID))

##what fraction of total individuals who met threshold did not have a single positive measurement
threshold_ids <- c()
threshold_neg_ids <- c()
for (id in IDs){
  IDsub <- subset(cohort1, cohort1$ID == id)
  if(nrow(IDsub) > 2 && (max(IDsub$Day - min(IDsub$Day))) > 14){
    threshold_ids <- append(threshold_ids, id)
    if(sum(IDsub$N_det, na.rm = T) == 0){
      threshold_neg_ids <- append(threshold_neg_ids, id)
    }
  }
}

threshold_pos_ids <- threshold_ids[! threshold_ids %in% threshold_neg_ids]

posprob <- 1- (length(threshold_neg_ids)/length(threshold_ids))
#_____________________________________________________

c1split <- cohort1 %>%
  split(~.$ID)

sarsthreshid <- c()
for (i in names(c1split)){
  idsub <- c1split[[as.character(i)]]
  
  if (nrow(idsub) >= 3 & max(idsub$Day) - min(idsub$Day) >= 10){
    sarsthreshid <- append(sarsthreshid,i)
  }
}
sarsthreshid

sarsthresh <- cohort1[which(cohort1$ID %in% sarsthreshid & is.na(cohort1$N_conc) == F),c('ID','Day','N_conc','N_det')]

sarsthresh.norm <- sarsthresh
##transform and normalize concentration data and day
sarsthresh.norm$N_conc <- log10(sarsthresh.norm$N_conc)

sarsthresh.norm$Day <- (sarsthresh.norm$Day-min(sarsthresh.norm$Day, na.rm = T))/
  max((sarsthresh.norm$Day-min(sarsthresh.norm$Day, na.rm = T)))

ggplot(sarsthresh.norm)+
  geom_point(aes(x = Day, y = N_conc, color = as.character(ID)))

cnormsplit <- sarsthresh.norm %>%
  split(~.$ID)

shedders <- c()
noshed <- c()
maxvec <- c()
for (i in names(cnormsplit)){
  
  maxvec <- append(maxvec, max(cnormsplit[[as.character(i)]]$N_conc))
  cnormsplit[[as.character(i)]]$N_conc <- cnormsplit[[as.character(i)]]$N_conc/max(cnormsplit[[as.character(i)]]$N_conc)
  cnormsplit[[as.character(i)]]$N_conc[cnormsplit[[as.character(i)]]$N_det == F] <- 0
  
  if(sum(cnormsplit[[as.character(i)]]$N_det) > 0 ){
    shedders <- append(shedders, i)
  } else {noshed <- append(noshed,i)}
}
names(maxvec) <- names(cnormsplit)

cnormsplit <- cnormsplit[names(cnormsplit) %in% shedders]

beta_cost_fun<- function(par, data){
  shape1 <- par[1]
  shape2 <- par[2]
  xscale <- par[3]
  yscale <- par[4]
  offset <- par[5]
  
  day <- append(0,data[,1]) ## Use 0,0 as data, how can I constrain this?
  conc <- append(0,data[,2])
  
  sim <- (dbeta(((day*xscale)+offset), shape1, shape2)/
            max(dbeta((day*xscale)+offset, shape1, shape2)))*abs(yscale)
  
  cost <- sqrt((sum((-sim+conc)^2))/length(conc))
  
}

ModIDs <- list('Beta' = c(4534,4547,4577,4579,4602,4520,4525,4532,4533,4538,4540,4541,4544,4545,4546,4549,4558,4576,4585,4586,4587,4591,4593,4596,4597,4598,4600))

beta_sims <- data.frame('ID' = as.character(ModIDs$Beta), 'Shape1' = NA, 'Shape2' = NA,'xScale' = NA, 'yScale' = NA, 'Offset' = NA,  'RMSE' = NA)
simdata <- list()

InitIDS <- list('A' = c('4525','4532','4533','4544','4545','4549','4584','4602'),
                'B' = c('4540','4541','4558','4597','4598','4600'),
                'C' = c('4544','4577','4585','4596','4593'),
                'D' = c('4591'),
                'E' = c('4546','4576','4587'),
                'F' = c('4520','4586'),
                'G' = c('4538'),
                'H' = c('4534'),
                'I' = c('4579'),
                'J' = c('4547'),
                'K' = c('4526'))

##Beta models
for (i in names(cnormsplit)){
  data <- cnormsplit[[as.character(i)]][,c(2,3)]
  
  if(i %in% names(cnormsplit)){
    for(j in names(InitIDS)){
      if (i %in% InitIDS[[j]]){ group <- j}
    }
    
    Inits <- list('A' = c(3,5,1,1,0.2-data$Day[data$N_conc == max(data$N_conc)]),
                  'B' = c(2,5,1,1,0.2-data$Day[data$N_conc == max(data$N_conc)]),
                  'C' = c(4,6,1,1,0.2-data$Day[data$N_conc == max(data$N_conc)]),
                  'D' = c(10,10,1,1,0.2-data$Day[data$N_conc == max(data$N_conc)]),
                  'E' = c(10,5,1,1,0.2-data$Day[data$N_conc == max(data$N_conc)]),
                  'F' = c(3,5,0.6,1,-0.1),
                  'G' = c(20,2,0.1,1,0),
                  'H' = c(7,2,2,1,-0.5),
                  'I' = c(2,5,0.3,1,0),
                  'J' = c(8,3,2,1,-0.4),
                  'K' = c(1,1,0.8,1,0.2-data$Day[data$N_conc == max(data$N_conc)]))
    
    par <- Inits[[group]]
    
    try(out <- optim(par,fn = beta_cost_fun, data = data,control=list(maxit=1E5)))
    
    beta_sims[beta_sims$ID == i,'Shape1'] <- out$par[1]
    beta_sims[beta_sims$ID == i,'Shape2'] <- out$par[2]
    beta_sims[beta_sims$ID == i,'xScale'] <- out$par[3]
    beta_sims[beta_sims$ID == i,'yScale'] <- out$par[4]
    beta_sims[beta_sims$ID == i, 'Offset'] <- out$par[5]
    beta_sims[beta_sims$ID == i,'RMSE'] <- out$value
    simdata[[as.character(i)]] <- out
  }
}
#renormalize data
predicts <- data.frame()
norm.predicts <- data.frame()

for (i in 1:nrow(beta_sims)){
  x.norm <- seq(0,5,length = 151)
  x <- ((x.norm)*30)-2
  
  c.norm <- dbeta((x.norm*beta_sims[i,'xScale'])+beta_sims[i,'Offset'], shape1 = beta_sims[i,'Shape1'], shape2 = beta_sims[i,'Shape2'])/
    max(dbeta((x.norm*beta_sims[i,'xScale'])+beta_sims[i,'Offset'], shape1 = beta_sims[i,'Shape1'], shape2 = beta_sims[i,'Shape2'])) *abs(beta_sims[i,'yScale'])
  
  c <- 10^((c.norm)*maxvec[beta_sims$ID[i]])
  
  predicts <- bind_rows(predicts, data.frame('ID' = beta_sims$ID[i],'Day'=x,'N_conc' = c))
  norm.predicts <- bind_rows(norm.predicts, data.frame('ID' = beta_sims$ID[i],'Day'=x.norm,'N_conc' = c.norm))
}

cnorm <- bind_rows(cnormsplit)

predicts <- predicts %>%
  split(~.$ID)
