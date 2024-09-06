library(tidyverse)
library(zoo)
library(patchwork)
library(scales)
library(deSolve)
library(fitdistrplus)

source('UMcolors.R')
source('SARS Fecal Shedding.R')
source('Biomarker Fecal Shedding.R')
source('Outbreak Model.R')
source('WW Modeling.R')
source('Shedding Functions.R')


N <- 10 #Population size

parms <- c(beta=0.2,gamma.1=0.1,gamma.e=0.2) #SEIR parameters

p = 1/N
init <- c(S = 1-p, E = 0, I = p, R = 0) #SEIR initial conditions

times <- 1:550 #How long to run the sim for, may take some trial and error

#Run SEIR Model and discretize
mod <- SEIR.model(init = init, beta.s = parms[1], gamma.e = parms[3], gamma.i = parms[2], times = times)
SEIR_out <- round(mod*N, digits = 0)

#Run shedding model
SARS_count <- create_infectioncounter(n = N)
for (day in 1:nrow(SEIR_out)){
  
  SARS_count <- day_infection_counts(I = SEIR_out[day,3],
                                     R = SEIR_out[day,4],
                                     shed_fits =  predicts, 
                                     n = N, 
                                     d = day, 
                                     infectioncounter = SARS_count,
                                     posprob = posprob)
  
}

SARS_tot_ww <- c()
for (jour in 5:ncol(SARS_count)){
  
  totdil <- 0
  
  for (ind in 1:nrow(SARS_count)){
    totdil <- totdil+ total_dil(Cf = SARS_count[ind,jour],Ptotal = N)
 }
  
  SARS_tot_ww <- append(SARS_tot_ww, totdil)
  }

ww_df <- data.frame('Day' = 1:length(SARS_tot_ww),'WW' = SARS_tot_ww)

ggplot(ww_df)+
  geom_line(data = SEIR_out, aes(x = 1:length(I),y = (I/max(I))*1.2*max(ww_df$WW)), linetype = 'dashed', color = 'red')+
  geom_line(aes(x = Day, y = WW))+
  guides(color = 'none')+
  labs(x = 'Day', y = 'gc/mL-total WW')+
  theme_bw()




