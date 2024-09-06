c1split <- cohort1 %>%
  split(~.$ID)

sample_size <- map_dbl(c1split,nrow)

idthresh <- names(sample_size)[sample_size > 2]

dayrange <- c()
for (i in unique(cohort1$ID)){
  idsub <- cohort1[which(cohort1$ID == i),]
  dayrange <- append(dayrange ,max(idsub$Day, na.rm = T)
                     - min(idsub$Day, na.rm = T))
}

names(dayrange) <- unique(cohort1$ID)

normtest <- data.frame()


CohThresh <- cohort1[cohort1$ID %in% names(sample_size >3),]

Cohnormp <- tibble('ID' = CohThresh$ID[which(CohThresh$PMMoV_det)], 
                   'Day' = CohThresh$Day[which(CohThresh$PMMoV_det)],
                   'PMMoV' = CohThresh$PMMoV_conc[which(CohThresh$PMMoV_det)])

Cohnormc <- tibble('ID' = CohThresh$ID[which(CohThresh$crAss_det)], 
                   'Day' = CohThresh$Day[which(CohThresh$crAss_det)],
                   'crAss' = CohThresh$crAss_conc[which(CohThresh$crAss_det)])


CohNorm <- full_join(Cohnormp, Cohnormc)
CohNorm$PMMoV <- log10(CohNorm$PMMoV)/log10(max(CohNorm$PMMoV, na.rm = T))
CohNorm$crAss <- log10(CohNorm$crAss)/log10(max(CohNorm$crAss, na.rm = T))

Cohsplit <- CohNorm %>%
  split(~.$ID)

pbeta <- data.frame()
cbeta <- data.frame()

pfit <- list()
cfit <- list()

pplot <- data.frame()
cplot <- data.frame()

CulFreyplotp <- ggplot()
CulFreyplotc <- ggplot()
for (i in idthresh){
  psub <- Cohsplit[[as.character(i)]]$PMMoV[which(is.na(Cohsplit[[i]]$PMMoV) == F)]
  csub <- Cohsplit[[as.character(i)]]$crAss[which(is.na(Cohsplit[[i]]$crAss) == F)]
  
  if (length(psub) >= 3){
    try(pfit <- fitdist(psub, distr = 'beta', method = 'mme'))
    pbeta <- bind_rows(pbeta, c('ID' = as.character(i),pfit$estimate[1], pfit$estimate[2]))
    dens <- dbeta(seq(0,1,length = 100), shape1 = pfit$estimate[1], shape2 = pfit$estimate[2])
    pplot <- bind_rows(pplot, data.frame('ID' = i, 
                                         'conc' = seq(0,1,length = 100),
                                         'density' = dens/max(dens)))
  }
  
  if (length(csub) >= 3){
    try(cfit <- fitdist(csub, distr = 'beta', method = 'mme'))
    cbeta <- bind_rows(cbeta, c('ID' = as.character(i), cfit$estimate[1], cfit$estimate[2])) 
    dens <- dbeta(seq(0,1,length = 100), shape1 = cfit$estimate[1], shape2 = cfit$estimate[2])
    cplot <- bind_rows(cplot, data.frame( 'ID' = i, 
                                          'conc' = seq(0,1,length = 100),
                                          'density' = dens/max(dens)))
  }
}

cplot$conc <- 10^(cplot$conc*log10(max(Cohnormc$crAss, na.rm = T)))

pplot$conc <- 10^(pplot$conc*log10(max(Cohnormp$PMMoV, na.rm = T)))

pplot <- pplot[which(is.infinite(pplot$density) == F),]
cplot <- cplot[which(is.infinite(cplot$density) == F),]


CohThresh <- cohort1[which(cohort1$ID %in% names(sample_size >3)),]

Cohnormp <- tibble('ID' = CohThresh$ID[CohThresh$PMMoV_det & is.na(CohThresh$PMMoV_det) == F], 
                   'Day' = CohThresh$Day[CohThresh$PMMoV_det & is.na(CohThresh$PMMoV_det) == F],
                   'PMMoV' = CohThresh$PMMoV_conc[CohThresh$PMMoV_det & is.na(CohThresh$PMMoV_det) == F])

Cohnormc <- tibble('ID' = CohThresh$ID[CohThresh$crAss_det & is.na(CohThresh$crAss_det) == F], 
                   'Day' = CohThresh$Day[CohThresh$crAss_det & is.na(CohThresh$crAss_det) == F],
                   'crAss' = CohThresh$crAss_conc[CohThresh$crAss_det & is.na(CohThresh$crAss_det) == F])


CohNorm <- full_join(Cohnormp, Cohnormc)
CohNorm$PMMoV <- log10(CohNorm$PMMoV)/log10(max(CohNorm$PMMoV, na.rm = T))
CohNorm$crAss <- log10(CohNorm$crAss)/log10(max(CohNorm$crAss, na.rm = T))

Cohsplit <- CohNorm %>%
  split(~.$ID)

pbeta <- data.frame()
cbeta <- data.frame()

pfit <- list()
cfit <- list()

pplot <- data.frame()
cplot <- data.frame()

CulFreyplotp <- ggplot()
CulFreyplotc <- ggplot()
for (i in unique(CohNorm$ID)){
  psub <- Cohsplit[[as.character(i)]]$PMMoV[is.na(Cohsplit[[as.character(i)]]$PMMoV) == F]
  csub <- Cohsplit[[as.character(i)]]$crAss[is.na(Cohsplit[[as.character(i)]]$crAss) == F]
  
  if (length(psub) >= 3){
    try(pfit <- fitdist(psub, distr = 'beta', method = 'mme'))
    pbeta <- bind_rows(pbeta, c('ID' = i,pfit$estimate[1], pfit$estimate[2]))
    dens <- dbeta(seq(0,1,length = 100), shape1 = pfit$estimate[1], shape2 = pfit$estimate[2])
    pplot <- bind_rows(pplot, data.frame('ID' = i, 
                                         'conc' = seq(0,1,length = 100),
                                         'density' = dens/max(dens)))
  }
  
  if (length(csub) >= 3){
    try(cfit <- fitdist(csub, distr = 'beta', method = 'mme'))
    cbeta <- bind_rows(cbeta, c('ID' = i, cfit$estimate[1], cfit$estimate[2])) 
    dens <- dbeta(seq(0,1,length = 100), shape1 = cfit$estimate[1], shape2 = cfit$estimate[2])
    cplot <- bind_rows(cplot, data.frame( 'ID' = i, 
                                          'conc' = seq(0,1,length = 100),
                                          'density' = dens/max(dens)))
  }
}

cplot$conc <- 10^(cplot$conc*log10(max(Cohnormc$crAss, na.rm = T)))

pplot$conc <- 10^(pplot$conc*log10(max(Cohnormp$PMMoV, na.rm = T)))

pplot <- pplot[which(is.infinite(pplot$density) == F),]
cplot <- cplot[which(is.infinite(cplot$density) == F),]



ggplot()+
  geom_histogram(data = Cohnormc,aes(x = crAss, y = after_stat(ndensity), color = 'crAss', fill = 'crAss'), alpha = 0.4)+
  geom_histogram(data = Cohnormp,aes(x = PMMoV, y = after_stat(ndensity), color = 'PMMoV', fill = 'PMMoV'), alpha = 0.5)+
  geom_line(data = cplot, aes(x = conc, y = density, color = 'crAss'))+
  geom_line(data = pplot, aes(x = conc, y = density, color = 'PMMoV'))+
  facet_wrap(vars(ID))+
  scale_x_log10(labels = trans_format('log10',math_format(10^.x)), breaks = 10^seq(2,8, by =2), minor_breaks = 10^c(1:9))+
  scale_color_manual(values = c('crAss' = '#00274C', 'PMMoV' = '#FFCB05'))+
  scale_fill_manual(values = c('crAss' = '#00274C', 'PMMoV' = '#FFCB05'))+
  labs(y = 'Normalized Density', x = 'Conc (copies/mL-WW)', fill = '')+
  guides(color = 'none')+
  theme_bw()


c1split <- cohort1 %>%
  split(~.$ID)

length(c1split[])/length(c1split)
Pfrac <- c()
Cfrac <- c()
for (i in names(c1split)){
  
  idsub <- c1split[[as.character(i)]]
  if(nrow(idsub) >=3){
    pp <- sum(idsub$PMMoV_det, na.rm = T) > 0
    cc <- sum(idsub$crAss_det,na.rm = T) > 0
    
    Pfrac <- append(Pfrac, pp)
    Cfrac <- append(Cfrac, cc)
  }
}

Pfrac <- sum(Pfrac)/length(Pfrac)
Cfrac <- sum(Cfrac)/length(Cfrac)



betabiomod <- function(pars, N, D, scale, fracshed){
  counter <- data.frame('Individual' = 1:N, 'ID' = 1:N)
  
  shedpop <- round(N*fracshed, digits = 0)
  
  counter$ID <- c(sample(1:(nrow(pars)), size = shedpop, replace = T), seq(0,0, length = N-shedpop))
  
  shed_vect <- c()
  for (day in 1:D){
    
    shed_vect <- rbeta(shedpop,
                       shape1 = as.numeric(pars[counter$ID[1:shedpop],2]), 
                       shape2 = as.numeric(pars[counter$ID[1:shedpop],3]))
    
    shed_vect <- c(10^(shed_vect*log10(scale)), seq(0,0, length = N-shedpop))
    
    counter <- cbind(counter,data.frame(shed_vect))
  }
  counter
}