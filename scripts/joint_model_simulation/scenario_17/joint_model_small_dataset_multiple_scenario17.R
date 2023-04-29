library(dplyr)
library(ggplot2)
library(magrittr)
library(ggpubr)


set.seed(4282432)
sample_size <- 300
end_time <- 26280

for(i in 1:2500){

  #Baseline characteristics simulation
  # set.seed(i)
  eso <- data.frame(ID = numeric(5000))
  eso$ID <- 1:5000
  eso$age <- rnorm(5000, 62.6, 8.23)
  summary(eso$age)

  #weight
  # set.seed(i+1)
  eso$wt <- rnorm(5000, 55.2, 10)
  summary(eso$wt)

  #renal function
  # set.seed(i+2)
  eso$egfr <- rnorm(5000, 86.8, 17.1)
  summary(eso$egfr)

  #albumin conc.
  # set.seed(i+3)
  eso$alb <- rnorm(5000, 3.94, 0.429)
  summary(eso$alb)

  #LDH
  # set.seed(i+4)
  u <- log(192)
  sigma <- 0.75
  eso$ldh <- rnorm(5000, u, sigma)
  eso$ldh <- exp(eso$ldh)
  summary(eso$ldh)

  #LDH uper limit of normal
  eso$uldh <- eso$ldh/250

  #sex
  # set.seed(i+5)
  eso$sex <- runif(5000,0,1)
  eso$sex <- ifelse(eso$sex > 0.15, 0, 1 )
  summary(eso$sex)

  #race
  # set.seed(i+6)
  eso$race <- runif(5000,0,1)
  eso$race <- ifelse(eso$sex > 0.05, 1, 2 )
  summary(eso$race)

  #performance status ECOG
  # set.seed(i+7)
  eso$ps <- runif(5000,0,1)
  eso$ps <- ifelse(eso$ps > 0.5, 0, 1 )
  summary(eso$ps)

  #tumor type
  eso$tt <- 1

  #Tumor size
  # set.seed(i+8)
  u <- log(4.027)
  sigma <- 0.6
  eso$TUS <- rnorm(5000, u, sigma)
  eso$TUS <- exp(eso$TUS)
  summary(eso$TUS)

  #pdl status
  # set.seed(i+9)
  eso$pdl <- runif(5000, 0,1)
  eso$pdl <- ifelse(eso$pdl > 0.5, 0, 1 )
  summary(eso$pdl)

  #treatment group
  # set.seed(i+10)
  eso$group <- runif(5000,0,1)
  eso$group <- ifelse(eso$group > 0.5, 1, 2 )
  summary(eso$group)

  #END OF BASELINE CHARACTERISTICS SIMULATION

  #Adding dosing records



  df_dose <- data.frame(TIME = c(0, end_time),
                        AMT = c(240,0),
                        ADDL = c(130,0),
                        II = c(336,0),
                        RATE = c(480,0),
                        DV = 0,
                        MDV = c(1,0),
                        EVID = c(1,0),
                        CMT = c(1,1),
                        FLAG = c(2,2))


  #multiplying the dataframe for 10000 subjects
  rows <- c(1:nrow(df_dose))
  times <- 5000
  df_dose<- df_dose[rep(rows, times),]

  #creating ids df to be merged with the df_bind dataset
  id <- seq(from = 1, to = 5000, by = 1)
  id <- as.data.frame(id)
  names(id)[1] <- 'ID'

  id <- id %>% slice(rep(1:n(), each =2))

  #merging the id and df_bind data set, so every subjects has a dosing record for three years

  df_full <- cbind(id, df_dose)

  #merging baseline covariates with the dosing records

  df_final <- merge(df_full, eso, by = 'ID')

  ###

  df_2drugs <- df_final

  df_2drugs$TRT <- ifelse(df_2drugs$group == 1, 1, 0)

  #### exclusion criteria wt < 34 and ldh <100
  df_2drugs_reduced <- subset(df_2drugs,
                              df_2drugs$wt >= 34)

  df_2drugs_reduced <- subset(df_2drugs_reduced,
                              df_2drugs_reduced$ldh >= 100)


  df_1drug_comp <- subset(df_2drugs_reduced, group == 1)

  df_1drug_nivo <- subset(df_2drugs_reduced, group == 2)


  selected_subjects_IDs_comp<- unique(df_1drug_comp$ID)

  selected_subjects_IDs_nivo<- unique(df_1drug_nivo$ID)


  # set.seed(i+11)
  enrolled_ids_comp<- sample(selected_subjects_IDs_comp,
                             sample_size/2, replace = FALSE,
                             prob = NULL)

  # set.seed(i+12)
  enrolled_ids_nivo<- sample(selected_subjects_IDs_nivo,
                             sample_size/2, replace = FALSE,
                             prob = NULL)


  enrolled_subjects_comp<- subset(df_2drugs_reduced, ID %in% enrolled_ids_comp)

  enrolled_subjects_nivo<- subset(df_2drugs_reduced, ID %in% enrolled_ids_nivo)

  simulation_dataset<- rbind(enrolled_subjects_nivo, enrolled_subjects_comp)

  ## writing to file
  write.table(simulation_dataset, file=paste('small_dataset',i,'.csv', sep=''),append =F, sep=",", na=".",row.names = FALSE,col.names = TRUE, quote = FALSE)
}
