#######################################################################
#####################Preparation before model building ################
#######################################################################

######load required packages######

library(JMbayes)
library(forcats)
library(magrittr)
library(tidyverse)

######load the data set######

Y <- read_csv("1-data/simdata.csv")


######Curate required variables######

D <- Y %>% 
  mutate(
    brainloss = recode_factor(brainloss, `0` = "No", `1` = "Yes"),
    dementia = recode_factor(dementia, `0` = "NoDementia", `1` = "Dementia"),
    male = recode_factor(male,  `0` = "Female", `1` = "Male"),
    time = years / 20,
    id = as.factor(id)
  )

######Keep observations where outcome variable globz are not missing######

D <- D %>% filter(!is.na(globz))


########################################################################
########################## model building ##############################
########################################################################

########################separate models#################################


########Step 1: Obtain initial LDA estimates: Fit the standard#########
########longitudinal mixed model separately and obtain initial#########
#######estimates for joint model########################################              

lmeFit <-
  lme(globz ~ brainloss * time + age0 + male ,
      random = ~ time | id,
      data = D,
      na.action = na.exclude)
summary(lmeFit)


######Step 2: Obtain initial EVENT estimates: Fit the Cox proportional#### 
#####hazards survival model to obtain initial estimates for joint model####

###transform dataset appropriately for running the survival model###
Dsurv <- D %>% group_by(id) %>% filter(row_number() == 1) %>% ungroup()

###The Cox proportional hazards survival model for dementia###

wfit<- coxph(Surv(demyears,as.numeric(dementia))~brainloss+age0+male,data = Dsurv,x=TRUE,na.action=na.exclude)

summary(wfit)


#############################joint model###################################

######Step 3: Obtain final SPM estimates by fitting the joint SPM ######

wfitJOINTBayes <- jointModelBayes(lmeFit, wfit,timeVar = "time", 
                                  param = "shared-RE") 

summary(wfitJOINTBayes)
