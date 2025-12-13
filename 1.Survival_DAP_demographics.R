library(lubridate)
library(survival)
library(flexsurv)
library(survminer)
library(dplyr)
datfolder <- "data"
resultsfolder <- "results"
figfolder <- "figures"
load(file.path(datfolder,"SurvivalData.RData"))
surv <- Surv(time=survivalData$first.age,
             time2=survivalData$last.age,
             event=survivalData$event,
             type='counting')
