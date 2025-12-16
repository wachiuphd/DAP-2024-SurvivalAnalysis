library(lubridate)
library(dplyr)
datfolder <- "data"
survivalfile <- "survivalData_2024_12_15a.csv" # "SurvivalAnalysisCuratedDogs_2024-12-31.csv"
# survivalData <- read.csv("survivalData_2024_12_15a.csv")
# overviewData <- read.csv(file.path(datfolder,"DAP_2024_DogOverview_v1.0.csv"))
# followupcols<-which(grepl("CSLB",names(overviewData)) |
#   grepl("AFUS",names(overviewData)) |
#   grepl("CognitiveGames",names(overviewData)) |
#   grepl("Morphometrics",names(overviewData)) |
#   grepl("DNA",names(overviewData)))
# overviewData.withfollowup <- data.frame(overviewData[,c("dog_id","Status")],
#                            has_followup = grepl("Y",apply(overviewData[,followupcols],1,paste0,collapse=""),
#                                                  ignore.case = TRUE))
# missing_ids <- setdiff(
#   subset(overviewData.withfollowup,Status=="Alive" & has_followup)$dog_id,
#   subset(survivalData,survival_status=="Alive")$dog_id
#   )
# 
## All survivalData in overviewData
## All alive in survivalData have follow-up in overviewData 
## All 


## In code, dap_pack_date vs. dd_form_date

survivalData <- read.csv(file.path(
  datfolder,survivalfile)) %>% 
  rename(dog_id=study_id) %>% left_join(
    read.csv(file.path(datfolder,"DAP_2024_DogOverview_v1.0.csv"))) 
datecols <- which((grepl("date",names(survivalData),ignore.case = TRUE) |
                     grepl("dob",names(survivalData),ignore.case = TRUE) |
                     grepl("dod",names(survivalData),ignore.case = TRUE)) &
                    !(grepl("method",names(survivalData),ignore.case = TRUE)) &
                    !(grepl("source",names(survivalData),ignore.case = TRUE)) &
                    !(grepl("certainty",names(survivalData),ignore.case = TRUE)))
survivalData <- survivalData %>%
  mutate(across(all_of(datecols), ~ as.Date(.,tryFormats =c("%Y-%m-%d","%m-%d-%Y", "%m/%d/%Y"))))

# Collapse male and female
survivalData$Sex <- "Female"
survivalData$Sex[grepl("Male",survivalData$Sex_Class_at_HLES)] <- "Male"
survivalData$Sex <- factor(survivalData$Sex,
                           levels=c("Female","Male"))

# Rename weight and breed classes
survivalData$Weight_Class_10KGBin_at_HLES <-
  paste(survivalData$Weight_Class_10KGBin_at_HLES,"kg")
survivalData$Weight_Class_10KGBin_at_HLES <- 
  factor(survivalData$Weight_Class_10KGBin_at_HLES,
         levels=sort(unique(survivalData$Weight_Class_10KGBin_at_HLES)))
survivalData$Breed_Class <- factor(survivalData$Breed_Class,
                                   levels=c("Non-AKC-Recognized or Mixed Breed",
                                            "AKC-Recognized Breed"))
survivalData$Size_Class_at_HLES <- survivalData$Breed_Size_Class_at_HLES
survivalData$Size_Class_at_HLES <- gsub("AKC-recognized purebred ","",
                                        survivalData$Size_Class_at_HLES)
survivalData$Size_Class_at_HLES <- gsub("non-AKC and mixed breed ","",
                                        survivalData$Size_Class_at_HLES)
survivalData$Size_Class_at_HLES <- factor(survivalData$Size_Class_at_HLES,
                                          levels=
                                            rev(c("Giant dogs",
                                                  "Large dogs",
                                                  "Standard dogs",
                                                  "Medium dogs",
                                                  "Toy and Small dogs")
                                            ))

survivalData$Breed_Size_Class_at_HLES <- factor(survivalData$Breed_Size_Class_at_HLES,
                                                levels=
                                                  c("Giant AKC-recognized purebred dogs",
                                                    "Giant non-AKC and mixed breed dogs",
                                                    "Large AKC-recognized purebred dogs",
                                                    "Large non-AKC and mixed breed dogs",
                                                    "Standard AKC-recognized purebred dogs",
                                                    "Standard non-AKC and mixed breed dogs",
                                                    "Medium AKC-recognized purebred dogs",
                                                    "Medium non-AKC and mixed breed dogs",
                                                    "Toy and Small AKC-recognized purebred dogs",
                                                    "Toy and Small non-AKC and mixed breed dogs"
                                                  ))

if ("Estimated_DOB" %in% names(survivalData)) {
  survivalData <- rename(survivalData,st_estimated_dob=Estimated_DOB)
}
# Should be able to use "baseline_age" for first.age
survivalData$first.age <- time_length(interval(survivalData$st_estimated_dob,survivalData$dap_pack_date),"years")
# Should be able to use "survival_age" for last.age
survivalData$last.age <- time_length(interval(survivalData$st_estimated_dob,survivalData$survival_date),"years")
survivalData$event <- as.numeric(survivalData$survival_status=="Dead")

survivalData <- subset(survivalData,last.age > first.age)

save(survivalData,file=file.path(datfolder,"SurvivalData.RData"))
write.csv(survivalData,file.path(datfolder,"SurvivalData.csv"),row.names = FALSE)

# 
# level.list <- list()
# level.list[[1]] <- c("Rural","Suburban","Urban")
# survivalData$home_area <- factor(
#   c("Urban","Suburban","Rural")[survivalData$de_home_area_type],
#   levels=level.list[[1]])
# level.list[[2]] <- c("New England","Middle Atlantic",
#                      "East North Central","West North Central",
#                      "South Atlantic","East South Central","West South Central",
#                      "Mountain","Pacific")
# survivalData$census_division <- factor( 
#   c("New England","Middle Atlantic",
#     "East North Central","West North Central",
#     "South Atlantic","East South Central","West South Central",
#     "Mountain","Pacific")[
#       survivalData$oc_primary_residence_census_division],
#   levels<-level.list[[2]])
# level.list[[3]] <- sort(unique(survivalData$oc_primary_residence_state))
# survivalData$state <- factor(
#   survivalData$oc_primary_residence_state,
#   levels=level.list[[3]])

