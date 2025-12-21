library(haven)
library(lubridate)
library(survival)
library(flexsurv)
library(survminer)
library(dplyr)
library(viridisLite)
library(usdata)
library(wCorr)
source("ggforest2.R")
datfolder <- "data"
resultsfolder <- "results"
figfolder <- "figures"

## Join geographic variables
load(file.path(datfolder,"SurvivalData.RData"))
load(file.path(datfolder,"DAP_2024_HLES_dog_owner_v1.0.RData"))
survivalData <- left_join(survivalData,
                          HLES_dog_owner[,c("dog_id",
                                            "de_home_area_type",
                                            "oc_primary_residence_census_division",
                                            "oc_primary_residence_state")])
remove(HLES_dog_owner)
level.list <- list()
# Largest N is reference for each one
level.list[[1]] <- c("Suburban","Rural","Urban") 
level.list[[2]] <- rev(c("New England","Middle Atlantic",
                         "East North Central","West North Central",
                         "South Atlantic","East South Central","West South Central",
                         "Mountain","Pacific")) 
level.list[[3]] <- sort(unique(survivalData$oc_primary_residence_state))
ref.state <- which(level.list[[3]]=="WA")
level.list[[3]] <- c(level.list[[3]][ref.state],level.list[[3]][-ref.state])

survivalData$de_home_area_type <- factor(
  c("Urban","Suburban","Rural")[survivalData$de_home_area_type],
  levels=level.list[[1]])

survivalData$oc_primary_residence_census_division <- factor( 
  c("New England","Middle Atlantic",
    "East North Central","West North Central",
    "South Atlantic","East South Central","West South Central",
    "Mountain","Pacific")[
      survivalData$oc_primary_residence_census_division],
  levels=level.list[[2]])
survivalData$oc_primary_residence_state <- factor(
  survivalData$oc_primary_residence_state,
  levels=level.list[[3]])

survivalData <- rename(survivalData,
                       home_area = de_home_area_type,
                       census_division = oc_primary_residence_census_division,
                       state = oc_primary_residence_state)
  
surv <- Surv(time=survivalData$first.age,
             time2=survivalData$last.age,
             event=survivalData$event,
             type='counting')



###### Unadjusted

fit.geo.unadj.ha <- survfit(surv ~ home_area,data=survivalData)
cox.geo.unadj.ha <- coxph(surv ~ home_area,data=survivalData)
fit.geo.unadj.cd <- survfit(surv ~ census_division,data=survivalData)
cox.geo.unadj.cd <- coxph(surv ~ census_division,data=survivalData)
fit.geo.unadj.st <- survfit(surv ~ state,data=survivalData)
cox.geo.unadj.st <- coxph(surv ~ state,data=survivalData)

cox.geo.unadj.pvalues <- data.frame(Analysis="Unadjusted",
                                    Variable=c("home_area","census_division","state"),
                                    log.rank.pvalue=c(summary(cox.geo.unadj.ha)$sctest["pvalue"],
                                                      summary(cox.geo.unadj.cd)$sctest["pvalue"],
                                                      summary(cox.geo.unadj.st)$sctest["pvalue"]),
                                    cox.zph.pvalue=c(
                                      cox.zph(cox.geo.unadj.ha)$table[1,"p"],
                                      cox.zph(cox.geo.unadj.cd)$table[1,"p"],
                                      cox.zph(cox.geo.unadj.st)$table[1,"p"])
)
print(cox.geo.unadj.pvalues)

###### Distribution of size, breed, sex

## Home area 

plt.size.by.area.ha <- ggplot(survivalData, 
                           aes(x = home_area, fill = Size_Class_at_HLES)) +
  scale_fill_viridis_d()+ggtitle("Demographics by Home Area Type")+
  scale_x_discrete(limits=c("Rural","Suburban","Urban"))+
  geom_bar(position="fill")+ylab("proportion")
plt.breed.by.area.ha <- ggplot(survivalData, 
                            aes(x = home_area, fill = Breed_Class)) +
  scale_fill_viridis_d(end=0.8)+
  scale_x_discrete(limits=c("Rural","Suburban","Urban"))+
  geom_bar(position="fill")+ylab("proportion")
plt.sex.by.area.ha <- ggplot(survivalData, 
                          aes(x = home_area, fill = Sex)) +
  scale_fill_viridis_d(end=0.8)+
  scale_x_discrete(limits=c("Rural","Suburban","Urban"))+
  geom_bar(position="fill")+ylab("proportion")

ggsave(file.path(figfolder,"Supp.Fig.Geoeffect_DemographicDist_homearea.pdf"),
       ggarrange(plt.size.by.area.ha,
                plt.breed.by.area.ha,
                plt.sex.by.area.ha,ncol=1),
       height=7,width=5,scale=1.2)

## Census Division
plt.size.by.area.cd <- ggplot(survivalData, 
                           aes(x = as.factor(as.integer(census_division)), fill = Size_Class_at_HLES)) +
  geom_bar(position="fill")+ggtitle("Demographics by Census Division")+
  scale_fill_viridis_d()+ylab("proportion")+xlab("Census Division #")
plt.breed.by.area.cd <- ggplot(survivalData, 
                            aes(x = as.factor(as.integer(census_division)), fill = Breed_Class)) +
  geom_bar(position="fill")+scale_fill_viridis_d(end=0.8)+ylab("proportion")+xlab("Census Division #")
plt.sex.by.area.cd <- ggplot(survivalData, 
                          aes(x = as.factor(as.integer(census_division)), fill = Sex)) +
  geom_bar(position="fill")+scale_fill_viridis_d(end=0.8)+ylab("proportion")+xlab("Census Division #")

ggsave(file.path(figfolder,"Supp.Fig.Geoeffect_DemographicDist_censusdivision.pdf"),
       ggarrange(plt.size.by.area.cd,
                 plt.breed.by.area.cd,
                 plt.sex.by.area.cd,ncol=1),
       height=7,width=5,scale=1.2)

## State

plt.size.by.area.st <- ggplot(survivalData, 
                           aes(y = state, fill = Size_Class_at_HLES)) +
  geom_bar(position="fill")+ggtitle("Demographics by State")+
  xlab("proportion")+ylab("state")+
  scale_y_discrete(limits=rev(sort(unique(as.character(survivalData$state)))))+
  scale_fill_viridis_d(guide = guide_legend(reverse = TRUE,position="bottom",nrow=2))


plt.breed.by.area.st <- ggplot(survivalData, 
                            aes(y = state, fill = Breed_Class)) +
  geom_bar(position="fill")+xlab("proportion")+ylab("state")+ggtitle("")+
  scale_y_discrete(limits=rev(sort(unique(as.character(survivalData$state)))))+
  scale_fill_viridis_d(end=0.8,guide = guide_legend(reverse = TRUE,position="bottom",nrow=2))

plt.sex.by.area.st <- ggplot(survivalData, 
                          aes(y = state, fill = Sex)) +
  geom_bar(position="fill")+xlab("proportion")+ylab("state")+ggtitle("")+
  scale_y_discrete(limits=rev(sort(unique(as.character(survivalData$state)))))+
  scale_fill_viridis_d(end=0.8,guide = guide_legend(reverse = TRUE,position="bottom",nrow=2))


ggsave(file.path(figfolder,"Supp.Fig.Geoeffect_DemographicDist_state.pdf"),
       ggarrange(plt.size.by.area.st,
                 plt.breed.by.area.st,
                 plt.sex.by.area.st,ncol=3,
                 widths=c(2,1,1)),
       height=7,width=8.5,scale=1.5)

###### Adjusted Analyses

## Common geo effect


fit.geo.adj.ha <- survfit(surv ~ home_area + 
                                 Size_Class_at_HLES + Breed_Class + Sex,
                               data=survivalData)
fit.geo.adj.cd <- survfit(surv ~ census_division + 
                                 Size_Class_at_HLES + Breed_Class + Sex,
                               data=survivalData)
fit.geo.adj.st <- survfit(surv ~ state + 
                                 Size_Class_at_HLES + Breed_Class +Sex,
                               data=survivalData)

cox.geo.adj.ha <- coxph(surv ~ home_area +
                                 strata(Size_Class_at_HLES,Breed_Class,Sex),
                             data=survivalData)
cox.geo.adj.cd <- coxph(surv ~ census_division +
                               strata(Size_Class_at_HLES,Breed_Class,Sex),
                             data=survivalData)
cox.geo.adj.st <- coxph(surv ~ state +
                               strata(Size_Class_at_HLES,Breed_Class,Sex),
                             data=survivalData)

cox.geo.adj.pvalues <- data.frame(Analysis="Adjusted",
                                  Variable=c("home_area","census_division","state"),
                                    log.rank.pvalue=c(summary(cox.geo.adj.ha)$sctest["pvalue"],
                                                      summary(cox.geo.adj.cd)$sctest["pvalue"],
                                                      summary(cox.geo.adj.st)$sctest["pvalue"]),
                                  cox.zph.pvalue=c(
                                    cox.zph(cox.geo.adj.ha)$table[1,"p"],
                                    cox.zph(cox.geo.adj.cd)$table[1,"p"],
                                    cox.zph(cox.geo.adj.st)$table[1,"p"]
                                  )
)
print(cox.geo.adj.pvalues)

hrplt.ha <- ggarrange(
  ggforest2(cox.geo.unadj.ha,data=survivalData,main="Unadjusted Hazard Ratio",xrange=1.5),
  ggforest2(cox.geo.adj.ha,data=survivalData,main="Adjusted Hazard Ratio",xrange=1.5),
  ncol=1)

hrplt.cd <- ggarrange(
  ggforest2(cox.geo.unadj.cd,data=survivalData,main="Unadjusted Hazard Ratio",xrange=1.5,
            noDigits = 3),
  ggforest2(cox.geo.adj.cd,data=survivalData,main="Adjusted Hazard Ratio",xrange=1.5,
            noDigits = 3),
  ncol=1)

hrplt.st <- ggarrange(
  ggforest2(cox.geo.unadj.st,data=survivalData,main="Unadjusted Hazard Ratio",xrange=2),
  ggforest2(cox.geo.adj.st,data=survivalData,main="Adjusted Hazard Ratio",xrange=2),
  ncol=2)

ggsave(file.path(figfolder,"Supp.Fig.Geoeffect_homearea.HR.pdf"),
       hrplt.ha,
       height=5,width=5,scale=1.25)
ggsave(file.path(figfolder,"Supp.Fig.Geoeffect_censusdistrict.HR.pdf"),
       hrplt.cd,
       height=7,width=5,scale=1.35)
ggsave(file.path(figfolder,"Supp.Fig.Geoeffect_state.HR.pdf"),
       hrplt.st,
       height=7,width=7,scale=2)
write.csv(rbind(cox.geo.unadj.pvalues,
                cox.geo.adj.pvalues),
          file=file.path(resultsfolder,
                         "Geoeffect-pvals.csv"),
          row.names = FALSE)

## Separate geo effect for each stratum

cox.geo.pvalues <- data.frame()
for (Size_Class_now in levels(survivalData$Size_Class_at_HLES)) {
  for (Breed_Class_now in levels(survivalData$Breed_Class)) {
    for (Sex_now in levels(survivalData$Sex)) {
      cox.geo.tmp.ha <- coxph(surv ~ home_area,data=survivalData,
                              subset=Size_Class_at_HLES==Size_Class_now &
                                Breed_Class == Breed_Class_now &
                                Sex == Sex_now)
      cox.geo.tmp.cd <- coxph(surv ~ census_division,data=survivalData,
                              subset=Size_Class_at_HLES==Size_Class_now &
                                Breed_Class == Breed_Class_now &
                                Sex == Sex_now)
      cox.geo.tmp.st <- coxph(surv ~ state,data=survivalData,
                              subset=Size_Class_at_HLES==Size_Class_now &
                                Breed_Class == Breed_Class_now &
                                Sex == Sex_now)
      zp.ha <- NA; zp.ha <- as.numeric(try(cox.zph(cox.geo.tmp.ha)$table[1,"p"]))
      zp.cd <- NA; zp.cd <- as.numeric(try(cox.zph(cox.geo.tmp.cd)$table[1,"p"]))
      zp.st <- NA; zp.st <- as.numeric(try(cox.zph(cox.geo.tmp.st)$table[1,"p"]))
      
      cox.geo.pvalues <- rbind(cox.geo.pvalues,
                               data.frame(Size_Class_at_HLES=Size_Class_now,
                                          Breed_Class=Breed_Class_now,
                                          Sex=Sex_now,
                                          log.rank.pvalue.ha=
                                            summary(cox.geo.tmp.ha)$sctest["pvalue"],
                                          log.rank.pvalue.cd=
                                            summary(cox.geo.tmp.cd)$sctest["pvalue"],
                                          log.rank.pvalue.st=
                                            summary(cox.geo.tmp.st)$sctest["pvalue"],
                                          cox.zph.pvalue.ha=zp.ha,
                                          cox.zph.pvalue.cd=zp.cd,
                                          cox.zph.pvalue.st=zp.st)
                               )
    }
  }
}
print(cox.geo.pvalues)

write.csv(cox.geo.pvalues,
          file=file.path(resultsfolder,
                         "Geoeffect-bystratum-pvals.csv"),
          row.names = FALSE)

###### Correlation between hazard ratios by states and US life expectancy?

hr.st <- exp(cbind(as.data.frame(cox.geo.adj.st$coefficients),
               as.data.frame(confint(cox.geo.adj.st))))
hr.st$Dog.SE <- log(hr.st$`97.5 %`/hr.st$`2.5 %`)/(2*qnorm(0.975)) # used for weights
names(hr.st)[1] <- "Dog.HR"
hr.st$state.abbr <- gsub("state","",rownames(hr.st))

le.st <- read.csv(file.path(datfolder,"U.S._State_Life_Expectancy_by_Sex__2021.csv"))
le.st$state.abbr <- state2abbr(le.st$State)
le.st <- subset(le.st,Sex=="Total")
le.st$LE.Ratio <- log(le.st$LE/subset(le.st,state.abbr=="WA")$LE)
hr.st <- left_join(hr.st,le.st)

lm.res <- lm(Dog.HR ~ LE,weights=1/Dog.SE^2,data=hr.st)
lm.res.coef <- coef(lm.res)
lm.res.sum <- summary(lm.res)
hr.le.pval <- signif(last(lm.res.sum$coefficients["LE",]),2)
hr.le.slope <- signif(first(lm.res.sum$coefficients["LE",]),2)
hr.le.r <- signif(weightedCorr(hr.st$Dog.HR, hr.st$LE, weights = 1/hr.st$Dog.SE^2, method = "Pearson"),2)
hr.le.rho <- signif(weightedCorr(hr.st$Dog.HR, hr.st$LE, weights = 1/hr.st$Dog.SE^2, method = "Spearman"),2)

plt.hr.le.state <-
ggplot(hr.st,aes(x=LE,y=Dog.HR))+
  geom_errorbar(aes(ymin=`2.5 %`,ymax=`97.5 %`),color="grey50")+
  geom_label(aes(label=state.abbr))+
  scale_x_continuous(breaks=seq(65,85))+
  theme_bw()+coord_trans(y="log10",ylim=c(0.4,2.1))+
  annotation_logticks(side="l",scaled=FALSE)+
  geom_smooth(aes(weight=1/Dog.SE^2),method="lm")+
  geom_abline(slope=-1/78.2,intercept=2,color="red",linetype="dashed")+
  annotate("text",x=71,y=2,label="Linear fit to state-level HRs",hjust=0,vjust=0)+
  annotate("text",x=71,y=2,label=bquote(italic(r) == .(hr.le.r)*","~italic(rho) == .(hr.le.rho)),hjust=0,vjust=1.5)+
  annotate("text",x=71,y=2,label=bquote(italic(p) == .(hr.le.pval)),hjust=0,vjust=3)+
  annotate("text",x=71,y=2,label=bquote(slope == .(hr.le.slope)),hjust=0,vjust=4.5)+
  xlab("Human Life Expectancy (yrs)\nby State")+
  ylab("Demographically Adjusted DAP Mortality Hazard Ratio (HR)")+
  theme(panel.grid.minor = element_blank())
print(plt.hr.le.state)

hr.le.slope.rel <- signif(78.2*first(lm.res.sum$coefficients["LE",]),2)
plt.hr.le.state.relative<-ggplot(hr.st,aes(x=LE/78.2,y=Dog.HR))+
  geom_hline(yintercept=1)+geom_vline(xintercept = 1)+
  geom_errorbar(aes(ymin=`2.5 %`,ymax=`97.5 %`),color="grey")+
  geom_point()+
  theme_bw()+theme(panel.grid = element_blank())+
  coord_cartesian(xlim=c(0.9,1.05),ylim=c(0.6,1.4))+
  geom_smooth(aes(weight=1/Dog.SE^2),method="lm")+
  geom_abline(slope=-1,intercept=2,color="red",linetype="dashed")+
  annotate("text",x=1.01,y=1.4,label=bquote(italic(r) == .(hr.le.r)),hjust=0,vjust=0)+
  annotate("text",x=1.01,y=1.4,label=bquote(italic(rho) == .(hr.le.rho)),hjust=0,vjust=1.5)+
  annotate("text",x=1.01,y=1.4,label=bquote(italic(p) == .(hr.le.pval)),hjust=0,vjust=3)+
  annotate("text",x=1.01,y=1.4,label=bquote(slope == .(hr.le.slope.rel)),hjust=0,vjust=4.5)+
  xlab("State Human Life Expectancy\nRelative to WA [78.2 years]")+
  ylab("Demographically-Adjusted DAP Mortality\nHazard Ratio by State (Reference=WA)")

print(plt.hr.le.state.relative)

#### Try using life expectancy as a continuous variable

rownames(le.st) <- le.st$state.abbr
survivalData$state.LE <- le.st[as.character(survivalData$state),"LE"]
cox.geo.adj.state.LE <- coxph(surv ~ state.LE +
                          strata(Size_Class_at_HLES,Breed_Class,Sex),
                        data=survivalData)
cox.geo.adj.state.LE.zph <- cox.zph(cox.geo.adj.state.LE)
survivalData$state.LE.5 <- survivalData$state.LE/5
cox.geo.adj.state.LE.5 <- coxph(surv ~ state.LE.5 +
                                strata(Size_Class_at_HLES,Breed_Class,Sex),
                              data=survivalData)
cox.geo.adj.state.LE.5.sum <- summary(cox.geo.adj.state.LE.5)
exp(c(coef(cox.geo.adj.state.LE.5),confint(cox.geo.adj.state.LE.5)))

cox.geo.adj.state.LE.5.coef <- 
  cbind(data.frame(Variable="State-Level\nHuman Life Expectancy",
                       logrankp=cox.geo.adj.state.LE.5.sum$sctest["pvalue"]),
            left_join(as.data.frame(cox.geo.adj.state.LE.5.sum$coefficients),
                      as.data.frame(cox.geo.adj.state.LE.5.sum$conf.int)),
        as.data.frame(t(cox.geo.adj.state.LE.zph$table[1,]))
  )
names(cox.geo.adj.state.LE.5.coef)[ncol(cox.geo.adj.state.LE.5.coef)-(2:0)] <- 
  paste0(names(cox.geo.adj.state.LE.5.coef)[ncol(cox.geo.adj.state.LE.5.coef)-(2:0)],".zph")
write.csv(cox.geo.adj.state.LE.5.coef,
          file=file.path(resultsfolder,
                         "Geoeffect-Cox-Human.LE.csv"),
                    row.names = FALSE)

hr.le.5.hr <- signif(cox.geo.adj.state.LE.5.coef$`exp(coef)`,2)
hr.le.5.hr.lcl <- signif(cox.geo.adj.state.LE.5.coef$`lower .95`,2)
hr.le.5.hr.ucl <- signif(cox.geo.adj.state.LE.5.coef$`upper .95`,2)
hr.le.5.pval <- as.numeric(signif(cox.geo.adj.state.LE.5.coef$logrankp,2))
plt.hr.le.5 <- 
ggplot(cox.geo.adj.state.LE.5.coef)+
  geom_hline(yintercept=1,linetype="dotted")+
  geom_point(aes(y=`exp(coef)`,x=Variable),shape=15,size=3,color="grey50")+
  geom_errorbar(aes(ymin=`lower .95`,ymax=`upper .95`,x=Variable),width=0)+
  annotate("text",x=0.5,y=2,label="Cox model for state",
           hjust=0,vjust=0)+
  annotate("text",x=0.5,y=2,label="human life expectancy",
           hjust=0,vjust=1.5)+
  annotate("text",x=0.5,y=2,label=
             bquote(HR == .(hr.le.5.hr)~"["*.(hr.le.5.hr.lcl) - .(hr.le.5.hr.ucl)*"]"),
           hjust=0,vjust=3)+
  annotate("text",x=0.5,y=2,label=
             bquote(italic(p) == .(hr.le.5.pval)),
           hjust=0,vjust=4.5)+
  ylab("")+
  xlab("per 5 yr increase in\nHuman Life Expectancy ")+
  scale_x_discrete(labels="")+
  theme_bw()+coord_trans(y="log10",ylim=c(0.4,2.1))+
  annotation_logticks(side="l",scaled=FALSE)+
  theme(panel.grid.minor = element_blank(),panel.grid.major.x=element_blank())
print(plt.hr.le.5)
fig4<-ggarrange(plt.hr.le.state,plt.hr.le.5,ncol=2,widths = c(2,1),labels = "AUTO")

ggsave(file.path(figfolder,"Fig.4.Geoeffect_State_DogvsHuman.pdf"),
       fig4,height=4,width=6,scale=1.25)


