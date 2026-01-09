library(lubridate)
library(survival)
library(flexsurv)
library(survminer)
library(dplyr)
library(viridisLite)
datfolder <- "data"
resultsfolder <- "results"
figfolder <- "figures"
load(file.path(datfolder,"SurvivalData.RData"))
surv <- Surv(time=survivalData$first.age,
             time2=survivalData$last.age,
             event=survivalData$event,
             type='counting')

fit.b.s <- survfit(surv ~ Size_Class_at_HLES + Breed_Class + Sex,data=survivalData)

fit.b.s.sum <- summary(fit.b.s)
fit.b.s.sum.df <- as.data.frame(fit.b.s.sum$table)
fit.b.s.quant <- quantile(fit.b.s)
fit.b.s.sum.df$q25 <- fit.b.s.quant$quantile[,1]
fit.b.s.sum.df$q25.lcl <- fit.b.s.quant$lower[,1]
fit.b.s.sum.df$q25.ucl <- fit.b.s.quant$upper[,1]
fit.b.s.sum.df$q75 <- fit.b.s.quant$quantile[,3]
fit.b.s.sum.df$q75.lcl <- fit.b.s.quant$lower[,3]
fit.b.s.sum.df$q75.ucl <- fit.b.s.quant$upper[,3]
fit.b.s.sum.df$stratum <- rownames(fit.b.s.sum$table)
strata<-trimws(unlist(strsplit(unlist(
  strsplit(rownames(fit.b.s.sum$table),",")),"=")))
strata <- strata[seq(2,length(strata),2)]
strata <- matrix(strata,ncol=3,byrow=TRUE)
fit.b.s.sum.df$Size_Class_at_HLES <-
  factor(strata[,1],levels=levels(survivalData$Size_Class_at_HLES))
fit.b.s.sum.df$Breed_Class <- 
  factor(strata[,2],levels=levels(survivalData$Breed_Class))
fit.b.s.sum.df$Sex <- 
  factor(strata[,3],levels=levels(survivalData$Sex))

write.csv(fit.b.s.sum.df,file.path(resultsfolder,"Survival_DAP_demographics.SumStats.csv"))

## Size effect by breed*sex strata

cox.b.s <- coxph(surv ~ Size_Class_at_HLES + Breed_Class + Sex,data=survivalData)
cox.b.s.pvalues <- data.frame()
for (Breed_Class_now in levels(survivalData$Breed_Class)) {
  for (Sex_now in levels(survivalData$Sex)) {
    cox.b.tmp <- coxph(surv ~ Size_Class_at_HLES,data=survivalData,
                       subset=Breed_Class == Breed_Class_now &
                         Sex == Sex_now)
    cox.b.s.pvalues <- rbind(cox.b.s.pvalues,
                             data.frame(Breed_Class=Breed_Class_now,
                                        Sex=Sex_now,
                                        log.rank.pvalue=summary(cox.b.tmp)$sctest["pvalue"]))
  }
}
# Make factors
cox.b.s.pvalues$Breed_Class <- factor(cox.b.s.pvalues$Breed_Class,
                                      levels=levels(survivalData$Breed_Class))
cox.b.s.pvalues$Sex <- factor(cox.b.s.pvalues$Sex,
                                      levels=levels(survivalData$Sex))
print(cox.b.s.pvalues)
write.csv(cox.b.s.pvalues,file.path(resultsfolder,"Survival_DAP_demographics_sizeeffect-pvals.csv"),
          row.names = FALSE)
# Reverse levels for figure
levels(survivalData$Size_Class_at_HLES) <- rev(levels(survivalData$Size_Class_at_HLES))
plt.b.km.s<-ggsurvplot_facet(fit.b.s,data=survivalData,
                             facet.by=c("Sex","Breed_Class"),
                             censor=FALSE,xlab="Age (yr)",
                             surv.median.line = "hv",short.panel.labs=TRUE
                             )+scale_color_viridis_d(option="turbo",end=0.8)+
  geom_text(aes(x=20,y=0.75,label=paste0("p=",signif(log.rank.pvalue,2))),
            data=cox.b.s.pvalues,hjust=0)+theme(legend.position = "bottom")

print(plt.b.km.s)

######

cox.b.breed.pvalues <- data.frame()
for (Size_Class_at_HLES_now in levels(survivalData$Size_Class_at_HLES)) {
  for (Sex_now in levels(survivalData$Sex)) {
    cox.b.tmp <- coxph(surv ~ Breed_Class,data=survivalData,
                       subset=Size_Class_at_HLES == Size_Class_at_HLES_now &
                         Sex == Sex_now)
    cox.b.breed.pvalues <- rbind(cox.b.breed.pvalues,
                                 data.frame(Size_Class_at_HLES=Size_Class_at_HLES_now,
                                            Sex=Sex_now,
                                            log.rank.pvalue=summary(cox.b.tmp)$sctest["pvalue"]))
  }
}
# Make factors
cox.b.breed.pvalues$Size_Class_at_HLES <- factor(cox.b.breed.pvalues$Size_Class_at_HLES,
                                      levels=levels(survivalData$Size_Class_at_HLES))
cox.b.breed.pvalues$Sex <- factor(cox.b.breed.pvalues$Sex,
                              levels=levels(survivalData$Sex))
print(cox.b.breed.pvalues)
write.csv(cox.b.breed.pvalues,
          file.path(resultsfolder,"Survival_DAP_demographics_breedeffect-pvals.csv"),
          row.names = FALSE)
plt.b.km.breed<-ggsurvplot_facet(fit.b.s,data=survivalData,
                                 facet.by=c("Sex","Size_Class_at_HLES"),
                                 censor=FALSE,xlab="Age (yr)",
                                 surv.median.line = "hv",short.panel.labs=TRUE
                                 )+scale_color_viridis_d(option="magma",end = 0.8)+
  geom_text(aes(x=15,y=0.75,label=paste0("p=",signif(log.rank.pvalue,2))),
            data=cox.b.breed.pvalues,hjust=0)
print(plt.b.km.breed)
ggsave(file.path(figfolder,"Supp.Fig.Survival_DAP_demographics.breed_effect.pdf"),
       plt.b.km.breed,
       height=4,width=7,scale=1.6)

######

cox.b.sex.pvalues <- data.frame()
for (Size_Class_at_HLES_now in levels(survivalData$Size_Class_at_HLES)) {
  for (Breed_Class_now in levels(survivalData$Breed_Class)) {
    cox.b.tmp <- coxph(surv ~ Sex,data=survivalData,
                       subset=Size_Class_at_HLES == Size_Class_at_HLES_now &
                         Breed_Class == Breed_Class_now)
    cox.b.sex.pvalues <- rbind(cox.b.sex.pvalues,
                               data.frame(Size_Class_at_HLES=Size_Class_at_HLES_now,
                                          Breed_Class=Breed_Class_now,
                                          log.rank.pvalue=summary(cox.b.tmp)$sctest["pvalue"]))
  }
}
cox.b.sex.pvalues$Size_Class_at_HLES <- factor(cox.b.sex.pvalues$Size_Class_at_HLES,
                                               levels=levels(survivalData$Size_Class_at_HLES))
cox.b.sex.pvalues$Breed_Class <- factor(cox.b.sex.pvalues$Breed_Class,
                                      levels=levels(survivalData$Breed_Class))

print(cox.b.sex.pvalues)
write.csv(cox.b.sex.pvalues,
          file.path(resultsfolder,"Survival_DAP_demographics_sexeffect-pvals.csv"),
          row.names = FALSE)
plt.b.km.sex<-ggsurvplot_facet(fit.b.s,data=survivalData,
                               facet.by=c("Breed_Class","Size_Class_at_HLES"),
                               censor=FALSE,xlab="Age (yr)",
                               surv.median.line = "hv",short.panel.labs=TRUE
                              )+scale_color_viridis_d(option="magma",end = 0.8)+
  geom_text(aes(x=15,y=0.75,label=paste0("p=",signif(log.rank.pvalue,2))),
            data=cox.b.sex.pvalues,hjust=0)

print(plt.b.km.sex)
ggsave(file.path(figfolder,"Supp.Fig.Survival_DAP_demographics.sex_effect.pdf"),
       plt.b.km.sex,
       height=4,width=7,scale=1.6)


#### 

cox.b.sex.breed.pvalues <- data.frame()
for (Size_Class_at_HLES_now in levels(survivalData$Size_Class_at_HLES)) {
    cox.b.tmp <- coxph(surv ~ Sex*Breed_Class,data=survivalData,
                       subset=Size_Class_at_HLES == Size_Class_at_HLES_now)
    cox.b.sex.breed.pvalues <- rbind(cox.b.sex.breed.pvalues,
                               data.frame(Size_Class_at_HLES=Size_Class_at_HLES_now,
                                          log.rank.pvalue=summary(cox.b.tmp)$sctest["pvalue"]))
}
cox.b.sex.breed.pvalues$Size_Class_at_HLES <- factor(cox.b.sex.breed.pvalues$Size_Class_at_HLES,
                                               levels=levels(survivalData$Size_Class_at_HLES))

print(cox.b.sex.breed.pvalues)
write.csv(cox.b.sex.breed.pvalues,
          file.path(resultsfolder,"Survival_DAP_demographics_sex.breed.effect-pvals.csv"),
          row.names = FALSE)

####

plt.forest.median.b.s<-ggplot(fit.b.s.sum.df)+
  geom_errorbarh(aes(xmin=`0.95LCL`,xmax=`0.95UCL`,
                     y=interaction(Sex,Breed_Class),color=Size_Class_at_HLES),height=0.05)+
  geom_errorbarh(aes(xmin=`q25.lcl`,xmax=`q25.ucl`,
                     y=interaction(Sex,Breed_Class),color=Size_Class_at_HLES),height=0)+
  geom_errorbarh(aes(xmin=`q75.lcl`,xmax=`q75.ucl`,
                     y=interaction(Sex,Breed_Class),color=Size_Class_at_HLES),height=0)+
  geom_point(aes(x=median,y=interaction(Sex,Breed_Class),color=Size_Class_at_HLES,shape="median"))+
  geom_point(aes(x=q25,y=interaction(Sex,Breed_Class),color=Size_Class_at_HLES,shape="IQR"))+
  geom_point(aes(x=q75,y=interaction(Sex,Breed_Class),color=Size_Class_at_HLES,shape="IQR"))+
  xlab("Median & IQR lifespan [CI]")+
  scale_shape_discrete("",limits=rev)+scale_color_viridis_d(option="turbo",end=0.8)+
  guides(color="none")+theme_bw()+theme(legend.position="bottom")+
  coord_cartesian(xlim=c(5,20))+scale_y_discrete(limits=rev)+
  scale_x_continuous(minor_breaks=seq(5,20))+
  facet_wrap(~Size_Class_at_HLES ,ncol=1)
print(plt.forest.median.b.s)

plt.forest.mean.b.s<-ggplot(fit.b.s.sum.df)+
  geom_errorbarh(aes(xmin=rmean-1.96*`se(rmean)`,
                     xmax=rmean+1.96*`se(rmean)`,
                     y=interaction(Sex,Breed_Class),color=Size_Class_at_HLES),height=0)+
  geom_point(aes(x=rmean,y=interaction(Sex,Breed_Class),color=Size_Class_at_HLES))+
  xlab("Mean lifespan [CI]")+scale_color_viridis_d(option="turbo",end=0.8)+
  guides(color="none")+theme_bw()+
  coord_cartesian(xlim=c(5,20))+scale_y_discrete(limits=rev)+
  scale_x_continuous(minor_breaks=seq(5,20))+
  facet_wrap(~Size_Class_at_HLES ,ncol=1)
print(plt.forest.mean.b.s)


ggsave(file.path(figfolder,"Fig.1.Survival_DAP_demographics.KM.forest.pdf"),
       ggarrange(plt.b.km.s,plt.forest.median.b.s,ncol=1,labels = "AUTO"),
       height=7,width=6,scale=1.3)
ggsave(file.path(figfolder,"Supp.Fig.Survival_DAP_demographics.forest.mean.pdf"),
       plt.forest.mean.b.s)
