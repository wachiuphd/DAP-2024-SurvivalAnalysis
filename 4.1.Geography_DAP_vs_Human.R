library(tidyverse)
library(survival)

datfolder <- "data"
resultsfolder <- "results"
figfolder <- "figures"

load(file.path(resultsfolder,"Geoeffect-Cox-results.Rdata"))
surv <- Surv(time=survivalData$first.age,
             time2=survivalData$last.age,
             event=survivalData$event,
             type='counting')

###### Correlation between hazard ratios by states and age-adjusted mortality rate

hr.st <- exp(cbind(as.data.frame(cox.geo.adj.st$coefficients),
                   as.data.frame(confint(cox.geo.adj.st))))
hr.st$Dog.SE <- log(hr.st$`97.5 %`/hr.st$`2.5 %`)/(2*qnorm(0.975)) # used for weights
names(hr.st)[1] <- "Dog.HR"
hr.st$state.abbr <- gsub("state","",rownames(hr.st))

mr.st <- subset(read.csv(file.path(datfolder,"HDPulse_data_export.csv"),skip=4),
                (FIPS > 0 ) & (!is.na(FIPS)))
mr.st$state.abbr <- state2abbr(mr.st$State)
mr.st <- subset(mr.st,state.abbr %in% c("WA",hr.st$state.abbr))
names(mr.st)[3] <- "MR"
mr.st$MR.Ratio <- mr.st$MR/subset(mr.st,state.abbr=="WA")$MR
hr.st <- left_join(hr.st,mr.st)

lm.res <- lm(log(Dog.HR) ~ log(MR.Ratio),weights=1/Dog.SE^2,data=hr.st)
lm.res.coef <- coef(lm.res)
lm.res.sum <- summary(lm.res)
hr.mr.pval <- signif(last(lm.res.sum$coefficients["log(MR.Ratio)",]),2)
hr.mr.slope <- signif(first(lm.res.sum$coefficients["log(MR.Ratio)",]),2)
hr.mr.r <- signif(weightedCorr(log(hr.st$Dog.HR), log(hr.st$MR), weights = 1/hr.st$Dog.SE^2, method = "Pearson"),2)
hr.mr.rho <- signif(weightedCorr(log(hr.st$Dog.HR), log(hr.st$MR), weights = 1/hr.st$Dog.SE^2, method = "Spearman"),2)

# plt.hr.mr.state <-
#   ggplot(hr.st,aes(x=MR.Ratio,y=Dog.HR))+
#   geom_hline(yintercept=1)+geom_vline(xintercept = 1)+
#   geom_errorbar(aes(ymin=`2.5 %`,ymax=`97.5 %`),color="grey50")+
#   geom_label(aes(label=state.abbr),alpha=0.7)+
#   scale_x_log10(breaks=seq(0.8,1.6,0.2))+
#   scale_y_log10(breaks=seq(0.4,2.0,0.2))+#coord_cartesian(xlim=c(0.7,1.7),ylim=c(0.5,2))+
#   theme_bw()+#coord_trans(x="log10",y="log10",xlim=c(0.4,2.1),ylim=c(0.4,2.1))+
#   annotation_logticks(side="bl")+
#   geom_smooth(aes(x=MR.Ratio,y=Dog.HR,weight=1/Dog.SE^2),method="lm")+
#   annotate("text",x=0.8,y=2,label=bquote(italic(r) == .(hr.mr.r)*","~italic(rho) == .(hr.mr.rho)),hjust=0,vjust=1.5)+
#   annotate("text",x=0.8,y=2,label=bquote(italic(p) == .(hr.mr.pval)),hjust=0,vjust=3)+
#   annotate("text",x=0.8,y=2,label=bquote(slope == .(hr.mr.slope)),hjust=0,vjust=4.5)+
#   xlab("Human Age-Adjusted Comparative Mortality Ratio (CMR) [2019-2023]")+
#   ylab("Canine Demographics-Adjusted Mortality Hazard Ratio (HR)")+
#   theme(panel.grid.minor = element_blank())

plt.hr.mr.state <- ggplot(hr.st,aes(x=MR.Ratio,y=Dog.HR))+
  geom_hline(yintercept=1)+geom_vline(xintercept = 1)+
  geom_point(aes(size=1/Dog.SE^2),shape=21,stroke=1.5,color="grey33")+
  geom_text(aes(label=state.abbr),hjust=1.3,vjust=-0.5,size=3)+
  scale_x_log10(breaks=seq(0.6,1.6,0.2))+
  scale_y_log10(breaks=seq(0.6,1.6,0.2))+coord_cartesian(xlim=c(0.75,1.5),ylim=c(0.6,1.4))+
  scale_size_continuous(range=c(0.01,6))+
  theme_bw()+#coord_trans(x="log10",y="log10",xlim=c(0.4,2.1),ylim=c(0.4,2.1))+
  annotation_logticks(side="bl")+
  geom_smooth(aes(x=MR.Ratio,y=Dog.HR,weight=1/Dog.SE^2),method="lm")+
  annotate("text",x=0.75,y=1.4,label=bquote(italic(r) == .(hr.mr.r)*","~italic(rho) == .(hr.mr.rho)),hjust=0,vjust=1.5)+
  annotate("text",x=0.75,y=1.4,label=bquote(italic(p) == .(hr.mr.pval)),hjust=0,vjust=3)+
  annotate("text",x=0.75,y=1.4,label=bquote(slope == .(hr.mr.slope)),hjust=0,vjust=4.5)+
  xlab("Human Age-Adjusted Comparative Mortality Ratio (CMR) [2019-2023]")+
  ylab("Canine Demographics-Adjusted Mortality Hazard Ratio (HR)")+
  theme(panel.grid.minor = element_blank(),legend.position = "none")
print(plt.hr.mr.state)

ggsave(file.path(figfolder,"Fig.2.Geoeffect_State_DogvsHumanMortRate.pdf"),
       plt.hr.mr.state,height=4,width=4,scale=1.33)

###### Sensitivity Analysis - Use only 2023 mortality rates

hr.st <- exp(cbind(as.data.frame(cox.geo.adj.st$coefficients),
                   as.data.frame(confint(cox.geo.adj.st))))
hr.st$Dog.SE <- log(hr.st$`97.5 %`/hr.st$`2.5 %`)/(2*qnorm(0.975)) # used for weights
names(hr.st)[1] <- "Dog.HR"
hr.st$state.abbr <- gsub("state","",rownames(hr.st))

mr23.st <- subset(read.csv(file.path(datfolder,"HDPulse_data_export-2023.csv"),skip=4),
                (FIPS > 0 ) & (!is.na(FIPS)))
mr23.st$state.abbr <- state2abbr(mr23.st$State)
mr23.st <- subset(mr23.st,state.abbr %in% c("WA",hr.st$state.abbr))
names(mr23.st)[3] <- "MR"
mr23.st$MR.Ratio <- mr23.st$MR/subset(mr23.st,state.abbr=="WA")$MR
hr.st <- left_join(hr.st,mr23.st)

lm.res <- lm(log(Dog.HR) ~ log(MR.Ratio),weights=1/Dog.SE^2,data=hr.st)
lm.res.coef <- coef(lm.res)
lm.res.sum <- summary(lm.res)
hr.mr23.pval <- signif(last(lm.res.sum$coefficients["log(MR.Ratio)",]),2)
hr.mr23.slope <- signif(first(lm.res.sum$coefficients["log(MR.Ratio)",]),2)
hr.mr23.r <- signif(weightedCorr(log(hr.st$Dog.HR), log(hr.st$MR), weights = 1/hr.st$Dog.SE^2, method = "Pearson"),2)
hr.mr23.rho <- signif(weightedCorr(log(hr.st$Dog.HR), log(hr.st$MR), weights = 1/hr.st$Dog.SE^2, method = "Spearman"),2)

plt.hr.mr23.state <-
  ggplot(hr.st,aes(x=MR.Ratio,y=Dog.HR))+
  geom_hline(yintercept=1)+geom_vline(xintercept = 1)+
  geom_point(aes(size=1/Dog.SE^2),shape=21,stroke=1.5,color="grey33")+
  geom_text(aes(label=state.abbr),hjust=1.3,vjust=-0.5,size=3)+
  scale_x_log10(breaks=seq(0.6,1.6,0.2))+
  scale_y_log10(breaks=seq(0.6,1.6,0.2))+coord_cartesian(xlim=c(0.75,1.5),ylim=c(0.6,1.4))+
  scale_size_continuous(range=c(0.01,6))+
  theme_bw()+
  annotation_logticks(side="bl")+
  geom_smooth(aes(x=MR.Ratio,y=Dog.HR,weight=1/Dog.SE^2),method="lm")+
  annotate("text",x=0.75,y=1.4,label=bquote(italic(r) == .(hr.mr23.r)*","~italic(rho) == .(hr.mr23.rho)),hjust=0,vjust=1.5)+
  annotate("text",x=0.75,y=1.4,label=bquote(italic(p) == .(hr.mr23.pval)),hjust=0,vjust=3)+
  annotate("text",x=0.75,y=1.4,label=bquote(slope == .(hr.mr23.slope)),hjust=0,vjust=4.5)+
  xlab("Human Age-Adjusted Comparative Mortality Ratio (CMR) [2023]")+
  ylab("Canine Demographics-Adjusted Mortality Hazard Ratio (HR)")+
  theme(panel.grid.minor = element_blank(),legend.position = "none")
print(plt.hr.mr23.state)

ggsave(file.path(figfolder,"SuppFig.Geoeffect_State_DogvsHuman2023MortRate.pdf"),
       plt.hr.mr23.state,height=4,width=4,scale=1.25)

###### Correlation between hazard ratios by states and US life expectancy

hr.st <- exp(cbind(as.data.frame(cox.geo.adj.st$coefficients),
                   as.data.frame(confint(cox.geo.adj.st))))
hr.st$Dog.SE <- log(hr.st$`97.5 %`/hr.st$`2.5 %`)/(2*qnorm(0.975)) # used for weights
names(hr.st)[1] <- "Dog.HR"
hr.st$state.abbr <- gsub("state","",rownames(hr.st))
# hr.st[nrow(hr.st)+1,] <- data.frame(1,NA,NA,NA,"WA")

le.st <- read.csv(file.path(datfolder,"U.S._State_Life_Expectancy_by_Sex__2021.csv"))
le.st$state.abbr <- state2abbr(le.st$State)
le.st <- subset(le.st,Sex=="Total")
le.st$LE.Ratio <- log(le.st$LE/subset(le.st,state.abbr=="WA")$LE)
hr.st <- left_join(hr.st,le.st)

lm.res <- lm(Dog.HR ~ LE,weights=1/Dog.SE^2,data=subset(hr.st,!is.na(Dog.SE)))
lm.res.coef <- coef(lm.res)
lm.res.sum <- summary(lm.res)
hr.le.pval <- signif(last(lm.res.sum$coefficients["LE",]),2)
hr.le.slope <- signif(first(lm.res.sum$coefficients["LE",]),2)
hr.le.r <- signif(weightedCorr(hr.st$Dog.HR, hr.st$LE, weights = 1/hr.st$Dog.SE^2, method = "Pearson"),2)
hr.le.rho <- signif(weightedCorr(hr.st$Dog.HR, hr.st$LE, weights = 1/hr.st$Dog.SE^2, method = "Spearman"),2)

plt.hr.le.state <-
  ggplot(hr.st,aes(x=LE,y=Dog.HR))+
  geom_point(aes(size=1/Dog.SE^2),shape=21,stroke=1.5,color="grey33")+
  geom_text(aes(label=state.abbr),hjust=1.3,vjust=-0.5,size=3)+
  scale_x_continuous(breaks=seq(65,85))+
  theme_bw()+coord_trans(y="log10",ylim=c(0.5,2.1))+
  annotation_logticks(side="l",scaled=FALSE)+
  geom_smooth(aes(weight=1/Dog.SE^2),method="lm")+
  geom_abline(slope=-1/78.2,intercept=2,color="red",linetype="dashed")+
  annotate("text",x=71,y=2,label="Linear fit to state-level HRs",hjust=0,vjust=0)+
  annotate("text",x=71,y=2,label=bquote(italic(r) == .(hr.le.r)*","~italic(rho) == .(hr.le.rho)),hjust=0,vjust=1.5)+
  annotate("text",x=71,y=2,label=bquote(italic(p) == .(hr.le.pval)),hjust=0,vjust=3)+
  annotate("text",x=71,y=2,label=bquote(slope == .(hr.le.slope)),hjust=0,vjust=4.5)+
  xlab("Human Life Expectancy (yrs)\nby State in 2021")+
  ylab("Demographically Adjusted DAP Mortality Hazard Ratio (HR)")+
  theme(panel.grid.minor = element_blank(),legend.position = "none")
print(plt.hr.le.state)

hr.le.slope.rel <- signif(78.2*first(lm.res.sum$coefficients["LE",]),2)
plt.hr.le.state.relative<-ggplot(hr.st,aes(x=LE/78.2,y=Dog.HR))+
  geom_hline(yintercept=1)+geom_vline(xintercept = 1)+
  geom_point(aes(size=1/Dog.SE^2),shape=21,stroke=1.5,color="grey33")+
  geom_text(aes(label=state.abbr),hjust=1.3,vjust=-0.5,size=3)+
  theme_bw()+theme(panel.grid = element_blank(),legend.position="none")+
  coord_cartesian(xlim=c(0.9,1.05),ylim=c(0.6,1.4))+
  geom_smooth(aes(weight=1/Dog.SE^2),method="lm")+
  geom_abline(slope=-1,intercept=2,color="red",linetype="dashed")+
  annotate("text",x=1.01,y=1.4,label=bquote(italic(r) == .(hr.le.r)),hjust=0,vjust=0)+
  annotate("text",x=1.01,y=1.4,label=bquote(italic(rho) == .(hr.le.rho)),hjust=0,vjust=1.5)+
  annotate("text",x=1.01,y=1.4,label=bquote(italic(p) == .(hr.le.pval)),hjust=0,vjust=3)+
  annotate("text",x=1.01,y=1.4,label=bquote(slope == .(hr.le.slope.rel)),hjust=0,vjust=4.5)+
  xlab("State Human Life Expectancy\nRelative to WA [78.2 years] in 2021")+
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
plt.hr.le.state.comb<-ggarrange(plt.hr.le.state,plt.hr.le.5,ncol=2,widths = c(2,1),labels = "AUTO")

ggsave(file.path(figfolder,"SuppFig.Geoeffect_State_DogvsHumanLE.pdf"),
       plt.hr.le.state.comb,height=4,width=6,scale=1.25)
