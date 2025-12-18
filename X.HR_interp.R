library(tidyverse)
library(survival)
library(ggplot2)
figfolder <- "Figures"
resultsfolder <- "Results"
## Interpreting HR in terms of years gained or lost
load(file.path(datfolder,"SurvivalData.RData"))
surv <- Surv(time=survivalData$first.age,
             time2=survivalData$last.age,
             event=survivalData$event,
             type='counting')
fit.b.s <- survfit(surv ~ Size_Class_at_HLES + Breed_Class + Sex,
                   data=survivalData)
fit.b.s.quant <- quantile(fit.b.s)
fit.b.s.sum <- summary(fit.b.s)
fit.b.s.sum.table <- fit.b.s.sum$table
n.strata<-nrow(fit.b.s.sum.table)
ifirst <- 1+c(0,1+which(diff(fit.b.s.sum$time)<0))
ilast <- c(which(diff(fit.b.s.sum$time)<0),length(fit.b.s.sum$time))
haz.b.s.list <- list()
gamma.df <- data.frame(stratum=names(fit.b.s$strata),gamma=NA,gamma.lcl=NA,gamma.ucl=NA)
# gamma.df.1 <- gamma.df
gamma.plt.list <- list()
for (j in 1:n.strata) {
  indx <- ifirst[j]:ilast[j]
  trange<-range(fit.b.s.quant$quantile[j,])
  x <- fit.b.s.sum$time[indx]
  y <- fit.b.s.sum$cumhaz[indx]
  model_formula <- y ~ exp(a + b * x) + c
  start_params <- list(a = 0, b = 0.4, c = 0.1) 
  dat<-data.frame(x, y)
  dat.iqr <- subset(dat,x > trange[1] & x < trange[2])
  # Fit the model
  fit <- nls(model_formula, 
             data = dat.iqr, 
             start = start_params)
  gamma.plt.list[[j]]<-ggplot(dat)+
    geom_point(aes(x=x,y=y))+
    geom_line(aes(x=x,y=predict(fit,newdata=dat)),col="red")+
    geom_line(aes(x=x,y=predict(fit)),data=dat.iqr,col="red",linewidth=2)+
    ylim(NA,max(y))+
    ggtitle(label="",subtitle = gsub(", ","\n",names(fit.b.s$strata)[j]))
  print(gamma.plt.list[[j]])
  gamma.df[j,2:4]<-c(coef(fit)["b"],confint(fit)["b",])
}

# Convert to fraction of lifespan
gamma.frac.df <- gamma.df
gamma.frac.df[,2:4] <- gamma.df[,2:4]*fit.b.s.quant$quantile[,"50"]


HR.interp <- gamma.df
HR.interp$delta.0.5yr <- exp(0.5*gamma.df$gamma)
HR.interp$delta.1.0yr <- exp(1*gamma.df$gamma)
HR.interp$delta.2.0yr <- exp(2*gamma.df$gamma)
HR.interp$delta.4.0yr <- exp(4*gamma.df$gamma)
HR.interp$delta.0.5yr.lcl <- exp(0.5*gamma.df$gamma.lcl)
HR.interp$delta.1.0yr.lcl <- exp(1*gamma.df$gamma.lcl)
HR.interp$delta.2.0yr.lcl <- exp(2*gamma.df$gamma.lcl)
HR.interp$delta.4.0yr.lcl <- exp(4*gamma.df$gamma.lcl)
HR.interp$delta.0.5yr.ucl <- exp(0.5*gamma.df$gamma.ucl)
HR.interp$delta.1.0yr.ucl <- exp(1*gamma.df$gamma.ucl)
HR.interp$delta.2.0yr.ucl <- exp(2*gamma.df$gamma.ucl)
HR.interp$delta.4.0yr.ucl <- exp(4*gamma.df$gamma.ucl)

write.csv(HR.interp,file.path(resultsfolder,"HR.interp.csv"),row.names = FALSE)
HR.interp.sum<-t(apply(HR.interp[,-1],2,quantile))
write.csv(HR.interp.sum,file.path(resultsfolder,"HR.interp.sum.csv"))
