library(survival)
library(coxme)

postscript("figure1.ps",horizontal=T,paper="letter")

################################################################################
# Frailty Cox PH models using time-to-event data
################################################################################

zlist <- list(time=0,event=0,age=0,female=0,arf=0,carddys=0,chf=0,crf=0,cvd=0,
  diabcomp=0,malig=0,pulmoned=0,shock=0,inst.volume=0,teaching.hosp=0,
  hosp.revasc=0,hosp.cath=0,inst="")

xdat <- data.frame(scan(r.level2.dat",zlist))

model1 <- coxph(Surv(time,event) ~ age + female + arf + carddys + chf + crf + 
  cvd + diabcomp + malig + pulmoned + shock + inst.volume + teaching.hosp +
  hosp.revasc + hosp.cath + frailty(inst,distribution="gaussian"),
  data=xdat)

frailty.var <- model1$history[[1]]$theta

sf1 <- survfit(model1)

surv.time <- sf1$time
H0 <- sf1$cumhaz

h0 <- c(H0[1],diff(H0))

frailty <- model1$frail
frailty.summary <- quantile(frailty)

h.min <- h0*exp(-2*sqrt(frailty.var))
h.q1 <- h0*exp(-sqrt(frailty.var))
h.q3 <- h0*exp(sqrt(frailty.var))
h.max <- h0*exp(2*sqrt(frailty.var))

H.min <- cumsum(h.min)
H.q1 <- cumsum(h.q1)
H.q3 <- cumsum(h.q3)
H.max <- cumsum(h.max)

S0 <- exp(-H0)
S.min <- exp(-H.min)
S.q1 <- exp(-H.q1)
S.q3 <- exp(-H.q3)
S.max <- exp(-H.max)

par(mfrow=c(1,2))

plot(surv.time,h.max,type="l",xlab="Time (days)",ylab="Hazard function",
  lty=1,col="red")
lines(surv.time,h.q3,type="l",lty=2,col="darkorange")
lines(surv.time,h0,type="l",lty=3,col="black")
lines(surv.time,h.q1,type="l",lty=4,col="darkgreen")
lines(surv.time,h.min,type="l",lty=5,col="blue")
title("Variation in hospital-specific hazard function")

legend("topright",
  legend = c("Hazard 2 SD above mean","Hazard 1 SD above mean",
   "Average hospital","Hazard 1 SD below mean","Hazard 2 SD below mean"),
  lty = 1:5,
  col = c("red","darkorange","black","darkgreen","blue"),
  bty = "n")

plot(surv.time,S.max,type="l",lty=1,col="red",
  ylim=c(min(S.max),max(S.min)),
  xlab="Time (days)",ylab="Survival probability")
lines(surv.time,S.q3,type="l",lty=2,col="darkorange")
lines(surv.time,S0,type="l",lty=3,col="black")
lines(surv.time,S.q1,type="l",lty=4,col="darkgreen")
lines(surv.time,S.min,type="l",lty=5,col="blue")
title("Variation in hospital-specific survival function")

legend("topright",
  legend = c("Hazard 2 SD above mean","Hazard 1 SD above mean",
   "Average hospital","Hazard 1 SD below mean","Hazard 2 SD below mean"),
  lty = 1:5,
  col = c("red","darkorange","black","darkgreen","blue"),
  bty = "n")

mtext("Figure 1. Variation in hospital-specific hazards and survival (frailty model)",
  side=3,outer=T,line=-1.25,font=2)
