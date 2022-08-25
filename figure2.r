# This code COMES WITH ABSOLUTELY NO WARRANTY and is provided for
# illustrative purposes only.

library(lme4)

postscript("figure2.ps",horizontal=T,paper="letter")

################################################################################
# Analyses where intervals were chosen using clinical judgment about hazards
# being approximately equal.
################################################################################

zlist <- list(time=0,event=0,age=0,female=0,arf=0,carddys=0,chf=0,crf=0,cvd=0,
  diabcomp=0,malig=0,pulmoned=0,shock=0,inst.volume=0,teaching.hosp=0,
  hosp.revasc=0,hosp.cath=0,interval="",logtime=0,inst="")

xdat <- data.frame(scan("r.level2.clin.dat",zlist))

################################################################################
# Poisson model with hospital random intercepts.
# This is equivalent to fitting a piecewise exponential survival model.
################################################################################

model3 <- glmer(event ~ -1 + age + female + arf + carddys + chf + crf + cvd +
  diabcomp + malig + pulmoned + shock + inst.volume + teaching.hosp +
  hosp.revasc + hosp.cath + interval + (1|inst),family=poisson,
  offset=logtime,data=xdat)
model3

log.hazard <- fixef(model3)[16:20]
hazard <- exp(log.hazard)
x.cut <- c(2,5,10,20)

var.re <- VarCorr(model3)$inst[[1]]
sd.re <- sqrt(var.re)

hazard.low2SD <- hazard*exp(-2*sd.re)
hazard.low1SD <- hazard*exp(-1*sd.re)
hazard.high2SD <- hazard*exp(2*sd.re)
hazard.high1SD <- hazard*exp(1*sd.re)

hf0 <- stepfun(x.cut,hazard,f=1,right=T)
hf1 <- stepfun(x.cut,hazard.low2SD,f=1,right=T)
hf2 <- stepfun(x.cut,hazard.low1SD,f=1,right=T)
hf3 <- stepfun(x.cut,hazard.high2SD,f=1,right=T)
hf4 <- stepfun(x.cut,hazard.high1SD,f=1,right=T)

plot(hf0,xlim=c(0,30),xlab="Time (days)",ylab="Hazard",lty=3,col="black",
  ylim=c(min(hazard.low2SD),max(hazard.high2SD)),
  main="Figure 2. Variation in hazard functions across hospitals (PWE model)")

lines(hf3,lty=1,col="red")
lines(hf4,lty=2,col="darkorange")
lines(hf2,lty=4,col="darkgreen")
lines(hf1,lty=5,col="blue")

legend("topright",
  legend = c("Hazard 2 SD above mean","Hazard 1 SD above mean",
   "Average hospital","Hazard 1 SD below mean","Hazard 2 SD below mean"),
  lty = 1:5,
  col = c("red","darkorange","black","darkgreen","blue"),
  bty = "n")

