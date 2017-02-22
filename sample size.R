
#importing historical data from the WHO
annual.tb <- read.csv('who tb.csv')
annual.tb <- annual.tb[annual.tb$country=='Kenya', ]
annual.tb <- data.frame(year=annual.tb$year, 
                        cases=rowSums(cbind(annual.tb$new_sp, annual.tb$new_sn, annual.tb$new_su, 
                                            annual.tb$new_ep, annual.tb$new_oth), na.rm=T))

#setting up the numerators and denominators
ke.pop <- 40e6
ke.tb <- 1e5 / ke.pop

kis.pop <- 1e6
kis.tb <- 3000 / kis.pop

cdr1 <- 0.75
cdr2 <- 0.90

#figuring the slope for Kenya
drop.years <- annual.tb[28:33, ]
drop.mod <- lm(cases / ke.pop ~ year, drop.years)
drop.slope <- coef(drop.mod)[2]
drop.se <- summary(drop.mod)$coefficients[2, 2]

#function for comparing slopes
slope.test.mod <- function(mod1, mod2){
  
  b1 <- coef(mod1)[2]
  b2 <- coef(mod2)[2]
  se1 <- summary(mod1)$coefficients[2, 2]
  se2 <- summary(mod2)$coefficients[2, 2]
  
  t <- as.numeric((b1 - b2) / sqrt(se1^2 + se2^2))
  df <- round((se1^2 + se2^2)^2 / (se1^4 / mod1$df + se2^4 /mod2$df))
  pval <- 1 - pt(t, df)
  
  return(list(stat=t, df=df, p.val=pval))
} 

#function for simulating the time series and estimating power
impact.sim <- function(pop=kis.pop, starting.count=3000, old.mod=drop.mod, 
                       slope.increase, nsim=1000, conf=.95, power.only=F){
  
  #calculating the intervention period slope based on the expected increase
  old.slope <- coef(old.mod)[2] * kis.pop
  new.slope <- old.slope + (slope.increase * old.slope)
  
  drop.vec <- 1:5 * new.slope
  p.vec <- (starting.count + drop.vec) / pop
  
  #simulating count data for 5 years
  sim.data <- data.frame(year=1:5, matrix(rbinom(nsim * 5, pop, p.vec), nrow=5) / pop)  
  
  #empty matrix to hold the stats from the t-test 
  slopes <- matrix(NA, nrow=1, ncol=numeric(nsim))
  stats <- matrix(NA, nrow=1, ncol=numeric(nsim))
  ses <- matrix(NA, nrow=1, ncol=numeric(nsim))
  pvals <- matrix(NA, nrow=1, ncol=numeric(nsim))
  
  #running the linear model loop
  for (i in 2:ncol(sim.data)){
    sim.mod <- lm(sim.data[, i] ~ sim.data$year)
    mod.test <- slope.test.mod(old.mod, sim.mod)
    stats[i] <- mod.test$stat
    pvals[i] <- mod.test$p.val
    slopes[i] <- coef(sim.mod)[2]
    ses[i] <- summary(sim.mod)$coefficients[2, 2]
  }
  
  if(!power.only){
    return(list(data=sim.data, slopes=slopes, ses=ses, stats=stats, pvals=pvals, power=sum(pvals <= 1 - conf, na.rm=T) / length(pvals)))   
  } else{
    return(power=sum(pvals <= 1 - conf, na.rm=T) / length(pvals))    
  }
}

slope.vec <- seq(.1, 2, by=.05)
# slope.vs.power <- mapply(impact.sim, slope.increase=slope.vec, power.only=T)

#plotting the lines
library(Hmisc, quietly=T)
plot(slope.vec, slope.vs.power, main=c('Slope vs. power'), xlab=c('Slope increase (proportion of original)'), 
     ylab=c('Power'))
lines(smooth.spline(slope.vec, slope.vs.power))
abline(h=.8, lty=3)
minor.tick(nx=5, ny=1, tick.ratio=.5)

#saving the predicted values from the spline
predicted.power <- predict(power.line, seq(.1, 2, by=.01))




