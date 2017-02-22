#script for analyzing BCBS MA opioid data for MMWR (6/2016)

#loading the packages
library(geepack)
library(car)

#function for comparing slopes
slope.test.t <- function(mod1, mod2){
  
  b1 <- coef(mod1)[2]
  b2 <- coef(mod2)[2]
  se1 <- summary(mod1)$coefficients[2, 2]
  se2 <- summary(mod2)$coefficients[2, 2]
  
  t <- as.numeric((b1 - b2) / sqrt(se1^2 + se2^2))
  df <- round((se1^2 + se2^2)^2 / (se1^4 / mod1$df.residual + se2^4/mod2$df.residual))
  pval <- 1 - pt(t, df)
  
  return(list(stat=t, df=df, pval=pval))
} 

#function for adding the coefficients of a regression model;
#takes the model and a vector of coefficient positions as input 
add_coeffs <- function(mod, pos, a=.05){
  
  b1 <- coef(mod)[pos[1]]
  b2 <- coef(mod)[pos[2]]
  se1 <- summary(mod)$coefficients[pos[1], 2]
  se2 <- summary(mod)$coefficients[pos[2], 2]
  covb1b2 <- vcov(mod)[pos[1], pos[2]]
  
  terms <- paste0(names(coef(mod))[pos[1]], ' + ', names(coef(mod))[pos[2]])
  bnew <- as.numeric(b1 + b2)
  senew <- as.numeric(sqrt(se1^2 + se2^2 + (2 * covb1b2)))
  crit <- qt(1 - a/2, nrow(mod$data) - 2)
  lower <- bnew - crit * senew
  upper <- bnew + crit * senew
  return(data.frame(terms=terms, beta=bnew, se=senew, lower=lower, upper=upper))
}

#function for normalizing counts and rates and stuff
normalize <- function(data){
  return((data - mean(data)) / sd(data))
}

#more stuff
bcbs <- read.csv('bcbsma.csv')
bcbs$month <- seq(1:nrow(bcbs))
bcbs$post_month <- rep(0, nrow(bcbs))

bcbs$post_month[bcbs$period=='post'] <- 1:36

bcbs$mem_short_perc <- as.numeric(bcbs$mem_short_perc)
bcbs_pre <- bcbs[bcbs$period == 'pre', ]
bcbs_post <- bcbs[bcbs$period == 'post', ]

onc <- bcbs[grepl('onc', colnames(bcbs)) & !grepl('ns', colnames(bcbs)) & !grepl('comb', colnames(bcbs))]
non_onc <- bcbs[!grepl('onc', colnames(bcbs)) & !grepl('ns', colnames(bcbs)) & 
                  !grepl('comb', colnames(bcbs)) & !grepl('tram', colnames(bcbs))]

#glm for ITS
shmem <- glm(members_short ~ month + post_month + (period == 'post'), family=poisson, offset=log(total_mem), data=bcbs)
shrx <- glm(rx_short ~ month + post_month + (period == 'post'), family=poisson, offset=log(total_mem), data=bcbs)
shoncmem <- glm(members_onc_short ~ month + post_month + (period == 'post'), family=poisson, offset=log(total_mem), data=bcbs)
shoncrx <- glm(rx_onc_short ~ month + post_month + (period == 'post'), family=poisson, offset=log(total_mem), data=bcbs)

lmem <- glm(members_long ~ month + post_month + (period == 'post'), family=poisson, offset=log(total_mem), data=bcbs)
lrx <- glm(rx_long ~ month + post_month + (period == 'post'), family=poisson, offset=log(total_mem), data=bcbs)
loncmem <- glm(members_onc_long ~ month + post_month + (period == 'post'), family=poisson, offset=log(total_mem), data=bcbs)
loncrx <- glm(rx_onc_long ~ month + post_month + (period == 'post'), family=poisson, offset=log(total_mem), data=bcbs)

cbmem <- glm(members_comb ~ month + post_month + (period == 'post'), family=poisson, offset=log(total_mem), data=bcbs)
cbrx <- glm(rx_comb ~ month + post_month + (period == 'post'), family=poisson, offset=log(total_mem), data=bcbs)
cboncmem <- glm(members_onc_comb ~ month + post_month + (period == 'post'), family=poisson, offset=log(total_mem), data=bcbs)
cboncrx <- glm(rx_onc_comb ~ month + post_month + (period == 'post'), family=poisson, offset=log(total_mem), data=bcbs)

#boxing up the models for easy summarization
modnames <- c('shmem', 'shrx', 'shoncmem', 'shoncrex', 'lmem', 'lrx', 'loncmem', 'loncrx', 
                'cbmem', 'cbrx', 'cboncmem', 'cboncrex')
modlist <- list(shmem, shrx, shoncmem, shoncrx, lmem, lrx, loncmem, loncrx, cbmem, cbrx, cboncmem, cboncrx)
coefs <- as.data.frame(do.call(rbind, lapply(modlist, coef)))
cis <- do.call(rbind, lapply(modlist, confint))
pre_cis <- data.frame(type=modnames, cis[rownames(cis) == 'month', ])

postslopes <- do.call(rbind, lapply(modlist, add_coeffs, pos=c(2,3)))
diff_cis <- data.frame(lower=cis[rownames(cis) == 'post_month', 1], upper=cis[rownames(cis) == 'post_month', 2])

#checking for autocorrelation
corstats <- data.frame(med_type=modnames, do.call(rbind, lapply(modlist, dwt)))

#combining the coefficients and 95% CIs, for annual percent change
modstats <- data.frame(pre_slope=coefs$month, post_slope=postslopes$beta, slope_diff=coefs$post_month, 
                       diff_lower=diff_cis$lower, diff_upper=diff_cis$upper)
modstats <- data.frame(med_type=modnames, round(modstats * 100 * 12, 4))

#pulling out the rows for the annual percent-changes
onlypost <- modstats[seq(3, nrow(modstats), by=4), ]

#id vector for GEE
id <- rep(1, 48)
id2 <- c(rep(1, 12), rep(2, 36))

#testing out GEE
test <- geeglm(members_onc_short ~ month + post_month + (period == 'post'), family=poisson(link='log'), 
               offset=log(total_mem), data=bcbs, id=id, corstr=c('ar1'))

#function for making plots
norm.plot <- function(x, y, color){
  y <- scale(y)
  points(x, y, col=color, cex=.75)
  lines(smooth.spline(x, y), col=color)
}

