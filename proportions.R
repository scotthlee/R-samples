###################################################################
#code for working with binomially-distributed data
###################################################################

#score test for the single proportion
score.test<-  function(k, n, null=.5){
  return(((k / n) - null) / sqrt(null * (1 - null) / n))
}

#wald test for single proportion
wald.test <-  function(n, k, null=.5){
  p <-  n / k
  return((p - null) / sqrt(p * (1 - p) / k))
}

#sample size to detect a difference in paired proportions
connor.size <-  function(s1, s2, d=.1, Za=1.96, Zb=.842){
  psi <-  (s1 + s2) - 2 * s1 * s2
  (Za * sqrt(psi) + Zb * sqrt(psi - d^2))^2 / (d^2)
}

#score test for differences in proportions, from Agresti et al. 2008
diff.score.test <-  function(k1, k2, n1, n2, th=0){
  p1 <-  k1 / n1
  p2 <-  k2 / n2
  pt <-  (k1 + k2) / (n1 + n2)
  se <-  sqrt(pt * (1 - pt) * (1 / n1 + 1 / n2))
  z <-  (p1 - p2) / se
  return(list(z=z, chi=z^2, pval=1 - pchisq(z^2, 1)))
}

#sums the chi scores from the diff.score.test 
score.sum.test <- function(k1, k2, n1, n2, th=0){
  chis <- diff.score.test(k1, k2, n1, n2, th)$chi
  sum <- sum(chis)
  pval <- 1 - pchisq(sum, length(chis))
  return(data.frame(df=length(chis), sum=sum, p.val=pval))
}

#sample size for score confidence interval
score.size <-  function(p, se=1, sp=1, a=.95, b=.8, null=.5, t=2){
  za <-  qtukey(a, t, Inf) / sqrt(2)
  zb <-  qnorm(b)
  z <-  za + zb
  return(z^2 * null * (1 - null) / (p - null)^2)
}

score.power <-  function(p, n, se=1, sp=1, a=.95, null=.5, t=2){
  za <-  qtukey(a, t, Inf) / sqrt(2)
  return(pnorm(score.test(p * n, n) - za))
}

#sample size calculation based on the score test, assuming
#population sizes are equal (n1=n2)
diff.score.size <-  function(p1, p2, se=1, sp=1, a=.95, b=.8, th=0){
  za <-  qnorm(1 - ((1 - a) / 2))
  zb <-  qnorm(b)
  z <-  za + zb
  p1 <-  p1 * se + (1 - sp) * (1 - p1)
  p2 <-  p2 * se + (1 - sp) * (1 - p2)
  pt <-  (p1+p2) / 2
  return(2 * z^2 * pt * (1 - pt) / (p1 - p2)^2)
}

#wald test for differences in proportions
diff.wald.test <-  function(n1, n2, k, th=0){
  p1 <-  n1 / k
  p2 <-  n2 / k
  return((p1 - p2 - th) / sqrt((p1 * (1 - p1) / k)+(p2 * (1 - p2) / k)))
}

#sample size calculation for Wald differences; may or may not be correct
diff.wald.size <-  function(p1, p2, se=1, sp=1, a=.95, b=.8){
  za <-  qnorm(1 - ((1 - a) / 2))
  zb <-  qnorm(b)
  z <-  za + zb
  p1 <-  p1 * se + (1 - sp) * (1 - p1)
  p2 <-  p2 * se + (1 - sp) * (1 - p2)
  return((z^2 * (p1 * (1 - p1)+p2 * (1 - p2)) / (p1 - p2)^2) / se)
}

################################################################################
#confidence intervals for single proportions and single pairwise comparisons
prop.sco.ci <-  function(num, den, conf=0.95, cont.corr=FALSE) {
  cc <-  as.numeric(cont.corr)
  z <-  qnorm((1+conf) / 2)
  zsq <-  z^2
  lower <-  pmax(0, (2 * num + zsq - cc - z * sqrt(zsq - 2 * cc  - cc / den + 4 * (num / den) * (den * (1 - num / den)+cc))) / (2 * (den+zsq)))
  upper <-  pmin(1, (2 * num + zsq + cc + z * sqrt(zsq + 2 * cc  - cc / den + 4 * (num / den) * (den * (1 - num / den) - cc))) / (2 * (den+zsq)))
  data.frame(num=num, den=den, conf=conf, cont.corr=cont.corr, prop=num / den, lower=lower, upper=upper)
}

prop.bin.ci <-  function(num, den, conf=0.95, cont.corr=FALSE) {
  cc <-  as.numeric(cont.corr)
  z <-  qnorm((1+conf) / 2)
  prop <-  num / den
  se <-  sqrt(prop * (1 - prop) / den + 0.5 * cc / den)
  data.frame(num=num, den=den, conf=conf, cont.corr=cont.corr, prop=prop, lower=prop - z * se, upper=prop+z * se)
}

prob.poi.ci <-  function(num, den, conf=0.95, cont.corr=FALSE) {
  cc <-  as.numeric(cont.corr)
  z <-  qnorm((1+conf) / 2)
  prop <-  num / den
  se <-  sqrt(prop / den + 0.5 * cc / den)
  data.frame(num=num, den=den, conf=conf, cont.corr=cont.corr, prop=prop, lower=prop - z * se, upper=prop+z * se)
}

prop.diff.ci.old <-  function(prop1, lower1, upper1, prop2, lower2, upper2) {
  prop.diff <-  prop2 - prop1
  lower.diff <-  prop.diff - sqrt((prop1 - lower1)^2 + (upper2 - prop2)^2)
  upper.diff <-  prop.diff + sqrt((upper1 - prop1)^2 + (prop2 - lower2)^2)
  data.frame(prop.diff=prop.diff, lower.diff=lower.diff, upper.diff=upper.diff)
}

prop.dd.wald.ci <-  function(prop1, lower1, upper1, prop2, lower2, upper2) {
  prop.diff <-  prop2 - prop1
  lower.diff <-  prop.diff - sqrt((prop1 - lower1)^2 + (upper2 - prop2)^2)
  upper.diff <-  prop.diff + sqrt((upper1 - prop1)^2 + (prop2 - lower2)^2)
  data.frame(prop.diff=prop.diff, lower.diff=lower.diff, upper.diff=upper.diff)
}

#formula for differences of proportions taken direct from Newcombe's paper;
#not continuity - corrected
prop.diff.ci <-  function(k1, k2, m, n, conf=.95){
  z <-  qnorm((1+conf) / 2)
  th <-  k1 / m - k2 / n
  a <-  k1
  c <-  m - a
  b <-  k2
  d <-  n - b
  l1 <-  ( -sqrt( -4 * a^2 * m * z^2 + 4 * a * m^2 * z^2 + m^2 * z^4) + 2 * a * m + m * z^2) / (2 * (m^2 + m * z^2))
  u1 <-  (sqrt( -4 * a^2 * m * z^2 + 4 * a * m^2 * z^2 + m^2 * z^4) + 2 * a * m + m * z^2) / (2 * (m^2 + m * z^2))
  l2 <-  ( - sqrt( -4 * b^2 * n * z^2 + 4 * b * n^2 * z^2 + n^2 * z^4) + 2 * b * n + n * z^2) / (2 * (n^2 + n * z^2))
  u2 <-  (sqrt( -4 * b^2 * n * z^2 + 4 * b * n^2 * z^2 + n^2 * z^4) + 2 * b * n + n * z^2) / (2 * (n^2 + n * z^2))
  delta <-  z * sqrt(l1 * (1 - l1) / m + u2 * (1 - u2) / n)
  eps <-  z * sqrt(u1 * (1 - u1) / m + l2 * (1 - l2) / n)
  data.frame(prop.diff=th, lower.diff=th-delta, upper.diff=th+eps)
}

paired.wald.ci <-  function(pairs, alpha=.05){
  n <-  nrow(pairs)
  z <-  qnorm(1 - (alpha / 2))
  e <-  sum(pairs[, 1]==1 & pairs[, 2]==1)
  f <-  sum(pairs[, 1]==1 & pairs[, 2]==0)
  g <-  sum(pairs[, 1]==0 & pairs[, 2]==1)
  h <-  sum(pairs[, 1]==0 & pairs[, 2]==0)
  theta <-  (f - g) / n
  se <-  sqrt(f+g - (f - g)^2 / n) / n
  return(data.frame(lower=theta - (z * se), upper=theta+(z * se)))
}

gen.pairs <-  function(total.pos, f.n, g.n, n){
  
  #making rows for the discordant pairs
  f.val <-  matrix(rep(c(0, 1), f.n), nrow=f.n, ncol=2, byrow=T)
  g.val <-  matrix(rep(c(1, 0), g.n), nrow=g.n, ncol=2,  byrow=T)
  
  #figuring out the e - cell based on the discordant pairs and the total n pos on either
  e.n <-  (total.pos * 2) - (f.n + g.n)
  e.val <-  matrix(rep(c(1, 1), e.n), nrow=e.n, ncol=2, byrow=T)
  
  #binding the rows for the e, f, and g cells
  some.pos <-  rbind(e.val, f.val, g.val)
  
  #figuring out how many rows are left
  h.n <-  n - nrow(some.pos)
  
  #binding all the rows together to make a n x 2 table of pairs
  some.pos <-  rbind(some.pos, matrix(rep(0, 2 * h.n), nrow=h.n, ncol=2, byrow=T))
  
  return(some.pos)
}

#Wilson score - based method without continuity correction
paired.sco.ci <-  function(pairs, alpha=.05, neg=2){
  z <-  qnorm(1 - (alpha / 2))
  e <-  as.numeric(sum(pairs[, 1]==1 & pairs[, 2]==1))
  f <-  as.numeric(sum(pairs[, 1]==1 & pairs[, 2]==0))
  g <-  as.numeric(sum(pairs[, 1]==0 & pairs[, 2]==1))
  h <-  as.numeric(sum(pairs[, 1]==0 & pairs[, 2]==0))
  n <-  sum(e, f, g, h)
  theta <-  (f - g) / n
  l2 <-  ( - sqrt(( - 2 * f * n - 2 * e * n - n * z^2)^2 - 4 * (e^2 + 2 * e * f + f^2) * (n^2 + n * z^2)) + 2 * e * n + 2 * f * n + n * z^2) / (2 * (n^2 + n * z^2))
  u2 <-  (sqrt(( - 2 * f * n - 2 * e * n - n * z^2)^2 - 4 * (e^2 + 2 * e * f + f^2) * (n^2 + n * z^2)) + 2 * e * n + 2 * f * n + n * z^2) / (2 * (n^2 + n * z^2))
  dl2 <-  ((e+f) / n) - l2
  du2 <-  u2 - ((e+f) / n)
  l3 <-  ( - sqrt(( - 2 * g * n - 2 * e * n - n * z^2)^2 - 4 * (e^2 + 2 * e * g + g^2) * (n^2 + n * z^2)) + 2 * e * n + 2 * g * n + n * z^2) / (2 * (n^2 + n * z^2))
  u3 <-  (sqrt(( - 2 * g * n - 2 * e * n - n * z^2)^2 - 4 * (e^2 + 2 * e * g + g^2) * (n^2 + n * z^2)) + 2 * e * n + 2 * g * n + n * z^2) / (2 * (n^2 + n * z^2))
  dl3 <-  ((e+g) / n) - l3
  du3 <-  u3 - ((e+g) / n)
  phi.denom <-  sqrt((e+f) * (g+h) * (e+g) * (f+h))
  if (phi.denom==0){
    phi <-  0
  }
  else {
    phi <-  (e * h - f * g) / phi.denom
  }
  delta <-  sqrt(dl2^2 - (2 * phi * dl2 * du3) + du3^2)
  epsilon <-  sqrt(du2^2 - (2 * phi * du2 * dl3) + dl3^2)
  return(data.frame(diff=theta, lower=theta - delta, upper=theta+epsilon))
}

#score CIs for paired data; takes a 2x2 contingency table as input
paired.tab.sco.ci <-  function(tab, alpha=.05){
  n <-  sum(tab)
  z <-  qnorm(1 - (alpha / 2))
  e <-  tab[1, 1]
  f <-  tab[1, 2]
  g <-  tab[2, 1]
  h <-  tab[2, 2]
  theta <-  (f - g) / n
  l2 <-  ( - sqrt(( - 2 * f * n - 2 * e * n - n * z^2)^2 - 4 * (e^2 + 2 * e * f + f^2) * (n^2 + n * z^2)) + 2 * e * n + 2 * f * n + n * z^2) / (2 * (n^2 + n * z^2))
  u2 <-  (sqrt(( - 2 * f * n - 2 * e * n - n * z^2)^2 - 4 * (e^2 + 2 * e * f + f^2) * (n^2 + n * z^2)) + 2 * e * n + 2 * f * n + n * z^2) / (2 * (n^2 + n * z^2))
  dl2 <-  ((e+f) / n) - l2
  du2 <-  u2 - ((e+f) / n)
  l3 <-  ( - sqrt(( - 2 * g * n - 2 * e * n - n * z^2)^2 - 4 * (e^2 + 2 * e * g + g^2) * (n^2 + n * z^2)) + 2 * e * n + 2 * g * n + n * z^2) / (2 * (n^2 + n * z^2))
  u3 <-  (sqrt(( - 2 * g * n - 2 * e * n - n * z^2)^2 - 4 * (e^2 + 2 * e * g + g^2) * (n^2 + n * z^2)) + 2 * e * n + 2 * g * n + n * z^2) / (2 * (n^2 + n * z^2))
  dl3 <-  ((e+g) / n) - l3
  du3 <-  u3 - ((e+g) / n)
  phi.denom <-  sqrt((e+f) * (g+h) * (e+g) * (f+h))
  if (phi.denom==0){
    phi <-  0
  }
  else {
    phi <-  (e * h - f * g) / phi.denom
  }
  delta <-  sqrt(dl2^2 - (2 * phi * dl2 * du3) + du3^2)
  epsilon <-  sqrt(du2^2 - (2 * phi * du2 * dl3) + dl3^2)
  return(data.frame(diff=theta, lower=theta - delta, upper=theta+epsilon))
}

#Wilson score - based method with continuity - corrected phi
paired.sco.cc.ci <-  function(pairs, alpha=.05){
  z <-  abs(qnorm(1 - (alpha / 2)))
  e <-  as.numeric(sum(pairs[, 1]==1 & pairs[, 2]==1))
  f <-  as.numeric(sum(pairs[, 1]==1 & pairs[, 2]==0))
  g <-  as.numeric(sum(pairs[, 1]==0 & pairs[, 2]==1))
  h <-  as.numeric(sum(pairs[, 1]==0 & pairs[, 2]==0))
  n <-  sum(e, f, g, h)
  theta <-  (f - g) / n
  l2 <-  ( - sqrt(( - 2 * f * n - 2 * e * n - n * z^2)^2 - 4 * (e^2 + 2 * e * f + f^2) * (n^2 + n * z^2)) + 2 * e * n + 2 * f * n + n * z^2) / (2 * (n^2 + n * z^2))
  u2 <-  (sqrt(( - 2 * f * n - 2 * e * n - n * z^2)^2 - 4 * (e^2 + 2 * e * f + f^2) * (n^2 + n * z^2)) + 2 * e * n + 2 * f * n + n * z^2) / (2 * (n^2 + n * z^2))
  dl2 <-  ((e+f) / n) - l2
  du2 <-  u2 - ((e+f) / n)
  l3 <-  ( - sqrt(( - 2 * g * n - 2 * e * n - n * z^2)^2 - 4 * (e^2 + 2 * e * g + g^2) * (n^2 + n * z^2)) + 2 * e * n + 2 * g * n + n * z^2) / (2 * (n^2 + n * z^2))
  u3 <-  (sqrt(( - 2 * g * n - 2 * e * n - n * z^2)^2 - 4 * (e^2 + 2 * e * g + g^2) * (n^2 + n * z^2)) + 2 * e * n + 2 * g * n + n * z^2) / (2 * (n^2 + n * z^2))
  dl3 <-  ((e+g) / n) - l3
  du3 <-  u3  - ((e+g) / n)
  if (e * h > f * g){
    phi <-  max(0, (e * h - f * g - n / 2)) / sqrt((e+f) * (g+h) * (e+g) * (f+h))
  }
  else {
    phi <-  (e * h - f * g) / sqrt((e+f) * (g+h) * (e+g) * (f+h))
  }
  delta <-  sqrt(dl2^2 - (2 * phi * dl2 * du3) + du3^2)
  epsilon <-  sqrt(du2^2 - (2 * phi * du2 * dl3) + dl3^2)
  return(data.frame(lower=theta - delta, upper=theta+epsilon))
}

#function for generating pairs of data based on the e, f, and g cells
gen.pairs <-  function(e.n, f.n, g.n, n){
  
  #making rows for the discordant pairs
  f.val <-  matrix(rep(c(0, 1), f.n), nrow=f.n, ncol=2, byrow=T)
  g.val <-  matrix(rep(c(1, 0), g.n), nrow=g.n, ncol=2,  byrow=T)
  e.val <-  matrix(rep(c(1, 1), e.n), nrow=e.n, ncol=2, byrow=T)
  
  #binding the rows for the e, f, and g cells
  some.pos <-  rbind(e.val, f.val, g.val)
  
  #figuring out how many rows are left
  h.n <-  n - nrow(some.pos)
  
  #binding all the rows together to make a n x 2 table of pairs
  pairs <-  rbind(some.pos, matrix(rep(0, 2 * h.n), nrow=h.n, ncol=2, byrow=T))
  
  return(pairs)
}
################################################################################
#Simulation code using the functions above

pair.sim <-  function(n){
  out <-  data.frame(matrix(NA, ncol=n+1, nrow=n+1))
  colnames(out) <-  0:n
  
  #running the loop, with i being the # of concordant positive calls
  for (i in 0:n){
    
    max.discord <-  n - i
    f.vec <-  0:max.discord
    g.vec <-  max.discord:0
  
    ci.mat <-  matrix(NA, nrow=max.discord + 1, ncol=1)
    rownames(ci.mat) <-  f.vec
    colnames(ci.mat) <-  i
  
    for (j in 0:max.discord + 1){
      ci <-  round(paired.sco.cc.ci(gen.pairs(i, f.vec[j], 0, n)), 3) 
      ci.mat[j, ] <-  paste0(ci[1], ', ', ci[2])
    }

   out[, i + 1] <-  c(ci.mat, rep(NA, (n+1) - nrow(ci.mat)))
    
  }
  
  rownames(out) <-  0:32
  return(out)
}


################################################################################
#confidence intervals for multiple comparisons of proportions, using the
#studentized range distribution instead of the standard normal

#this function assumes n to be equal for all groups
prop.diff.multi <-  function(k1, k2, n1, n2, t=2, conf=.95, cont.corr=F){
  a <-  1 - (1 - pnorm(qtukey(conf, t, Inf) / sqrt(2))) * 2
  return(prop.diff.ci(k1, k2, n1, n2, conf=a))
}

# prop.ci.multi <- function(num, den, t=2, conf=.95, cont.corr=F){
#   a <-  1 - (1 - pnorm(qtukey(conf, t, Inf) / sqrt(2))) * 2
#   return(prop.sco.ci(num, den, conf=a))
# }

paired.diff.multi <-  function(pairs, t=2, conf=.95){
  a <-  (1 - pnorm(qtukey(conf, t, Inf) / sqrt(2))) * 2
  return(paired.sco.ci(pairs, alpha=a))
}

################################################################################
#code for working with multinomial data
################################################################################

#calculating observed proportions positive for multivariate categorical data
prop.pos <-  function(p, n, se, sp){
 
  p <-  p / sum(p) 
  x <-  round(n * p) 
  k <-  length(x) 
  mat <-  matrix(0, ncol=k, nrow=k)
  pos <-  round(x * se)
  diag(mat) <-  pos
  fpr <-  1 - sp
  
  #filling in the false - positive counts
  for (i in 1:k){
    nfp <-  round(fpr[i] * x[i])
    for(j in 1:nfp){
      row <-  sample(c(1:k)[ - i], 1)
      mat[row, i] <-  mat[row, i] + 1
    }
  }
  
  mat <-  rbind(round(mat), numeric(k))
  
  #filling in the "other" category at the bottom of the table
  for (i in 1:k){
    mat[k+1, i] <-  x[i] - sum(mat[, i])
  }
  
  counts <-  c(rowSums(mat))
  
  return(list(se=se, fpr=fpr, expected=x, observed=counts, tp=pos, props=counts / sum(x), 
              cross.tab=as.table(mat)))
  
}

#returns sample size for multinomial proportions using two statistics:
#Thompson's (1987) normal approximation, and Pearson's chi - squared
multinom.size <-  function(p, a=.95, b=.8, d=.05){
  
  #code based on Thompson's method
  m <-  length(p)
  alpha <-  (1 - a)
  z <-  qnorm(1 -  alpha / (2 * m))
  th.n <-  z^2 * (1 / m) * (1 - (1 / m)) / d^2
  
  #code based on Pearson's chi - squared statistic
  p0 <-  1 / length(p)
  p <-  p / sum(p)
  df <-  length(p) - 1
  
  xa <-  qchisq(a, m - 1)
  xb <-  qchisq(b, m - 1)
  x <-  xa + xb
  prs.n <-  x / sum(m * (p - p0)^2)
  
  return(list(thompson=th.n, pearson=prs.n))
}


#simulates draws from a multinomial distribution and returns summary stats, 
#including chi - square p - values and alpha - level quantile CIs
multinom.sim <-  function(p, n, nsim, a=.05, alt=FALSE){
  
  #setting null and alternative probabilities
  p <-  p / sum(p)
  p0 <-  rep(1 / length(p), length(p))
  if (alt){
    p0 <-  p
  }
  
  #generating random samples from the alternative distribution
  stats <-  matrix(nrow=nsim)
  observed <-  t(rmultinom(nsim, n, p))
  props <-  observed / n
  quants <-  t(apply(props, quantile, MARGIN=2, probs=c(a / 2, 1 - (a / 2))))
  
  return(list(counts=observed, props=props, quantiles=quants))
}

multinom.loop <-  function(p, nvec, nsim, a=.05, alt=FALSE){
  iqr <-  matrix(NA, nrow=length(p), ncol=length(nvec))
  quantiles <-  matrix(NA, nrow=length(p), ncol=2 * length(nvec))
  
  for (i in 1:length(nvec)){
    run <-  multinom.sim(p, nvec[i], nsim, a, alt)
    quantiles[, 2 * i - 1] <-  run$quantiles[, 1]
    quantiles[, 2 * i] <-  run$quantiles[, 2]
    iqr[, i] <-  run$quantiles[, 2] - run$quantiles[, 1]
  }
  
  return(list(n=nvec, quantiles=quantiles, iqr=iqr, error=iqr / 2))
}

#gets permutation-based p-value for sum of scores (permutation omnibus for difference in proportions)
score.sum.permtest <- function(k1, k2, n1, n2, nset=1000){
  
  obs.scores <- diff.score.test(k1, k2, n1, n2)$chi
  props <- data.frame(p1=k1/n1, p2=k2/n2)
  stats <- matrix(nrow=length(k1), ncol=nset)
  stat.sums <- matrix(nrow=length(k1), ncol=1)
  
  for (i in 1:nset){
    
    data <- data.frame(k1=rbinom(length(k1), n1, props$p1),
                       k2=rbinom(length(k1), n2, props$p1))
    
    stats[, i] <- diff.score.test(data$k1, data$k2, n1, n2)$chi
    stat.sums[i] <- sum(stats[, i])
    
  }
  
  p.val <- sum(stat.sums >= sum(obs.scores))
  
  return(list(data=data.frame(k1, k2, n1, n2), props=props, obs.scores=obs.scores,
              stats=stats, stat.sums=stat.sums, p.val=p.val))
}

