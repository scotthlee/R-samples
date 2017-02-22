
#set the working directory and load the WHO height/weight for age scripts
source('WHO/igrowup_restricted.R')
source('WHO/igrowup_standard.R')

#reading in vitals and screening tables
vitals <- read.csv('vitals.csv')
screening <- read.csv('screening.csv')

#import the lookup tables from the WHO
who.ht <- as.data.frame(read.table('WHO/lenanthro1.txt', col.names=c('sex', 'age', 'l', 'm', 's', 'loh'), 
                     colClasses=c(rep('numeric', 5), 'character')))
who.wt <- as.data.frame(read.table('WHO/weianthro1.txt', col.names=c('sex', 'age', 'l', 'm', 's'), 
                                   colClasses=c(rep('numeric', 5))))
who.wfh <- as.data.frame(read.table('WHO/wfhanthro1.txt', col.names=c('sex', 'height', 'l', 'm', 's', 'lorh'), 
                                colClasses=c(rep('numeric', 5), 'character')))
who.wfl <- as.data.frame(read.table('WHO/wflanthro1.txt', col.names=c('sex', 'length', 'l', 'm', 's', 'lorh'), 
                                colClasses=c(rep('numeric', 5), 'character')))
who.ac <- as.data.frame(read.table('WHO/acanthro1.txt', col.names=c('sex', 'age', 'l', 'm', 's'), 
                                   colClasses=rep('numeric', 5)))

#function that retrieves height and weight LMS for a given age
lms <- function(sex, val, table){
  return(cbind(l=table[table[, 2]==val & table[, 1]==sex, 'l'], 
               m=table[table[, 2]==val & table[, 1]==sex, 'm'], 
               s=table[table[, 2]==val & table[, 1]==sex, 's']))
}

#generic function that calculates z-score based on LMS lookup
z <- function(sex, val1, val2, table){
  v.lms <- lms(sex, val2, table)
  v.z <- (val1-v.lms[2])/(v.lms[2] * v.lms[3])
  v.sd3pos <- v.lms[2] * ((1  +  v.lms[1] * v.lms[3] * 3)^(1/v.lms[1]))
  v.sd2pos <- v.lms[2] * ((1  +  v.lms[1] * v.lms[3] * 2)^(1/v.lms[1]))
  v.sd3neg <- v.lms[2] * ((1  +  v.lms[1] * v.lms[3] * (-3))^(1/v.lms[1]))
  v.sd2neg <- v.lms[2] * ((1  +  v.lms[1] * v.lms[3] * (-2))^(1/v.lms[1]))
  if(v.z > 3){
    v.z <- 3  +  ((val1 - v.sd3pos)/(v.sd3pos-v.sd2pos))
  }
  if(v.z < -3){
    v.z <- -3  +  ((val1 - v.sd3neg)/(v.sd2neg-v.sd3neg))
  }
  else{
    v.z <- v.z
  }
  return(v.z)
}

#specific function that takes anthropometric data and does the z-scores
all.z <-Vectorize(function(sex, age, ht, wt, ac){
  h.lms <- lms(sex, age, who.ht)
  hfa <- (ht-h.lms[2])/(h.lms[3] * h.lms[2])
  wfa <- z(sex, wt, age, who.wt)
  if(ht >= min(who.wfh$height)){
    wfh <- z(sex, wt, round(ht, 1), who.wfh)
  }
  else{
    wfh <- z(sex, wt, round(ht, 1), who.wfl)
  }
  if (age >= 91){
    ac.lms <- lms(sex, age, who.ac)
    ac.z <- (ac-ac.lms[2])/(ac.lms[3] * ac.lms[2])
  }
  else {
    ac.z <- NA
  }
  return(data.frame(hfa.z=hfa, hfa.cent=pnorm(hfa) * 100, 
                    wfa.z=wfa, wfa.cent=pnorm(wfa) * 100, 
                    wfh.z=wfh, wfh.cent=pnorm(wfh) * 100, 
                    ac.z=ac.z, ac.cent=pnorm(ac.z) * 100))
})


#pulling out study ID,  sex,  and DOB for the merge
age <- data.frame(StudyID=screening$StudyID, sex=screening$Sex, DOB=screening$DOB)

#converting ages to the same format
age$DOB <- as.Date(strptime(as.character(age$DOB),  "%d-%b-%Y"))
vitals$DOV <- as.Date(strptime(vitals$DateOfVitals,  "%d-%b-%Y"))

#merges tables to get vitals plus age
vitals <- merge(vitals, age, by='StudyID', all.x=F)
age.vit <- as.numeric(vitals$DOV-vitals$DOB)
hw <- data.frame(StudyID=vitals$StudyID, sex=vitals$sex, age=age.vit, 
                  height=vitals$Height, weight=vitals$Weight, muac=vitals$MUAC, 
                 DOV=vitals$DOV, VisitType=vitals$VisitType)

#limiting table to kids that the WHO formula will process 
hw <- hw[hw$age>=0 & hw$age<=1856, ]

#running the kids through the WHO formula to get z-scores and centiles
hw <- data.frame(hw, t(all.z(hw$sex, hw$age, hw$height, hw$weight, hw$muac)))

#collapsing the lists with z-scores and centiles
hw$hfa.z <- as.numeric(do.call(rbind, hw$hfa.z))
hw$hfa.cent <- as.numeric(do.call(rbind, hw$hfa.cent))
hw$wfa.z <- as.numeric(do.call(rbind, hw$wfa.z))
hw$wfa.cent <- as.numeric(do.call(rbind, hw$wfa.cent))
hw$wfh.z <- as.numeric(do.call(rbind, hw$wfh.z))
hw$wfh.cent <- as.numeric(do.call(rbind, hw$wfh.cent))
hw$ac.z <- as.numeric(do.call(rbind, hw$ac.z))
hw$ac.cent <- as.numeric(do.call(rbind, hw$ac.cent))

#removing errors in MUAC
hw$ac.z[hw$ac.z > 20] <- NA

