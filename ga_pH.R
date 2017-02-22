
#loading the required packages
require(nnet)
require(survival)
require(geepack)
require(multgee)

#import the patient-level lab results
mtball <- read.csv('data/rdata/ByPtNo_dot.csv')
mtball <- mtball[grepl('strng_xp', colnames(mtball)) | grepl('strng_mg', colnames(mtball)) |
                   grepl('gsasp_xp', colnames(mtball)) | grepl('gsasp_mg', colnames(mtball)) |
                   grepl('StudyID', colnames(mtball)) | grepl('anypos', colnames(mtball))]
one.mtb <- data.frame(StudyID=mtball$StudyID, mtball[grepl('1', colnames(mtball))])
two.mtb <- data.frame(StudyID=mtball$StudyID, mtball[grepl('2', colnames(mtball))])
colnames(one.mtb)[-1] <- c('ga_xp', 'str_xp', 'ga_mg', 'str_mg')
colnames(two.mtb)[-1] <- c('ga_xp', 'str_xp', 'ga_mg', 'str_mg')

pos_names <- mtball$StudyID[mtball$anypos == 1]

#import the lab sammple reception data
reception <- read.csv('data/rdata/TBLabSampleReception.csv')
reception <- reception[reception$SampleType == 1,]

##########################################################
#import the ga and xpert results data
ga.csv <- read.csv('data/rdata/GastricAspirate.csv')

#fixing bad feed dates
ga.csv$DateCollected <- as.character(ga.csv$DateCollected)
ga.csv$DateLastfeed <- as.character(ga.csv$DateLastfeed)
ga.csv$DateCollected[527] <- '7/3/2014'
ga.csv$DateLastfeed[459] <- '7/11/2015'
ga.csv$DateCollected[234] <- '9/30/2014'
ga.csv$DateCollected[500] <- '3/26/2015'
ga.csv$DateLastfeed[500] <- '3/25/2015'
ga.csv$DateCollected[543] <- '6/4/2015'
ga.csv$DateLastfeed[227] <- '9/24/2014'

#fixing up the times
CollDates <- as.character(ga.csv$DateCollected)
FeedDates <- as.character(ga.csv$DateLastfeed)

CollTimes <- as.character(ga.csv$TimeCollected)
CollTimes <- gsub(';', ':', CollTimes)
CollTimes <- gsub('03:', '3:', CollTimes)
CollTimes <- gsub('04:', '4:', CollTimes)
CollTimes <- gsub('05:', '5:', CollTimes)
CollTimes <- gsub('06:', '6:', CollTimes)

FeedTimes <- as.character(ga.csv$Timelastfeed)
FeedTimes <- gsub(';', ':', FeedTimes)
FeedTimes <- gsub('03:', '3:', FeedTimes)
FeedTimes <- gsub('04:', '4:', FeedTimes)
FeedTimes <- gsub('05:', '5:', FeedTimes)
FeedTimes <- gsub('06:', '6:', FeedTimes)

CollTimes <- strptime(paste(CollDates, CollTimes), format='%m/%d/%Y %H:%M')
FeedTimes <- strptime(paste(FeedDates, FeedTimes), format='%m/%d/%Y %H:%M')

#taking the difference between collection time and last feeding time
FastTimes <- difftime(CollTimes, FeedTimes, units='hours')

#making a new variable in the ga.csv table
ga.csv$FastTimeHours <- as.numeric(FastTimes)
ga.csv$FastTimeHours[ga.csv$FastTimeHours < 0] <- ga.csv$FastTimeHours[ga.csv$FastTimeHours < 0] + 24
#########################################

#divide ga results by sample number
one.clin <- ga.csv[ga.csv$SampleNo == 1, c(1, 2, 5, 13, 14, 30)]
two.clin <- ga.csv[ga.csv$SampleNo == 2, c(1, 2, 5, 13, 14, 30)]

#merging the data for each sample number
one <- merge(one.clin, one.mtb, by='StudyID', all.y=T)
two <- merge(two.clin, two.mtb, by='StudyID', all.y=T)

#putting things back together
all <- rbind(one, two) 
all <- all[order(all$StudyID),] #ordering by StudyID
all <- all[!is.na(all$SampleNo), ] #leaving out missing samples
all$PHafterBicarb[all$PHafterBicarb == 999] <- all$PHPriorBicarb[all$PHafterBicarb == 999]
colnames(all)[4:5] <- c('ph1', 'ph2')

#making categories for the ph levels
all$ph_cat <- 'normal'
all$ph_cat[all$ph1 < 7] <- 'low'
#all$ph_cat[all$ph1 >= 8] <- 'high'

#pulling out the tb-positives
all_pos <- all[all$StudyID %in% pos_names, ]

#removing children with missing ga culture results
all_pos <- all_pos[!is.na(all_pos$ga_mg), ]
all <- all[!is.na(all_pos$ga_mg), ]

#makes a table for each of the tests
test_tab <- apply(all_pos[,6:9], MARGIN=2, FUN=table, all_pos$ph_cat)

########################################################
#running some binary models
geepos <- geeglm(ga_mg==1 ~ ph_cat + VolumeMls + FastTimeHours, family=binomial, data=all_pos, id=StudyID) 
geecon <- geeglm(ga_mg==3 ~ ph_cat + VolumeMls, family=binomial, data=all, id=StudyID)

#running some multinomial models
multi_unadj <- multinom(relevel(as.factor(ga_mg), 2) ~ ph_cat, data=all_pos)
multimod <- multinom(relevel(as.factor(ga_mg), 2) ~ ph_cat + VolumeMls + FastTimeHours, data=all_pos)

#piecing together a multinomial gee
gee12 <- geeglm(ga_mg==1 ~ ph1 + VolumeMls + FastTimeHours, 
                family=binomial, id=StudyID, data=all_pos[all_pos$ga_mg %in% c(1, 2), ], corstr='exchangeable')
gee32 <- geeglm(ga_mg==3 ~ ph_cat + VolumeMls + FastTimeHours,
                family=binomial, id=StudyID, data=all_pos[all_pos$ga_mg %in% c(2, 3), ], corstr='exchangeable')

#######################################################
#test code for getting bootstrap quantile CIs
resamp_kids <- function(data){
  data$StudyID <- as.character(data$StudyID)
  old_ids <- unique(data$StudyID)
  n_kids <- length(unique(data$StudyID))
  new_ids <- sample(old_ids, n_kids, replace=TRUE)
  new_rows <- 0
  for (i in 1:length(new_ids)){
    new_rows <- new_rows + sum(data$StudyID == new_ids[i])
  }
  
  out <- data.frame(matrix(nrow=new_rows, ncol=ncol(data)))
  colnames(out) <- colnames(data)
  k <- 1
  for (i in 1:length(new_ids)){
    base <- data[data$StudyID == new_ids[i], ]
    for (j in 1:nrow(base)){
      out[k, ] <-  base[sample(1:nrow(base), 1, replace=TRUE), ]
      k <- k + 1
    }
  }
  return(out)
}

boot_loop <- function(data, n_sim=1000){
  coef <- matrix(nrow=n_sim, ncol=1)
  for (i in 1:n_sim){
    if (i == 100){
      print(i)
    }
    boot_set <- resamp_kids(data)
    coef[i] <- coef(multinom(ga_mg==1 ~ ph_cat + VolumeMls + FastTimeHours, data=boot_set))[2]
  }
  return(coef)
}

#making some tables
ph_tab <- table(pH=all_pos$ph_cat, mg=all_pos$ga_mg)
