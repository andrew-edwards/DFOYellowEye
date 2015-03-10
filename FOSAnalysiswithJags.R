### read Data 
#mainDir <- '/media/metienne_h/Consulting/DFO/2015-YelloweyeSurvey'
#mainDir <- '/home/metienne/bianca/Consulting/DFO/2015-YelloweyeSurvey/'
#mainDir <- '/home/metienne/Consulting/DFO/2015-YelloweyeSurvey/'
mainDir <- '/home/metienne/EnCours/2015-YelloweyeSurvey'
dataDir <- 'Data'
resDir <- 'Results/FOS'
codeDir <- 'DFOYellowEye'


nIter <- 10000
library('rjags')
library('stringr')
library('scales')
library('ggplot2')

source(file.path(mainDir, codeDir,"glmFunction.R"))
yearBreaks=c(1986,1989 , 1994, 1996, 2002, 2006,2016)
depthBreaks = c(0, 101, 201, 401, 801 )
source(file.path(mainDir, codeDir, 'formattingFunction.R'))

if(!file.exists(file.path(mainDir, resDir)))
  dir.create(file.path(mainDir, resDir))
### 
FOSData <- read.table(file.path(mainDir, dataDir, 'FOS_non-trawl_sets_catch_countsMPE_FEB15.csv'),
                      sep=",", header=T, nrows=58171)

# PreProcessing of the Data -----------------------------------------------
### some graphics
hist(FOSData$YEYE_count[FOSData$YEYE_count>0], freq=FALSE)
mean(FOSData$YEYE_count[FOSData$YEYE_count>0])

## replace all na in YEYE_counts by 0
FOSData$YEYE_count[which(is.na(FOSData$YEYE_count))]=0

### Standardizing columns name and units
FOSData$Depth <- FOSData$BEST_DEPTH_FM*2
FOSData$Year <- FOSData$YEAR

## exclude all depth >600 (m)
FOSData <- FOSData[ which(FOSData$Depth<=600),]

#remove records for 2001, 2004 and 2005 and 2015
table(FOSData$YEAR)
FOSData <- subset(FOSData, YEAR>=2006 & YEAR < 2015) 
FOSData$NT <- FOSData$YEYE_count

## create classes for years and depth
FOSData$TimePeriod <- splitInClasses(FOSData$Year, yearBreaks)
FOSData$Fact1 <- splitInClasses(FOSData$Depth, depthBreaks) #depth
FOSData$Fact2 <- FOSData$FISHERY_SECTOR #fishery
#FOSData$Fact3 <- FOSData$GEAR #gear

data <- FOSData[,c("Year","TimePeriod", "Fact1", "Fact2", "NT")]


# Running CPUE analysis -----------

Fact=names(data)[which(!names(data) %in% c("Year", "TimePeriod", "NT"))]
##
TimePeriod = sort(unique(data$TimePeriod))
lapply(levels(TimePeriod), function(d){
  dataPeriod <- data[data$TimePeriod==d,]
  analysis <- glmAnalysis(dataPeriod, Fact, fig.path=file.path(mainDir, resDir), nIter=nIter, pref="FOS")
  })       

  


