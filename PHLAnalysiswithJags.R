### read Data 
#mainDir <- '/media/metienne_h/Consulting/DFO/2015-YelloweyeSurvey'
#mainDir <- '/home/metienne/bianca/Consulting/DFO/2015-YelloweyeSurvey/'
#mainDir <- '/home/metienne/Consulting/DFO/2015-YelloweyeSurvey/'
mainDir <- '/home/metienne/EnCours/2015-YelloweyeSurvey'
dataDir <- 'Data'
resDir <- 'Results/PHL'
codeDir <- 'DFOYellowEye'

library('rjags')
library('stringr')
library('ggplot2')
source(file.path(mainDir, codeDir,"glmFunction.R"))
source(file.path(mainDir, codeDir, 'formattingFunction.R'))

yearBreaks=c(1986,1989 , 1994, 1996, 2002, 2006,2016)
depthBreaks = c(0, 101, 201, 401, 801 )
nIter <- 10000

if(!file.exists(file.path(mainDir, resDir)))
  dir.create(file.path(mainDir, resDir))

### 
PHLData <- read.table(file.path(mainDir, dataDir, 'PHL_line_gear_sets_catch_counts_MPE_Feb15.csv'),
                      sep=",", header=T)


# PreProcessing of the Data -----------------------------------------------
### some graphics
hist(PHLData$Count_YE[PHLData$Count_YE>0], freq=FALSE)
mean(PHLData$Count_YE[PHLData$Count_YE>0])

## replace all na in YEYE_counts by 0
PHLData$NT <- PHLData$Count_YE
PHLData$NT[which(is.na(PHLData$NT))]=0

### Standardizing columns name and units
PHLData$Depth <- PHLData$BEST_DEPTH
PHLData$Year <- PHLData$YEAR

## exclude all depth >600 (m)
PHLData <- PHLData[ which(PHLData$Depth<=600),]
##dealing with unspecified gear
## following Lynn recommandations if Depth <50 HANDLINE, else 
PHLData$GEAR[PHLData$GEAR=="Unspecified line gear" & PHLData$Depth<50] <- "Handline"
PHLData$GEAR[PHLData$GEAR=="Unspecified line gear" & PHLData$Depth>=50] <- "Longline"
PHLData$GEAR <- droplevels(PHLData$GEAR)


#remove records 1985
table(PHLData$Year)
PHLData <- subset(PHLData, Year>=1986) 

## create classes for years and depth
PHLData$TimePeriod <- splitInClasses(PHLData$Year, yearBreaks)
PHLData$Fact1 <- splitInClasses(PHLData$Depth, depthBreaks) #depth
PHLData$Fact2 <- PHLData$FISHERY #fishery
PHLData$Fact3 <- PHLData$GEAR #gear

data <- PHLData[,c("Year","TimePeriod", "Fact1",  "Fact2", "Fact3", "NT")]
summary(PHLData)

# Running CPUE analysis -----------

Fact=names(data)[which(!names(data) %in% c("Year", "TimePeriod", "NT"))]

TimePeriod = sort(unique(data$TimePeriod))
res<-lapply(levels(TimePeriod), function(d){
  dataPeriod <- data[data$TimePeriod==d,]
  analysis <- glmAnalysis(dataPeriod, Fact, fig.path=file.path(mainDir, resDir), nIter=nIter, pref=paste0("PHL-",d))
  })  
resFinal <- Reduce(rbind, res)
write.table(x =resFinal,file = file.path(mainDir,resDir,"PHL-Indices.csv"), row.names=F)
