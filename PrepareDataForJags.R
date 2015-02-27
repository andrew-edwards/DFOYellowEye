### read Data 
mainDir <- '/media/metienne_h/Consulting/DFO/2015-YelloweyeSurvey'
dataDir <- 'Data'
resDir <- 'Results'
codeDir <- '2015YellowEye'

breaks=c(1986,1989 , 1994, 1996, 2002, 2006,2016)
### 
FOSData <- read.table(file.path(mainDir, dataDir, 'FOS_non-trawl_sets_catch_countsMPE_FEB15.csv'),
                      sep=",", header=T, nrows=58171)



# PreProcessing of the Data -----------------------------------------------
## replace all na in YEYE_counts by 0
FOSData$YEYE_count[which(is.na(FOSData$YEYE_count))]=0

## exclude all depth >300
FOSData <- FOSData[ which(FOSData$BEST_DEPTH_FM<=300),]
FOSData$AreaCat <- as.numeric(FOSData$AREA)




## create numeric time period
FOSData$TimePeriod <- createTimePeriod(FOSData, breaks=breaks)
FOSData$GearCat <- createGearCat(FOSData, cat=c("LONGLINE", "HOOK AND LINE"))
FOSData$FishSectorCat <- as.numeric(FOSData$FISHERY_SECTOR)

FOSData$Pres <- 1*(FOSData$YEYE_count>0)

### some graphics
hist(FOSData$YEYE_count[FOSData$YEYE_count>0], freq=FALSE)
mean(FOSData$YEYE_count[FOSData$YEYE_count>0])


#remove records for 2001, 2004 and 2005 and 2015
table(FOSData$YEAR)
ind <- which(FOSData$YEAR>=2006 & FOSData$YEAR<2015)


## by periods
by(FOSData, FOSData$TimePeriod, summary)
histTimePeriod <- lapply(sort(unique(FOSData$TimePeriod)), function(d){
  hist(FOSData$YEYE_count[FOSData$TimePeriod==d],freq=F, 
       main = paste0("From ", breaks[d], " to ", breaks[d+1]))
})
                         
FOSDataToAnalyze <- subset(FOSData, YEAR>=2006 & YEAR < 2015) 


# Preparing Files for jags ------------------------------------------------
FOSDataToAnalyze$YEARFact <-  as.factor(FOSDataToAnalyze$YEAR)
FOSDataToAnalyze$FISHERY_SECTOR <- droplevels(FOSDataToAnalyze$FISHERY_SECTOR)
nYear<- nlevels(FOSDataToAnalyze$YEARFact)
nGear <- nlevels(FOSDataToAnalyze$GEAR)
nFish <- nlevels(FOSDataToAnalyze$FISHERY_SECTOR)


## Presence data

presence <- by(FOSDataToAnalyze$YEYE_count, INDICES=list(FOSDataToAnalyze$YEARFact,FOSDataToAnalyze$GEAR, FOSDataToAnalyze$FISHERY_SECTOR), function(d) sum(d==0), simplify=F)
ysets <- by(FOSDataToAnalyze$YEYE_count, INDICES=list(FOSDataToAnalyze$YEARFact,FOSDataToAnalyze$GEAR, FOSDataToAnalyze$FISHERY_SECTOR), function(d) length(d))
ysets <- as.numeric(unlist(ysets))
possiblepresYear <- rep(attr(presence, "dimnames")[[1]], prod(attr(presence, "dim")[2:3]))
possiblepresGear <- rep(unlist(lapply(attr(presence, "dimnames")[[2]], 
                                                 function(d) {rep(d, attr(presence, "dim")[1])})), 
                                   attr(presence, "dim")[3])
possiblepresFish = unlist(lapply(attr(presence, "dimnames")[[3]], 
                         function(d) {rep(d, prod(attr(presence, "dim")[1:2]))}))


ex <- unlist(lapply(presence,function(d){!is.null(d)}))
formatPres <- data.frame(pres=as.numeric(unlist(presence)),
                         ysets=ysets[ex],
           presYear=  possiblepresYear[ex], 
           presGear=possiblepresGear[ex],
           presFish=possiblepresFish[ex]
        )

nPres <- nrow(formatPres)
pres <- formatPres$pres
presYear <- as.numeric(formatPres$presYear)
presGear <- as.numeric(formatPres$presGear)
presFish <- as.numeric(formatPres$presFish)
presYsets <- formatPres$ysets
formatPres$p = formatPres$pres/formatPres$ysets
plot(by(formatPres$p, formatPres$presYear, mean))

## Count Data
FOSDataToAnalyzePos <- FOSDataToAnalyze[which(FOSDataToAnalyze$YEYE_count>0),]
nCpue <- nrow(FOSDataToAnalyzePos)
cpue <- FOSDataToAnalyzePos$YEYE_count-1
cpueYear <- as.numeric(as.factor(FOSDataToAnalyzePos$YEAR))
cpueGear <- as.numeric(FOSDataToAnalyzePos$GEAR)
cpueFish <- as.numeric(FOSDataToAnalyzePos$FISHERY_SECTOR)




data.list <- list(nYear=nYear, nGear=nGear, nFish=nFish, 
                  nPres=nPres, pres=pres, presYear=presYear, 
                   presGear=presGear, presFish=presFish, presYsets=presYsets,
                   nCpue=nCpue, cpue=cpue, cpueYear=cpueYear, cpueGear=cpueGear, cpueFish=cpueFish,
                  betay=c(0, rep(NA, nYear-1)), betag=c(0, rep(NA, nGear-1)), betaf=c(0, rep(NA, nFish-1)),
                   alphay=c(0, rep(NA, nYear-1)), alphag=c(0, rep(NA, nGear-1)), alphaf=c(0, rep(NA, nFish-1))
)
attach(data.list)
dump(list= names(data.list), file="dataJags.txt")
detach(data.list)
# ### init
# indRef <- which(as.numeric(formatPres$presYear)==1 & as.numeric(formatPres$presGear)==1 & as.numeric(formatPres$presFish)==1)
# p0 <- sum(formatPres$pres[indRef])/sum(formatPres$ysets[indRef])
# beta0 <- log(p0/(1-p0))
# betay
# betag
# betaf
# alpha0
# alphay
# alphag
# alphaf


# Calling Jags ------------------------------------------------------------

library('rjags')
FOSmodel <- jags.model(file.path(mainDir, codeDir, 'PoissonModel.txt'), data=data.list,    n.chains = 1, n.adapt=1000, quiet=FALSE)
FOS.samples <- jags.samples(FOSmodel, c("alphay","betay","r"), n.iter=1000)
postBetay <- matrix(FOS.samples$betay, nrow=9, ncol=1000)
plot(postBetay[2,])
