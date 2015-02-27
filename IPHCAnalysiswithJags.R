### read Data 
#mainDir <- '/media/metienne_h/Consulting/DFO/2015-YelloweyeSurvey'
mainDir <- '/home/metienne/bianca/Consulting/DFO/2015-YelloweyeSurvey/'
#mainDir <- '/home/metienne/EnCours/2015-YelloweyeSurvey'
dataDir <- 'Data'
resDir <- 'Results'
codeDir <- '2015YellowEye'


library('rjags')
library('stringr')
source(file.path(mainDir, codeDir,"glmFunction.R"))
source(file.path(mainDir, codeDir, 'formattingFunction.R'))
nIter <- 10

breaks=c(1986,1989 , 1994, 1996, 2002, 2006,2016)
### 
IPHCData <- read.table(file.path(mainDir, dataDir, 'IPHC2003_2014.csv'),
                      sep=",", header=T)


summary(IPHCData)
# PreProcessing of the Data -----------------------------------------------
source(file.path(mainDir, codeDir, 'formattingFunction.R'))
## replace all na in YEYE_counts by 0
if(sum(is.na(IPHCData$NT))>0){
IPHCData$NT[which(is.na(IPHCData$NT))]=0
}

## exclude all depth >300
length(which(IPHCData$Depth>300))
nrow(IPHCData)
IPHCData <- IPHCData[ which(IPHCData$Depth<=300),]


IPHCData$YEAR <- IPHCData$Year
## create numeric time period
IPHCData$TimePeriod <- createTimePeriod(IPHCData, breaks=breaks)

IPHCData$Pres <- 1*(IPHCData$NT>0)

### some graphics
hist(IPHCData$NT[IPHCData$NT>0], freq=FALSE)
mean(IPHCData$NT[IPHCData$NT>0])

## which year represented
table(IPHCData$YEAR)

## for each time period
timePeriod <- unique(sort(IPHCData$TimePeriod))
res <- lapply(timePeriod, function(d){
  ## by periods
  Data <- IPHCData[which(IPHCData$TimePeriod==d),]
  by(Data, Data$TimePeriod, summary)
  pdf(file=file.path(mainDir, resDir,paste0("Hist",d,".pdf")))
  hist(Data$NT,freq=F, 
       main = paste0("From ", breaks[d], " to ", breaks[d+1]-1))
  dev.off()
  # Preparing Files for jags ------------------------------------------------
  Data$YEAR <-  as.factor(Data$YEAR)
  YEAR=levels(Data$YEAR)
  Data$AREA <- as.factor(Data$StatArea)
  nArea <- nlevels(Data$AREA)
  nYear<- nlevels(Data$YEAR)
  cat(paste(YEAR), "\n")
  ## Presence data
  
  presence <- by(Data$NT, INDICES=list(Data$YEAR, Data$AREA), function(d) sum(d==0), simplify=F)
  ysets <- by(Data$NT, INDICES=list(Data$YEAR, Data$AREA), function(d) length(d))
  ysets <- as.numeric(unlist(ysets))
  possiblepresYear <- rep(attr(presence, "dimnames")[[1]], prod(attr(presence, "dim")[2]))
    possiblepresArea = unlist(lapply(attr(presence, "dimnames")[[2]], 
                                   function(d) {rep(d, prod(attr(presence, "dim")[1]))}))
  
  
  ex <- unlist(lapply(presence,function(d){!is.null(d)}))
  formatPres <- data.frame(pres=as.numeric(unlist(presence)),
                           ysets=ysets[ex],
                           presYear=  possiblepresYear[ex], 
                           presArea=possiblepresArea[ex]
  )
  
  nPres <- nrow(formatPres)
  pres <- formatPres$pres
  presYear <- as.numeric(formatPres$presYear)
  presArea <- as.numeric(formatPres$presArea)
  presYsets <- formatPres$ysets
  formatPres$p = formatPres$pres/formatPres$ysets
  plot(by(formatPres$p, formatPres$presYear, mean))
  
  ## Count Data
  DataPos <- Data[which(Data$NT>0),]
  nCpue <- nrow(DataPos)
  cpue <- DataPos$NT-1
  cpueYear <- as.numeric(as.factor(DataPos$YEAR))
  cpueArea <- as.numeric(DataPos$AREA)
  
  
  
  
  data.list <- list(nYear=nYear,  nArea=nArea, 
                    nPres=nPres, pres=pres, presYear=presYear, 
                    presArea=presArea, presYsets=presYsets,
                    nCpue=nCpue, cpue=cpue, cpueYear=cpueYear, 
                    cpueArea=cpueArea,
                    #betay=c(0, rep(NA, nYear-1)),
                     betaf=c(0, rep(NA, nArea-1)),
                    #alphay=c(0, rep(NA, nYear-1)), 
                     alphaf=c(0, rep(NA, nArea-1))
  )
  
  # ### init
  if(nYear==1)
  {
    glm1 <- glm(cbind(formatPres$pres, formatPres$ysets-formatPres$pres)~formatPres$presGear+formatPres$presArea, family=binomial)
    }else{
    glm1 <- glm(cbind(formatPres$pres, formatPres$ysets-formatPres$pres)~formatPres$presYear+formatPres$presArea, family=binomial)
  }
  prescoef <- coef(glm1)
    betaf <- rep(NA, nArea)
    betay<- prescoef[1:nYear]
    betaf[2:nArea]<- prescoef[(nYear+1):(nYear+nArea-1)]
 
  r <- 1
  if(nYear==1)
  {
    glm2 <- glm(cpue~ as.factor(cpueArea), family=poisson)
  }else{
    glm2 <- glm(cpue~as.factor(cpueYear)+ as.factor(cpueArea), family=poisson)
  }
  cpuecoef=coef(glm2)
   alphay <- cpuecoef[1:nYear]
  alphaf <- rep(NA, nArea)
  alphaf[2:nArea]<- cpuecoef[(nYear+1):(nYear+nArea-1)]
  
  init.list <- list(r=r,  betay =betay,  betaf=betaf,
                     alphay =alphay, alphaf=alphaf
  )
  # Calling Jags ------------------------------------------------------------
  
  
  IPHCmodel <- jags.model(file.path(mainDir, codeDir, 'PoissonModel.txt'), 
                         data=data.list,    n.chains = 1, n.adapt=1000, 
                         quiet=FALSE, inits=init.list)
  IPHC.samples <- jags.samples(IPHCmodel, c("alphay","betay","r", "predccpue"), n.iter=nIter)
  save(list=c("YEAR", "IPHC.samples"),file=file.path(mainDir, resDir, paste0(d,"IPHCsamples.Rd")))
  
  
  # Producing Results -------------------------------------------------------
  
  load(file.path(mainDir, resDir, paste0(d,"IPHCsamples.Rd")))
  nIter <- dim(IPHC.samples$betay)[2]
  postBetay <- data.frame(index=as.numeric(IPHC.samples$betay), year=rep(YEAR, each=nIter))
  postPredCPUE <- data.frame(index=as.numeric(IPHC.samples$predccpue), year=rep(YEAR, nIter))
  write.csv(summary(as.mcmc.list(IPHC.samples$predccpue))[[1]], 
            file=file.path(mainDir, resDir, paste0(d, "IPHC-Abundance.csv")),
  )
  pdf(file.path(mainDir, resDir, paste0(d, "IPHC-Abundance.pdf")))
  plot(postPredCPUE$index~as.factor(postPredCPUE$year), ylim=c(0, 0.05), xlab="Year", ylab="Relative Abundance")
 dev.off()
  return(summary(as.mcmc.list(IPHC.samples$predccpue))[[1]])
}
)
resFinal <- Reduce(rbind, res)