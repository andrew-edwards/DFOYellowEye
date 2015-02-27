require('rjags')
require('coda')
### function to perform Bayesian tends analysis for cpue Using a binomial Negative  model
## dataPeriod contains the data with at least 
##          a column NT containiung the Catch
##          a column Year containing Year of each record
##          one column for every factor to be accounted for (names Fact1, Fact2, Fact3)
## Fact a character vector containing the factor of interest
## fig.path an optional path to store the produced graphs
## nIter the number of iteration for the MCMC approach
## pref a prefix for the name of the output files

glmAnalysis <- function(dataPeriod, Fact, fig.path="", nIter=10, pref="")
{
  per <- unique(dataPeriod$TimePeriod)
  ms <- setModel(dataPeriod, Fact = Fact)
  model <- ms$model
  if( length(ms$Fact)<length(Fact)){
    Fact <- ms$Fact
    dataPeriod <- ms$data
  }
    yearLim <- levels(dataPeriod$Year)[c(1,nlevels(dataPeriod$Year))]
    pdf(file.path(fig.path,paste0( "CountHist-", per,".pdf")))  
    hist(dataPeriod$NT,freq=F, 
         main = paste0("From ", yearLim[1], " to ", yearLim[2]))
    dev.off()
  
  # Preparing Files for jags ------------------------------------------------
  dataPeriod$Year <-  as.factor(dataPeriod$Year)
  YEAR=levels(dataPeriod$Year)
  
  ## reordering level to have the most present as references
  invisible(lapply(Fact, function(d){
    dataPeriod[,d] <<- droplevels(dataPeriod[,d])
    ref <- which.max(table(dataPeriod[,d]))
    dataPeriod[,d] <-relevel(dataPeriod[,d], ref=ref)
  }))
  
  nYear<- nlevels(dataPeriod$Year)
  
  ## Presence data
  presenceData <- as.list(by(dataPeriod, 
                         INDICES=lapply(c("Year", Fact), function(j) {dataPeriod[,j]}), 
                         function(d){ if(!is.null(d)) return(as.numeric(d[1,c("Year", Fact), ]))}))
  presenceData <- as.data.frame(Reduce(rbind, presenceData))
  names(presenceData) <- paste0("pres",c("Year", Fact))
  presenceData$pres <- unlist(by(dataPeriod$NT, 
                 INDICES=lapply(c("Year", Fact), function(j) {dataPeriod[,j]}), 
                 function(d) sum(d>0), simplify=F))
  presenceData$presYsets <- unlist(by(dataPeriod$NT, 
              INDICES=lapply(c("Year", Fact), function(j) {dataPeriod[,j]}), 
              function(d) {if(!is.null(d))  sum(d>=0)}, simplify=F))
  nPres <- nrow(presenceData)
  
  
  
  ## Count Data
  cpueData <- dataPeriod[which(dataPeriod$NT>0),c("Year", Fact)]
  cpueData <- as.data.frame(sapply(names(cpueData), function(d) as.numeric(cpueData[,d])))
  names(cpueData)<- paste0("cpue", c("Year", Fact))
  nCpue <- nrow(cpueData)
  cpueData$cpue <- dataPeriod[which(dataPeriod$NT>0),"NT"]-1
  invisible(lapply(paste0("cpue", c("Year", Fact)), function(d){ cpueData[,d]<<-as.factor(cpueData[,d])}))
  invisible(lapply(paste0("pres", c("Year", Fact)), function(d){ presenceData[,d]<<-as.factor(presenceData[,d])}))
  
  l1 <- lapply(1:length(Fact), function(d){ 
    c(0, rep(NA, nlevels(dataPeriod[,paste0("Fact",d)])-1))})
  l2 <- l1
  l3 <- lapply(1:length(Fact), function(d){  nlevels(dataPeriod[,paste0("Fact",d)])})
  names(l1) <- paste0("beta", 1:length(Fact))
  names(l2) <- paste0("alpha", 1:length(Fact))
  names(l3) <- paste0("n", Fact)
  
  data.list<- c(list(nYear=nYear, nPres=nPres, nCpue=nCpue), cpueData, presenceData,
                l1,l2,l3)
  # ### init presence
  if(nYear==1)
  {
    formula1 <- paste0( "cbind(pres,presYsets-pres)~ ", paste(paste0("pres",Fact), collapse = " +"))
    formula2 <- paste0( "cpue~ ", paste(paste0("cpue",Fact), collapse = " +"))
  }else{
    formula1 <- paste0( "cbind(pres,presYsets-pres)~ presYear+", paste(paste0("pres",Fact), collapse = " +"))
    formula2 <- paste0( "cpue~ cpueYear+", paste(paste0("cpue",Fact), collapse = " +"))
  }
  
  glm1 <- glm(formula1, data=presenceData, family=binomial)
  glm2 <- glm(formula2, data=cpueData, family=poisson)
  
  prescoef <- coef(glm1)
  cpuecoef <- coef(glm2)
  betay<-prescoef[1:nYear]
  alphay<-cpuecoef[1:nYear]
  
  initl1 <- lapply(Fact, function(d){ 
    n <- nlevels(dataPeriod[,d])
    res<- c(NA, rep(1, n-1))
    names(res) <- paste0("pres", d, 1:n)
    ind1 <- which(!is.na(str_locate(paste0(names(prescoef)), pattern = d)[,1]))
    ind2 <- sapply(1:length(ind1), function(i){
      which(!is.na(str_locate(names(res), pattern = names(prescoef)[ind1[i]])[,1]))})
    res[ind2] <- prescoef[ind1]
    return(res)
  })
  names(initl1)<- paste0("beta", 1:length(Fact))
  
  initl2 <- lapply(Fact, function(d){ 
    n <- nlevels(dataPeriod[,d])
    res<- c(NA, rep(1, n-1))
    names(res) <- paste0("cpue", d, 1:n)
    ind1 <- which(!is.na(str_locate(paste0(names(cpuecoef)), pattern = d)[,1]))
    ind2 <- sapply(1:length(ind1), function(i){
      which(!is.na(str_locate(names(res), pattern = names(cpuecoef)[ind1[i]])[,1]))})
    res[ind2] <- cpuecoef[ind1]
    return(res)
  })
  names(initl2)<- paste0("alpha", 1:length(Fact))
  r <- 1
  
  init.list <- c(list(r=r,  betay =betay, alphay=alphay), initl1, initl2)
  
  # Calling Jags ------------------------------------------------------------
  
  jagsRun <- jags.model(file.path(mainDir, codeDir, paste0('Poisson', model)), 
                         data=data.list,    n.chains = 1, n.adapt=500, 
                         quiet=FALSE, inits=init.list)
  cat("Launch estimation\n")
  out.samples <- jags.samples(jagsRun, c("alphay","betay","r", "predccpue"), n.iter=nIter)
  save(list=c("YEAR","out.samples"),file=file.path(mainDir, resDir, paste0(pref,"-outSamples.Rd")))
  
  
  # Producing Results -------------------------------------------------------

load(file.path(mainDir, resDir, paste0(pref, "-outSamples.Rd")))
nIter <- dim(out.samples$betay)[2]
postBetay <- data.frame(index=as.numeric(out.samples$betay), 
                        year=rep(YEAR, nIter))
postPredCPUE <- data.frame(index=as.numeric(out.samples$predccpue), 
                           year=rep(YEAR, nIter))
summaryCpue <- summary(as.mcmc.list(out.samples$predccpue))[[1]]
write.csv(summaryCpue , 
          file=file.path(mainDir, resDir,  paste0(pref,"-Abundance.csv")))

pdf(file.path(mainDir, resDir, paste0(pref ,"-Abundance.pdf")))
plot(postPredCPUE$index~as.factor(postPredCPUE$year), xlab="Year", ylab="Relative Abundance")
dev.off()

lapply(YEAR, function(y){
  pdf(file.path(fig.path, paste0(pref ,"-History-beta-",y,".pdf")))
  plot(1:nIter,postBetay[postBetay$year==y,1],"l", xlab="iterations", ylab="Year effect", main=y)
  dev.off()
  pdf(file.path(fig.path, paste0(pref ,"-History-cpue-",y,".pdf")))
    plot(1:nIter,postPredCPUE[postPredCPUE$year==y,1],"l", xlab="iterations", ylab="Relative Abundance Index", main=y)
  dev.off()
})
return(summaryCpue)
}
