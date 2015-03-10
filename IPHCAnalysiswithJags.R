### read Data 
#mainDir <- '/media/metienne_h/Consulting/DFO/2015-YelloweyeSurvey'
#mainDir <- '/home/metienne/bianca/Consulting/DFO/2015-YelloweyeSurvey/'
mainDir <- '/home/metienne/EnCours/2015-YelloweyeSurvey'
dataDir <- 'Data'
resDir <- 'Results/IPHC'
codeDir <- 'DFOYellowEye'


library('rjags')
library('longline')
library('stringr')
library('ggplot2')
library('scales')

source(file.path(mainDir, codeDir,"glmFunction.R"))
source(file.path(mainDir, codeDir, 'formattingFunction.R'))
nIter <- 10

if(!file.exists(file.path(mainDir, resDir)))
  dir.create(file.path(mainDir, resDir))

depthBreaks = c(0, 101, 201,  801 )
### 
IPHCData1 <- read.table(file.path(mainDir, dataDir, 'IPHC_1996-2002.csv'),
                        sep=",", header=T)
IPHCData1 <- IPHCData1[,which(names(IPHCData1)%in% c("Year", "Depth", "soaktime", "NT", "NNT", 
                                                    "Nb", "Ne"))]
IPHCData2 <- read.table(file.path(mainDir, dataDir, 'IPHC2003_2014.csv'),
                        sep=",", header=T)
IPHCData2 <- IPHCData2[,which(names(IPHCData2)%in% c("Year", "Depth", "soaktime", "NT", "NNT", 
                                                     "Nb", "Ne"))]


IPHCData <- rbind(IPHCData1, IPHCData2)
IPHCData$Fact2 <- splitInClasses(data = IPHCData$Depth, depthBreaks)
IPHCData$Fact1 <- IPHCData$Year

table(IPHCData$Fact1)
summary(IPHCData)
longline.IPHC <-  longline(fact1 = IPHCData$Fact1,
                           fact2 = IPHCData$Fact2,
                           nt = IPHCData$NT, nnt = IPHCData$NNT, nb = IPHCData$Nb, ne = IPHCData$Ne, s = IPHCData$soaktime)



CPUE <- ComputeCPUE.longline(longline.IPHC)

summary(longline.IPHC)

### 

for(m in 1:2)
{
  posIndices <- RunBayesEstimation(llData = longline.IPHC, MEM = m, peConst = F, filePath =   file.path(mainDir,resDir), nIter = 1000, burnin = 100)
}

