### function to define classes from a quantitative variable
### data contain the quantitative variable
### returns a a numerical variabl with the class number
splitInClasses <- function(data, breaks){
  res <- sapply(data, function(d){
    sum(d>=breaks)
  })
  return(as.factor(res))
}  


## Fact is a character vector containing the list of factor to be accounted (except year which is always in the model)
setModel <- function(data, Fact)
{
  nl <- sapply(Fact, function(d) nlevels(droplevels(data[[d]])))
  ## if only one level in  the factor
  pl <- which(nl==1)
  if(length(pl)>=1) {
    data <- data[,which(!names(data)%in%Fact[pl])]
    col <- which(!is.na(str_locate(names(data), pattern="Fact")[,1]))
    Fact <- paste0("Fact", 1:length(col))
    names(data)[col]<- Fact
    print(names(data))
    nl <- nl[-pl]
  }
  nFact <- sum(nl>1)
  res<- switch(as.character(nFact), 
               "0"= 'Model.txt',
               "1" = 'Model1Fact.txt',
               "2" = 'Model2Fact.txt',
               "3" = 'Model3Fact.txt')
  return(list(model=res, Fact=Fact, data=data))
}


