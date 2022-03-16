itemp <- read.delim(file=paste(wd,flexname,"-inf.txt",sep=""), header=FALSE, sep = "\t")
nitems <- nrow(itemp)-1

iteminfo <- as.matrix(itemp[1:nitems,3:83],dimnames=NULL)
testinfo <- itemp[(nitems+1),3:83]

xtemp <- scan(file=paste(wd,flexname,"-irt.txt",sep=""), what="real")

for (i in 1:length(xtemp)) {
  if (xtemp[i] == "Categories") newx <- xtemp[(i-1):((i+(nitems+1)*5)-2)]
}
newx <- matrix(newx,(nitems+1),5,byrow=TRUE)
newx <- matrix(as.numeric(newx[2:(nitems+1),1:2]),nitems,2)
rnames <- as.list(newx[,1])
mincat <- min(newx[,2])
maxcat <- max(newx[,2])
for(k in mincat:maxcat){
 cat <- matrix(0,nitems,1)
for(j in 1:nitems){
  if(newx[j,2]==k) cat[j]<-newx[j,1]
  else cat[j] <- NA
}
cat<-as.list(na.omit(cat))
assign(paste("cat",k,sep=""),cat)
}
for (i in 1:length(xtemp)) {
  if (xtemp[i] == "(Bock,") newz <- xtemp[(i-1):(i+17+(nitems*2)*maxcat*2)]
}

theta <- seq(-4.0,4.0,0.1)
cpar <- matrix(0,nitems,maxcat)
apar <- matrix(0,nitems,maxcat)
CBD <- matrix(0,nitems,(maxcat-1))
int <- matrix(0,nitems,(maxcat-1))
P <- NULL
info <- NULL
relinfo <- NULL

for(q in mincat:maxcat){
for (k in get(paste("cat",q,sep=""))) {
  for (i in 1:length(newz)) {
    if ((newz[i] == k) && (newz[i+2]=="a"))  {
      apar[k,1:q] <- newz[(i+3):(i+2+q)]
      cpar[k,1:q] <- newz[(i+4+q):(i+3+q*2)] 
      rnames[[k]] <- as.numeric(gsub("[^[:digit:]]", "",newz[i+1]))
      break
}}}}

apar <- matrix(as.numeric(apar),nitems,maxcat)
cpar <- matrix(as.numeric(cpar),nitems,maxcat)

logit <- NULL
catP <- matrix(0,length(theta),maxcat)
catinfo <- matrix(0,maxcat,ncol(iteminfo))
for(q in mincat:maxcat){
for(k in get(paste("cat",q,sep=""))) {
for(j in 1:q){
  logit[[j]]<-exp(apar[k,j] * theta + cpar[k,j])}
total <- Reduce("+",logit)
for(j in 1:q){
  catP[,j] <- logit[[j]]/total}
  P[[k]] <- catP
for(j in 1:q){
    catinfo[j,] <- iteminfo[k,] * catP[,j]
  }
  info[[k]] <- rbind(iteminfo[k,],catinfo)
  relinfo[[k]] <- info[[k]]/q
}}

for(k in mincat:maxcat){
for(m in 1:nitems){
  if(apar[m,k]==0){apar[m,k]<-NA}
  if(cpar[m,k]==0){cpar[m,k]<-NA}
}}

for(q in 1:(maxcat-1)){
CBD[,q]<-(apar[,(q+1)] - apar[,q])
int[,q] <- (cpar[,q] - cpar[,(q+1)]) / (CBD[,q])
}
rownames(CBD) <- unlist(rnames)
colnames(CBD) <- c(paste("CBD", 1:(maxcat-1)))
print(CBD)
rownames(int) <- unlist(rnames)
colnames(int) <- c(paste("Int", 1:(maxcat-1)))
print(round(int,2))


setwd(wd)
for(m in 1:(nitems)){
  bmp(paste(flexname,"CRC Item",rnames[[m]],".bmp"))
  plotP <- P[[m]]
  matplot(theta,plotP[,1:newx[m,2]],ylim=c(0,1),xlim=c(-4,4),xlab=expression(theta),ylab="Probability", type="l",lty=1,lwd=3,col=c(2:(maxcat+1)),main=paste("Category Response Curves
Item",rnames[[m]])) 
}
graphics.off()
for(m in 1:(nitems)){
  plotinfo <- t(relinfo[[m]])
  bmp(paste(flexname,"Info - Item",rnames[[m]],".bmp"))
  matplot(theta,plotinfo[,1:(newx[m,2]+1)],ylim=c(0,1),xlim=c(-4,4),xlab=expression(theta),ylab="Information", type="l",lty=1,lwd=3,col=c(1:(maxcat+1)),main=paste("Item and Category Information
Item", rnames[[m]]))
}
graphics.off()
bmp(paste(flexname,"Test Info.bmp",sep=""))
totalpars <- sum(newx[,2])-nitems
reltest <- testinfo/totalpars
matplot(theta,t(reltest),ylim=c(0,1),xlim=c(-4,4),xlab=expression(theta),ylab="Relative Information", type="l",lty=c(1),lwd=c(3),col=c(1),main=paste("Relative Test Information"))
graphics.off()

  

