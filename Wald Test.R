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

for(k in mincat:maxcat){
  for(m in 1:nitems){
    if(apar[m,k]==0){apar[m,k]<-NA}
    if(cpar[m,k]==0){cpar[m,k]<-NA}
  }}

for(q in 1:(maxcat-1)){
  CBD[,q]<-(apar[,(q+1)] - apar[,q])
  int[,q] <- (cpar[,q] - cpar[,(q+1)]) / (CBD[,q])
}
pnum <- matrix(0,nitems,1)
for (k in 1:nitems) {
for(n in 1:length(xtemp)){
  if (xtemp[n] == paste("v",rnames[[k]], sep="")) { 
    pnum[k,] <- as.numeric(xtemp[n+1])
    break
}}}

for(q in mincat:maxcat){
  for (k in get(paste("cat",q,sep=""))) {
    pnum[k,] <- pnum[k,]-(q-2)
  }}

sigma <- read.csv(file=paste(wd,flexname,"-cov.txt",sep=""), header=FALSE,sep=",", dec=".")

ahat <- matrix(0,(maxcat-1),1)
sigmahat <- matrix(0,(maxcat-1),(maxcat-1))
Q <- matrix(0,nitems,1)
pro <- matrix(0,nitems,1)
df <- matrix(0,nitems,1)

for(q in mincat:maxcat){
  for (k in get(paste("cat",q,sep=""))) {
    ahat <- matrix(0,(q-1),1)
    sigmahat <- matrix(0,(q-1),(q-1))
    ahat <- as.matrix(CBD[k,1:(q-1)])
    sigmahat <- as.matrix(sigma[pnum[k]:(pnum[k]+(q-2)),pnum[k]:(pnum[k]+(q-2))])
    
  l <- t(contr.poly(q-1))
  lahat <- l %*% ahat
  lsigmalt <- l %*% sigmahat %*% t(l)
  
  Q[k,] <- t(lahat) %*% solve(lsigmalt) %*% lahat
  pro[k,] <- pchisq(Q[k,], df= (q-2), lower.tail=FALSE)
  df[k,] <- (q-2)
}}
# 2 because constraint on both a and c
(  out <- round(cbind(Q,df,pro),digits=3))

rownames(out) <- unlist(rnames)
colnames(out) <- c("Q",'df',"p-val")
print(out)
