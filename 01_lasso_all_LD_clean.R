### This script capture the lasso program for the neighbors in metabolite network.
### metabolite network based on stoichiometric matrix. 
### For furher information: tonghao0605@gmail.com

dir <- "D:/01_neigh/metaneighbor"
setwd(dir)

library(lattice)
library(Matrix)
library(foreach)
library(glmnet)
library(genetics)
library(rrBLUP)

## genotype data
xx <- read.table("IL_numeric_final.csv",sep="\t",header=F)
## metabolite data
pheno <- read.table("phenoblup.csv",sep=",",header=T)[,2:39]
pheno <- t(pheno)
metaid <- read.table("metaid_final.csv",sep=",",header=T) #match ID in Smatrix
## GPR path data
neigh <- read.table("Smatrix.csv",sep=",",header=F)
react <- read.table("reactgene_final.csv",sep=",",header=F,stringsAsFactors=F)[,2:69]
genesnp <- read.table("snpgene.txt",sep="\t",header=T,stringsAsFactors=F)[,c(1,4)]
snpid <- read.table("5SNP_info.csv",sep=",",header=F)

####################################################################################
#Function Sections
####################################################################################

###################################################################
## Lasso function  #####elastic net model
## The basic lasso model  #####elastic net model
lasso<-function(x,y,geno,phe,alpha){

fold <- 10

fit<-cv.glmnet(x=x,y=y,alpha=alpha)
#plot(fit)
theta<-coef(fit,s="lambda.min")
beta<-theta[1]
yhat<-predict(fit,newx=x,s="lambda.min")
goodness<-drop(cor(y,yhat)^2)
lambda<-fit$lambda.min
index<-match(lambda,fit$lambda)
#mse.cv<-fit$cvm[index]
nonzero<-fit$nzero[index]
parm<-data.frame(beta,alpha,lambda,nonzero,goodness)
#write.csv(x=parm,file=paste("parms_nsnp",geno,"_meta",phe,".csv",sep=""),row.names=T)

return(c(parm))

}


## The classic lasso model   #####elastic net model
lassof<-function(x,y,geno,phe,lambda,alpha,sid){

fit<-glmnet(x=x,y=y,alpha=alpha,lambda=lambda)
theta<-coef(fit)
beta<-theta[1]
yhat<-predict(fit,newx=x,s=lambda)
expdev<-fit$dev.ratio
nulldev <- fit$nulldev
nonzero<-fit$df
goodness<-drop(cor(y,yhat)^2)

n <- fit$nobs

# definition one
lnL0 <- fit$nulldev - deviance(fit)

# definition two
resid <- y - yhat
rss <- sum(resid^2) 
lnL <- n*log(rss/n)

# definition one
k <- fit$df

# definition two
#ld <- lambda * diag(ncol(x))
#H <- x %*% solve(t(x) %*% x + ld) %*% t(x)
#df <- psych::tr(H)

aic <- 2*k + lnL

aicc <- aic + 2*k*(k+1)/(n-k-1)

bic <- log(n)*k + lnL

#tss <- sum((y-mean(y))^2)
#r2 <- 1-(rss/tss)

parm<-data.frame(beta,nonzero,goodness,expdev,lnL0,lnL,aic,aicc,bic,nulldev)

gamma<-theta[-1]
sigma2<-sum((y-yhat)^2)/nrow(y)

m<-ncol(x)
effectall <- NULL
for(k in 1:m){
if (gamma[k] != 0){
  effect <- NULL
  t1<-t(x[,k])%*%x[,k]
  t2<-gamma[k]^2*t1
  t3<-sqrt(t2^2+4*t2*sigma2)
  sigma2k<-(t2+t3)/(2*t1)
  var<-sigma2*sigma2k/(sigma2+sigma2k*t1)
  stderr<-sqrt(var)
  if(gamma[k]==0) {
     wald<-0
  }  else {
     wald<-gamma[k]^2/var
  }
  lod<-wald/log(100)
  pvalue<-pchisq(wald,1,lower.tail=F)
  
  sidk <- sid[k] #genome SNP IDs
  effect <- cbind(k,sidk,gamma[k],stderr,wald,lod,pvalue)
  effectall <- rbind(effectall,effect)
  }
}
colnames(effectall) <- c("msnpID","allsnpID","gamma","stderr","wald","LOD","pvalue")

write.csv(x=effectall,file=paste("eneteffect",geno,phe,".csv",sep=""),row.names=F)

return(c(parm))

}


###################################################################
## Extract SNPs within gene set
snpfromgene<-function(gene){
sidd <- NULL
ngene <- length(gene)
	for (j in 1:ngene){
		sid <- which(genesnp[,2]==gene[j])
		sname <- genesnp[sid,1]
		sidnn <- NULL
		for (k in 1:length(sname)){
			sidn <- which(snpid[,1]==sname[k])
			sidnn <- c(sidnn,sidn)##gene SNP
		}
		sidd <- c(sidd, sidnn)
	}
	
siddu <- unique(sidd)
return(siddu)		
}


####################################################################################
#STEP 0 (extract SNPs of each reaction)
####################################################################################

nr <- ncol(neigh)  #reaction number

nsnp <- NULL
metallsnp <- NULL
ngene <- NULL
metallgene <- NULL

for (n in 1:nr){

	print(paste("Reaction",n,"!",sep=""))
	gidd <- NULL #final gene id for a reaction
	sidd <- NULL #final snp id for a reaction

	gid <- as.character(react[n,])	##gene name
	gidd <- gid[which(gid != 0)]
	sidd <- NULL
	for (j in 1:68){
	if (gid[j] != "0"){
		sid <- which(genesnp[,2]==gid[j])
		sname <- genesnp[sid,1]
		sidnn <- NULL
		for (k in 1:length(sname)){
		sidn <- which(snpid[,1]==sname[k])
		sidnn <- c(sidnn,sidn)##a gene SNP
		}
		
		sidd <- c(sidd, sidnn)
	}}
	#sidf <- c(sidf, sidd)
	#gidf <- c(gidf, gidd)
	
	#sidf <- unique(sidf)
	nsnpm <- length(sidd)
	#gidf <- unique(gidf)
	ngenem <- length(gidd)
	
	ngene <- c(ngene,ngenem)
	metallgene <- c(metallgene,gidd)
	nsnp <- c(nsnp,nsnpm)
	metallsnp <- c(metallsnp,sidd)
}

write.csv(nsnp,"nsnpall.csv",row.names=F)
write.csv(metallsnp,"metallsnpall.csv",row.names=F)
write.csv(ngene,"ngeneall.csv",row.names=F)
write.csv(metallgene,"metallgeneall.csv",row.names=F)


####################################################################################
#STEP 1 (elastic net, 1st neighborhood SNPs, eatimate lambda)
####################################################################################
########## gene SNPs and neiborhood SNPs ##########
mt <- c(1:28)
nmeta <- NULL
metallmeta <- NULL
nsnp <- NULL
metallsnp <- NULL
ngene <- NULL
metallgene <- NULL
lall <- NULL
#lallall <- NULL

for (n in mt){
print(paste("Metabolite",n,"!",sep=""))
sidf <- NULL #final snp id for a metabolite
gidf <- NULL #final gene id for a metabolite
metaidreact <- metaid[n,1]
nid <- neigh[metaidreact ,]
nidd <- which(nid != 0)
nnidd <- length(nidd)
metallmeta <- c(metallmeta,nidd)
nmeta <- c(nmeta,nnidd)

phe <- metaid[n,2]
meta <- as.numeric(pheno[metaid[n,2],])

for (i in nidd){
	#neighbor metabolite
	gid <- as.character(react[i,])	##gene name
	gidd <- gid[which(gid != 0)]
	sidd <- NULL
	for (j in 1:68){
	if (gid[j] != "0"){
		sid <- which(genesnp[,2]==gid[j])
		sname <- genesnp[sid,1]
		sidnn <- NULL
		for (k in 1:length(sname)){
		sidn <- which(snpid[,1]==sname[k])
		sidnn <- c(sidnn,sidn)##gene SNP
		}
		
		sidd <- c(sidd, sidnn)
	}}
	sidf <- c(sidf, sidd)
	gidf <- c(gidf, gidd)
	}
	
	sidf <- unique(sidf)
	nsnpm <- length(sidf)
	gidf <- unique(gidf)
	ngenem <- length(gidf)
	
	#### Lasso model selection
	x <- xx[sidf,]
	x <- t(as.matrix(x))

	phe <- paste("meta",n,sep="")
	meta <- as.numeric(pheno[metaid[n,2],])
	meta <- meta - mean(meta)
	meta <- matrix(meta,513,1)
	geno <- "nb1"
	
	### diffenent fold to select best lambda
	nn <- 50 #number of alpha
	parmrept <- NULL
	for (i in 1:nn){
	print(paste("-----",i,"out of",nn,"DONE! -----",sep=" "))	
	alpha <- i/nn
	
	rept = 100 #number of corss validation

	
	for (i in 1:rept){
		parm <- NULL
		parm <- lasso(x,meta,geno,phe,alpha)
		parmrept <- rbind(parmrept, parm)
	}
	### repeat end.
	
	}

	parmreptbest <- parmrept[which(parmrept[,5]==max(as.numeric(parmrept[,5]),na.rm=T)),]
	write.csv(parmrept,paste("enealltresults",geno,phe,".csv",sep=""),row.names=F)
	
	nsnp <- c(nsnp,nsnpm)
	metallsnp <- c(metallsnp,sidf)
	#metallsnp <- unique(metallsnp)
	#### metabolites SNPs including overlaps!!!!
	
	ngene <- c(ngene,ngenem)
	metallgene <- c(metallgene,gidf)
	
	#lallall <- rbind(lallall,parmrept)
	lall <- rbind(lall,parmreptbest)
	
}

write.csv(lall,paste("enetresults",n,".csv",sep=""),row.names=F)
write.csv(nmeta,paste("nmeta",n,".csv",sep=""),row.names=F)
write.csv(metallmeta,paste("metallmeta",n,".csv",sep=""),row.names=F)
write.csv(nsnp,paste("nsnp",n,".csv",sep=""),row.names=F)
write.csv(metallsnp,paste("metallsnp",n,".csv",sep=""),row.names=F)
write.csv(ngene,paste("ngene",n,".csv",sep=""),row.names=F)
write.csv(metallgene,paste("metallgene",n,".csv",sep=""),row.names=F)
#write.csv(lallall,"enealltresults.csv",row.names=F)


####################################################################################
#STEP 2 (elastic net, 1st neighborhood SNPs, controlled lambda)
####################################################################################

###read neighborhood snp id dataset
nsnp <- read.table("nsnp.csv",sep=",",header=T)
metallsnp <- read.table("metallsnp.csv",sep=",",header=T)
nsnp <- c(as.matrix(nsnp))
metallsnp <- c(as.matrix(metallsnp))

#### extract the best lambda ####
parmall <- read.table("enetresultsall.csv",sep=",",header=T)
alphall <- as.vector(parmall[,2])
lambdall <- as.vector(parmall[,3])

########## 1st neiborhood SNPs final model ##########
mt <- c(1:24,26:27)

parmall <- NULL
for (n in mt){
	
	#print(paste("metabolite",n,"!",sep=""))
	
	#####################################################
	#extract 1st neighborhood dataset#
	nsnpn <- nsnp[n]
	if (n != 1) {sidfn <- metallsnp[(sum(nsnp[1:(n-1)])+1):(sum(nsnp[1:n]))]}
	if (n ==1) {sidfn <- metallsnp[1:nsnp[n]]}
	x <- xx[sidfn,]
	x <- t(as.matrix(x))
	
	meta <- as.numeric(pheno[metaid[n,2],])
	meta <- meta - mean(meta)
	meta <- matrix(meta,513,1)
	
	geno <- "n1"
	phe <- paste("meta",n,sep="")
	
	###report significant SNPs in the best model 
	alphabest <- as.numeric(alphall[n])
	lambdabest <- as.numeric(lambdall[n])
	#lassoeffect(x,meta,geno,phe,lambdabest,alphabest,sidfn)
	
	parm <- lassof(x,meta,geno,phe,lambdabest,alphabest,sidfn)
	parm <- c(parm,alphabest,lambdabest)
	parmall <- rbind(parmall,parm)
	colnames(parmall) <- c("intercept","nonzero","goodness","expdev","lnL0","lnL","aic","aicc","bic","nulldev","alpha","lambda")

}

write.csv(parmall,"n1SNPcontrolled.csv",row.names=F)

####################################################
#### cross-validation R2 ####

idsall <- read.table("04_rrBLUP/foldid.csv",sep=",",header=F)
rr <- 5 #replicate number
f <- 5 #fold number
nn <- 513
		
corall <- NULL
mt <- c(1:24,26:27)
for (n in mt){

	print(paste("Metabolite ",n," Start!",sep=""))
	
	nsnpn <- nsnp[n]
	if (n != 1) {sidfn <- metallsnp[(sum(nsnp[1:(n-1)])+1):(sum(nsnp[1:n]))]}
	if (n ==1) {sidfn <- metallsnp[1:nsnp[n]]}
	x <- xx[sidfn,]
	x <- t(as.matrix(x))
	
	meta <- as.numeric(pheno[metaid[n,2],])
	meta <- meta - mean(meta)
	meta <- matrix(meta,513,1)
	
	#alphabest <- as.numeric(alphall[n])
	#lambdabest <- as.numeric(lambdall[n])

## cross-validation

corr <- NULL
for (r in 1:rr){

#print(paste("Repetition ",r," Start!",sep=""))
ids <- idsall[,r]

Y <- meta	

for (p in 1:f){
  trset <- which(ids!=p)
  teset <- which(ids==p)
  xxtr <- x[trset,]
  xxte <- x[teset,]
  y <- Y[trset,]
  sol <- mixed.solve(y,Z=xxtr,K=NULL,SE=F)
  ucoef <- as.matrix(sol$u)
  ypredx <- rep(sol$beta,length(teset)) + as.numeric(xxte %*% ucoef)
  yp <- Y[teset,]
  cor <- cor(yp,ypredx,method="pearson")
  corr <- rbind(corr,cor)
}


}

corm <- mean(corr[,1])
corrm <- rbind(corr,corm)
corall <- cbind(corall,corrm)

write.csv(corall,"04_rrBLUP/n1SNP_rrBLIP_all.csv",row.names=F)

}


####################################################################################
#STEP 3 (elastic net, all metabolite SNPs, controlled lambda)
####################################################################################

metallsnpall <- read.table("metallsnpall.csv",sep=",",header=T)
metallsnpall <- c(as.matrix(metallsnpall))
index <- unique(metallsnpall)

#### extract the best lambda ####
parmall <- read.table("enetresultsall.csv",sep=",",header=T)
alphall <- as.vector(parmall[,2])
lambdall <- as.vector(parmall[,3])

xm <- xx[index,]
xm <- t(as.matrix(xm))
		
mt <- c(1:24,26:27)
parmall <- NULL
for (n in mt){
	parm <- NULL
	geno <- "metasnp"
	phe <- paste("meta",n,sep="")
	meta <- as.numeric(pheno[metaid[n,2],])
	meta <- meta - mean(meta)
	meta <- matrix(meta,513,1)
	alphabest <- as.numeric(alphall[n])
	lambdabest <- as.numeric(lambdall[n])
	parm <- lassof(xm,meta,geno,phe,lambdabest,alphabest,index)
	parm <- c(parm,alphabest,lambdabest)
	parmall <- rbind(parmall,parm)
	colnames(parmall) <- c("intercept","nonzero","goodness","expdev","lnL0","lnL","aic","aicc","bic","nulldev","alpha","lambda")
}
write.csv(parmall,"mSNPcontrolled.csv",row.names=F)

####################################################
#### cross-validation R2 ####

metallsnpall <- read.table("metallsnpall.csv",sep=",",header=T)
metallsnpall <- c(as.matrix(metallsnpall))
index <- unique(metallsnpall)

#### extract the best lambda ####
parmall <- read.table("enetresultsall.csv",sep=",",header=T)
alphall <- as.vector(parmall[,2])
lambdall <- as.vector(parmall[,3])

xm <- xx[index,]
xm <- t(as.matrix(xm))

idsall <- read.table("04_cv/foldid.csv",sep=",",header=F)
rr <- 50 #replicate number
f <- 5 #fold number
nn <- 513
		
corall <- NULL
mt <- c(1:24,26:27)
for (n in mt){

	print(paste("Metabolite ",n," Start!",sep=""))
	
	#enetresults <- read.table(paste("eneteffectmetasnpmeta",n,".csv",sep=""),sep=",",header=T)
	#idnz <- enetresults[,2]
	#xm <- xx[idnz,]
	#xm <- t(as.matrix(xm))
	
	meta <- as.numeric(pheno[metaid[n,2],])
	meta <- meta - mean(meta)
	meta <- matrix(meta,513,1)
	
	alphabest <- as.numeric(alphall[n])
	lambdabest <- as.numeric(lambdall[n])

## cross-validation

corr <- NULL
for (r in 1:rr){

#print(paste("Repetition ",r," Start!",sep=""))
ids <- idsall[,r]

Y <- meta	

for (p in 1:f){
	trset <- which(ids!=p)
	teset <- which(ids==p)
	xxtr <- xm[trset,]
	xxte <- xm[teset,]
	y <- matrix(Y[trset,],length(trset),1)
	
	fit<-glmnet(x=xxtr,y=y,alpha=alphabest,lambda=lambdabest)
	#theta<-coef(fit)
	#beta<-theta[1]
	#gamma<-theta[-1]
	yhat<-predict(fit,newx=xxte,s=lambdabest)
	#expdev<-fit$dev.ratio
	#nonzero<-fit$df
	#goodness<-cor(y,yhat,method="pearson")
	#parm<-data.frame(beta,nonzero,goodness,expdev,yhat)
	#parmall <- rbind(parmall,parm)
	#colnames(parmall) <- c("intercept","nonzero","goodness","expdev","lnL0","lnL","aic","aicc","bic","nulldev","alpha","lambda")

	yp <- Y[teset,]
	cor <- cor(yp,yhat,method="pearson")
	corr <- rbind(corr,cor)
}


}

corm <- mean(corr[,1])
corrm <- rbind(corr,corm)
corall <- cbind(corall,corrm)

write.csv(corall,"04_cv/mSNPcv_all.csv",row.names=F)

}


###intersect significant 1st neighbor SNP and significant metabolite SNP
mt <- c(1:24,26:27)
nsnpmgall <- NULL
for (n in mt){
	sign1snp <- read.table(file=paste("eneteffectn1effectmeta",n,".csv",sep=""),sep=",",header=T)
	sigmsnp <- read.table(file=paste("eneteffectmetasnpmeta",n,".csv",sep=""),sep=",",header=T)
	snpmg <- intersect(sign1snp[,2],sigmsnp[,2])
	nsnpmg <- length(snpmg)
	nsnpmgall  <- rbind(nsnpmgall, nsnpmg)
}

write.csv(nsnpmgall,"mSNPandn1SNP.csv",row.names=F)

####################################################################################
#STEP 4 (elastic net, all genome SNPs, controlled lambda)
####################################################################################

#### extract the best lambda ####
parmall <- read.table("enetresultsall.csv",sep=",",header=T)
alphall <- as.vector(parmall[,2])
lambdall <- as.vector(parmall[,3])

index <- c(1:551616)
xa <- t(as.matrix(xx))
		
mt <- c(1:24,26:27)
parmall <- NULL
for (n in mt){
	print(paste("Metabolite",n,"!",sep=""))
	parm <- NULL
	geno <- "genomesnp"
	phe <- paste("meta",n,sep="")
	meta <- as.numeric(pheno[metaid[n,2],])
	meta <- meta - mean(meta)
	meta <- matrix(meta,513,1)
	alphabest <- as.numeric(alphall[n])
	lambdabest <- as.numeric(lambdall[n])
	parm <- lassof(xa,meta,geno,phe,lambdabest,alphabest,index)
	parm <- c(parm,alphabest,lambdabest)
	parmall <- rbind(parmall,parm)
	
	colnames(parmall) <- c("intercept","nonzero","goodness","expdev","lnL0","lnL","aic","aicc","bic","nulldev","alpha","lambda")

	write.csv(parmall,"gSNPcontrolled.csv",row.names=F)
}

####################################################
#### cross-validation R2 ####

#### extract the best lambda ####
parmall <- read.table("enetresultsall.csv",sep=",",header=T)
alphall <- as.vector(parmall[,2])
lambdall <- as.vector(parmall[,3])

index <- c(1:551616)
xa <- t(as.matrix(xx))

idsall <- read.table("04_cv/foldid.csv",sep=",",header=F)
rr <- 50 #replicate number
f <- 5 #fold number
nn <- 513
		
corall <- NULL
mt <- c(1:24,26:27)
for (n in mt){

	print(paste("Metabolite ",n," Start!",sep=""))
	
	#enetresults <- read.table(paste("eneteffectmetasnpmeta",n,".csv",sep=""),sep=",",header=T)
	#idnz <- enetresults[,2]
	#xm <- xx[idnz,]
	#xm <- t(as.matrix(xm))
	
	meta <- as.numeric(pheno[metaid[n,2],])
	meta <- meta - mean(meta)
	meta <- matrix(meta,513,1)
	
	alphabest <- as.numeric(alphall[n])
	lambdabest <- as.numeric(lambdall[n])

## cross-validation

corr <- NULL
for (r in 1:rr){

print(paste("Repetition ",r," Start!",sep=""))
ids <- idsall[,r]

Y <- meta	

for (p in 1:f){
	trset <- which(ids!=p)
	teset <- which(ids==p)
	xxtr <- xa[trset,]
	xxte <- xa[teset,]
	y <- matrix(Y[trset,],length(trset),1)
	
	fit<-glmnet(x=xxtr,y=y,alpha=alphabest,lambda=lambdabest)
	#theta<-coef(fit)
	#beta<-theta[1]
	#gamma<-theta[-1]
	yhat<-predict(fit,newx=xxte,s=lambdabest)
	#expdev<-fit$dev.ratio
	#nonzero<-fit$df
	#goodness<-cor(y,yhat,method="pearson")
	#parm<-data.frame(beta,nonzero,goodness,expdev,yhat)
	#parmall <- rbind(parmall,parm)
	#colnames(parmall) <- c("intercept","nonzero","goodness","expdev","lnL0","lnL","aic","aicc","bic","nulldev","alpha","lambda")

	yp <- Y[teset,]
	cor <- cor(yp,yhat,method="pearson")
	corr <- rbind(corr,cor)
}


}

corm <- mean(corr[,1])
corrm <- rbind(corr,corm)
corall <- cbind(corall,corrm)

write.csv(corall,"04_cv/gSNPcv_all.csv",row.names=F)

}

###intersect significant metabolite SNP and significant genome SNP
mt <- c(1:24,26:27)
nsnpmgall <- NULL
for (n in mt){
	sigmsnp <- read.table(file=paste("eneteffectmetasnpmeta",n,".csv",sep=""),sep=",",header=T)
	siggsnp <- read.table(file=paste("eneteffectgenomesnpmeta",n,".csv",sep=""),sep=",",header=T)
	snpmg <- intersect(siggsnp[,2],sigmsnp[,2])
	nsnpmg <- length(snpmg)
	nsnpmgall  <- rbind(nsnpmgall, nsnpmg)
}

write.csv(nsnpmgall,"gSNPandmSNP.csv",row.names=F)

####################################################################################
### END OF SCRIPTS ###
####################################################################################


