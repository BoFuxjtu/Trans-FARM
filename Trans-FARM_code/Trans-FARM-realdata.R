rm(list=ls())
library("rstudioapi")
setwd(dirname(getActiveDocumentContext()$path))
###Read in Data####
#D = read.csv('realdata_189_126.csv')#using part of 1992-2007 as target
D = read.csv('realdata_115_126.csv')#using part of 2010-2020 as target
D = D[,-1]
D<-as.matrix(data.frame(D))
D=D[,-c(81,83)]
D2 = read.csv('realdata_189_126.csv')#source data
D2 = D2[,-1]
D2<-as.matrix(data.frame(D2))
D2=D2[,-c(81,83)]
##-----------------------------------------------------------------
#Loading required packages
##------------------------------------------------------------------
library(glmnet)
library(glmtrans)
set.seed(100)
Tn = nrow(D)
N <- ncol(D)
j = 81 #set j=81 is the prediction result for GS5 / j=49 for HOUSTNE
X = D[,-j] #Use rest of the data to be covariate
Y = D[,j] #Response variable
#val.idx = (sample(1:Tn, size = round(0.1*Tn)))
#train.idx = setdiff(1:Tn, val.idx)
#X.val = X[val.idx,]
#Y.val = Y[val.idx]
#X.train.full = X[train.idx,]
#Y.train.full = Y[train.idx]
#source.ratio = 0.5
#n.train.full = length(train.idx)
#n.source = round(source.ratio * n.train.full)
#source.idx = sample(train.idx, n.source)

T = 90  ## moving window approach, window size                         
M = Tn-T  ## predict sample size 
n2 = nrow(D2)
#---------------------------------------------------------------------------------------
#########################################################################################################
#### The following values are initialized for storing prediction difference, R^2, predicted values, obtained from different methods
## like FARM, FLasso, Lasso, Sample Mean, PCR, Ridge Elastic Net, Random Forest, respectively. Details are in Section 5.4 of the paper.##########
#########################################################################################################
#-----------------------------------
#Store out-of-sample R^2 for several methods mentioned above
#------------------------------------
R2.FARM = numeric(N)           #FARM  
R2.LASSO = numeric(N)          #LASSO
R2.Trans_FARM = numeric(N)           #Trans_FARM  
R2.Trans_LASSO = numeric(N)          #Trans_LASSO

Pred.FARM = matrix(0,M,N)                       ## FARM for estimation and prediction
Pred.Trans_FARM = matrix(0,M,N)                  
Pred.Lasso = matrix(0,M,N)                      ## Lasso for prediction
Pred.Trans_Lasso = matrix(0,M,N)
Pred.MEAN = matrix(0,M,N)
#----------------------------------------------
#Predictions is conducted Next
#-----------------------------------------------
for(i in 1:M){
  idx = i:(i+T-1) #Moving window prediction, every window has length T
  x = X[idx,] #Training Data Covariate
  y = Y[idx]  #Training Data Response
  x.new = X[i+T,] #Prediction Covariate
  y.new = Y[i+T] #Prediction Response
X.mu = colMeans(x)
X.sd = as.numeric(apply(x,2,sd))
Y.mu = mean(y)
Y.sd = sd(y)
X.target = t((t(x)-X.mu)/X.sd)                               ## Data normalization, we standardize every column of X to be mean zero and sd 1
Y.target = (y-Y.mu)/Y.sd
X.new = (x.new-X.mu)/X.sd
Y.new = (y.new-Y.mu)/Y.sd

source.idx = setdiff((1:Tn), (i:i+T))
X.source1 = X[source.idx,]
Y.source1 = Y[source.idx]
X.source2 = D2[ ,-j]
Y.source2 = D2[ ,j]
#half = floor(n2/2)
#X.source2 = D2[1:half,-j]
#Y.source2 = D2[1:half,j]
#X.source3 = D2[(half+1):n2,-j]
#Y.source3 = D2[(half+1):n2,j]
X.source1 = t((t(X.source1)-X.mu)/X.sd)                               ## Data normalization, we standardize every column of X to be mean zero and sd 1
Y.source1 = (Y.source1-Y.mu)/Y.sd
X.source2 = t((t(X.source2)-colMeans(X.source2))/X.sd)                               ## Data normalization, we standardize every column of X to be mean zero and sd 1
Y.source2= (Y.source2-mean(Y.source2))/Y.sd
#X.source3 = t((t(X.source3)-colMeans(X.source3))/X.sd)                               ## Data normalization, we standardize every column of X to be mean zero and sd 1
#Y.source3 = (Y.source3-mean(Y.source2))/Y.sd
##-------------------------------------------------------------------
###In the following, we conduct several prediction comparison using sample Mean, Lasso, Ridge, Elastic Net, FARM, PCR, FLasso
##--------------------------------------------------------------------
Pred.MEAN[i,j] = (Y.new* Y.sd)^2   #Sample Mean Prediction
##--------------------------------------------------------------------
##Lasso
##--------------------------------------------------------------------
cv.fit.x = cv.glmnet(X.target,Y.target,intercept=FALSE)	
lambda.fit.x = cv.fit.x$lambda.min
fit.x = glmnet(X.target,Y.target,intercept=FALSE,lambda=lambda.fit.x)  ## Lasso Estimation
beta.hat.x = as.vector(fit.x$beta)
Pred.Lasso[i,j] = (Y.new - X.new %*% beta.hat.x)^2 * Y.sd^2
##--------------------------------------------------------------------
##Trans-Lasso
##--------------------------------------------------------------------
source.Trans_Lasso <- lapply(1:2, function(k) {
  list(
    x = get(paste0("X.source", k)),
    y = get(paste0("Y.source", k))
  )
})
target.Trans_Lasso <- list(x = X.target,y = Y.target)
fit.Trans_Lasso=glmtrans(target=target.Trans_Lasso,source=source.Trans_Lasso,detection.info=FALSE)
beta.Trans_Lasso<-as.vector(fit.Trans_Lasso$beta)[-1]
Pred.Trans_Lasso[i,j] = (Y.new - X.new %*% beta.Trans_Lasso)^2 * Y.sd^2
## --------------------------------------------------------------------
##Factor Estimation
## --------------------------------------------------------------------
Sigma.x = tcrossprod(X.target)/nrow(X.target)      #covariance matrix
eigenx = eigen(Sigma.x)        
eigvec = eigenx$vectors        #Eigen vector
eigvalue = eigenx$values       #Eigen values
K.hat = max(1,which.min(diff(log(eigvalue[1:10]))))      ## Factor estimation
F.hat = eigvec[,1:K.hat]*sqrt(nrow(X.target))        #Estimated Factor
B.hat = t(t(F.hat)%*%X.target)/nrow(X.target)          #Estimated Factor Loading
U.hat = X.target-F.hat%*%t(B.hat)        #Estimated idiosyncratic component
##--------------------------------------------------------------------
#FARM
##--------------------------------------------------------------------
lmY.F = lm(Y.target~F.hat-1)           
gamma.hat = coef(lmY.F)    #Estimate \gamma vector in the FARM model and PCR
Y.tilde = resid(lmY.F)                             
cv.fit.U = cv.glmnet(U.hat,Y.tilde,intercept=FALSE)	 
lambda.fit.U = cv.fit.U$lambda.min #tuning parameter selection
fit.U = glmnet(U.hat,Y.tilde,intercept=FALSE,lambda=lambda.fit.U) ##FARM Estimation
beta.hat.U = as.vector(fit.U$beta) #Obtain the fitted beta for FARM
lmx.B = lm(X.new~B.hat-1)
F.new = coef(lmx.B)
U.new = resid(lmx.B)
Pred.FARM[i,j] = (Y.new - F.new %*% gamma.hat - U.new %*% beta.hat.U)^2 * Y.sd^2
##--------------------------------------------------------------------
#Trans-FARM
##--------------------------------------------------------------------
for (k in 1:2){
  x = get(paste0("X.source", k))
  y = get(paste0("Y.source", k))
  Sigma.x = tcrossprod(x)/nrow(x)      #covariance matrix
  eigenx = eigen(Sigma.x)        
  eigvec = eigenx$vectors        #Eigen vector
  eigvalue = eigenx$values       #Eigen values
  K.hatk = max(1,which.min(diff(log(eigvalue[1:10]))))      ## Factor estimation
  F.hatk = eigvec[,1:K.hatk]*sqrt(nrow(x))        #Estimated Factor
  B.hatk = t(t(F.hatk)%*%x)/nrow(x)          #Estimated Factor Loading
  U.hatk = x-F.hatk%*%t(B.hatk)        #Estimated idiosyncratic component
  Y.hatk = (diag(nrow(x))-F.hatk%*%t(F.hatk)/nrow(x))%*%y
  assign(paste0("U.hat", k), U.hatk)
  assign(paste0("Y.hat", k), Y.hatk)
}
source.Trans_FARM <- lapply(1:2, function(k) {
  list(
    x = get(paste0("U.hat", k)),
    y = get(paste0("Y.hat", k))
  )
})
target.Trans_FARM <- list(x = U.hat,y = Y.tilde)
fit.Trans_FARM=glmtrans(target=target.Trans_FARM,source=source.Trans_FARM,detection.info=FALSE)
beta.Trans_FARM<-(as.vector(fit.Trans_FARM$beta))[-1]
Pred.Trans_FARM[i,j] = (Y.new - F.new %*% gamma.hat - U.new %*% beta.Trans_FARM)^2 * Y.sd^2
}
#————————————————————————————————————————————————————————————
##Output R^2 for Table 3
#------------------------------------------------------------
R2.LASSO[j] = 1-sum(Pred.Lasso[,j])/sum(Pred.MEAN[,j])
print(R2.LASSO[j])                 #Output R^2 value for LASSO
R2.Trans_LASSO[j] = 1-sum(Pred.Trans_Lasso[,j])/sum(Pred.MEAN[,j])
print(R2.Trans_LASSO[j]) 
R2.FARM[j] = 1-sum(Pred.FARM[,j])/sum(Pred.MEAN[,j])	
print(R2.FARM[j])                    #Output R^2 value for FARM for PCR
R2.Trans_FARM[j] = 1-sum(Pred.Trans_FARM[,j])/sum(Pred.MEAN[,j])	
print(R2.Trans_FARM[j])    