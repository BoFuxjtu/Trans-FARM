rm(list=ls())
library("rstudioapi")
setwd(dirname(getActiveDocumentContext()$path))

library(glmtrans)
library(glmnet)
library(mvtnorm)
library(doParallel)
library(foreach)
library(doSNOW)
#library(clime)
library(SIS)


rcv<-function(Y,X,n){
  Idx = 1:(n/2)
  X0 = X[Idx,]; X2 = X[-Idx,]
  Sigma0 = tcrossprod(X0)*2/n; Sigma2 = tcrossprod(X2)*2/n
  eigen0 = eigen(Sigma0); eigen2 = eigen(Sigma2)
  eigvec0 = eigen0$vectors; eigvalue0 = eigen0$values
  eigvec2 = eigen2$vectors; eigvalue2 = eigen2$values
  K0 = which.min(diff(log(eigvalue0[1:100])))
  K2 = which.min(diff(log(eigvalue2[1:100])))
  F.hat0 = eigvec0[,1:K0]*sqrt(n/2); F.hat2 = eigvec2[,1:K2]*sqrt(n/2)
  B.hat0.T = t(F.hat0)%*%X0*2/n; B.hat2.T = t(F.hat2)%*%X2*2/n
  U.hat0 = X0 - F.hat0%*%B.hat0.T; U.hat2 = X2 - F.hat2%*%B.hat2.T
  Y0 = Y[Idx]; Y2 = Y[-Idx]
  tmp0 = tcrossprod(F.hat0); tmp2 = tcrossprod(F.hat2)
  Y0.new = Y0 - tmp0%*%Y0*2/n
  Y2.new = Y2 - tmp2%*%Y2*2/n
  SIS0 = SIS(U.hat0,Y0.new)
  SIS2 = SIS(U.hat2,Y2.new)
  S.hat0.ISIS = SIS0$ix; S.hat0.SIS = SIS0$ix0
  S.hat2.ISIS = SIS2$ix; S.hat2.SIS = SIS2$ix0 
  X0=cbind(F.hat0,U.hat0)
  X2=cbind(F.hat2,U.hat2)
  X0.ISIS = X0[,c(1:K0,S.hat2.ISIS+K0)]
  X0.SIS = X0[,c(1:K0,K0+S.hat2.SIS)]
  X2.ISIS = X2[,c(1:K2,K2+S.hat0.ISIS)]; X2.SIS = X2[,c(1:K2,K2+S.hat0.SIS)]
  lm0.ISIS = lm(Y0~X0.ISIS-1) 
  lm0.SIS = lm(Y0~X0.SIS-1)
  lm2.ISIS = lm(Y2~X2.ISIS-1) 
  lm2.SIS = lm(Y2~X2.SIS-1)
  Q0.ISIS.H1 = sum((resid(lm0.ISIS))^2)
  sigma.hat0.ISIS = Q0.ISIS.H1/(n/2 - K0 - length(S.hat2.ISIS))
  Q0.SIS.H1 = sum((resid(lm0.SIS))^2)
  sigma.hat0.SIS = Q0.SIS.H1/(n/2 - K0 - length(S.hat2.SIS))
  Q2.ISIS.H1 = sum((resid(lm2.ISIS))^2)
  sigma.hat2.ISIS = Q2.ISIS.H1/(n/2 - K2 - length(S.hat0.ISIS))
  Q2.SIS.H1 = sum((resid(lm2.SIS))^2)
  sigma.hat2.SIS = Q2.SIS.H1/(n/2 - K2 - length(S.hat0.SIS))
  sigma.hat.ISIS = mean(c(sigma.hat0.ISIS, sigma.hat2.ISIS))
  sigma.hat.SIS = mean(c(sigma.hat0.SIS, sigma.hat2.SIS))
  c(sigma.hat.ISIS,sigma.hat.SIS)}

Trans_FARM = function(n0,nk,p,r,K,K_A,s,q,input_type="varphi",w,eta,alpha=0.05){
  k_CI <- p/4
  beta0 = c(rep(q,s),rep(0,p-s))
  #beta0 = c(runif(s,-1,1),rep(0,p-s))
  B0 = matrix(runif(p*r,-1,1),nrow=p)
  F0 = matrix(rnorm(n0*r),nrow=n0)
  meanx=rep(0,p)
  sigma0=diag(p)
  for (i in 1:p){
    for (j in 1:p){sigma0[i,j]=0.5^(abs(i-j))}}
  #U0 = rmvnorm(n0,mean=meanx,sigma=sigma0)
  U0 = matrix(rnorm(n0*p),nrow=n0)
  #U0 = rmvt(n0,delta=meanx,sigma=sigma0,df=10)
  X0 = F0%*%t(B0)+U0
  #E0=rnorm(n0,mean=0,sd=1)
  E0=rt(n0,df=5)
  if(input_type=="gamma"){
    gamma0 = rep(w,r)
    Y0 = U0%*%beta0+F0%*%gamma0+E0
  }else if(input_type=="varphi"){
    varphi0 = rep(w,r)
    Y0 = X0%*%beta0+F0%*%varphi0+E0
  }
  
  Sigma0 = tcrossprod(X0)/n0
  Eig0 = eigen(Sigma0,only.values=TRUE)
  Eigval0 = Eig0$values
  r.est0 = which.min(diff(log(Eigval0[1:20])))
  Svd0 = svd(X0,nu=r.est0,nv=0)
  Eigvec0 = Svd0$u
  F.hat0 = sqrt(n0)*Eigvec0[,1:r.est0]
  mol0 = lm(X0~F.hat0-1)
  U.hat0 = resid(mol0)
  Y.hat0 = (diag(n0)-F.hat0%*%t(F.hat0)/n0)%*%Y0
  
  A_eta <- sample(1:K,size=K_A)
  for (k in 1:K){
    Bk = matrix(runif(p*r,-1,1),nrow=p)
    Fk = matrix(rnorm(nk*r),nrow=nk)
    epsk = rmvnorm(nk,mean=meanx,sigma=0.2^2*diag(p))
    sigmak = sigma0+t(epsk)%*%epsk
    Uk = rmvnorm(nk,mean=meanx,sigma=sigmak)
    #Uk = rmvt(nk,delta=meanx,sigma=sigmak,df=10)
    Xk = Fk%*%t(Bk)+Uk
    Ek=rt(nk,df=5)
    #Ek=rnorm(nk,mean=0,sd=1)
    Rade1 = 2*rbinom(p,1,0.5)-1
    Rade2 = 2*rbinom(r,1,0.5)-1
    betak = beta0
    if(k %in% A_eta){
      betak = beta0+(eta/p)*Rade1
    }else{
      betak = beta0+(2*eta/p)*Rade1
      S_k <- sample((2*s+1):p, s)
      S_k <- c(S_k,(s+1):(2*s))
      betak[S_k] = betak[S_k]+0.5
    }
    if(input_type=="gamma"){
      if(k %in% A_eta){
        gammak = gamma0+0.1*Rade2
      }else{
        gammak = gamma0+0.5*Rade2
      }
      Yk = Uk%*%betak+Fk%*%gammak+Ek
    }else if(input_type=="varphi"){
      if(k %in% A_eta){
        varphik = varphi0+0.1*Rade2
      }else{
        varphik = varphi0+0.5*Rade2
      }
      Yk = Xk%*%betak+Fk%*%varphik+Ek
    }
    
    Sigmak = tcrossprod(Xk)/nk
    Eigk = eigen(Sigmak,only.values=TRUE)
    Eigvalk = Eigk$values
    r.estk = which.min(diff(log(Eigvalk[1:20])))
    Svdk = svd(Xk,nu=r.estk,nv=0)
    Eigveck = Svdk$u
    F.hatk = sqrt(nk)*Eigveck[,1:r.estk]
    molk = lm(Xk~F.hatk-1)
    U.hatk = resid(molk)
    Y.hatk = (diag(nk)-F.hatk%*%t(F.hatk)/nk)%*%Yk
    
    assign(paste0("U.hat", k), U.hatk)
    assign(paste0("Y.hat", k), Y.hatk)
    assign(paste0("X", k), Xk)
    assign(paste0("Y", k), Yk)
  }
  sigma_hat <- rcv(Y0,X0,n0)[1]
  sigma_hat <- sqrt(sigma_hat)
  
  #-------------------Use Nodewise regression to estimate Theta, the invese matrix of the population covariance------
  C<-as.matrix(diag(rep(1,p)))
  C <- C[(1:k_CI),]
  T<-c()
  for( j in 1:k_CI){
    fit_u = glmnet(U.hat0[,-j], U.hat0[,j], intercept=FALSE,
                   lambda=cv.glmnet(U.hat0[,-j], U.hat0[,j],intercept=FALSE)$lambda.1se)
    beta<-as.vector(fit_u$beta)
    C[j,-j]<--beta
    T<-c(T,1/n0*sum((U.hat0[,j]-U.hat0[,-j]%*%beta)^2)+fit_u$lambda/2*sum(abs(beta)))
    #if (j%%10==0){
    #  print(j)
    #}
  }
  T1<-diag(1/T)
  Theta_U_hat<-T1%*%C #estimated Theta
  std = sqrt(diag(Theta_U_hat))
  
  fit.only_FARM=glmnet(U.hat0, Y.hat0,
                       lambda=cv.glmnet(U.hat0, Y.hat0)$lambda.1se)  #Estimate Lasso model
  beta.only_FARM <- as.vector(c(fit.only_FARM$beta[1:p],fit.only_FARM$a0))
  REG.only_FARM <- cbind(U.hat0,1)
  beta_debiased_only_FARM <- beta.only_FARM[1:k_CI]+Theta_U_hat%*%t(U.hat0)%*%(Y.hat0-REG.only_FARM%*%beta.only_FARM)/n0
  
  c_n.only_FARM = sapply(1:500,function(xxx,n0,sigma_hat,Theta_U_hat,U.hat0,Y.hat0){
    ei = rnorm(n0)
    Q_e.hat = Theta_U_hat%*% apply(U.hat0*ei, 2, mean)*sigma_hat
    boots=max(abs(Q_e.hat))
    return(boots)
  },n0,sigma_hat,Theta_U_hat,U.hat0,Y.hat0)
  quantiles.only_FARM <- quantile(c_n.only_FARM, probs = 1 - alpha)
  cp.only_FARM <- max(abs(beta0[1:k_CI]-beta_debiased_only_FARM))<=quantiles.only_FARM
  il.only_FARM <- 2*quantiles.only_FARM
  
  
  source.Trans_FARM <- lapply(1:K, function(k) {
    list(
      x = get(paste0("U.hat", k)),
      y = get(paste0("Y.hat", k))
    )
  })
  target.Trans_FARM <- list(x = get(paste0("U.hat", 0)),y = get(paste0("Y.hat", 0)))
  fit.Trans_FARM=glmtrans(target=target.Trans_FARM,source=source.Trans_FARM)
  beta.Trans_FARM<-(as.vector(fit.Trans_FARM$beta))
  REG.Trans_FARM <- cbind(1,U.hat0)
  beta_debiased_Trans_FARM <- beta.Trans_FARM[2:(k_CI+1)]+Theta_U_hat%*%t(U.hat0)%*%(Y.hat0-REG.Trans_FARM%*%beta.Trans_FARM)/n0
  
  
  c_n.Trans_FARM = sapply(1:500,function(xxx,n0,Theta_U_hat,U.hat0,sigma_hat,std){
    ei = rnorm(n0)
    Q_ne.hat = diag(std^(-1)) %*% Theta_U_hat %*% apply(U.hat0*ei, 2, mean) * sigma_hat
    boots=max(abs(Q_ne.hat))
    return(boots)
  },n0,Theta_U_hat,U.hat0,sigma_hat,std)
  quantiles.Trans_FARM <- quantile(c_n.Trans_FARM, probs = 1 - alpha)*std
  cp.Trans_FARM <- all(abs(beta0[1:k_CI]-beta_debiased_Trans_FARM)<=quantiles.Trans_FARM)
  il.Trans_FARM <- 2*mean(quantiles.Trans_FARM)
  
  list(cp1=cp.only_FARM,il1=il.only_FARM,cp2=cp.Trans_FARM,il2=il.Trans_FARM)
}
#-----------------------------------------

#-----------------------------------------
run_Trans_FARM_parallel <- function(n0=500,nk=500,p=1000,r=2,K=10,K_A,s=3,q=q,input_type="varphi",w=0,eta=40,itr.max=500,ncores=16,alpha=0.05){
  cl <- makeCluster(ncores)
  registerDoSNOW(cl)
  pb <- txtProgressBar(max = itr.max, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  results <- foreach(itr=1:itr.max, .combine=rbind,
                     .packages=c("glmnet","conquer","mvtnorm","glmtrans","clime","SIS"),
                     .export=c("Trans_FARM","rcv"),
                     .options.snow = opts) %dopar% {
                       Trans_FARM.est <- Trans_FARM(n0=n0,nk=nk,p=p,r=r,K=K,K_A=K_A,s=s,q=q,input_type=input_type,w=w,eta=eta,alpha=alpha)
                       cp1 <- Trans_FARM.est$cp1
                       cp2 <- Trans_FARM.est$cp2
                       il1 <- Trans_FARM.est$il1
                       il2 <- Trans_FARM.est$il2
                       c(cp1,cp2,il1,il2)
                     }
  
  stopCluster(cl)
  close(pb)
  
  cp1 <- mean(results[,1])
  cp2 <- mean(results[,2])
  il1 <- mean(results[,3])
  il2 <- mean(results[,4])
  
  cat("n0=", n0, "nk=", nk, "p=", p, "s=", s, ",w=", w, "r=", r, "K=", K, "K_A=", K_A, "eta=", eta, "q=", q, "Normal U; t5 error", "\n",
      "cp1 =", cp1, "\n",
      "cp2 =", cp2, "\n",
      "il1 =", il1, "\n",
      "il2 =", il2, "\n")
  
  return(list(cp1=cp1,cp2=cp2,il1=il1,il2=il2))
}

#-----------------------------------------
# run
#-----------------------------------------
n0 <- 500
nk <- 500
p <- 500
r <- 2
K <- 10
K_A <- 5
s <- 5
itr.max <- 500
ncores <- 20
#K_A_values <- 0:10
#eta_values <- c(5, 10)
eta <- 5
input_type <- "gamma"
w <- 0.5
q <- 1
results_eta5_n500_p500 <- run_Trans_FARM_parallel(n0 = n0, nk = nk, p = p, r = r, K = K, K_A = K_A,
                                                            s = s, q=q, input_type = input_type, w = w, eta = eta,
                                                            itr.max = itr.max, ncores = ncores)

results_eta5_n200_p500 <- run_Trans_FARM_parallel(n0 = 200, nk = nk, p = p, r = r, K = K, K_A = K_A,
                                                        s = s, q=q, input_type = input_type, w = w, eta = eta,
                                                        itr.max = itr.max, ncores = ncores)

results_eta5_n200_p200 <- run_Trans_FARM_parallel(n0 = 200, nk = nk, p = 200, r = r, K = K, K_A = K_A,
                                                  s = s, q=q, input_type = input_type, w = w, eta = eta,
                                                  itr.max = itr.max, ncores = ncores)

results_eta5_n500_p200 <- run_Trans_FARM_parallel(n0 = n0, nk = nk, p = 200, r = r, K = K, K_A = K_A,
                                                  s = s, q=q, input_type = input_type, w = w, eta = eta,
                                                  itr.max = itr.max, ncores = ncores)
