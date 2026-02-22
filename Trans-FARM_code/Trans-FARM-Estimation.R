rm(list=ls())
library("rstudioapi")
setwd(dirname(getActiveDocumentContext()$path))

library(glmtrans)
library(glmnet)
library(mvtnorm)
library(doParallel)
library(foreach)
library(doSNOW)


Trans_FARM = function(n0,nk,p,r,K,K_A,s,input_type="varphi",w,eta){
  beta0 = c(rep(0.5,s),rep(0,p-s))
  #beta0 = c(runif(s,-1,1),rep(0,p-s))
  B0 = matrix(runif(p*r,-1,1),nrow=p)
  F0 = matrix(rnorm(n0*r),nrow=n0)
  meanx=rep(0,p)
  sigma0=diag(p)
  for (i in 1:p){
    for (j in 1:p){sigma0[i,j]=0.5^(abs(i-j))}}
  #U0 = rmvnorm(n0,mean=meanx,sigma=sigma0)
  U0 = rmvt(n0,delta=meanx,sigma=sigma0,df=10)
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
    #Uk = rmvnorm(nk,mean=meanx,sigma=sigmak)
    Uk = rmvt(nk,delta=meanx,sigma=sigmak,df=10)
    Xk = Fk%*%t(Bk)+Uk
    #Ek=rnorm(nk,mean=0,sd=1)
    Ek=rt(nk,df=5)
    Rade1 = 2*rbinom(p,1,0.5)-1
    Rade2 = 2*rbinom(r,1,0.5)-1
    if(k %in% A_eta){
      betak = beta0+(eta/p)*Rade1
    }else{
      betak = beta0+(2*eta/p)*Rade1
      S_k <- sample((2*s+1):p, s)
      S_k <- c(S_k,(s+1):2*s)
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
        varphik = varphi0+0*Rade2
      }else{
        varphik = varphi0+0*Rade2
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
  
  fit.only_FARM=glmnet(U.hat0, Y.hat0,intercept=FALSE,
                       lambda=cv.glmnet(U.hat0, Y.hat0,intercept=FALSE)$lambda.1se)  #Estimate Lasso model
  beta.only_FARM <- (as.vector(fit.only_FARM$beta))
  
  fit.only_Lasso=glmnet(X0, Y0,intercept=FALSE,lambda=cv.glmnet(X0, Y0,intercept=FALSE)$lambda.1se)  #Estimate Lasso model
  beta.only_Lasso<-as.vector(fit.only_Lasso$beta)
  
  source.Trans_FARM <- lapply(1:K, function(k) {
    list(
      x = get(paste0("U.hat", k)),
      y = get(paste0("Y.hat", k))
    )
  })
  target.Trans_FARM <- list(x = get(paste0("U.hat", 0)),y = get(paste0("Y.hat", 0)))
  fit.Trans_FARM=glmtrans(target=target.Trans_FARM,source=source.Trans_FARM)
  beta.Trans_FARM<-(as.vector(fit.Trans_FARM$beta))[-1]
  
  source.Trans_Lasso <- lapply(1:K, function(k) {
    list(
      x = get(paste0("X", k)),
      y = get(paste0("Y", k))
    )
  })
  target.Trans_Lasso <- list(x = get(paste0("X", 0)),y = get(paste0("Y", 0)))
  fit.Trans_Lasso=glmtrans(target=target.Trans_Lasso,source=source.Trans_Lasso)
  beta.Trans_Lasso<-(as.vector(fit.Trans_Lasso$beta))[-1]
  
  fit.Oracle_Trans_FARM=glmtrans(target=target.Trans_FARM,source=source.Trans_FARM,transfer.source.id = A_eta)  #Estimate Lasso model
  beta.Oracle_Trans_FARM<-(as.vector(fit.Oracle_Trans_FARM$beta))[-1]
  
  
  fit.Oracle_Trans_Lasso=glmtrans(target=target.Trans_Lasso,source=source.Trans_Lasso,transfer.source.id = A_eta)  #Estimate Lasso model
  beta.Oracle_Trans_Lasso<-(as.vector(fit.Oracle_Trans_Lasso$beta))[-1]
  
  
  fit.Pooled_Trans_FARM=glmtrans(target=target.Trans_FARM,source=source.Trans_FARM,transfer.source.id = "all")  #Estimate Lasso model
  beta.Pooled_Trans_FARM<-(as.vector(fit.Pooled_Trans_FARM$beta))[-1]
  
  
  fit.Pooled_Trans_Lasso=glmtrans(target=target.Trans_Lasso,source=source.Trans_Lasso,transfer.source.id = "all")  #Estimate Lasso model
  beta.Pooled_Trans_Lasso<-(as.vector(fit.Pooled_Trans_Lasso$beta))[-1]
  
  
  error.only_FARM = norm(beta.only_FARM-beta0,"2")
  error.only_Lasso = norm(beta.only_Lasso-beta0,"2")
  error.Trans_FARM = norm(beta.Trans_FARM-beta0,"2")
  error.Trans_Lasso = norm(beta.Trans_Lasso-beta0,"2")
  error.Oracle_Trans_FARM = norm(beta.Oracle_Trans_FARM-beta0,"2")
  error.Oracle_Trans_Lasso = norm(beta.Oracle_Trans_Lasso-beta0,"2")
  error.Pooled_Trans_FARM = norm(beta.Pooled_Trans_FARM-beta0,"2")
  error.Pooled_Trans_Lasso = norm(beta.Pooled_Trans_Lasso-beta0,"2")
  
  
  
  list(error.only_FARM=error.only_FARM,error.only_Lasso=error.only_Lasso,error.Trans_FARM=error.Trans_FARM,
       error.Trans_Lasso=error.Trans_Lasso,error.Oracle_Trans_FARM=error.Oracle_Trans_FARM,
       error.Oracle_Trans_Lasso=error.Oracle_Trans_Lasso,error.Pooled_Trans_FARM=error.Pooled_Trans_FARM,
       error.Pooled_Trans_Lasso=error.Pooled_Trans_Lasso)
}

#-----------------------------------------
# parallel
#-----------------------------------------
run_Trans_FARM_parallel <- function(n0=500,nk=500,p=1000,r=2,K=10,K_A,s=3,input_type="varphi",w=0,eta=40,itr.max=500,ncores=16){
  cl <- makeCluster(ncores)
  registerDoSNOW(cl)
  pb <- txtProgressBar(max = itr.max, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  results <- foreach(itr=1:itr.max, .combine=rbind,
                     .packages=c("glmnet","conquer","mvtnorm","glmtrans"),
                     .export=c("Trans_FARM"),
                     .options.snow = opts) %dopar% {
                       cat("Iteration:", itr, "\n")
                       Trans_FARM.est <- Trans_FARM(n0=n0,nk=nk,p=p,r=r,K=K,K_A=K_A,s=s,input_type=input_type,w=0,eta=eta)
                       error.only_FARM <- Trans_FARM.est$error.only_FARM
                       error.only_Lasso <- Trans_FARM.est$error.only_Lasso
                       error.Trans_FARM <- Trans_FARM.est$error.Trans_FARM
                       error.Trans_Lasso <- Trans_FARM.est$error.Trans_Lasso
                       error.Oracle_Trans_FARM <- Trans_FARM.est$error.Oracle_Trans_FARM
                       error.Oracle_Trans_Lasso <- Trans_FARM.est$error.Oracle_Trans_Lasso
                       error.Pooled_Trans_FARM <- Trans_FARM.est$error.Pooled_Trans_FARM
                       error.Pooled_Trans_Lasso <- Trans_FARM.est$error.Pooled_Trans_Lasso
                       c(error.only_FARM,error.only_Lasso,error.Trans_FARM,error.Trans_Lasso,
                         error.Oracle_Trans_FARM,error.Oracle_Trans_Lasso,error.Pooled_Trans_FARM,
                         error.Pooled_Trans_Lasso)
                     }
  
  stopCluster(cl)
  close(pb)
  
  error.only_FARM <- mean(results[,1])
  error.only_Lasso <- mean(results[,2])
  error.Trans_FARM <- mean(results[,3])
  error.Trans_Lasso <- mean(results[,4])
  error.Oracle_Trans_FARM <- mean(results[,5])
  error.Oracle_Trans_Lasso <- mean(results[,6])
  error.Pooled_Trans_FARM <- mean(results[,7])
  error.Pooled_Trans_Lasso <- mean(results[,8])
  
  cat("n0=", n0, "nk=", nk, "p=", p, "s=", s, ",w=", w, "r=", r, "K=", K, "K_A=", K_A, "eta=", eta, "\n",
      "error.only_FARM =", error.only_FARM, "\n",
      "error.only_Lasso =", error.only_Lasso, "\n",
      "error.Trans_FARM =", error.Trans_FARM, "\n",
      "error.Trans_Lasso =", error.Trans_Lasso, "\n",
      "error.Oracle_Trans_FARM =", error.Oracle_Trans_FARM, "\n",
      "error.Oracle_Trans_Lasso =", error.Oracle_Trans_Lasso, "\n",
      "error.Pooled_Trans_FARM =", error.Pooled_Trans_FARM, "\n",
      "error.Pooled_Trans_Lasso =", error.Pooled_Trans_Lasso, "\n")
  
  return(list(error.only_FARM=error.only_FARM, error.only_Lasso=error.only_Lasso, error.Trans_FARM=error.Trans_FARM,error.Trans_Lasso=error.Trans_Lasso,
              error.Oracle_Trans_FARM=error.Oracle_Trans_FARM,error.Oracle_Trans_Lasso=error.Oracle_Trans_Lasso,
              error.Pooled_Trans_FARM=error.Pooled_Trans_FARM,error.Pooled_Trans_Lasso=error.Pooled_Trans_Lasso))
}

#-----------------------------------------
# run
#-----------------------------------------
n0 <- 300
nk <- 300
p <- 500
r <- 2
K <- 10
s <- 20
itr.max <- 200
ncores <- 24
K_A_values <- 0:10
eta_values <- c(5, 10)
input_type <- "gamma"
w <- 0.5
results_l2 <- list()

for (eta in eta_values) {
  results_l2[[paste0("eta", eta)]] <- lapply(K_A_values, function(KA) {
    run_Trans_FARM_parallel(n0 = n0, nk = nk, p = p, r = r, K = K, K_A = KA,
                            s = s, input_type = input_type, w = w, eta = eta,
                            itr.max = itr.max, ncores = ncores)
  })
}


error.only_FARM_l2_eta5 <- sapply(results_l2[["eta5"]], function(res) res$error.only_FARM)
error.only_Lasso_l2_eta5 <- sapply(results_l2[["eta5"]], function(res) res$error.only_Lasso)
error.Trans_FARM_l2_eta5 <- sapply(results_l2[["eta5"]], function(res) res$error.Trans_FARM)
error.Trans_Lasso_l2_eta5 <- sapply(results_l2[["eta5"]], function(res) res$error.Trans_Lasso)
error.Oracle_Trans_FARM_l2_eta5 <- sapply(results_l2[["eta5"]], function(res) res$error.Oracle_Trans_FARM)
error.Oracle_Trans_Lasso_l2_eta5 <- sapply(results_l2[["eta5"]], function(res) res$error.Oracle_Trans_Lasso)
error.Pooled_Trans_FARM_l2_eta5 <- sapply(results_l2[["eta5"]], function(res) res$error.Pooled_Trans_FARM)
error.Pooled_Trans_Lasso_l2_eta5 <- sapply(results_l2[["eta5"]], function(res) res$error.Pooled_Trans_Lasso)

error.only_FARM_l2_eta10 <- sapply(results_l2[["eta10"]], function(res) res$error.only_FARM)
error.only_Lasso_l2_eta10 <- sapply(results_l2[["eta10"]], function(res) res$error.only_Lasso)
error.Trans_FARM_l2_eta10 <- sapply(results_l2[["eta10"]], function(res) res$error.Trans_FARM)
error.Trans_Lasso_l2_eta10 <- sapply(results_l2[["eta10"]], function(res) res$error.Trans_Lasso)
error.Oracle_Trans_FARM_l2_eta10 <- sapply(results_l2[["eta10"]], function(res) res$error.Oracle_Trans_FARM)
error.Oracle_Trans_Lasso_l2_eta10 <- sapply(results_l2[["eta10"]], function(res) res$error.Oracle_Trans_Lasso)
error.Pooled_Trans_FARM_l2_eta10 <- sapply(results_l2[["eta10"]], function(res) res$error.Pooled_Trans_FARM)
error.Pooled_Trans_Lasso_l2_eta10 <- sapply(results_l2[["eta10"]], function(res) res$error.Pooled_Trans_Lasso)

par(mar = c(5, 4, 3, 2))
labs <- c("only_FARM","only_Lasso","Trans_FARM","Trans_Lasso",
          "Oracle_Trans_FARM","Oracle_Trans_Lasso",
          "Pooled_Trans_FARM","Pooled_Trans_Lasso")

plot(K_A_values, error.only_FARM_l2_eta5, type="o", col="blue", pch=16,
     ylim = range(c(error.only_FARM_l2_eta5, error.only_Lasso_l2_eta5,
                    error.Trans_FARM_l2_eta5, error.Trans_Lasso_l2_eta5,
                    error.Oracle_Trans_FARM_l2_eta5, error.Oracle_Trans_Lasso_l2_eta5,
                    error.Pooled_Trans_FARM_l2_eta5, error.Pooled_Trans_Lasso_l2_eta5)),
     xlab = "",  
     ylab = expression(paste(l[2], " error")),
     #main = expression(eta == 5 ~ "," ~U[k] %~% N(0[p], Sigma[k]) ~ "," ~ E[k] %~% t[5]),mgp = c(2, 0.5, 0))
     main = expression(eta == 5 ~ "," ~U[k] %~% t[10](0[p], Sigma[k]) ~ "," ~ E[k] %~% N(0,1)),mgp = c(2, 0.5, 0))
     #main = expression(eta == 5 ~ "," ~U[k] %~% N(0[p], Sigma[k]) ~ "," ~ E[k] %~% N(0,1)),mgp = c(2, 0.5, 0))
     #main = expression(eta == 5 ~ "," ~U[k] %~% t[10](0[p], Sigma[k]) ~ "," ~ E[k] %~% t[5]),mgp = c(2, 0.5, 0))
mtext(expression("|" * A[eta] * "|"), side = 1, line = 1.5, cex = 1)
lines(K_A_values,error.only_Lasso_l2_eta5, type="o", col="blue", pch=12)
lines(K_A_values,error.Trans_FARM_l2_eta5, type="o", col="red", pch=16)
lines(K_A_values,error.Trans_Lasso_l2_eta5, type="o",  col="red", pch=12)
lines(K_A_values,error.Oracle_Trans_FARM_l2_eta5, type="o",  col="purple", pch=16)
lines(K_A_values,error.Oracle_Trans_Lasso_l2_eta5, type="o",  col="purple", pch=12)
lines(K_A_values,error.Pooled_Trans_FARM_l2_eta5, type="o",  col="orange", pch=16)
lines(K_A_values,error.Pooled_Trans_Lasso_l2_eta5, type="o",  col="orange", pch=12)
usr <- par("usr") 
x_left <- usr[1] - 0.10 * diff(usr[1:2])
y_bottom <- usr[3] - 0.18 * diff(usr[3:4])
text(x_left, y_bottom + 0.00 * diff(usr[3:4]), "Method:", adj = 0, cex = 0.9, xpd = TRUE)
maxwidth <- max(strwidth(labs, units = "user"))
legend("bottomleft",
       inset = c(-0.01, -0.25), 
       legend = labs,
       text.width = maxwidth * 1,
       col = c("blue","blue","red","red","purple","purple","orange","orange"),
       pch = c(16,12,16,12,16,12,16,12),
       lty = 1,
       ncol = 4,              
       horiz = FALSE,        
       xpd = TRUE,
       cex = 0.8,
       x.intersp = 0.6,
       bty = "n")
labs <- c("only_FARM","only_Lasso","Trans_FARM","Trans_Lasso",
          "Oracle_Trans_FARM","Oracle_Trans_Lasso",
          "Pooled_Trans_FARM","Pooled_Trans_Lasso")


plot(K_A_values, error.only_FARM_l2_eta10, type="o", col="blue", pch=16, ylim=range(c(error.only_FARM_l2_eta10,error.only_Lasso_l2_eta10,
                                                                                      error.Trans_FARM_l2_eta10,error.Trans_Lasso_l2_eta10,
                                                                                      error.Oracle_Trans_FARM_l2_eta10,error.Oracle_Trans_Lasso_l2_eta10,
                                                                                      error.Pooled_Trans_FARM_l2_eta10,error.Pooled_Trans_Lasso_l2_eta10)),
     xlab = "",  
     ylab = expression(paste(l[2], " error")),
     #main = expression(eta == 10 ~ "," ~U[k] %~% N(0[p], Sigma[k]) ~ "," ~ E[k] %~% t[5]),mgp = c(2, 0.5, 0))
     main = expression(eta == 10 ~ "," ~U[k] %~% t[10](0[p], Sigma[k]) ~ "," ~ E[k] %~% N(0,1)),mgp = c(2, 0.5, 0))
     #main = expression(eta == 10 ~ "," ~U[k] %~% N(0[p], Sigma[k]) ~ "," ~ E[k] %~% N(0,1)),mgp = c(2, 0.5, 0))
     #main = expression(eta == 10 ~ "," ~U[k] %~% t[10](0[p], Sigma[k]) ~ "," ~ E[k] %~% t[5]),mgp = c(2, 0.5, 0))
mtext(expression("|" * A[eta] * "|"), side = 1, line = 1.5, cex = 1)
lines(K_A_values,error.only_Lasso_l2_eta10, type="o", col="blue", pch=12)
lines(K_A_values,error.Trans_FARM_l2_eta10, type="o", col="red", pch=16)
lines(K_A_values,error.Trans_Lasso_l2_eta10, type="o",  col="red", pch=12)
lines(K_A_values,error.Oracle_Trans_FARM_l2_eta10, type="o",  col="purple", pch=16)
lines(K_A_values,error.Oracle_Trans_Lasso_l2_eta10, type="o",  col="purple", pch=12)
lines(K_A_values,error.Pooled_Trans_FARM_l2_eta10, type="o",  col="orange", pch=16)
lines(K_A_values,error.Pooled_Trans_Lasso_l2_eta10, type="o",  col="orange", pch=12)

text(x_left, y_bottom + 0.03 * diff(usr[3:4]), "Method:", adj = 0, cex = 0.9, xpd = TRUE)
legend("bottomleft",
       inset = c(-0.01, -0.25), 
       legend = labs,
       text.width = maxwidth * 1,
       col = c("blue","blue","red","red","purple","purple","orange","orange"),
       pch = c(16,12,16,12,16,12,16,12),
       lty = 1,
       ncol = 4,              
       horiz = FALSE,        
       xpd = TRUE,
       cex = 0.8,
       x.intersp = 0.6,
       bty = "n")
     par(mar = c(5, 4, 3, 2))
labs <- c("only_FARM","only_Lasso","Trans_FARM","Trans_Lasso",
          "Oracle_Trans_FARM","Oracle_Trans_Lasso",
          "Pooled_Trans_FARM","Pooled_Trans_Lasso")
