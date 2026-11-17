# Table 1: 95% nonstudentized simultaneous intervals for Trans-FARM.
# Each setting name gives the target sample size n0 followed by the dimension p.
settings_to_run <- c("n200_p200", "n200_p500", "n500_p500")
cases_to_run <- c("a", "b", "c", "d")
normal_target_toeplitz <- TRUE  # For cases a/b. Set FALSE for an identity target matrix.

itr.max <- 500L
B <- 500L
ncores <- 24L
alpha <- 0.05
output_dir <- file.path(getwd(), "CI_nonstudentized")

library(glmtrans)
library(glmnet)
library(mvtnorm)
library(foreach)
library(doSNOW)
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

Trans_FARM_nonstudentized = function(n0,nk,p,r,K,K_A,s,q,input_type,w,eta,
                                      alpha,B,ci_case,normal_target_toeplitz){
  normal_U <- ci_case %in% c("a", "b")
  normal_E <- ci_case %in% c("a", "c")
  k_CI <- p/4
  beta0 = c(rep(q,s),rep(0,p-s))
  #beta0 = c(runif(s,-1,1),rep(0,p-s))
  B0 = matrix(runif(p*r,-1,1),nrow=p)
  F0 = matrix(rnorm(n0*r),nrow=n0)
  meanx=rep(0,p)
  sigma0=diag(p)
  for (i in 1:p){
    for (j in 1:p){sigma0[i,j]=0.5^(abs(i-j))}}
  U0 = if (!normal_U) {
    rmvt(n0,delta=meanx,sigma=sigma0,df=10)
  } else if (normal_target_toeplitz) {
    rmvnorm(n0,mean=meanx,sigma=sigma0)
  } else {
    matrix(rnorm(n0*p),nrow=n0)
  }
  X0 = F0%*%t(B0)+U0
  E0 = if (normal_E) rnorm(n0,mean=0,sd=1) else rt(n0,df=5)
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
    epsk = rmvnorm(1,mean=meanx,sigma=0.2^2*diag(p))
    sigmak = sigma0+t(epsk)%*%epsk
    Uk = if (normal_U) {
      rmvnorm(nk,mean=meanx,sigma=sigmak)
    } else {
      rmvt(nk,delta=meanx,sigma=sigmak,df=10)
    }
    Xk = Fk%*%t(Bk)+Uk
    Ek = if (normal_E) rnorm(nk,mean=0,sd=1) else rt(nk,df=5)
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

  # A single unstandardized critical value is already on the coefficient scale.
  boot_max <- replicate(B, {
    ei <- rnorm(n0)
    Q_e_hat <- sigma_hat * (Theta_U_hat %*% colMeans(U.hat0 * ei))
    max(abs(Q_e_hat))
  })
  half_width <- unname(quantile(boot_max, probs = 1 - alpha))
  covered <- max(abs(beta0[seq_len(k_CI)] - beta_debiased_Trans_FARM)) <= half_width
  c(covered = as.numeric(covered), interval_length = 2 * half_width)
}

run_Trans_FARM_nonstudentized <- function(n0, nk, p, r, K, K_A, s, q, input_type,
                                          w, eta, alpha, B, ci_case, normal_target_toeplitz,
                                          itr.max, ncores) {
  cl <- parallel::makeCluster(min(ncores, itr.max))
  on.exit(parallel::stopCluster(cl), add = TRUE)
  registerDoSNOW(cl)
  pb <- txtProgressBar(max = itr.max, style = 3)
  on.exit(close(pb), add = TRUE)
  opts <- list(progress = function(n) setTxtProgressBar(pb, n))

  reps <- foreach(itr = seq_len(itr.max),
                  .packages = c("glmnet", "glmtrans", "mvtnorm", "SIS"),
                  .export = c("rcv", "Trans_FARM_nonstudentized"),
                  .options.snow = opts) %dopar% {
    Trans_FARM_nonstudentized(
      n0 = n0, nk = nk, p = p, r = r, K = K, K_A = K_A, s = s,
      q = q, input_type = input_type, w = w, eta = eta, alpha = alpha, B = B,
      ci_case = ci_case, normal_target_toeplitz = normal_target_toeplitz
    )
  }
  colMeans(do.call(rbind, reps))
}

# Section 5.2 settings, in Table 1 column order.
settings <- data.frame(
  n0 = c(200L, 200L, 500L, 500L),
  p  = c(200L, 500L, 200L, 500L)
)
rownames(settings) <- sprintf("n%d_p%d", settings$n0, settings$p)
if (identical(settings_to_run, "all")) settings_to_run <- rownames(settings)
if (!length(settings_to_run) || !all(settings_to_run %in% rownames(settings))) {
  stop("Choose a setting from n200_p200, n200_p500, n500_p200, n500_p500, or use 'all'.")
}
if (!length(cases_to_run) || !all(cases_to_run %in% c("a", "b", "c", "d"))) {
  stop("Choose one or more cases from a, b, c, d.")
}
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
message("Results directory: ", normalizePath(output_dir, mustWork = TRUE))

for (ci_case in unique(cases_to_run)) {
  for (setting in unique(settings_to_run)) {
    n0_i <- settings[setting, "n0"]
    p_i <- settings[setting, "p"]
    message("Running nonstudentized Trans-FARM: case ", ci_case, ", n0 = ", n0_i, ", p = ", p_i)
    result <- run_Trans_FARM_nonstudentized(
      n0 = n0_i, nk = 500L, p = p_i, r = 2L, K = 10L, K_A = 5L, s = 5L,
      q = 1, input_type = "gamma", w = 0.5, eta = 5, alpha = alpha, B = B,
      ci_case = ci_case, normal_target_toeplitz = normal_target_toeplitz,
      itr.max = itr.max, ncores = ncores
    )
    target_matrix <- if (ci_case %in% c("a", "b") && !normal_target_toeplitz) "identity" else "Toeplitz"
    row <- data.frame(
      case = ci_case, n0 = n0_i, p = p_i, nk = 500L, target_U_matrix = target_matrix,
      replications = itr.max, bootstrap_draws = B, alpha = alpha,
      ACP_Trans_FARM_nonstudentized = unname(result["covered"]),
      AL_Trans_FARM_nonstudentized = unname(result["interval_length"])
    )
    filename <- sprintf("CI_nonstudentized_case_%s_%s_%s.csv", ci_case, setting, target_matrix)
    write.csv(row, file.path(output_dir, filename), row.names = FALSE)
    print(row, row.names = FALSE)
  }
}
