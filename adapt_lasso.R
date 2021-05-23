library(MASS)
library(lars)
library(tsDyn)
library(ncvreg)
library(leaps)
library(MASS)
library(glmnet)

##表1

#产生模型0数据
data_xy_model0=function(i,n,sigma_eps){
  set.seed(i*10)
  p=4;p0=3;rho1=-0.39;rho2=0.23
  C=matrix(,p,p)
  C[1:p0,1:p0]=(1-rho1)*diag(1,p0)+rho1*matrix(1,p0,p0)
  C[1:p0,p]=rho2*matrix(1,p0,1)
  C[p,p]=1
  C[p,1:p0]=t(rho2*matrix(1,p0,1))
  mu=rep(0,p)
  X=mvrnorm(n,mu=mu,Sigma = C)
  eps=rnorm(n,mean=0,sd=sigma_eps)
  beta=c(5.6,5.6,5.6,0)
  Y=X%*%beta+eps
  return(list(Y,X))
}

#lasso_beta函数
lasso_beta_f=function(X,Y){
  lar1=lars(X,Y,type="lasso")
  lasso_beta=lar1$beta[which.min(lar1$Cp),]
  return(lasso_beta)
}

#adaptlasso_beta,gamma
adalasso_beta_f=function(X,Y,gamma){
  ols_coef=lm(Y~X+0)$coef
  adalasso_cv=cv.glmnet(X,Y,alpha=1,penalty.factor=abs(ols_coef)^(-gamma))
  adalasso_lambda=adalasso_cv$lambda[which.min(adalasso_cv$cvm)]
  adalasso=glmnet(X,Y,alpha=1,penalty.factor=abs(ols_coef)^(-gamma),lambda=adalasso_lambda)
  adalasso_beta=adalasso$beta
  return(adalasso_beta)
}


#cv,five-fold
cv_f=function(X,Y,n,i,gamma,sigma_eps){
  n1=n/5
  set.seed(i*5)
  idx=seq(n)
  k1=sample(idx,n1);k2=sample(idx[-k1],n1);k3=sample(idx[-c(k1,k2)],n1)
  k4=sample(idx[-c(k1,k2,k3)],n1);k5=idx[-c(k1,k2,k3,k4)]
  k=matrix(c(k1,k2,k3,k4,k5),5,n1,T)
  rpe=rep(0,n1)
  for (j in 1:5) {
    y_test=Y[k[j,]];x_test=X[k[j,],]
    y_train=Y[-k[j,]];x_train=X[-k[j,],]
    ols_coef=lm(y_train~x_train+0)$coef
    adalasso_cv=cv.glmnet(x_train,y_train,alpha=1,penalty.factor=abs(ols_coef)^(-gamma))
    adalasso_lambda=adalasso_cv$lambda[which.min(adalasso_cv$cvm)]
    adalasso=glmnet(x_train,y_train,alpha=1,penalty.factor=abs(ols_coef)^(-gamma),lambda=adalasso_lambda)
    adalasso_beta=adalasso$beta
    adalasso_hat_y=x_test%*%adalasso_beta
    rpe[j]=sum((adalasso_hat_y-y_test)^2)/n1/sigma_eps^2
  }
  return(mean(rpe))
}

#产生表1结果
table1=function(n,times,sigma_eps){
  lasso_true=0
  adalasso_0.5_true=0
  adalasso_1_true=0
  adalasso_2_true=0
  adalasso_cv_true=0
  for (i in 1:times) {
    Y=data_xy_model0(i=i,n=n,sigma_eps=sigma_eps)[[1]]
    X=data_xy_model0(i=i,n=n,sigma_eps=sigma_eps)[[2]]
    #Lasso
    lasso_beta=lasso_beta_f(X=X,Y=Y)
    if(lasso_beta[1]!=0&lasso_beta[2]!=0&lasso_beta[3]!=0&lasso_beta[4]==0){
      lasso_true=lasso_true+1
    }
    
    #adaptlasso,gamma=0.5
    adalasso_0.5_beta=adalasso_beta_f(X=X,Y=Y,gamma=0.5)
    if(adalasso_0.5_beta[1]!=0 & adalasso_0.5_beta[2]!=0
       &adalasso_0.5_beta[3]!=0&adalasso_0.5_beta[4]==0){
      adalasso_0.5_true=adalasso_0.5_true+1
    }

    #adaptlasso,gamma=1
    adalasso_1_beta=adalasso_beta_f(X=X,Y=Y,gamma=1)
    if(adalasso_1_beta[1]!=0&adalasso_1_beta[2]!=0&
       adalasso_1_beta[3]!=0&adalasso_1_beta[4]==0){
      adalasso_1_true=adalasso_1_true+1
    }
    
    #adaptlasso,gamma=2
    adalasso_2_beta=adalasso_beta_f(X=X,Y=Y,gamma=2)
    if(adalasso_2_beta[1]!=0&adalasso_2_beta[2]!=0&
       adalasso_2_beta[3]!=0&adalasso_2_beta[4]==0){
      adalasso_2_true=adalasso_2_true+1
    }
    
    #adaptlasso_cv,five-fold=5
    cv=matrix(c(0.5,1,2),3,2)
    cv[1,1]=cv_f(X=X,Y=Y,n=n,i=i,gamma=0.5,sigma_eps=sigma_eps)
    cv[2,1]=cv_f(X=X,Y=Y,n=n,i=i,gamma=1,sigma_eps=sigma_eps)
    cv[3,1]=cv_f(X=X,Y=Y,n=n,i=i,gamma=2,sigma_eps=sigma_eps)
    gamma=cv[which.min(cv[1:3,1]),2]
    adalasso_beta_cv=adalasso_beta_f(X=X,Y=Y,gamma=gamma)
    if(adalasso_beta_cv[1]!=0&adalasso_beta_cv[2]!=0&
       adalasso_beta_cv[3]!=0&adalasso_beta_cv[4]==0){
      adalasso_cv_true=adalasso_cv_true+1
    }
  }
  result=list(lasso=lasso_true/times,adalasso_0.5=adalasso_0.5_true/times,
              adalasso_1=adalasso_1_true/times,adalasso_2=adalasso_2_true/times,
              adalasso_cv=adalasso_cv_true/times)
  return(result)
}

list(col1=table1(n=60,sigma_eps=9,times=100),
col2=table1(n=120,sigma_eps=5,times=100),
col3=table1(n=300,sigma_eps=3,times=100))


##表2

#产生模型1数据
data_xy_model1=function(i,n,sigma_eps){
  p=8
  rho=0.5
  sigma_x=matrix(0,p,p)
  for (t in 1:p) {
    for (j in 1:p) {
      sigma_x[t,j]=rho^(abs(t-j)) 
    }
  }
  mu=rep(0,p)
  set.seed(i*10)
  X=mvrnorm(n,mu=mu,Sigma=sigma_x)
  eps=rnorm(n,mean=0,sd=sigma_eps)
  beta=c(3,1.5,0,0,2,0,0,0)
  Y=X%*%beta+eps
  return(list(X,Y))
}

#相对预测误差
RPE=function(sigma,hat_y,Y){
  rpe=sum((hat_y-Y)^2)/n/sigma^2
  return(rpe)
}

#产生表2结果
table1=function(n,sigma_eps,times){
  lasso_RPE=matrix(,times,1)
  lasso_RPE_sd=matrix(,times,1)
  adaptlasso_RPE=matrix(,times,1)
  adaptlasso_RPE_sd=matrix(,times,1)
  SCAD_RPE=matrix(,times,1)
  SCAD_RPE_sd=matrix(,times,1)
  for (i in 1:times) {
    set.seed(i)
    X=data_xy_model1(i,n,sigma_eps)[[1]]
    Y=data_xy_model1(i,n,sigma_eps)[[2]]
    #Lasso
    lar1=lars(X,Y,type="lasso")
    lasso_beta=lar1$beta[which.min(lar1$Cp),]
    #lasso_p=length(which(lasso_beta!=0))
    lasso_hat_y=X%*%lasso_beta
    lasso_RPE[i]=RPE(sigma_eps,lasso_hat_y,Y)
    
    #SCAD
    scad=ncvreg(X, Y,penalty="SCAD",gamma=3.7)
    scad_cv=cv.ncvreg(X, Y,penalty="SCAD")
    scad_lambda=scad_cv$lambda[which.min(scad_cv$cve)]
    SCAD_beta=coef(scad, lambda=scad_lambda)
    #SCAD_p=length(which(SCAD_beta!=0))
    SCAD_hat_y=X%*%SCAD_beta[-1]
    SCAD_RPE[i]=RPE(sigma_eps,SCAD_hat_y,Y)
    
    #adapt lasso
    ols_coef=lm(Y~X+0)$coef
    adaptlasso_cv=cv.glmnet(X,Y,alpha=1,penalty.factor=1/ols_coef)
    adaptlasso_lambda=adaptlasso_cv$lambda[which.min(adaptlasso_cv$cvm)]
    adaptlasso=glmnet(X,Y,alpha=1,penalty.factor=1/ols_coef,lambda=adaptlasso_lambda)
    adaptlasso_beta=adaptlasso$beta
    adaptlasso_hat_y=X%*%adaptlasso_beta
    adaptlasso_RPE[i]=RPE(sigma_eps,adaptlasso_hat_y,Y)
  }
  adaptlasso_RPE_median=median(adaptlasso_RPE)
  lasso_RPE_median=median(lasso_RPE)
  SCAD_RPE_median=median(SCAD_RPE)
  for (j in 1:times) {
    lasso_RPE_sd[j]=sd(sample(lasso_RPE,times,replace=TRUE))
    adaptlasso_RPE_sd[j]=sd(sample(adaptlasso_RPE,times,replace = TRUE))
    SCAD_RPE_sd[j]=sd(sample(SCAD_RPE,times,replace=TRUE)) 
  }
  lasso_RPE_sd_median=median(lasso_RPE_sd)
  SCAD_RPE_sd_median=median(SCAD_RPE_sd)
  adaptlasso_RPE_sd_median=median(adaptlasso_RPE_sd)
  result=list(lasso=c(lasso_RPE_median,lasso_RPE_sd_median),
              adaptlasso=c(adaptlasso_RPE_median,adaptlasso_RPE_sd_median),
              SCAD=c(SCAD_RPE_median,SCAD_RPE_sd_median))
  return(result)
}

times=10
n=20
g1=data.frame("sigma=1"=c(table1(n,1,times)[[1]],table1(n,1,times)[[2]],table1(n,1,times)[[3]]),
           "sigma=3"=c(table1(n,3,times)[[1]],table1(n,3,times)[[2]],table1(n,3,times)[[3]]),
           "sigma=6"=c(table1(n,6,times)[[1]],table1(n,6,times)[[2]],table1(n,6,times)[[3]])
           ,row.names=c("Lasso","Lasso_sd","adaptlasso","adaptlasso_sd","SCAD","SCAD_sd"))
n=60
g2=data.frame("sigma=1"=c(table1(n,1,times)[[1]],table1(n,1,times)[[2]],table1(n,1,times)[[3]]),
              "sigma=3"=c(table1(n,3,times)[[1]],table1(n,3,times)[[2]],table1(n,3,times)[[3]]),
              "sigma=6"=c(table1(n,6,times)[[1]],table1(n,6,times)[[2]],table1(n,6,times)[[3]])
              ,row.names=c("Lasso","Lasso_sd","adaptlasso","adaptlasso_sd","SCAD","SCAD_sd"))
list("Model1(n=20)"=g1,"Model1(n=60)"=g2)

