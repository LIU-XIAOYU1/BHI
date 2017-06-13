data<-generate(n,intercept,a1,a2,b1,b2,p,v)

library(HI)
library(mvtnorm)
function(xcov,Z,L,R,status,k,cof_range,niter=10,burn=2500)
{
  # help function
  Ispline <- function(x, order, knots) {
    k = order + 1
    m = length(knots)
    n = m - 2 + k
    t = c(rep(1, k) * knots[1], knots[2:(m - 1)], rep(1, 
                                                      k) * knots[m])
    yy1 = array(rep(0, (n + k - 1) * length(x)), dim = c(n + 
                                                           k - 1, length(x)))
    for (l in k:n) {
      yy1[l, ] = (x >= t[l] & x < t[l + 1])/(t[l + 1] - 
                                               t[l])
    }
    yytem1 = yy1
    for (ii in 1:order) {
      yytem2 = array(rep(0, (n + k - 1 - ii) * length(x)), 
                     dim = c(n + k - 1 - ii, length(x)))
      for (i in (k - ii):n) {
        yytem2[i, ] = (ii + 1) * ((x - t[i]) * yytem1[i, 
                                                      ] + (t[i + ii + 1] - x) * yytem1[i + 1, ])/(t[i + 
                                                                                                      ii + 1] - t[i])/ii
      }
      yytem1 = yytem2
    }
    index = rep(0, length(x))
    for (i in 1:length(x)) {
      index[i] = sum(t <= x[i])
    }
    yy = array(rep(0, (n - 1) * length(x)), dim = c(n - 1, 
                                                    length(x)))
    if (order == 1) {
      for (i in 2:n) {
        yy[i - 1, ] = (i < index - order + 1) + (i == 
                                                   index) * (t[i + order + 1] - t[i]) * yytem2[i, 
                                                                                               ]/(order + 1)
      }
    }
    else {
      for (j in 1:length(x)) {
        for (i in 2:n) {
          if (i < (index[j] - order + 1)) {
            yy[i - 1, j] = 1
          }
          else if ((i <= index[j]) && (i >= (index[j] - 
                                             order + 1))) {
            yy[i - 1, j] = (t[(i + order + 1):(index[j] + 
                                                 order + 1)] - t[i:index[j]]) %*% yytem2[i:index[j], 
                                                                                         j]/(order + 1)
          }
          else {
            yy[i - 1, j] = 0
          }
        }
      }
    }
    return(yy)
  }
  positivepoissonrnd <- function(lambda) {
    samp = rpois(1, lambda)
    while (samp == 0) {
      samp = rpois(1, lambda)
    }
    return(samp)
  }
  beta_fun<-function(x,beta,xx,sigma.beta,te1,te2,j)
  {
    beta[j]<-x
    tt<-sum(xx[,j]*x*te1)-sum(exp(xx%*%beta)*te2)-x^2/sigma.beta^2/2
    return(tt)
  }
  ind_fun<-function(x,beta,xx,sigma.beta,te1,te2,j)
    (x > -cof_range)*(x < cof_range)
  
  
  #metropolis-hastings for cure
  # t distribution
  alpha_fun<-function(X,u,alpha,alpha_0,sigma.alpha,covG)
  {
    Z<-cbind(1,X)
    l<-ncol(Z)
    alpha.new<-c(alpha+rmvt(1,sigma =covG,l))
    
    #log prior, multivariate normal distribution
    lpold<-dmvnorm(alpha,mean=rep(alpha_0,l),sigma = diag(sigma.alpha,l),log = T)
    lpnew<-dmvnorm(alpha.new,mean = rep(alpha_0,l),sigma = diag(sigma.alpha,l),log = T)
    
    #cure rate
    
    piold<-exp(Z%*%alpha)/(1+exp(Z%*%alpha))
    pinew<-exp(Z%*%alpha.new)/(1+exp(Z%*%alpha.new))
    
    lold<-sum(u*piold+(1-u)*(1-piold))
    lnew<-sum(u*pinew+(1-u)*(1-pinew))
    
    ratio<-(lnew+lpnew)-(lold+lpold)   #random walk m-h algorithm is symestic
    U<-runif(1)
    
    accept<-0
    if(is.na(ratio)==FALSE)
    {
      if(log(U)<ratio)
      {
        alpha<-alpha.new
        accept<-accept+1
      }
    }
     list(alpha=alpha,accept=accept)
    }
 
  
  #  hyper parameter
  
  cof_range=2
  alpha_0=0                 #prior mean of alpha
  sigma.alpha=100             #prior variance of alpha
  covG= diag(0.1,ncol(xcov)+1)                    #proposal setting of alpha in t distribution
  sigma.beta=100
  a_lam=1
  b_lam=1
  
  #  data agumentation process 
    
   L=matrix(data$L,ncol=1)
   R=matrix(ifelse(is.na(data$R),0,data$R),ncol = 1)
   status=matrix(data$status,ncol = 1)
   xcov=as.matrix(cbind(data$z1,data$z2),ncol=2)
   p=ncol(xcov)
   n=nrow(L)
   t1=R*as.numeric(status==0)+L*as.numeric(status==1)
   t2=R*as.numeric(status==1)+L*as.numeric(status==2)
   obs=cbind(t1,t2)
   order=3
   knots=seq(min(obs),max(obs),length=10)
   grids=seq(min(obs),max(obs),length=100)
   k=length(knots)-2+order
   bis1=Ispline(t(obs[, 1]),order,knots)  # k*n matrix
   bis2=Ispline(t(obs[, 2]),order,knots)
   bisg=Ispline(grids,order,knots)
   lambda = rgamma(1, a_lam, rate = b_lam)
   gamcoef=matrix(rgamma(k,1,1),ncol = k)
   #beta=matrix(rep(0,p),ncol = 1)
   Lamb1=t(gamcoef%*%bis1)
   Lamb2=t(gamcoef%*%bis2)
   Lambg=t(gamcoef%*%bisg)
   
   
   
   # build the track of parameter 
  
   n=length(L)
   alpha.track=matrix(0,nrow = niter+1,ncol = ncol(xcov)+1)
   beta.track=matrix(0,nrow=niter+1,ncol=ncol(xcov))
   gamma.track=matrix(1,nrow = niter+1,ncol=k)
   u.track=matrix(1,nrow = niter+1,ncol = n)
   S.track=matrix(1,nrow = niter+1,ncol = n)
   accepta<-rep(0,niter)
   
   # intial parameters 
   
   alpha.track[1,]=c(1,1,-1)
   beta.track[1,]=c(1,-1)
   gamma.track[1,]=c(rep(1,k))
   u.track[1,]=c(rep(1,n))
  
   
   
   iter=1
   for(iter in 1:niter)
   {
    
     choose<-u.track[iter,]==1  #index of u=1
     x_new<-xcov[choose,]
     L_new<-L[choose,]
     R_new<-R[choose,]
     status_new<-status[choose,]
     m<-nrow(x_new)
     t1_new<-t1[choose,]
     t2_new<-t2[choose,]
     obs_new=cbind(t1_new,t2_new)
     order=3
     knots=seq(min(obs_new),max(obs_new),length=10)
     k=length(knots)-2+order
     bis1_new=Ispline(t(obs_new[, 1]),order,knots)  # k*m matrix
     bis2_new=Ispline(t(obs_new[, 2]),order,knots)
     gamcoef=matrix(rgamma(k,1,1),ncol = k)
     beta=matrix(rep(0,p),ncol = 1)
     Lamb1_new=t(gamcoef%*%bis1_new)
     Lamb2_new=t(gamcoef%*%bis2_new)
     
     
     
     
     
     
     
     #sample possion process
     
     v=array(rep(0,m),dim = c(m,1))
     w=v
     vv=array(rep(0,m*k),dim = c(m,k))
     ww=vv
     
     for (i in 1:m)
     {
       if(status_new[i]==0){
         templam1=Lamb1_new[i] * exp(x_new[i,]%*%beta)
         v[i]=positivepoissonrnd(templam1)
         vv[i,]=rmultinom(1,v[i], gamcoef* t(bis1_new[,i]))
       }
      else if(status_new[i]==1)
       {
        templam2=(Lamb2_new[i]-Lamb1_new[i]) * exp(x_new[i,]%*%beta)
        w[i] = positivepoissonrnd(templam2)
        ww[i,]=rmultinom(1,w[i],gamcoef*t(bis2_new[,i]-bis1_new[,i]))
       }
     }
   
     
      te1 = v * as.numeric(status_new == 0) + w * as.numeric(status_new == 
                                                           1)
      te2 = (Lamb1_new * as.numeric(status_new == 0) + Lamb2_new * 
               as.numeric(status_new == 1) + Lamb2_new * as.numeric(status_new == 
                                                                 2))
      for(j in 1:p)
      {
        beta[j]=arms(beta[j],beta_fun,ind_fun,1,
                     xx=matrix(x_new,ncol = ncol(xcov)),beta=beta,j=j,
                     sigma.beta=sigma.beta,te1=te1,te2=te2)

      }
      
      beta.track[iter,]<-beta
     
      
      for(l in 1:k)
      {
        tempa = 1 + sum(vv[, l] * as.numeric(status_new == 0) * 
                          (bis1_new[l, ] > 0) + ww[, l] * as.numeric(status_new == 
                                                                   1) * ((bis2_new[l, ] - bis1_new[l, ]) > 0))
        tempb = lambda + sum((bis1_new[l, ] * as.numeric(status_new == 
                                                    0) + bis2_new[l, ] * as.numeric(status_new == 1) + bis2_new[l, 
                                                                                                    ] * as.numeric(status_new == 2)) * exp(x_new %*% beta))
        gamcoef[l] = rgamma(1, tempa, rate = tempb)
      }
     
     gamma.track[iter,]<-gamcoef
     
     lam<-rgamma(1,a_lam+k, b_lam+sum(gamcoef))
   
     
     
     # sample for alpha #
     if (iter >= 2*burn) covbG<-cov(alpha.track[(burn+1):(2*burn),])  #????
     f<-alpha_fun(X=x_new,u=u.track[iter,],alpha=alpha.track[iter,],alpha_0,sigma.alpha,covG)
     accepta[iter]<-f$accept
     alpha.track[iter,]<-alpha<-f$alpha
    
     # sample for u 
     
     Lamb1=t(gamcoef%*%bis1)
     Lamb2=t(gamcoef%*%bis2) # updata the whole Lamb1
     Lambg=t(gamcoef%*%bisg)
     
     S.track[iter,]<-S<-exp(-Lamb1*exp(X%*%beta))
     cure<-exp(Z%*%alpha)/(1+exp(Z%*%alpha))
     p_u<-cure*S/(1-cure+cure*S)
     u.temp<-rbinom(rep(1,n),rep(1,n),p_u)
     u.tem1<-ifelse(status==2,u.temp,1)
     R_max<-max(obs*(1-as.numeric(status==2))) 
     u.new<-ifelse(L>R_max,0,u.tem1)
     u.track[iter,]<-u.new
     
    
    if(iter%%1000==0)
    {
      cat(paste(iter, " iterations", " have finished.\n", sep=""))
    }
 
   }
  list(alpha.track=alpha.track,beta.track=beta.track,gamma.track=gamma.track,
       u.track=u.track,S.track=S.track,accept=accepta)  
    
    
}

### help function ####
 