library(HI)
function(xcov,Z,L,R,status,k,niter=10,burn=)
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
  beta_fun<-function(x,beta,xx,sigma.beta,te1,te2)
  {
    beta[j]<-x
    tt<-sum(xx[,j]*x*te1)-sum(exp(xx%*%beta)*te2)-x^2/sigma.beta^2/2
    return(tt)
  }
  ind_fun<-function()
  alpha_fun<-function()
  
  
  
  
  
  # build the track of parameter 
  
  
  n=length(T)
  alpha.track=matrix(0,nrow = niter+1,ncol = ncol(Z))
  beta.track=matrix(0,nrow=niter+1,ncol=ncol(X))
  gamma.track=matrix(1,nrow = niter+1,ncol=k)
  u.track=matrix(1,nrow = niter+1,ncol = n)
  lam.track=matrix(1,nrow = niter+1,ncol=n)
  
  
  # intial parameters 
  
  alpha.track[1,]=c(1,1,-1)
  beta.track[1,]=c(1,-1)
  gamma.track[1,]=c(rep(1,k))
  u.track[1,]=c(rep(1,n))
  sigma.alpha=
  sigma.beta=
  lambda=
    
    
  #  data agumentation process 
    
   L=matrix(L,ncol=1)
   R=matrix(ifelse(is.na(R),0,R),ncol = 1)
   status=matrix(status,ncol = 1)
   xcov=as.matrix(xcov)
   p=ncol(xcov)
   n=nrow(L)
   t1=R*as.numeric(status==0)+L*as.numeric(status==1)
   t2=R*as.numeric(status==1)+L*as.numeric(status==2)
   obs=cbind(t1,t2)
   order=3
   knots=seq(min(obs),max(obs),length=10)
   k=length(knots)-2+order
   bis1=Ispline(t(obs[, 1]),order,knots)  # k*n matrix
   bis2=Ispline(t(obs[, 2]),order,knots)
   gamcoef=matrix(rgamma(k,1,1),ncol = k)
   beta=matrix(rep(0,p),ncol = 1)
   Lamb1=t(gamcoef%*%bis1)
   Lamb2=t(gammcoef%*%bis2)
   
   
   
   # build the track of parameter 
  
   n=length(T)
   alpha.track=matrix(0,nrow = niter+1,ncol = ncol(Z))
   beta.track=matrix(0,nrow=niter+1,ncol=ncol(X))
   gamma.track=matrix(1,nrow = niter+1,ncol=k)
   u.track=matrix(1,nrow = niter+1,ncol = n)
   lam.track=matrix(1,nrow = niter+1,ncol=n)
   
   iter=1
   while(iter<niter+1)
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
     Lamb2_new=t(gammcoef%*%bis2_new)
     
     
     
     
     
     
     
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
         vv[i,]=rmultinom(1,v[i], gamcoef%*% t(bis1_new[,i]))
       }
      else if(status_new[i]==1)
       {
        templam2=(Lamb2_new[i]-Lamb1_new[i]) * exp(x_new[i,]%*%beta)
        w[i] = positivepoissonrnd(templam1)
        ww[i,]=rmultinom(1,w[i],gamcoef%*% t(bis2_new[,i]-bis1_new[,i]))
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
                     sigma.beta=sigma.beta,te1=te1_new,te2=te2_new)

      }
     
      for(l in 1:k)
      {
        
      }
     
     
     
     
     #update Lamb1 and Lamb2
   }
    
    
    
}

### help function ####
 