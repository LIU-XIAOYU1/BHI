##################################################
######## MFVB FOR BHI ############################
##################################################

## Generate data of BHI
n<-200
MaxIter<-100
intercept<-1
a1<-1
a2<--1
b1<-1
b2<--1

p<-5
v<-0.5 # decide the left point
order<-3 #decide the order of i spline
cof_range=2


generate1<-function(n,intercept,a1,a2,b1,b2,p,v)
{
  #n sample size
  #intercept cofficients of cure rate
  #a1,a2     cofficients of cure rate
  #b1,b2     cofficients of PH model
  #p         examination times
  #v         first generate uniform range
  #status    censording indicator, 0-left censored,1-interval censored, 2-right censored
  #L         left point of all objects
  #R         right point of all objects
  #Y         nopossible examination time
  #cure rate indicator
  z1<-rbinom(n,1,0.5) 
  z2<-rnorm(n,0,0.5^2) 
  linpred<-cbind(1,z1,z2) %*% c(intercept,a1,a2)
  prob<-exp(linpred)/(1+exp(linpred))
  Y<-rbinom(n=n, size=1, prob=prob)
  #survival probabilities exponentional(1) and time calculation
  u<-runif(n,0,1)
  x<-cbind(z1,z2)
  time<--log(1-u)*exp(-x%*%c(b1,b2))
  ###generate examination times####
  
  C<-matrix(rep(0,n),ncol = 1)   # all right point of subject
  status<-matrix(rep(0,n),ncol=1)
  L<-matrix(rep(0,n),ncol=1)
  R<-matrix(rep(0,n),ncol=1)
  for(j in 1:n)
  { 
    P<-rpois(1,p)+2    #  the number of examination times 
    obs<-matrix(rep(0,P),ncol =1)
    obs[1]<-runif(1,0,v)
    len<-runif(1,1,3) # decide the length of interval 
    for (i in 2:P) 
    {
      obs[i]<-obs[1]+len*(i-1)  
    }               
    C[j]<-obs[P]
    T<-time*Y+C*(1-Y) # generate exact failure data
    if(Y[j]==0) #cure fraction
    {
      status[j]<-2
      L[j]<-obs[P]
      R[j]<-NA    #(L_i,R_i)=(obs_p, infinite)
    } 
    if(Y[j]==1) #uncure fraction
    {
      if(T[j]<obs[1]) 
      {
        status[j]<-0
        L[j]<-0
        R[j]<-obs[1]
      } #left censored, (L_j,R_j)=(0, Y_1)
      
      if(T[j]>obs[P])
      {
        status[j]<-2
        L[j]<-obs[P]
        R[j]<-NA
      }  #right censored,(L_j,R_j)=(obs_p,NA)
      if(T[j]<obs[P]&T[j]>obs[1])
      { status[j]<-1 
      i<-1
      repeat
      {
        i<-i+1
        if(T[j]<obs[i])
        {break}
      }
      L[j]<-obs[i-1]
      R[j]<-obs[i]  # find the interval contain T,(L_j,R_j)=(obs_i-1,obs_i)
      }
    }
  }
  data<-data.frame(T,L,R,status,z1,z2)
  return(data)
}

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



#### Main routine ###
MFVB<-function(data,)
  
  
L=matrix(data$L,ncol=1)
R=matrix(data$R,ncol=1)
R2=ifelse(is.na(R),0,R)
xcov=matrix(c(data$z1,data$z2),ncol=2)
status=data$status

X=xcov
Z=cbind(1,X)
nb<-ncol(xcov)
nu<-nrow(xcov)
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





#### buliding tracks 

mu.beta.track<-matrix(rep(0,nb*MaxIter),ncol = nb)
extr.track<-matrix(rep(0,k*MaxIter),ncol = k)
mu.gamma.track<-matrix(rep(1,k*MaxIter),ncol = k)
mu.u.track<-matrix(rep(1,nu*MaxIter),ncol = nu)

Sigma.q.alpha.track<-vector("list",length = MaxIter)
Sigma.q.beta.track<-vector("list",length= MaxIter)
Sigma.q.gamma.track<-vector("list",length= MaxIter)


#### intialize parameter ####

mu.u.track[1,]<-rep(1,nu)
mu.gamma.track[1,]<-rep(1,k)
gamcoef<-mu.gamma.track[1,]
Lamb1=t(gamcoef%*%bis1)
Lamb2=t(gamcoef%*%bis2)

###convergence trackers

log.p.track = rep(NA,MaxIter)
StepConv = rep(NA,MaxIter)
StepConvBound = 0.01 #relconv 1% relative convergence for now



####obejective function 
h_beta<-function(u,X,status,beta,Lamb1,Lamb2,gamma,v,w)
{
  D=(status==1)+(status==2)
  logp_v=crossprod(X,u*Lamb1*exp(X%*%beta))+crossprod(X,u*v)
  logp_w=crossprod(X,u*)
}
