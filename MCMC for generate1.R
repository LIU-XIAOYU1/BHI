################## modify the proposal distribution of alpha ########
################## covG =0.6,0.1,0.1= cov(alpha.track[burn+1:2*burn,])###########
##########################################################
library(HI)
library(mvtnorm)
##### generate data ####
#set.seed(2)

n<-200

repl<-1
intercept<-1
a1<-1
a2<--1
b1<-1
b2<--1

p<-5
v<-0.5 # decide the left point
order<-3 #decide the order of i spline
cof_range=2

#  hyper parameter

cof_range=2
alpha_0=0                 #prior mean of alpha
sigma.alpha=10            #prior variance of alpha
covG= diag(0.1,3)                    #proposal setting of alpha in t distribution
#covG=diag(c(0.5,0.1,0.1))
sigma.beta=10
a_lam=1
b_lam=1
lambda = rgamma(1, a_lam, rate = b_lam)


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
  
  lold<-sum(u*log(piold)+(1-u)*log(1-piold))
  lnew<-sum(u*log(pinew)+(1-u)*log(1-pinew))
  
  ratio<-(lnew+lpnew)-(lold+lpold)   #random walk m-h algorithm is symestic
  U<-runif(1)
  
  accept<-0
  if(is.na(ratio)==FALSE)
  {
    if(log(U)<ratio)
    {
      alpha<-alpha.new
      accept<-1
    }
  }
  list(alpha=alpha,accept=accept)
}



BHIsur<-function(xcov,Z,L,R,status,cof_range=2,niter=10000,burn=niter/4)
{
  
  #  data agumentation process 
  
  # L=matrix(data$L,ncol=1)
  # R=matrix(ifelse(is.na(data$R),0,data$R),ncol = 1)
  # status=matrix(data$status,ncol = 1)
  # xcov=as.matrix(cbind(data$z1,data$z2),ncol=2)
  p=ncol(xcov)
  n=nrow(L)
  X=xcov
  Z=cbind(1,X)
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
  
  
  
  for (iter in 1:niter)
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
    gamcoef=gamma.track[iter,]
    beta=beta.track[iter,]
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
        templam1=(Lamb2_new[i]-Lamb1_new[i]) * exp(x_new[i,]%*%beta)
        w[i] = positivepoissonrnd(templam1)
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
    
    beta.track[iter+1,]<-beta
    
    
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
    
    gamma.track[iter+1,]<-gamcoef
    
    lam<-rgamma(1,a_lam+k, b_lam+sum(gamcoef))
    #### update of  Lamb1_new ?
    
    
    # sample for alpha #
    if (iter >= 2*burn) covG<-cov(alpha.track[(burn+1):(2*burn),])  # modify the covariance of alpha
    f<-alpha_fun(X=xcov,u=u.track[iter,],alpha=alpha.track[iter,],alpha_0,sigma.alpha,covG)
    accepta[iter+1]<-f$accept
    alpha.track[iter+1,]<-alpha<-f$alpha
    
    # sample for u 
    
    Lamb1=t(gamcoef%*%bis1)
    Lamb2=t(gamcoef%*%bis2) # updata the whole Lamb1
    Lambg=t(gamcoef%*%bisg)
    
    S.track[iter+1,]<-S<-exp(-Lamb1*exp(X%*%beta))
    cure<-exp(Z%*%alpha)/(1+exp(Z%*%alpha))
    p_u<- cure*S/(1-cure+cure*S)
    u.temp<-rbinom(rep(1,n),rep(1,n),p_u)
    u.tem1<-ifelse(status==2,u.temp,1)
    R_max<-max(obs*(1-as.numeric(status==2))) 
    u.new<-ifelse(L>R_max,0,u.tem1)
    u.track[iter+1,]<-u.new
    
    
    if(iter%%1000==0)
    {
      cat(paste(iter, " iterations", " have finished.\n", sep=""))
    }
  }
  
  f.alpha<-apply(alpha.track[burn:niter,],2,mean)
  f.beta<-apply(beta.track[burn:niter,],2,mean)
  f.gamma<-apply(gamma.track[burn:niter,],2,mean)
  list(alpha=f.alpha,beta=f.beta,gamma=f.gamma,accept=sum(accepta)/niter,a_lam=a_lam,b_lam=b_lam)  
  
}


##  replication #

alpha.rep=matrix(0,nrow = repl,ncol = 3)
beta.rep=matrix(0,nrow=repl,ncol=2)
accept.rep=matrix(0,nrow = repl,ncol = 3)
gamma.rep=matrix(1,nrow =repl,ncol=11)
S.rep=matrix(1,nrow = repl,ncol = n)

for(k in 1:repl)
{
  data<-generate1(n,intercept,a1,a2,b1,b2,p,v)
  L=matrix(data$L,ncol=1)
  R=matrix(ifelse(is.na(data$R),0,data$R),ncol = 1)
  status=matrix(data$status,ncol = 1)
  xcov=as.matrix(cbind(data$z1,data$z2),ncol=2)
  result<-BHIsur(xcov=xcov,L=L, R=R, status=status)
  cat(paste(k, " replication", " have finished.\n", sep=""))
  print(result$accept)
  alpha.rep[k,]<-result$alpha
  beta.rep[k,]<-result$beta
  gamma.rep[k,]<-result$gamma
  accept.rep[k,]<-result$accept
  
}
##### plot cunmulative hazard function


data<-generate1(n,intercept,a1,a2,b1,b2,p,v)
L=matrix(data$L,ncol=1)
R=matrix(ifelse(is.na(data$R),0,data$R),ncol = 1)
xcov=as.matrix(cbind(data$z1,data$z2),ncol=2)
obser<-(cbind(L,R))
timepoint<-seq(min(obser),max(obser),length.out = nrow(xcov))
knots=seq(min(obser),max(obser),length=10)
bist<-Ispline(timepoint,order,knots)


# plot cumulative hazard funtion
plot(timepoint,exp(-timepoint),type = "l",xlab = "T",ylab = "Lambda[0](T)",lwd=3, col="black")
for(k in 1:q)
{
  Lam=gamma.rep[k,]%*%bist
  lines(timepoint,Lam,col="red")
}


## plot survival function

plot(timepoint,exp(-timepoint),type="l",col="black",lwd=3,xlab="T", 
     ylab=expression(S(T)))
q<-nrow(gamma.rep)
for(l in 1:q)  
{
  Lambt<-t(gamma.rep[l,]%*%bist)
  #Sur<-exp(-Lambt*exp(xcov%*%beta.rep[l,]))
  Sur<-exp(-Lambt)
  lines(timepoint,Sur,col = "red")
  
}

apply(alpha.rep,2,mean)
apply(beta.rep,2,mean)
mydata<-data.frame(alpha.rep,beta.rep,gamma.rep,accept.rep)
write.csv(mydata,file = paste("covG =0.6,0.1,0.1 ",repl,".csv",sep=""))






