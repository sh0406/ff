

CPmodelChibReg			= function(
  data		= ySIM,
  
  K 		= K,
  
  start 		= list(
    mu			= rnorm(K,0,1),
    mu2			= rnorm(K,0,1),
    sigma2		= abs(rnorm(K,0,1)),
    #sigma2		= c(sigma2SIM,rep(20,Kmax)),
    #st			= rep(1:2,each=1000)[1:1500],
    st			= stSIM,
    pi			= 	rdirichlet(1,rep(1,K))								
  ),
  
  prior 		= list(
    mu 			= c(0,1000),
    mu2 			= c(0,1000),
    sigma2 		= c(1,1),
    beta		= c(0,1000) 
  ),
  
  MCMC_par 	= list(
    
    iter		= 15000, 
    burnin		= 11000, 
    thin		= 10
    
  ),
  
  
  
  sd_prop			= 0.5
  
)
{
  
  rmnorm=function(n = 1, mean = rep(0, d), varcov)  # mean, var를 가진거의 rmnorm
  {
    d <- if (is.matrix(varcov)) 
      ncol(varcov)
    else 1
    z <- matrix(rnorm(n * d), n, d) %*% chol(varcov) #임의의 x행렬을 정규화  
    y <- t(mean + t(z)) #mean[dX1]
    return(y)
  }
  ##### DATA
  nT 									= nrow(data)
  
  
  #################
  
  MCMCpar = c(
    MCMC_par$iter, 
    MCMC_par$burnin, 
    MCMC_par$thin, 
    round((MCMC_par$iter-MCMC_par$burnin)/MCMC_par$thin)
  )
  
  
  
  
  sdprop 				= sd_prop
  
  indexMCMC 			= MCMCpar[2]
  
  stMCMC 				= start$st
  TT							= table(stMCMC)
  nvec						= rep(0,K )
  nvec[as.numeric(names(TT))] = TT
  
  
  
  muMCMC				= rep(0,K)
  muMCMC[1:K]			= start$mu[1:K]
  
  mu2MCMC				= rep(0,K)
  mu2MCMC[1:K]			= start$mu2[1:K]
  
  sigma2MCMC			= rep(1,K)
  sigma2MCMC[1:K]		= start$sigma2[1:K]
  
  piMCMC			= rep(1,K)
  piMCMC[1:K]		= start$pi[1:K]
  piMCMC[K]		= 1
  
  betaMCMC			= start$beta
  
  
  
  
  
  muOUT 				= matrix(NA, nrow=MCMCpar[4], ncol=K)
  mu2OUT 				= matrix(NA, nrow=MCMCpar[4], ncol=K)
  sigma2OUT 			= matrix(NA, nrow=MCMCpar[4], ncol=K)
  betaOUT 			= matrix(NA, nrow=MCMCpar[4], ncol=1)
  stOUT 				= matrix(NA, nrow=MCMCpar[4], ncol=nT)
  piOUT 				= matrix(NA, nrow=MCMCpar[4], ncol=K)
  nvecOUT 				= matrix(NA, nrow=MCMCpar[4], ncol=K)
  
  appst					= matrix(0, nrow=nT,ncol=K)
  MatrixP 				= matrix(0,nrow=K,ncol=K)
  MatrixB 				= matrix(0,nrow=K,ncol=K)
  iterations			= 0
  
  BIC = c()
  AIC = c()
  DIC = c()
  
  
  #dataT		= c(1:nT)
  dataT=data$x
  data=data$y
  
  for(iMCMC in 1:MCMCpar[4])
  {
    #browser()
    for(jMCMC in 1:indexMCMC)
    {
      
      iterations = iterations+1
      
      ##### beta MCMC
      
      
      betaProp = exp(rnorm(1,log(betaMCMC), sdprop))
      
      
      MHratio		= -0.5*(betaProp-prior$beta[1])^2/prior$beta[2]
      MHratio		= MHratio-(-0.5*(betaMCMC-prior$beta[1])^2/prior$beta[2])
      MHratio 		= MHratio+log(betaProp)-log(betaMCMC)
      
      MHratio		= MHratio+(K-1)*log(betaProp)-(K-1)*log(betaMCMC)
      
      for(k in 1:K)
      {
        MHratio 	=  MHratio+(lgamma(betaProp+1)-lgamma(nvec[k]+betaProp+1-(k==K)) )-(lgamma(betaMCMC+1)-lgamma(nvec[k]+betaMCMC+1-(k==K)))
      }
      
      
      if(runif(1,0,1)<exp(MHratio))
      {
        #browser()
        betaMCMC = betaProp
      }
      
      #print(iterations)
      
      ## 1. parameter<THETA: mu, sigma>
      ##### mu
      
      for(k in 1:K)
      {
        # browser()
        w				= stMCMC==k
        
        # V				= matrix(ncol=2, nrow=2)
        # V[1,1]	= sum(w) # 상태가 k인 애들의 개수
        # V[1,2]	= sum(dataT[w]) #dataT = 1:nT(데이터 개수)_time series의 t 
        # V[2,1]	= sum(dataT[w])
        # V[2,2]	= sum(dataT[w]^2)
        bigX    = cbind(rep(1,sum(w)),dataT[w])
        
        V       = t(bigX)%*%bigX
        
        V				= solve(V/sigma2MCMC[k]+diag(c(1/prior$mu[2],1/prior$mu2[2])))
        
        M				= matrix(ncol=1, nrow=2)
        
        # M[1,1]	= sum(data[w])
        # M[2,1]	= sum(dataT[w]*data[w])
        M       = t(bigX)%*%data[w]
        
        M				= V%*%(M/sigma2MCMC[k]+  c(prior$mu[1]/prior$mu[2],prior$mu2[1]/prior$mu2[2])   )
        
        
        # meanMU	= (prior$mu[2]*sum(data[xiMCMC==k])+sigma2MCMC[k]*prior$mu[1])/( nvec[k]*prior$mu[2]+sigma2MCMC[k])
        #
        # varMU		= sigma2MCMC[k]*prior$mu[2]/( nvec[k]*prior$mu[2]+sigma2MCMC[k])
        # muMCMC[k] = rnorm(1, meanMU,varMU^0.5)
        
        simm = rmnorm(1,M,V)
        
        muMCMC[k] = simm[1]
        mu2MCMC[k] = simm[2]
        
      }
      
      ####### sigma2
      
      
      for(k in 1:K)
      {
        
        Asigma2	= prior$sigma2[1]+nvec[k]/2
        
        w				= stMCMC==k
        Bsigma2 = 1/prior$sigma2[2]+ sum( (data[w]-muMCMC[k]-mu2MCMC[k]*dataT[w])^2 )/2
        
        sigma2MCMC[k] = 1.0/rgamma(1,shape=Asigma2,rate=Bsigma2)
      }
      
      
      ######### ######### ######### ######### 
      ######### sampling st
      ######### ######### ######### ######### k
      if(K>1)
      {
        #browser()
        appst					= matrix(0, nrow=nT,ncol=K)
        diag(MatrixP) = piMCMC
        for(k in 1:(K-1))
        {
          MatrixP[k,k+1] = 1-MatrixP[k,k]
        }
        
        #browser()
        MatrixP[K,K] = 1
        appst[1,1] = 1
        if(K>1)
        {
          if(K>2)
          {
            for(t in 2:(K-1))
            {
              appst[t,] =  log((appst[t-1,]%*%MatrixP )  )+dnorm(data[t], muMCMC+mu2MCMC*dataT[t], sigma2MCMC^0.5, log=T) 
              appst[t,] =   exp(appst[t,])
              appst[t,(t+1):K] =  0
              appst[t,] = appst[t,]/sum(appst[t,]) 
              
            }
          }
          
          for(t in (K):nT)
          {
            appst[t,] =  log((appst[t-1,]%*%MatrixP )  )+dnorm(data[t], muMCMC+mu2MCMC*dataT[t], sigma2MCMC^0.5, log=T) 
            appst[t,] =   exp(appst[t,])
            appst[t,] = appst[t,]/sum(appst[t,])
          }
          
        }else{
          appst[,1] = 1
        }
        
        
        stMCMC[nT] = K
        #stMCMC[nT] = max(which(appst[nT,]>0)) #regime K 찾는거인듯...
        
        # 역순 
        for(t in (nT-1):1)
        {
          # browser()
          k = stMCMC[t+1]
          if(k>1)
          {
            stMCMC[t] = sample(c(k-1,k),1, prob = c( (1-piMCMC[k-1])*appst[t,k-1], piMCMC[k]*appst[t,k]  ))
          }else{
            stMCMC[t] = 1 # (t+1)시점의st가 1이면 t시점은 당연히 st=1임
          }
          
        }
        
        #browser()
        ### nvec contains the number of observations in the regimes
        TT							= table(stMCMC)
        nvec						= rep(0,K )
        nvec[as.numeric(names(TT))] = TT
        
        #### pi
        for(k in 1:(K-1))
        {
          #piMCMC[k]	= rbeta(1,1+nvec[k], betaMCMC+1)
          #piMCMC[k]	= rbeta(1,prior$beta[1]+nvec[k], prior$beta[2]+1)
          piMCMC[k]	= rbeta(1,prior$beta[1]+(nvec[k]-1), prior$beta[2]+1)
        }
      }
      
      
      
      
    }
    indexMCMC = MCMCpar[3]
    
    
    
    ###### SALVO PARAMETRI
    
    muOUT[iMCMC,1:K] 			= muMCMC[1:K]
    mu2OUT[iMCMC,1:K] 			= mu2MCMC[1:K]
    sigma2OUT[iMCMC,1:K] 	= sigma2MCMC[1:K]
    piOUT[iMCMC,1:K] 		= piMCMC[1:K]
    
    #betaOUT[iMCMC,1] 		= betaMCMC
    stOUT[iMCMC,] 				= stMCMC
    nvecOUT[iMCMC,]     =nvec
    
    #### Bayes factor
    
    appst					= matrix(0, nrow=nT,ncol=K)
    diag(MatrixP) = piMCMC
    for(k in 1:(K-1))
    {
      MatrixP[k,k+1] = 1-MatrixP[k,k]
    }
    MatrixP[K,K] = 1
    appst[1,1] = 1
    
    SS = c(1)
    
    if(K>1)
    {
      #browser()
      if(K>2){
        for(t in 2:(K-1))
        {
          appst[t,] =  log((appst[t-1,]%*%MatrixP )  )+dnorm(data[t], muMCMC+mu2MCMC*dataT[t], sigma2MCMC^0.5, log=T)
          appst[t,] =   exp(appst[t,])
          appst[t,(t+1):K] =  0
          SS[t]			= sum(appst[t,])
          appst[t,]= appst[t,]/sum(appst[t,])
        }
      }
      for(t in (K):nT)
      {
        appst[t,] =  log((appst[t-1,]%*%MatrixP )  )+dnorm(data[t], muMCMC+mu2MCMC*dataT[t], sigma2MCMC^0.5, log=T)
        appst[t,] =   exp(appst[t,])
        SS[t]			= sum(appst[t,])
        appst[t,]= appst[t,]/sum(appst[t,])
      }
      
    }else{
      for(t in 1:nT)
      {
        #browser()
        #SS[t] = exp(dpois(data[t], lambdaMCMC[1], log=T))
        SS[t] = exp(dnorm(data[t], muMCMC+mu2MCMC*dataT[t], sigma2MCMC^0.5, log=T))
      }
    }
    
    
    
    
    BIC[iMCMC] = 0
    BIC[iMCMC] = BIC[iMCMC]+sum(log(SS))
    
    AIC[iMCMC] = -2*BIC[iMCMC]+2*K+(K-1)
    BIC[iMCMC] = -2*BIC[iMCMC]+(2*K+(K-1))*log(nT)
    
    
    
    
  }	
  
  
  
  #####
  out = list(muOUT=muOUT,mu2OUT=mu2OUT,sigma2OUT=sigma2OUT,
             betaOUT=betaOUT,
             stOUT=stOUT, 
             BIC=max(BIC),
             piOUT=piOUT,nvec=nvec,nvecOUT=nvecOUT)
  
  
  
  return(out)
  
}





