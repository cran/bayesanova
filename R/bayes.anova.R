bayes.anova <- function(n=10000,first,second,third,fourth=NULL,fifth=NULL,sixth=NULL,hyperpars="custom",burnin=n/2,sd="sd",q=0.1,ci=0.95){
  
  # utility function for computing the mode
  getmode <- function(v) {
    uniqv <- unique(v)
    uniqv[which.max(tabulate(match(v, uniqv)))]
  }
  
  Nsim=n # number of steps performed by the Gibbs sampler
  if(length(first)>1){
    firstComponent = first # set first component
  } else {
    warning("Please provide at least data for three groups. Data for first group is missing.")
  }
  
  if(length(second)>1){
    secondComponent = second # set second component
  } else {
    warning("Please provide at least data for three groups. Data for second group is missing.")
  }
  
  if(length(third)>1){
    thirdComponent = third # set second component
  } else {
    warning("Please provide at least data for three groups. Data for third group is missing.")
  }
  
  if(length(fourth)>1){
    fourthComponent = fourth # set second component
  }
  if(length(fifth)>1){
    fifthComponent = fifth # set second component
  }
  if(length(sixth)>1){
    sixthComponent = sixth # set second component
  }
  
  mu1 <- numeric(Nsim) # arrays for parameters in first chain
  mu2 <- numeric(Nsim)
  mu3 <- numeric(Nsim)
  if(length(fourth)>1){
    mu4 <- numeric(Nsim)
  }
  if(length(fifth)>1){
    mu5 <- numeric(Nsim)
  }
  if(length(sixth)>1){
    mu6 <- numeric(Nsim)
  }
  sigma1Sq <- numeric(Nsim)
  sigma2Sq <- numeric(Nsim)
  sigma3Sq <- numeric(Nsim)
  if(length(fourth)>1){
    sigma4Sq <- numeric(Nsim)
  }
  if(length(fifth)>1){
    sigma5Sq <- numeric(Nsim)
  }
  if(length(sixth)>1){
    sigma6Sq <- numeric(Nsim)
  }
  
  mu1SecC <- numeric(Nsim) # arrays for parameters in first chain
  mu2SecC <- numeric(Nsim)
  mu3SecC <- numeric(Nsim)
  if(length(fourth)>1){
    mu4SecC <- numeric(Nsim)
  }
  if(length(fifth)>1){
    mu5SecC <- numeric(Nsim)
  }
  if(length(sixth)>1){
    mu6SecC <- numeric(Nsim)
  }
  sigma1SqSecC <- numeric(Nsim)
  sigma2SqSecC <- numeric(Nsim)
  sigma3SqSecC <- numeric(Nsim)
  if(length(fourth)>1){
    sigma4SqSecC <- numeric(Nsim)
  }
  if(length(fifth)>1){
    sigma5SqSecC <- numeric(Nsim)
  }
  if(length(sixth)>1){
    sigma6SqSecC <- numeric(Nsim)
  }
  
  mu1ThdC <- numeric(Nsim) # arrays for parameters in first chain
  mu2ThdC <- numeric(Nsim)
  mu3ThdC <- numeric(Nsim)
  if(length(fourth)>1){
    mu4ThdC <- numeric(Nsim)
  }
  if(length(fifth)>1){
    mu5ThdC <- numeric(Nsim)
  }
  if(length(sixth)>1){
    mu6ThdC <- numeric(Nsim)
  }
  sigma1SqThdC <- numeric(Nsim)
  sigma2SqThdC <- numeric(Nsim)
  sigma3SqThdC <- numeric(Nsim)
  if(length(fourth)>1){
    sigma4SqThdC <- numeric(Nsim)
  }
  if(length(fifth)>1){
    sigma5SqThdC <- numeric(Nsim)
  }
  if(length(sixth)>1){
    sigma6SqThdC <- numeric(Nsim)
  }
  
  mu1FrtC <- numeric(Nsim) # arrays for parameters in first chain
  mu2FrtC <- numeric(Nsim)
  mu3FrtC <- numeric(Nsim)
  if(length(fourth)>1){
    mu4FrtC <- numeric(Nsim)
  }
  if(length(fifth)>1){
    mu5FrtC <- numeric(Nsim)
  }
  if(length(sixth)>1){
    mu6FrtC <- numeric(Nsim)
  }
  sigma1SqFrtC <- numeric(Nsim)
  sigma2SqFrtC <- numeric(Nsim)
  sigma3SqFrtC <- numeric(Nsim)
  if(length(fourth)>1){
    sigma4SqFrtC <- numeric(Nsim)
  }
  if(length(fifth)>1){
    sigma5SqFrtC <- numeric(Nsim)
  }
  if(length(sixth)>1){
    sigma6SqFrtC <- numeric(Nsim)
  }
  
  mu1[1]=mean(firstComponent) # initialize parameters from priors of mu_k and sigma_k^2
  mu2[1]=mean(secondComponent)
  mu3[1]=mean(thirdComponent)
  if(length(fourth)>1){
    mu4[1]=mean(fourthComponent)
  }
  if(length(fifth)>1){
    mu5[1]=mean(fifthComponent)
  }
  if(length(sixth)>1){
    mu6[1]=mean(sixthComponent)
  }
  
  sigma1Sq[1]=var(firstComponent)
  sigma2Sq[1]=var(secondComponent)
  sigma3Sq[1]=var(thirdComponent)
  if(length(fourth)>1){
    sigma4Sq[1]=var(fourthComponent)
  }
  if(length(fifth)>1){
    sigma5Sq[1]=var(fifthComponent)
  }
  if(length(sixth)>1){
    sigma6Sq[1]=var(sixthComponent)
  }
  
  mu1SecC[1]=20*mean(firstComponent) # initialize parameters from priors of mu_k and sigma_k^2
  mu2SecC[1]=20*mean(secondComponent)
  mu3SecC[1]=20*mean(thirdComponent)
  if(length(fourth)>1){
    mu4SecC[1]=20*mean(fourthComponent)
  }
  if(length(fifth)>1){
    mu5SecC[1]=20*mean(fifthComponent)
  }
  if(length(sixth)>1){
    mu6SecC[1]=20*mean(sixthComponent)
  }
  sigma1SqSecC[1]=20*var(firstComponent)
  sigma2SqSecC[1]=20*var(secondComponent)
  sigma3SqSecC[1]=20*var(thirdComponent)
  if(length(fourth)>1){
    sigma4SqSecC[1]=20*var(fourthComponent)
  }
  if(length(fifth)>1){
    sigma5SqSecC[1]=20*var(fifthComponent)
  }
  if(length(sixth)>1){
    sigma6SqSecC[1]=20*var(sixthComponent)
  }
  
  mu1ThdC[1]=1/20*mean(firstComponent) # initialize parameters from priors of mu_k and sigma_k^2
  mu2ThdC[1]=1/20*mean(secondComponent)
  mu3ThdC[1]=1/20*mean(thirdComponent)
  if(length(fourth)>1){
    mu4ThdC[1]=1/20*mean(fourthComponent)
  }
  if(length(fifth)>1){
    mu5ThdC[1]=1/20*mean(fifthComponent)
  }
  if(length(sixth)>1){
    mu6ThdC[1]=1/20*mean(sixthComponent)
  }
  sigma1SqThdC[1]=1/20*var(firstComponent)
  sigma2SqThdC[1]=1/20*var(secondComponent)
  sigma3SqThdC[1]=1/20*var(thirdComponent)
  if(length(fourth)>1){
    sigma4SqThdC[1]=1/20*var(fourthComponent)
  }
  if(length(fifth)>1){
    sigma5SqThdC[1]=1/20*var(fifthComponent)
  }
  if(length(sixth)>1){
    sigma6SqThdC[1]=1/20*var(sixthComponent)
  }
  
  
  mu1FrtC[1]=10*mean(firstComponent) # initialize parameters from priors of mu_k and sigma_k^2
  mu2FrtC[1]=10*mean(secondComponent)
  mu3FrtC[1]=10*mean(thirdComponent)
  if(length(fourth)>1){
    mu4FrtC[1]=10*mean(fourthComponent)
  }
  if(length(fifth)>1){
    mu5FrtC[1]=10*mean(fifthComponent)
  }
  if(length(sixth)>1){
    mu6FrtC[1]=10*mean(sixthComponent)
  }
  sigma1SqFrtC[1]=10*var(firstComponent)
  sigma2SqFrtC[1]=10*var(secondComponent)
  sigma3SqFrtC[1]=10*var(thirdComponent)
  if(length(fourth)>1){
    sigma4SqFrtC[1]=10*var(fourthComponent)
  }
  if(length(fifth)>1){
    sigma5SqFrtC[1]=10*var(fifthComponent)
  }
  if(length(sixth)>1){
    sigma6SqFrtC[1]=10*var(sixthComponent)
  }
  
  # Formula parts
  N_1_S=length(firstComponent)
  N_2_S=length(secondComponent)
  N_3_S=length(thirdComponent)
  if(length(fourth)>1){
    N_4_S=length(fourthComponent)
  }
  if(length(fifth)>1){
    N_5_S=length(fifthComponent)
  }
  if(length(sixth)>1){
    N_6_S=length(sixthComponent)
  }
  vary1=var(firstComponent)
  vary2=var(secondComponent)
  vary3=var(thirdComponent)
  if(length(fourth)>1){
    vary4=var(fourthComponent)
  }
  if(length(fifth)>1){
    vary5=var(fifthComponent)
  }
  if(length(sixth)>1){
    vary6=var(sixthComponent)
  }
  bary1=mean(firstComponent)
  bary2=mean(secondComponent)
  bary3=mean(thirdComponent)
  if(length(fourth)>1){
    bary4=mean(fourthComponent)
  }
  if(length(fifth)>1){
    bary5=mean(fifthComponent)
  }
  if(length(sixth)>1){
    bary6=mean(sixthComponent)
  }
  
  # Set hyperparameters
  if(hyperpars=="rafterys"){
    if(length(first)>1 && length(second)>1 && length(third)>1){
      smple=c(firstComponent,secondComponent,thirdComponent)
    }
    if(length(first)>1 && length(second)>1 && length(third)>1 && length(fourth)>1){
      smple=c(firstComponent,secondComponent,thirdComponent,fourthComponent)
    }
    if(length(first)>1 && length(second)>1 && length(third)>1 && length(fourth)>1 && length(fifth)>1){
      smple=c(firstComponent,secondComponent,thirdComponent,fourthComponent,fifthComponent)
    }
    if(length(first)>1 && length(second)>1 && length(third)>1 && length(fourth)>1 && length(fifth)>1 && length(sixth)>1){
      smple=c(firstComponent,secondComponent,thirdComponent,fourthComponent,fifthComponent,sixthComponent)
    }
    b0=mean(smple) # Raftery's hyperparameters
    B0=var(smple)/2.6*(max(smple)-min(smple))
    #B0=sd(smple)*10
    c0=1.28 # Raftery's hyperparameters
    C0=0.36*var(smple) # Raftery's hyperparameters
    #C0=10*sd(smple)
  }
  if(hyperpars=="custom"){
    #b0=var(smple)*100*((mean(smple)/var(smple))*4) # Raftery's hyperparameters
    #B0=var(smple)*100
    #c0=1.28 # Raftery's hyperparameters
    #C0=var(smple)*2.0
    
    #b0=mean(smple) # Raftery's hyperparameters
    #B0=var(smple)/0.8*(max(smple)-min(smple))
    #c0=30*var(smple) # Raftery's hyperparameters
    #C0=900*var(smple) # Raftery's hyperparameters
    
    #b0=mean(smple) # Raftery's hyperparameters
    #B0=100*var(smple)
    #c0=50*var(smple) # Raftery's hyperparameters
    #C0=900*var(smple) # Raftery's hyperparameters
    
    # Funktionierte bisher ganz gut
    if(length(first)>1 && length(second)>1 && length(third)>1){
      smple=c(firstComponent,secondComponent,thirdComponent)
    }
    if(length(first)>1 && length(second)>1 && length(third)>1 && length(fourth)>1){
      smple=c(firstComponent,secondComponent,thirdComponent,fourthComponent)
    }
    if(length(first)>1 && length(second)>1 && length(third)>1 && length(fourth)>1 && length(fifth)>1){
      smple=c(firstComponent,secondComponent,thirdComponent,fourthComponent,fifthComponent)
    }
    if(length(first)>1 && length(second)>1 && length(third)>1 && length(fourth)>1 && length(fifth)>1 && length(sixth)>1){
      smple=c(firstComponent,secondComponent,thirdComponent,fourthComponent,fifthComponent,sixthComponent)
    }
    b0=mean(smple) # Raftery's hyperparameters
    B0=250*var(smple)
    c0=1.28 # Raftery's hyperparameters
    #C0=-0.5*var(smple)+(c0+0.5*(length(firstComp)+length(secondComp)))*5 # Raftery's hyperparameters
    if(length(first)>1 && length(second)>1 && length(third)>1){
      C0=-0.5*var(smple)+(c0+1/3*(length(firstComponent)+length(secondComponent)+length(thirdComponent)))*q
    }
    if(length(first)>1 && length(second)>1 && length(third)>1 && length(fourth)>1){
      C0=-0.5*var(smple)+(c0+1/3*(length(firstComponent)+length(secondComponent)+length(thirdComponent)+length(fourthComponent)))*q
    }
    if(length(first)>1 && length(second)>1 && length(third)>1 && length(fourth)>1 && length(fifth)>1){
      C0=-0.5*var(smple)+(c0+1/3*(length(firstComponent)+length(secondComponent)+length(thirdComponent)+length(fourthComponent)+length(fifthComponent)))*q
    }
    if(length(first)>1 && length(second)>1 && length(third)>1 && length(fourth)>1 && length(fifth)>1 && length(sixth)>1){
      C0=-0.5*var(smple)+(c0+1/3*(length(firstComponent)+length(secondComponent)+length(thirdComponent)+length(fourthComponent)+length(fifthComponent)+length(sixthComponent)))*q
    }
  }
  
  
  # Gibbs sampling via full conditionals
  for(t in 2:Nsim){
    # FIRST CHAIN
    # sample sigma_1^2|mu1,mu2,mu3,sigma_2^2,sigma_3^2,S,y
    c_1_S=c0+0.5*N_1_S
    C_1_S=C0+0.5*sum((firstComponent-mu1[t-1])^2)
    sigma1Sq[t]=MCMCpack::rinvgamma(1,shape=c_1_S,scale=C_1_S)
    
    # sample sigma_2^2|mu1,mu2,mu3,sigma_1^2,sigma_3^2,S,y
    c_2_S=c0+0.5*N_2_S
    C_2_S=C0+0.5*sum((secondComponent-mu2[t-1])^2)
    sigma2Sq[t]=MCMCpack::rinvgamma(1,shape=c_2_S,scale=C_2_S)
    
    # sample sigma_3^2|mu1,mu2,mu3,sigma_1^2,sigma_2^2,S,y
    c_3_S=c0+0.5*N_3_S
    C_3_S=C0+0.5*sum((thirdComponent-mu3[t-1])^2)
    sigma3Sq[t]=MCMCpack::rinvgamma(1,shape=c_3_S,scale=C_3_S)
    
    if(length(fourth>1)){
      # sample sigma_4^2|mu1,mu2,mu3,mu4,sigma_1^2,sigma_2^2,sigma_3^2,S,y
      c_4_S=c0+0.5*N_4_S
      C_4_S=C0+0.5*sum((fourthComponent-mu4[t-1])^2)
      sigma4Sq[t]=MCMCpack::rinvgamma(1,shape=c_4_S,scale=C_4_S)
    }
    if(length(fifth>1)){
      # sample sigma_5^2|mu1,mu2,mu3,mu4,mu5,sigma_1^2,sigma_2^2,sigma_3^2,sigma_4^2,S,y
      c_5_S=c0+0.5*N_5_S
      C_5_S=C0+0.5*sum((fifthComponent-mu5[t-1])^2)
      sigma5Sq[t]=MCMCpack::rinvgamma(1,shape=c_5_S,scale=C_5_S)
    }
    if(length(sixth>1)){
      # sample sigma_6^2|mu1,mu2,mu3,mu4,mu5,mu6,sigma_1^2,sigma_2^2,sigma_3^2,sigma_4^2,sigma_5^2,S,y
      c_6_S=c0+0.5*N_6_S
      C_6_S=C0+0.5*sum((sixthComponent-mu6[t-1])^2)
      sigma6Sq[t]=MCMCpack::rinvgamma(1,shape=c_6_S,scale=C_6_S)
    }
    
    # sample mu1|mu2,mu3,sigma_1^2,sigma_2^2,sigma_3^2,S,y
    B_1_S=1/((1/B0)+(1/sigma1Sq[t])*N_1_S) # use updated sigma1Sq[t] here for sigma_1^2
    b_1_S=B_1_S*((1/sigma1Sq[t])*N_1_S*bary1+(1/B0)*b0)
    mu1[t]=rnorm(1,mean=b_1_S,sd=B_1_S)
    
    # sample mu2|mu1,mu3,sigma_1^2,sigma_2^2,sigma_3^2,S,y
    B_2_S=1/((1/B0)+(1/sigma2Sq[t])*N_2_S) # use updated sigma2Sq[t] here for sigma_2^2
    b_2_S=B_2_S*((1/sigma2Sq[t])*N_2_S*bary2+(1/B0)*b0)
    mu2[t]=rnorm(1,mean=b_2_S,sd=B_2_S)
    
    # sample mu3|mu1,mu2,sigma_1^2,sigma_2^2,sigma_3^2,S,y
    B_3_S=1/((1/B0)+(1/sigma3Sq[t])*N_3_S) # use updated sigma3Sq[t] here for sigma_3^2
    b_3_S=B_3_S*((1/sigma3Sq[t])*N_3_S*bary3+(1/B0)*b0)
    mu3[t]=rnorm(1,mean=b_3_S,sd=B_3_S)
    
    if(length(fourth>1)){
      # sample mu4|mu1,mu2,mu3,sigma_1^2,sigma_2^2,sigma_3^2,sigma_4^2,S,y
      B_4_S=1/((1/B0)+(1/sigma4Sq[t])*N_4_S) # use updated sigma4Sq[t] here for sigma_4^2
      b_4_S=B_4_S*((1/sigma4Sq[t])*N_4_S*bary4+(1/B0)*b0)
      mu4[t]=rnorm(1,mean=b_4_S,sd=B_4_S)
    }
    if(length(fifth>1)){
      # sample mu5|mu1,mu2,mu3,mu4,sigma_1^2,sigma_2^2,sigma_3^2,sigma_4^2,sigma_5^2,S,y
      B_5_S=1/((1/B0)+(1/sigma5Sq[t])*N_5_S) # use updated sigma5Sq[t] here for sigma_5^2
      b_5_S=B_5_S*((1/sigma5Sq[t])*N_5_S*bary5+(1/B0)*b0)
      mu5[t]=rnorm(1,mean=b_5_S,sd=B_5_S)
    }
    if(length(sixth>1)){
      # sample mu6|mu1,mu2,mu3,mu4,mu5,sigma_1^2,sigma_2^2,sigma_3^2,sigma_4^2,sigma_5^2,sigma_6^2,S,y
      B_6_S=1/((1/B0)+(1/sigma6Sq[t])*N_6_S) # use updated sigma6Sq[t] here for sigma_6^2
      b_6_S=B_6_S*((1/sigma6Sq[t])*N_6_S*bary6+(1/B0)*b0)
      mu6[t]=rnorm(1,mean=b_6_S,sd=B_6_S)
    }
    
    ###############################################################
    # SECOND CHAIN
    # sample sigma_1^2|mu1,mu2,mu3,sigma_2^2,sigma_3^2,S,y
    C_1_S=C0+0.5*sum((firstComponent-mu1SecC[t-1])^2)
    sigma1SqSecC[t]=MCMCpack::rinvgamma(1,shape=c_1_S,scale=C_1_S)
    
    # sample sigma_2^2|mu1,mu2,mu3,sigma_1^2,sigma_3^2,S,y
    C_2_S=C0+0.5*sum((secondComponent-mu2SecC[t-1])^2)
    sigma2SqSecC[t]=MCMCpack::rinvgamma(1,shape=c_2_S,scale=C_2_S)
    
    # sample sigma_3^2|mu1,mu2,mu3,sigma_1^2,sigma_2^2,S,y
    C_3_S=C0+0.5*sum((thirdComponent-mu3SecC[t-1])^2)
    sigma3SqSecC[t]=MCMCpack::rinvgamma(1,shape=c_3_S,scale=C_3_S)
    
    if(length(fourth>1)){
      # sample sigma_4^2|mu1,mu2,mu3,mu4,sigma_1^2,sigma_2^2,sigma_3^2,S,y
      C_4_S=C0+0.5*sum((fourthComponent-mu4SecC[t-1])^2)
      sigma4SqSecC[t]=MCMCpack::rinvgamma(1,shape=c_4_S,scale=C_4_S)
    }
    if(length(fifth>1)){
      # sample sigma_5^2|mu1,mu2,mu3,mu4,mu5,sigma_1^2,sigma_2^2,sigma_3^2,sigma_4^2,S,y
      C_5_S=C0+0.5*sum((fifthComponent-mu5SecC[t-1])^2)
      sigma5SqSecC[t]=MCMCpack::rinvgamma(1,shape=c_5_S,scale=C_5_S)
    }
    if(length(sixth>1)){
      # sample sigma_6^2|mu1,mu2,mu3,mu4,mu5,mu6,sigma_1^2,sigma_2^2,sigma_3^2,sigma_4^2,sigma_5^2,S,y
      C_6_S=C0+0.5*sum((sixthComponent-mu6SecC[t-1])^2)
      sigma6SqSecC[t]=MCMCpack::rinvgamma(1,shape=c_6_S,scale=C_6_S)
    }
    
    
    
    # sample mu1|mu2,mu3,sigma_1^2,sigma_2^2,sigma_3^2,S,y
    B_1_S=1/((1/B0)+(1/sigma1SqSecC[t])*N_1_S) # use updated sigma1Sq[t] here for sigma_1^2
    b_1_S=B_1_S*((1/sigma1SqSecC[t])*N_1_S*bary1+(1/B0)*b0)
    mu1SecC[t]=rnorm(1,mean=b_1_S,sd=B_1_S)
    
    # sample mu2|mu1,mu3,sigma_1^2,sigma_2^2,sigma_3^2,S,y
    B_2_S=1/((1/B0)+(1/sigma2SqSecC[t])*N_2_S) # use updated sigma2Sq[t] here for sigma_2^2
    b_2_S=B_2_S*((1/sigma2SqSecC[t])*N_2_S*bary2+(1/B0)*b0)
    mu2SecC[t]=rnorm(1,mean=b_2_S,sd=B_2_S)
    
    # sample mu3|mu1,mu2,sigma_1^2,sigma_2^2,sigma_3^2,S,y
    B_3_S=1/((1/B0)+(1/sigma3SqSecC[t])*N_3_S) # use updated sigma3Sq[t] here for sigma_3^2
    b_3_S=B_3_S*((1/sigma3SqSecC[t])*N_3_S*bary3+(1/B0)*b0)
    mu3SecC[t]=rnorm(1,mean=b_3_S,sd=B_3_S)
    
    if(length(fourth>1)){
      # sample mu4|mu1,mu2,mu3,sigma_1^2,sigma_2^2,sigma_3^2,sigma_4^2,S,y
      B_4_S=1/((1/B0)+(1/sigma4SqSecC[t])*N_4_S) # use updated sigma4Sq[t] here for sigma_4^2
      b_4_S=B_4_S*((1/sigma4SqSecC[t])*N_4_S*bary4+(1/B0)*b0)
      mu4SecC[t]=rnorm(1,mean=b_4_S,sd=B_4_S)
    }
    if(length(fifth>1)){
      # sample mu5|mu1,mu2,mu3,mu4,sigma_1^2,sigma_2^2,sigma_3^2,sigma_4^2,sigma_5^2,S,y
      B_5_S=1/((1/B0)+(1/sigma5SqSecC[t])*N_5_S) # use updated sigma5Sq[t] here for sigma_5^2
      b_5_S=B_5_S*((1/sigma5SqSecC[t])*N_5_S*bary5+(1/B0)*b0)
      mu5SecC[t]=rnorm(1,mean=b_5_S,sd=B_5_S)
    }
    if(length(sixth>1)){
      # sample mu6|mu1,mu2,mu3,mu4,mu5,sigma_1^2,sigma_2^2,sigma_3^2,sigma_4^2,sigma_5^2,sigma_6^2,S,y
      B_6_S=1/((1/B0)+(1/sigma6SqSecC[t])*N_6_S) # use updated sigma6Sq[t] here for sigma_6^2
      b_6_S=B_6_S*((1/sigma6SqSecC[t])*N_6_S*bary6+(1/B0)*b0)
      mu6SecC[t]=rnorm(1,mean=b_6_S,sd=B_6_S)
    }
    
    ###############################################################
    # THIRD CHAIN
    # sample sigma_1^2|mu1,mu2,mu3,sigma_2^2,sigma_3^2,S,y
    C_1_S=C0+0.5*sum((firstComponent-mu1ThdC[t-1])^2)
    sigma1SqThdC[t]=MCMCpack::rinvgamma(1,shape=c_1_S,scale=C_1_S)
    
    # sample sigma_2^2|mu1,mu2,mu3,sigma_1^2,sigma_3^2,S,y
    C_2_S=C0+0.5*sum((secondComponent-mu2ThdC[t-1])^2)
    sigma2SqThdC[t]=MCMCpack::rinvgamma(1,shape=c_2_S,scale=C_2_S)
    
    # sample sigma_3^2|mu1,mu2,mu3,sigma_1^2,sigma_2^2,S,y
    C_3_S=C0+0.5*sum((thirdComponent-mu3ThdC[t-1])^2)
    sigma3SqThdC[t]=MCMCpack::rinvgamma(1,shape=c_3_S,scale=C_3_S)
    
    if(length(fourth>1)){
      # sample sigma_4^2|mu1,mu2,mu3,mu4,sigma_1^2,sigma_2^2,sigma_3^2,S,y
      C_4_S=C0+0.5*sum((fourthComponent-mu4ThdC[t-1])^2)
      sigma4SqThdC[t]=MCMCpack::rinvgamma(1,shape=c_4_S,scale=C_4_S)
    }
    if(length(fifth>1)){
      # sample sigma_5^2|mu1,mu2,mu3,mu4,mu5,sigma_1^2,sigma_2^2,sigma_3^2,sigma_4^2,S,y
      C_5_S=C0+0.5*sum((fifthComponent-mu5ThdC[t-1])^2)
      sigma5SqThdC[t]=MCMCpack::rinvgamma(1,shape=c_5_S,scale=C_5_S)
    }
    if(length(sixth>1)){
      # sample sigma_6^2|mu1,mu2,mu3,mu4,mu5,mu6,sigma_1^2,sigma_2^2,sigma_3^2,sigma_4^2,sigma_5^2,S,y
      C_6_S=C0+0.5*sum((sixthComponent-mu6ThdC[t-1])^2)
      sigma6SqThdC[t]=MCMCpack::rinvgamma(1,shape=c_6_S,scale=C_6_S)
    }
    
    # sample mu1|mu2,mu3,sigma_1^2,sigma_2^2,sigma_3^2,S,y
    B_1_S=1/((1/B0)+(1/sigma1SqThdC[t])*N_1_S) # use updated sigma1Sq[t] here for sigma_1^2
    b_1_S=B_1_S*((1/sigma1SqThdC[t])*N_1_S*bary1+(1/B0)*b0)
    mu1ThdC[t]=rnorm(1,mean=b_1_S,sd=B_1_S)
    
    # sample mu2|mu1,mu3,sigma_1^2,sigma_2^2,sigma_3^2,S,y
    B_2_S=1/((1/B0)+(1/sigma2SqThdC[t])*N_2_S) # use updated sigma2Sq[t] here for sigma_2^2
    b_2_S=B_2_S*((1/sigma2SqThdC[t])*N_2_S*bary2+(1/B0)*b0)
    mu2ThdC[t]=rnorm(1,mean=b_2_S,sd=B_2_S)
    
    # sample mu3|mu1,mu2,sigma_1^2,sigma_2^2,sigma_3^2,S,y
    B_3_S=1/((1/B0)+(1/sigma3SqThdC[t])*N_3_S) # use updated sigma3Sq[t] here for sigma_3^2
    b_3_S=B_3_S*((1/sigma3SqThdC[t])*N_3_S*bary3+(1/B0)*b0)
    mu3ThdC[t]=rnorm(1,mean=b_3_S,sd=B_3_S)
    
    if(length(fourth>1)){
      # sample mu4|mu1,mu2,mu3,sigma_1^2,sigma_2^2,sigma_3^2,sigma_4^2,S,y
      B_4_S=1/((1/B0)+(1/sigma4SqThdC[t])*N_4_S) # use updated sigma4Sq[t] here for sigma_4^2
      b_4_S=B_4_S*((1/sigma4SqThdC[t])*N_4_S*bary4+(1/B0)*b0)
      mu4ThdC[t]=rnorm(1,mean=b_4_S,sd=B_4_S)
    }
    if(length(fifth>1)){
      # sample mu5|mu1,mu2,mu3,mu4,sigma_1^2,sigma_2^2,sigma_3^2,sigma_4^2,sigma_5^2,S,y
      B_5_S=1/((1/B0)+(1/sigma5SqThdC[t])*N_5_S) # use updated sigma5Sq[t] here for sigma_5^2
      b_5_S=B_5_S*((1/sigma5SqThdC[t])*N_5_S*bary5+(1/B0)*b0)
      mu5ThdC[t]=rnorm(1,mean=b_5_S,sd=B_5_S)
    }
    if(length(sixth>1)){
      # sample mu6|mu1,mu2,mu3,mu4,mu5,sigma_1^2,sigma_2^2,sigma_3^2,sigma_4^2,sigma_5^2,sigma_6^2,S,y
      B_6_S=1/((1/B0)+(1/sigma6SqThdC[t])*N_6_S) # use updated sigma6Sq[t] here for sigma_6^2
      b_6_S=B_6_S*((1/sigma6SqThdC[t])*N_6_S*bary6+(1/B0)*b0)
      mu6ThdC[t]=rnorm(1,mean=b_6_S,sd=B_6_S)
    }
    
    ###############################################################
    # FOURTH CHAIN
    # sample sigma_1^2|mu1,mu2,mu3,sigma_2^2,sigma_3^2,S,y
    C_1_S=C0+0.5*sum((firstComponent-mu1FrtC[t-1])^2)
    sigma1SqFrtC[t]=MCMCpack::rinvgamma(1,shape=c_1_S,scale=C_1_S)
    
    # sample sigma_2^2|mu1,mu2,mu3,sigma_1^2,sigma_3^2,S,y
    C_2_S=C0+0.5*sum((secondComponent-mu2FrtC[t-1])^2)
    sigma2SqFrtC[t]=MCMCpack::rinvgamma(1,shape=c_2_S,scale=C_2_S)
    
    # sample sigma_3^2|mu1,mu2,mu3,sigma_1^2,sigma_2^2,S,y
    C_3_S=C0+0.5*sum((thirdComponent-mu3FrtC[t-1])^2)
    sigma3SqFrtC[t]=MCMCpack::rinvgamma(1,shape=c_3_S,scale=C_3_S)
    
    if(length(fourth>1)){
      # sample sigma_4^2|mu1,mu2,mu3,mu4,sigma_1^2,sigma_2^2,sigma_3^2,S,y
      C_4_S=C0+0.5*sum((fourthComponent-mu4FrtC[t-1])^2)
      sigma4SqFrtC[t]=MCMCpack::rinvgamma(1,shape=c_4_S,scale=C_4_S)
    }
    if(length(fifth>1)){
      # sample sigma_5^2|mu1,mu2,mu3,mu4,mu5,sigma_1^2,sigma_2^2,sigma_3^2,sigma_4^2,S,y
      C_5_S=C0+0.5*sum((fifthComponent-mu5FrtC[t-1])^2)
      sigma5SqFrtC[t]=MCMCpack::rinvgamma(1,shape=c_5_S,scale=C_5_S)
    }
    if(length(sixth>1)){
      # sample sigma_6^2|mu1,mu2,mu3,mu4,mu5,mu6,sigma_1^2,sigma_2^2,sigma_3^2,sigma_4^2,sigma_5^2,S,y
      C_6_S=C0+0.5*sum((sixthComponent-mu6FrtC[t-1])^2)
      sigma6SqFrtC[t]=MCMCpack::rinvgamma(1,shape=c_6_S,scale=C_6_S)
    }
    
    # sample mu1|mu2,mu3,sigma_1^2,sigma_2^2,sigma_3^2,S,y
    B_1_S=1/((1/B0)+(1/sigma1SqFrtC[t])*N_1_S) # use updated sigma1Sq[t] here for sigma_1^2
    b_1_S=B_1_S*((1/sigma1SqFrtC[t])*N_1_S*bary1+(1/B0)*b0)
    mu1FrtC[t]=rnorm(1,mean=b_1_S,sd=B_1_S)
    
    # sample mu2|mu1,mu3,sigma_1^2,sigma_2^2,sigma_3^2,S,y
    B_2_S=1/((1/B0)+(1/sigma2SqFrtC[t])*N_2_S) # use updated sigma2Sq[t] here for sigma_2^2
    b_2_S=B_2_S*((1/sigma2SqFrtC[t])*N_2_S*bary2+(1/B0)*b0)
    mu2FrtC[t]=rnorm(1,mean=b_2_S,sd=B_2_S)
    
    # sample mu3|mu1,mu2,sigma_1^2,sigma_2^2,sigma_3^2,S,y
    B_3_S=1/((1/B0)+(1/sigma3SqFrtC[t])*N_3_S) # use updated sigma3Sq[t] here for sigma_3^2
    b_3_S=B_3_S*((1/sigma3SqFrtC[t])*N_3_S*bary3+(1/B0)*b0)
    mu3FrtC[t]=rnorm(1,mean=b_3_S,sd=B_3_S)
    
    if(length(fourth>1)){
      # sample mu4|mu1,mu2,mu3,sigma_1^2,sigma_2^2,sigma_3^2,sigma_4^2,S,y
      B_4_S=1/((1/B0)+(1/sigma4SqFrtC[t])*N_4_S) # use updated sigma4Sq[t] here for sigma_4^2
      b_4_S=B_4_S*((1/sigma4SqFrtC[t])*N_4_S*bary4+(1/B0)*b0)
      mu4FrtC[t]=rnorm(1,mean=b_4_S,sd=B_4_S)
    }
    if(length(fifth>1)){
      # sample mu5|mu1,mu2,mu3,mu4,sigma_1^2,sigma_2^2,sigma_3^2,sigma_4^2,sigma_5^2,S,y
      B_5_S=1/((1/B0)+(1/sigma5SqFrtC[t])*N_5_S) # use updated sigma5Sq[t] here for sigma_5^2
      b_5_S=B_5_S*((1/sigma5SqFrtC[t])*N_5_S*bary5+(1/B0)*b0)
      mu5FrtC[t]=rnorm(1,mean=b_5_S,sd=B_5_S)
    }
    if(length(sixth>1)){
      # sample mu6|mu1,mu2,mu3,mu4,mu5,sigma_1^2,sigma_2^2,sigma_3^2,sigma_4^2,sigma_5^2,sigma_6^2,S,y
      B_6_S=1/((1/B0)+(1/sigma6SqFrtC[t])*N_6_S) # use updated sigma6Sq[t] here for sigma_6^2
      b_6_S=B_6_S*((1/sigma6SqFrtC[t])*N_6_S*bary6+(1/B0)*b0)
      mu6FrtC[t]=rnorm(1,mean=b_6_S,sd=B_6_S)
    }
    
    ##############################################################
    # All chains sample, reiterate
  }
  
  # Adapt chains for selected burnin
  mu1=mu1[burnin:Nsim]
  mu1SecC=mu1SecC[burnin:Nsim]
  mu1ThdC=mu1ThdC[burnin:Nsim]
  mu1FrtC=mu1FrtC[burnin:Nsim]
  
  mu2=mu2[burnin:Nsim]
  mu2SecC=mu2SecC[burnin:Nsim]
  mu2ThdC=mu2ThdC[burnin:Nsim]
  mu2FrtC=mu2FrtC[burnin:Nsim]
  
  mu3=mu3[burnin:Nsim]
  mu3SecC=mu3SecC[burnin:Nsim]
  mu3ThdC=mu3ThdC[burnin:Nsim]
  mu3FrtC=mu3FrtC[burnin:Nsim]
  
  if(length(fourth)>1){
    mu4=mu4[burnin:Nsim]
    mu4SecC=mu4SecC[burnin:Nsim]
    mu4ThdC=mu4ThdC[burnin:Nsim]
    mu4FrtC=mu4FrtC[burnin:Nsim]
  }
  
  if(length(fifth)>1){
    mu5=mu5[burnin:Nsim]
    mu5SecC=mu5SecC[burnin:Nsim]
    mu5ThdC=mu5ThdC[burnin:Nsim]
    mu5FrtC=mu5FrtC[burnin:Nsim]
  }
  
  if(length(sixth)>1){
    mu6=mu6[burnin:Nsim]
    mu6SecC=mu6SecC[burnin:Nsim]
    mu6ThdC=mu6ThdC[burnin:Nsim]
    mu6FrtC=mu6FrtC[burnin:Nsim]
  }
  
  sigma1Sq=sigma1Sq[burnin:Nsim]
  sigma1SqSecC=sigma1SqSecC[burnin:Nsim]
  sigma1SqThdC=sigma1SqThdC[burnin:Nsim]
  sigma1SqFrtC=sigma1SqFrtC[burnin:Nsim]
  
  sigma2Sq=sigma2Sq[burnin:Nsim]
  sigma2SqSecC=sigma2SqSecC[burnin:Nsim]
  sigma2SqThdC=sigma2SqThdC[burnin:Nsim]
  sigma2SqFrtC=sigma2SqFrtC[burnin:Nsim]
  
  sigma3Sq=sigma3Sq[burnin:Nsim]
  sigma3SqSecC=sigma3SqSecC[burnin:Nsim]
  sigma3SqThdC=sigma3SqThdC[burnin:Nsim]
  sigma3SqFrtC=sigma3SqFrtC[burnin:Nsim]
  
  if(length(fourth)>1){
    sigma4Sq=sigma4Sq[burnin:Nsim]
    sigma4SqSecC=sigma4SqSecC[burnin:Nsim]
    sigma4SqThdC=sigma4SqThdC[burnin:Nsim]
    sigma4SqFrtC=sigma4SqFrtC[burnin:Nsim]
  }
  
  if(length(fifth)>1){
    sigma5Sq=sigma5Sq[burnin:Nsim]
    sigma5SqSecC=sigma5SqSecC[burnin:Nsim]
    sigma5SqThdC=sigma5SqThdC[burnin:Nsim]
    sigma5SqFrtC=sigma5SqFrtC[burnin:Nsim]
  }
  
  if(length(sixth)>1){
    sigma6Sq=sigma6Sq[burnin:Nsim]
    sigma6SqSecC=sigma6SqSecC[burnin:Nsim]
    sigma6SqThdC=sigma6SqThdC[burnin:Nsim]
    sigma6SqFrtC=sigma6SqFrtC[burnin:Nsim]
  }
  
  if(sd=="sd"){
    sigma1Sq=sqrt(sigma1Sq)
    sigma1SqSecC=sqrt(sigma1SqSecC)
    sigma1SqThdC=sqrt(sigma1SqThdC)
    sigma1SqFrtC=sqrt(sigma1SqFrtC)
    
    sigma2Sq=sqrt(sigma2Sq)
    sigma2SqSecC=sqrt(sigma2SqSecC)
    sigma2SqThdC=sqrt(sigma2SqThdC)
    sigma2SqFrtC=sqrt(sigma2SqFrtC)
    
    sigma3Sq=sqrt(sigma3Sq)
    sigma3SqSecC=sqrt(sigma3SqSecC)
    sigma3SqThdC=sqrt(sigma3SqThdC)
    sigma3SqFrtC=sqrt(sigma3SqFrtC)
    
    if(length(fourth)>1){
      sigma4Sq=sqrt(sigma4Sq)
      sigma4SqSecC=sqrt(sigma4SqSecC)
      sigma4SqThdC=sqrt(sigma4SqThdC)
      sigma4SqFrtC=sqrt(sigma4SqFrtC)
    }
    
    if(length(fifth)>1){
      sigma5Sq=sqrt(sigma5Sq)
      sigma5SqSecC=sqrt(sigma5SqSecC)
      sigma5SqThdC=sqrt(sigma5SqThdC)
      sigma5SqFrtC=sqrt(sigma5SqFrtC)
    }
    
    if(length(sixth)>1){
      sigma6Sq=sqrt(sigma6Sq)
      sigma6SqSecC=sqrt(sigma6SqSecC)
      sigma6SqThdC=sqrt(sigma6SqThdC)
      sigma6SqFrtC=sqrt(sigma6SqFrtC)
    }
  }
  
  diffOfMeansMu2MinusMu1=mu2-mu1
  diffOfMeansMu2MinusMu1SecC=mu2SecC-mu1SecC
  diffOfMeansMu2MinusMu1ThdC=mu2ThdC-mu1ThdC
  diffOfMeansMu2MinusMu1FrtC=mu2FrtC-mu1FrtC
  
  diffOfMeansMu3MinusMu1=mu3-mu1
  diffOfMeansMu3MinusMu1SecC=mu3SecC-mu1SecC
  diffOfMeansMu3MinusMu1ThdC=mu3ThdC-mu1ThdC
  diffOfMeansMu3MinusMu1FrtC=mu3FrtC-mu1FrtC
  
  diffOfMeansMu3MinusMu2=mu3-mu2
  diffOfMeansMu3MinusMu2SecC=mu3SecC-mu2SecC
  diffOfMeansMu3MinusMu2ThdC=mu3ThdC-mu2ThdC
  diffOfMeansMu3MinusMu2FrtC=mu3FrtC-mu2FrtC
  
  if(length(fourth)>1){
    diffOfMeansMu4MinusMu1=mu4-mu1
    diffOfMeansMu4MinusMu1SecC=mu4SecC-mu1SecC
    diffOfMeansMu4MinusMu1ThdC=mu4ThdC-mu1ThdC
    diffOfMeansMu4MinusMu1FrtC=mu4FrtC-mu1FrtC
    
    diffOfMeansMu4MinusMu2=mu4-mu2
    diffOfMeansMu4MinusMu2SecC=mu4SecC-mu2SecC
    diffOfMeansMu4MinusMu2ThdC=mu4ThdC-mu2ThdC
    diffOfMeansMu4MinusMu2FrtC=mu4FrtC-mu2FrtC
    
    diffOfMeansMu4MinusMu3=mu4-mu3
    diffOfMeansMu4MinusMu3SecC=mu4SecC-mu3SecC
    diffOfMeansMu4MinusMu3ThdC=mu4ThdC-mu3ThdC
    diffOfMeansMu4MinusMu3FrtC=mu4FrtC-mu3FrtC
  }
  
  if(length(fourth)>1 && length(fifth)>1){
    diffOfMeansMu5MinusMu1=mu5-mu1
    diffOfMeansMu5MinusMu1SecC=mu5SecC-mu1SecC
    diffOfMeansMu5MinusMu1ThdC=mu5ThdC-mu1ThdC
    diffOfMeansMu5MinusMu1FrtC=mu5FrtC-mu1FrtC
    
    diffOfMeansMu5MinusMu2=mu5-mu2
    diffOfMeansMu5MinusMu2SecC=mu5SecC-mu2SecC
    diffOfMeansMu5MinusMu2ThdC=mu5ThdC-mu2ThdC
    diffOfMeansMu5MinusMu2FrtC=mu5FrtC-mu2FrtC
    
    diffOfMeansMu5MinusMu3=mu5-mu3
    diffOfMeansMu5MinusMu3SecC=mu5SecC-mu3SecC
    diffOfMeansMu5MinusMu3ThdC=mu5ThdC-mu3ThdC
    diffOfMeansMu5MinusMu3FrtC=mu5FrtC-mu3FrtC
    
    diffOfMeansMu5MinusMu4=mu5-mu4
    diffOfMeansMu5MinusMu4SecC=mu5SecC-mu4SecC
    diffOfMeansMu5MinusMu4ThdC=mu5ThdC-mu4ThdC
    diffOfMeansMu5MinusMu4FrtC=mu5FrtC-mu4FrtC
  }
  
  if(length(fourth)>1 && length(fifth)>1 && length(sixth)>1){
    diffOfMeansMu6MinusMu1=mu6-mu1
    diffOfMeansMu6MinusMu1SecC=mu6SecC-mu1SecC
    diffOfMeansMu6MinusMu1ThdC=mu6ThdC-mu1ThdC
    diffOfMeansMu6MinusMu1FrtC=mu6FrtC-mu1FrtC
    
    diffOfMeansMu6MinusMu2=mu6-mu2
    diffOfMeansMu6MinusMu2SecC=mu6SecC-mu2SecC
    diffOfMeansMu6MinusMu2ThdC=mu6ThdC-mu2ThdC
    diffOfMeansMu6MinusMu2FrtC=mu6FrtC-mu2FrtC
    
    diffOfMeansMu6MinusMu3=mu6-mu3
    diffOfMeansMu6MinusMu3SecC=mu6SecC-mu3SecC
    diffOfMeansMu6MinusMu3ThdC=mu6ThdC-mu3ThdC
    diffOfMeansMu6MinusMu3FrtC=mu6FrtC-mu3FrtC
    
    diffOfMeansMu6MinusMu4=mu6-mu4
    diffOfMeansMu6MinusMu4SecC=mu6SecC-mu4SecC
    diffOfMeansMu6MinusMu4ThdC=mu6ThdC-mu4ThdC
    diffOfMeansMu6MinusMu4FrtC=mu6FrtC-mu4FrtC
    
    diffOfMeansMu6MinusMu5=mu6-mu5
    diffOfMeansMu6MinusMu5SecC=mu6SecC-mu5SecC
    diffOfMeansMu6MinusMu5ThdC=mu6ThdC-mu5ThdC
    diffOfMeansMu6MinusMu5FrtC=mu6FrtC-mu5FrtC
  }
  
  diffOfVariancesS2MinusS1=sigma2Sq-sigma1Sq
  diffOfVariancesS2MinusS1SecC=sigma2SqSecC-sigma1SqSecC
  diffOfVariancesS2MinusS1ThdC=sigma2SqThdC-sigma1SqThdC
  diffOfVariancesS2MinusS1FrtC=sigma2SqFrtC-sigma1SqFrtC
  
  diffOfVariancesS3MinusS1=sigma3Sq-sigma1Sq
  diffOfVariancesS3MinusS1SecC=sigma3SqSecC-sigma1SqSecC
  diffOfVariancesS3MinusS1ThdC=sigma3SqThdC-sigma1SqThdC
  diffOfVariancesS3MinusS1FrtC=sigma3SqFrtC-sigma1SqFrtC
  
  diffOfVariancesS3MinusS2=sigma3Sq-sigma2Sq
  diffOfVariancesS3MinusS2SecC=sigma3SqSecC-sigma2SqSecC
  diffOfVariancesS3MinusS2ThdC=sigma3SqThdC-sigma2SqThdC
  diffOfVariancesS3MinusS2FrtC=sigma3SqFrtC-sigma2SqFrtC
  
  if(length(fourth)>1){
    diffOfVariancesS4MinusS1=sigma4Sq-sigma1Sq
    diffOfVariancesS4MinusS1SecC=sigma4SqSecC-sigma1SqSecC
    diffOfVariancesS4MinusS1ThdC=sigma4SqThdC-sigma1SqThdC
    diffOfVariancesS4MinusS1FrtC=sigma4SqFrtC-sigma1SqFrtC
    
    diffOfVariancesS4MinusS2=sigma4Sq-sigma2Sq
    diffOfVariancesS4MinusS2SecC=sigma4SqSecC-sigma2SqSecC
    diffOfVariancesS4MinusS2ThdC=sigma4SqThdC-sigma2SqThdC
    diffOfVariancesS4MinusS2FrtC=sigma4SqFrtC-sigma2SqFrtC
    
    diffOfVariancesS4MinusS3=sigma4Sq-sigma3Sq
    diffOfVariancesS4MinusS3SecC=sigma4SqSecC-sigma3SqSecC
    diffOfVariancesS4MinusS3ThdC=sigma4SqThdC-sigma3SqThdC
    diffOfVariancesS4MinusS3FrtC=sigma4SqFrtC-sigma3SqFrtC
  }
  
  if(length(fourth)>1 && length(fifth)>1){
    diffOfVariancesS5MinusS1=sigma5Sq-sigma1Sq
    diffOfVariancesS5MinusS1SecC=sigma5SqSecC-sigma1SqSecC
    diffOfVariancesS5MinusS1ThdC=sigma5SqThdC-sigma1SqThdC
    diffOfVariancesS5MinusS1FrtC=sigma5SqFrtC-sigma1SqFrtC
    
    diffOfVariancesS5MinusS2=sigma5Sq-sigma2Sq
    diffOfVariancesS5MinusS2SecC=sigma5SqSecC-sigma2SqSecC
    diffOfVariancesS5MinusS2ThdC=sigma5SqThdC-sigma2SqThdC
    diffOfVariancesS5MinusS2FrtC=sigma5SqFrtC-sigma2SqFrtC
    
    diffOfVariancesS5MinusS3=sigma5Sq-sigma3Sq
    diffOfVariancesS5MinusS3SecC=sigma5SqSecC-sigma3SqSecC
    diffOfVariancesS5MinusS3ThdC=sigma5SqThdC-sigma3SqThdC
    diffOfVariancesS5MinusS3FrtC=sigma5SqFrtC-sigma3SqFrtC
    
    diffOfVariancesS5MinusS4=sigma5Sq-sigma4Sq
    diffOfVariancesS5MinusS4SecC=sigma5SqSecC-sigma4SqSecC
    diffOfVariancesS5MinusS4ThdC=sigma5SqThdC-sigma4SqThdC
    diffOfVariancesS5MinusS4FrtC=sigma5SqFrtC-sigma4SqFrtC
  }
  
  if(length(fourth)>1 && length(fifth)>1 && length(sixth)>1){
    diffOfVariancesS6MinusS1=sigma6Sq-sigma1Sq
    diffOfVariancesS6MinusS1SecC=sigma6SqSecC-sigma1SqSecC
    diffOfVariancesS6MinusS1ThdC=sigma6SqThdC-sigma1SqThdC
    diffOfVariancesS6MinusS1FrtC=sigma6SqFrtC-sigma1SqFrtC
    
    diffOfVariancesS6MinusS2=sigma6Sq-sigma2Sq
    diffOfVariancesS6MinusS2SecC=sigma6SqSecC-sigma2SqSecC
    diffOfVariancesS6MinusS2ThdC=sigma6SqThdC-sigma2SqThdC
    diffOfVariancesS6MinusS2FrtC=sigma6SqFrtC-sigma2SqFrtC
    
    diffOfVariancesS6MinusS3=sigma6Sq-sigma3Sq
    diffOfVariancesS6MinusS3SecC=sigma6SqSecC-sigma3SqSecC
    diffOfVariancesS6MinusS3ThdC=sigma6SqThdC-sigma3SqThdC
    diffOfVariancesS6MinusS3FrtC=sigma6SqFrtC-sigma3SqFrtC
    
    diffOfVariancesS6MinusS4=sigma6Sq-sigma4Sq
    diffOfVariancesS6MinusS4SecC=sigma6SqSecC-sigma4SqSecC
    diffOfVariancesS6MinusS4ThdC=sigma6SqThdC-sigma4SqThdC
    diffOfVariancesS6MinusS4FrtC=sigma6SqFrtC-sigma4SqFrtC
    
    diffOfVariancesS6MinusS5=sigma6Sq-sigma5Sq
    diffOfVariancesS6MinusS5SecC=sigma6SqSecC-sigma5SqSecC
    diffOfVariancesS6MinusS5ThdC=sigma6SqThdC-sigma5SqThdC
    diffOfVariancesS6MinusS5FrtC=sigma6SqFrtC-sigma5SqFrtC
  }
  
  # Effect sizes
  if(sd=="var"){
    effectSize21=(mu1-mu2)/(sqrt(((N_1_S-1)*sigma1Sq+(N_2_S-1)*sigma2Sq)/(N_1_S+N_2_S-2)))
    effectSize21SecC=(mu1SecC-mu2SecC)/(sqrt(((N_1_S-1)*sigma1SqSecC+(N_2_S-1)*sigma2SqSecC)/(N_1_S+N_2_S-2)))
    effectSize21ThdC=(mu1ThdC-mu2ThdC)/(sqrt(((N_1_S-1)*sigma1SqThdC+(N_2_S-1)*sigma2SqThdC)/(N_1_S+N_2_S-2)))
    effectSize21FrtC=(mu1FrtC-mu2FrtC)/(sqrt(((N_1_S-1)*sigma1SqFrtC+(N_2_S-1)*sigma2SqFrtC)/(N_1_S+N_2_S-2)))
    
    effectSize13=(mu1-mu3)/(sqrt(((N_1_S-1)*sigma1Sq+(N_3_S-1)*sigma3Sq)/(N_1_S+N_3_S-2)))
    effectSize13SecC=(mu1SecC-mu3SecC)/(sqrt(((N_1_S-1)*sigma1SqSecC+(N_3_S-1)*sigma3SqSecC)/(N_1_S+N_3_S-2)))
    effectSize13ThdC=(mu1ThdC-mu3ThdC)/(sqrt(((N_1_S-1)*sigma1SqThdC+(N_3_S-1)*sigma3SqThdC)/(N_1_S+N_3_S-2)))
    effectSize13FrtC=(mu1FrtC-mu3FrtC)/(sqrt(((N_1_S-1)*sigma1SqFrtC+(N_3_S-1)*sigma3SqFrtC)/(N_1_S+N_3_S-2)))
    
    effectSize23=(mu2-mu3)/(sqrt(((N_2_S-1)*sigma2Sq+(N_3_S-1)*sigma3Sq)/(N_2_S+N_3_S-2)))
    effectSize23SecC=(mu2SecC-mu3SecC)/(sqrt(((N_2_S-1)*sigma2SqSecC+(N_3_S-1)*sigma3SqSecC)/(N_2_S+N_3_S-2)))
    effectSize23ThdC=(mu2ThdC-mu3ThdC)/(sqrt(((N_2_S-1)*sigma2SqThdC+(N_3_S-1)*sigma3SqThdC)/(N_2_S+N_3_S-2)))
    effectSize23FrtC=(mu2FrtC-mu3FrtC)/(sqrt(((N_2_S-1)*sigma2SqFrtC+(N_3_S-1)*sigma3SqFrtC)/(N_2_S+N_3_S-2)))
    
    if(length(fourth>1)){
      effectSize14=(mu1-mu4)/(sqrt(((N_1_S-1)*sigma1Sq+(N_4_S-1)*sigma4Sq)/(N_1_S+N_4_S-2)))
      effectSize14SecC=(mu1SecC-mu4SecC)/(sqrt(((N_1_S-1)*sigma1SqSecC+(N_4_S-1)*sigma4SqSecC)/(N_1_S+N_4_S-2)))
      effectSize14ThdC=(mu1ThdC-mu4ThdC)/(sqrt(((N_1_S-1)*sigma1SqThdC+(N_4_S-1)*sigma4SqThdC)/(N_1_S+N_4_S-2)))
      effectSize14FrtC=(mu1FrtC-mu4FrtC)/(sqrt(((N_1_S-1)*sigma1SqFrtC+(N_4_S-1)*sigma4SqFrtC)/(N_1_S+N_4_S-2)))
      
      effectSize24=(mu2-mu4)/(sqrt(((N_2_S-1)*sigma2Sq+(N_4_S-1)*sigma4Sq)/(N_2_S+N_4_S-2)))
      effectSize24SecC=(mu2SecC-mu4SecC)/(sqrt(((N_2_S-1)*sigma2SqSecC+(N_4_S-1)*sigma4SqSecC)/(N_2_S+N_4_S-2)))
      effectSize24ThdC=(mu2ThdC-mu4ThdC)/(sqrt(((N_2_S-1)*sigma2SqThdC+(N_4_S-1)*sigma4SqThdC)/(N_2_S+N_4_S-2)))
      effectSize24FrtC=(mu2FrtC-mu4FrtC)/(sqrt(((N_2_S-1)*sigma2SqFrtC+(N_4_S-1)*sigma4SqFrtC)/(N_2_S+N_4_S-2)))
      
      effectSize34=(mu3-mu4)/(sqrt(((N_3_S-1)*sigma3Sq+(N_4_S-1)*sigma4Sq)/(N_3_S+N_4_S-2)))
      effectSize34SecC=(mu3SecC-mu4SecC)/(sqrt(((N_3_S-1)*sigma3SqSecC+(N_4_S-1)*sigma4SqSecC)/(N_3_S+N_4_S-2)))
      effectSize34ThdC=(mu3ThdC-mu4ThdC)/(sqrt(((N_3_S-1)*sigma3SqThdC+(N_4_S-1)*sigma4SqThdC)/(N_3_S+N_4_S-2)))
      effectSize34FrtC=(mu3FrtC-mu4FrtC)/(sqrt(((N_3_S-1)*sigma3SqFrtC+(N_4_S-1)*sigma4SqFrtC)/(N_3_S+N_4_S-2)))
    }
    if(length(fourth>1) && length(fifth)>1){
      effectSize15=(mu1-mu5)/(sqrt(((N_1_S-1)*sigma1Sq+(N_5_S-1)*sigma5Sq)/(N_1_S+N_5_S-2)))
      effectSize15SecC=(mu1SecC-mu5SecC)/(sqrt(((N_1_S-1)*sigma1SqSecC+(N_5_S-1)*sigma5SqSecC)/(N_1_S+N_5_S-2)))
      effectSize15ThdC=(mu1ThdC-mu5ThdC)/(sqrt(((N_1_S-1)*sigma1SqThdC+(N_5_S-1)*sigma5SqThdC)/(N_1_S+N_5_S-2)))
      effectSize15FrtC=(mu1FrtC-mu5FrtC)/(sqrt(((N_1_S-1)*sigma1SqFrtC+(N_5_S-1)*sigma5SqFrtC)/(N_1_S+N_5_S-2)))
      
      effectSize25=(mu2-mu5)/(sqrt(((N_2_S-1)*sigma2Sq+(N_5_S-1)*sigma5Sq)/(N_2_S+N_5_S-2)))
      effectSize25SecC=(mu2SecC-mu5SecC)/(sqrt(((N_2_S-1)*sigma2SqSecC+(N_5_S-1)*sigma5SqSecC)/(N_2_S+N_5_S-2)))
      effectSize25ThdC=(mu2ThdC-mu5ThdC)/(sqrt(((N_2_S-1)*sigma2SqThdC+(N_5_S-1)*sigma5SqThdC)/(N_2_S+N_5_S-2)))
      effectSize25FrtC=(mu2FrtC-mu5FrtC)/(sqrt(((N_2_S-1)*sigma2SqFrtC+(N_5_S-1)*sigma5SqFrtC)/(N_2_S+N_5_S-2)))
      
      effectSize35=(mu3-mu5)/(sqrt(((N_3_S-1)*sigma3Sq+(N_5_S-1)*sigma5Sq)/(N_3_S+N_5_S-2)))
      effectSize35SecC=(mu3SecC-mu5SecC)/(sqrt(((N_3_S-1)*sigma3SqSecC+(N_5_S-1)*sigma5SqSecC)/(N_3_S+N_5_S-2)))
      effectSize35ThdC=(mu3ThdC-mu5ThdC)/(sqrt(((N_3_S-1)*sigma3SqThdC+(N_5_S-1)*sigma5SqThdC)/(N_3_S+N_5_S-2)))
      effectSize35FrtC=(mu3FrtC-mu5FrtC)/(sqrt(((N_3_S-1)*sigma3SqFrtC+(N_5_S-1)*sigma5SqFrtC)/(N_3_S+N_5_S-2)))
      
      effectSize45=(mu4-mu5)/(sqrt(((N_4_S-1)*sigma4Sq+(N_5_S-1)*sigma5Sq)/(N_4_S+N_5_S-2)))
      effectSize45SecC=(mu4SecC-mu5SecC)/(sqrt(((N_4_S-1)*sigma4SqSecC+(N_5_S-1)*sigma5SqSecC)/(N_4_S+N_5_S-2)))
      effectSize45ThdC=(mu4ThdC-mu5ThdC)/(sqrt(((N_4_S-1)*sigma4SqThdC+(N_5_S-1)*sigma5SqThdC)/(N_4_S+N_5_S-2)))
      effectSize45FrtC=(mu4FrtC-mu5FrtC)/(sqrt(((N_4_S-1)*sigma4SqFrtC+(N_5_S-1)*sigma5SqFrtC)/(N_4_S+N_5_S-2)))
    }
    if(length(fourth>1)&& length(fifth)>1 && length(sixth)>1){
      effectSize16=(mu1-mu6)/(sqrt(((N_1_S-1)*sigma1Sq+(N_6_S-1)*sigma6Sq)/(N_1_S+N_6_S-2)))
      effectSize16SecC=(mu1SecC-mu6SecC)/(sqrt(((N_1_S-1)*sigma1SqSecC+(N_6_S-1)*sigma6SqSecC)/(N_1_S+N_6_S-2)))
      effectSize16ThdC=(mu1ThdC-mu6ThdC)/(sqrt(((N_1_S-1)*sigma1SqThdC+(N_6_S-1)*sigma6SqThdC)/(N_1_S+N_6_S-2)))
      effectSize16FrtC=(mu1FrtC-mu6FrtC)/(sqrt(((N_1_S-1)*sigma1SqFrtC+(N_6_S-1)*sigma6SqFrtC)/(N_1_S+N_6_S-2)))
      
      effectSize26=(mu2-mu6)/(sqrt(((N_2_S-1)*sigma2Sq+(N_6_S-1)*sigma6Sq)/(N_2_S+N_6_S-2)))
      effectSize26SecC=(mu2SecC-mu6SecC)/(sqrt(((N_2_S-1)*sigma2SqSecC+(N_6_S-1)*sigma6SqSecC)/(N_2_S+N_6_S-2)))
      effectSize26ThdC=(mu2ThdC-mu6ThdC)/(sqrt(((N_2_S-1)*sigma2SqThdC+(N_6_S-1)*sigma6SqThdC)/(N_2_S+N_6_S-2)))
      effectSize26FrtC=(mu2FrtC-mu6FrtC)/(sqrt(((N_2_S-1)*sigma2SqFrtC+(N_6_S-1)*sigma6SqFrtC)/(N_2_S+N_6_S-2)))
      
      effectSize36=(mu3-mu6)/(sqrt(((N_3_S-1)*sigma3Sq+(N_6_S-1)*sigma6Sq)/(N_3_S+N_6_S-2)))
      effectSize36SecC=(mu3SecC-mu6SecC)/(sqrt(((N_3_S-1)*sigma3SqSecC+(N_6_S-1)*sigma6SqSecC)/(N_3_S+N_6_S-2)))
      effectSize36ThdC=(mu3ThdC-mu6ThdC)/(sqrt(((N_3_S-1)*sigma3SqThdC+(N_6_S-1)*sigma6SqThdC)/(N_3_S+N_6_S-2)))
      effectSize36FrtC=(mu3FrtC-mu6FrtC)/(sqrt(((N_3_S-1)*sigma3SqFrtC+(N_6_S-1)*sigma6SqFrtC)/(N_3_S+N_6_S-2)))
      
      effectSize46=(mu4-mu6)/(sqrt(((N_4_S-1)*sigma4Sq+(N_6_S-1)*sigma6Sq)/(N_4_S+N_6_S-2)))
      effectSize46SecC=(mu4SecC-mu6SecC)/(sqrt(((N_4_S-1)*sigma4SqSecC+(N_6_S-1)*sigma6SqSecC)/(N_4_S+N_6_S-2)))
      effectSize46ThdC=(mu4ThdC-mu6ThdC)/(sqrt(((N_4_S-1)*sigma4SqThdC+(N_6_S-1)*sigma6SqThdC)/(N_4_S+N_6_S-2)))
      effectSize46FrtC=(mu4FrtC-mu6FrtC)/(sqrt(((N_4_S-1)*sigma4SqFrtC+(N_6_S-1)*sigma6SqFrtC)/(N_4_S+N_6_S-2)))
      
      effectSize56=(mu5-mu6)/(sqrt(((N_5_S-1)*sigma5Sq+(N_6_S-1)*sigma6Sq)/(N_5_S+N_6_S-2)))
      effectSize56SecC=(mu5SecC-mu6SecC)/(sqrt(((N_5_S-1)*sigma5SqSecC+(N_6_S-1)*sigma6SqSecC)/(N_5_S+N_6_S-2)))
      effectSize56ThdC=(mu5ThdC-mu6ThdC)/(sqrt(((N_5_S-1)*sigma5SqThdC+(N_6_S-1)*sigma6SqThdC)/(N_5_S+N_6_S-2)))
      effectSize56FrtC=(mu5FrtC-mu6FrtC)/(sqrt(((N_5_S-1)*sigma5SqFrtC+(N_6_S-1)*sigma6SqFrtC)/(N_5_S+N_6_S-2)))
    }
  }
  if(sd=="sd"){
    effectSize12=(mu1-mu2)/(sqrt(((N_1_S-1)*sigma1Sq+(N_2_S-1)*sigma2Sq)/(N_1_S+N_2_S-2)))
    effectSize12SecC=(mu1SecC-mu2SecC)/(sqrt(((N_1_S-1)*sigma1SqSecC+(N_2_S-1)*sigma2SqSecC)/(N_1_S+N_2_S-2)))
    effectSize12ThdC=(mu1ThdC-mu2ThdC)/(sqrt(((N_1_S-1)*sigma1SqThdC+(N_2_S-1)*sigma2SqThdC)/(N_1_S+N_2_S-2)))
    effectSize12FrtC=(mu1FrtC-mu2FrtC)/(sqrt(((N_1_S-1)*sigma1SqFrtC+(N_2_S-1)*sigma2SqFrtC)/(N_1_S+N_2_S-2)))
    
    effectSize13=(mu1-mu3)/(sqrt(((N_1_S-1)*sigma1Sq+(N_3_S-1)*sigma3Sq)/(N_1_S+N_3_S-2)))
    effectSize13SecC=(mu1SecC-mu3SecC)/(sqrt(((N_1_S-1)*sigma1SqSecC+(N_3_S-1)*sigma3SqSecC)/(N_1_S+N_3_S-2)))
    effectSize13ThdC=(mu1ThdC-mu3ThdC)/(sqrt(((N_1_S-1)*sigma1SqThdC+(N_3_S-1)*sigma3SqThdC)/(N_1_S+N_3_S-2)))
    effectSize13FrtC=(mu1FrtC-mu3FrtC)/(sqrt(((N_1_S-1)*sigma1SqFrtC+(N_3_S-1)*sigma3SqFrtC)/(N_1_S+N_3_S-2)))
    
    effectSize23=(mu2-mu3)/(sqrt(((N_2_S-1)*sigma2Sq+(N_3_S-1)*sigma3Sq)/(N_2_S+N_3_S-2)))
    effectSize23SecC=(mu2SecC-mu3SecC)/(sqrt(((N_2_S-1)*sigma2SqSecC+(N_3_S-1)*sigma3SqSecC)/(N_2_S+N_3_S-2)))
    effectSize23ThdC=(mu2ThdC-mu3ThdC)/(sqrt(((N_2_S-1)*sigma2SqThdC+(N_3_S-1)*sigma3SqThdC)/(N_2_S+N_3_S-2)))
    effectSize23FrtC=(mu2FrtC-mu3FrtC)/(sqrt(((N_2_S-1)*sigma2SqFrtC+(N_3_S-1)*sigma3SqFrtC)/(N_2_S+N_3_S-2)))
    
    if(length(fourth>1)){
      effectSize14=(mu1-mu4)/(sqrt(((N_1_S-1)*sigma1Sq+(N_4_S-1)*sigma4Sq)/(N_1_S+N_4_S-2)))
      effectSize14SecC=(mu1SecC-mu4SecC)/(sqrt(((N_1_S-1)*sigma1SqSecC+(N_4_S-1)*sigma4SqSecC)/(N_1_S+N_4_S-2)))
      effectSize14ThdC=(mu1ThdC-mu4ThdC)/(sqrt(((N_1_S-1)*sigma1SqThdC+(N_4_S-1)*sigma4SqThdC)/(N_1_S+N_4_S-2)))
      effectSize14FrtC=(mu1FrtC-mu4FrtC)/(sqrt(((N_1_S-1)*sigma1SqFrtC+(N_4_S-1)*sigma4SqFrtC)/(N_1_S+N_4_S-2)))
      
      effectSize24=(mu2-mu4)/(sqrt(((N_2_S-1)*sigma2Sq+(N_4_S-1)*sigma4Sq)/(N_2_S+N_4_S-2)))
      effectSize24SecC=(mu2SecC-mu4SecC)/(sqrt(((N_2_S-1)*sigma2SqSecC+(N_4_S-1)*sigma4SqSecC)/(N_2_S+N_4_S-2)))
      effectSize24ThdC=(mu2ThdC-mu4ThdC)/(sqrt(((N_2_S-1)*sigma2SqThdC+(N_4_S-1)*sigma4SqThdC)/(N_2_S+N_4_S-2)))
      effectSize24FrtC=(mu2FrtC-mu4FrtC)/(sqrt(((N_2_S-1)*sigma2SqFrtC+(N_4_S-1)*sigma4SqFrtC)/(N_2_S+N_4_S-2)))
      
      effectSize34=(mu3-mu4)/(sqrt(((N_3_S-1)*sigma3Sq+(N_4_S-1)*sigma4Sq)/(N_3_S+N_4_S-2)))
      effectSize34SecC=(mu3SecC-mu4SecC)/(sqrt(((N_3_S-1)*sigma3SqSecC+(N_4_S-1)*sigma4SqSecC)/(N_3_S+N_4_S-2)))
      effectSize34ThdC=(mu3ThdC-mu4ThdC)/(sqrt(((N_3_S-1)*sigma3SqThdC+(N_4_S-1)*sigma4SqThdC)/(N_3_S+N_4_S-2)))
      effectSize34FrtC=(mu3FrtC-mu4FrtC)/(sqrt(((N_3_S-1)*sigma3SqFrtC+(N_4_S-1)*sigma4SqFrtC)/(N_3_S+N_4_S-2)))
    }
    if(length(fourth>1) && length(fifth)>1){
      effectSize15=(mu1-mu5)/(sqrt(((N_1_S-1)*sigma1Sq^2+(N_5_S-1)*sigma5Sq^2)/(N_1_S+N_5_S-2)))
      effectSize15SecC=(mu1SecC-mu5SecC)/(sqrt(((N_1_S-1)*sigma1SqSecC^2+(N_5_S-1)*sigma5SqSecC^2)/(N_1_S+N_5_S-2)))
      effectSize15ThdC=(mu1ThdC-mu5ThdC)/(sqrt(((N_1_S-1)*sigma1SqThdC^2+(N_5_S-1)*sigma5SqThdC^2)/(N_1_S+N_5_S-2)))
      effectSize15FrtC=(mu1FrtC-mu5FrtC)/(sqrt(((N_1_S-1)*sigma1SqFrtC^2+(N_5_S-1)*sigma5SqFrtC^2)/(N_1_S+N_5_S-2)))
      
      effectSize25=(mu2-mu5)/(sqrt(((N_2_S-1)*sigma2Sq^2+(N_5_S-1)*sigma5Sq^2)/(N_2_S+N_5_S-2)))
      effectSize25SecC=(mu2SecC-mu5SecC)/(sqrt(((N_2_S-1)*sigma2SqSecC^2+(N_5_S-1)*sigma5SqSecC^2)/(N_2_S+N_5_S-2)))
      effectSize25ThdC=(mu2ThdC-mu5ThdC)/(sqrt(((N_2_S-1)*sigma2SqThdC^2+(N_5_S-1)*sigma5SqThdC^2)/(N_2_S+N_5_S-2)))
      effectSize25FrtC=(mu2FrtC-mu5FrtC)/(sqrt(((N_2_S-1)*sigma2SqFrtC^2+(N_5_S-1)*sigma5SqFrtC^2)/(N_2_S+N_5_S-2)))
      
      effectSize35=(mu3-mu5)/(sqrt(((N_3_S-1)*sigma3Sq^2+(N_5_S-1)*sigma5Sq^2)/(N_3_S+N_5_S-2)))
      effectSize35SecC=(mu3SecC-mu5SecC)/(sqrt(((N_3_S-1)*sigma3SqSecC^2+(N_5_S-1)*sigma5SqSecC^2)/(N_3_S+N_5_S-2)))
      effectSize35ThdC=(mu3ThdC-mu5ThdC)/(sqrt(((N_3_S-1)*sigma3SqThdC^2+(N_5_S-1)*sigma5SqThdC^2)/(N_3_S+N_5_S-2)))
      effectSize35FrtC=(mu3FrtC-mu5FrtC)/(sqrt(((N_3_S-1)*sigma3SqFrtC^2+(N_5_S-1)*sigma5SqFrtC^2)/(N_3_S+N_5_S-2)))
      
      effectSize45=(mu4-mu5)/(sqrt(((N_4_S-1)*sigma4Sq^2+(N_5_S-1)*sigma5Sq^2)/(N_4_S+N_5_S-2)))
      effectSize45SecC=(mu4SecC-mu5SecC)/(sqrt(((N_4_S-1)*sigma4SqSecC^2+(N_5_S-1)*sigma5SqSecC^2)/(N_4_S+N_5_S-2)))
      effectSize45ThdC=(mu4ThdC-mu5ThdC)/(sqrt(((N_4_S-1)*sigma4SqThdC^2+(N_5_S-1)*sigma5SqThdC^2)/(N_4_S+N_5_S-2)))
      effectSize45FrtC=(mu4FrtC-mu5FrtC)/(sqrt(((N_4_S-1)*sigma4SqFrtC^2+(N_5_S-1)*sigma5SqFrtC^2)/(N_4_S+N_5_S-2)))
    }
    if(length(fourth>1)&& length(fifth)>1 && length(sixth)>1){
      effectSize16=(mu1-mu6)/(sqrt(((N_1_S-1)*sigma1Sq^2+(N_6_S-1)*sigma6Sq^2)/(N_1_S+N_6_S-2)))
      effectSize16SecC=(mu1SecC-mu6SecC)/(sqrt(((N_1_S-1)*sigma1SqSecC^2+(N_6_S-1)*sigma6SqSecC^2)/(N_1_S+N_6_S-2)))
      effectSize16ThdC=(mu1ThdC-mu6ThdC)/(sqrt(((N_1_S-1)*sigma1SqThdC^2+(N_6_S-1)*sigma6SqThdC^2)/(N_1_S+N_6_S-2)))
      effectSize16FrtC=(mu1FrtC-mu6FrtC)/(sqrt(((N_1_S-1)*sigma1SqFrtC^2+(N_6_S-1)*sigma6SqFrtC^2)/(N_1_S+N_6_S-2)))
      
      effectSize26=(mu2-mu6)/(sqrt(((N_2_S-1)*sigma2Sq^2+(N_6_S-1)*sigma6Sq^2)/(N_2_S+N_6_S-2)))
      effectSize26SecC=(mu2SecC-mu6SecC)/(sqrt(((N_2_S-1)*sigma2SqSecC^2+(N_6_S-1)*sigma6SqSecC^2)/(N_2_S+N_6_S-2)))
      effectSize26ThdC=(mu2ThdC-mu6ThdC)/(sqrt(((N_2_S-1)*sigma2SqThdC^2+(N_6_S-1)*sigma6SqThdC^2)/(N_2_S+N_6_S-2)))
      effectSize26FrtC=(mu2FrtC-mu6FrtC)/(sqrt(((N_2_S-1)*sigma2SqFrtC^2+(N_6_S-1)*sigma6SqFrtC^2)/(N_2_S+N_6_S-2)))
      
      effectSize36=(mu3-mu6)/(sqrt(((N_3_S-1)*sigma3Sq^2+(N_6_S-1)*sigma6Sq^2)/(N_3_S+N_6_S-2)))
      effectSize36SecC=(mu3SecC-mu6SecC)/(sqrt(((N_3_S-1)*sigma3SqSecC^2+(N_6_S-1)*sigma6SqSecC^2)/(N_3_S+N_6_S-2)))
      effectSize36ThdC=(mu3ThdC-mu6ThdC)/(sqrt(((N_3_S-1)*sigma3SqThdC^2+(N_6_S-1)*sigma6SqThdC^2)/(N_3_S+N_6_S-2)))
      effectSize36FrtC=(mu3FrtC-mu6FrtC)/(sqrt(((N_3_S-1)*sigma3SqFrtC^2+(N_6_S-1)*sigma6SqFrtC^2)/(N_3_S+N_6_S-2)))
      
      effectSize46=(mu4-mu6)/(sqrt(((N_4_S-1)*sigma4Sq^2+(N_6_S-1)*sigma6Sq^2)/(N_4_S+N_6_S-2)))
      effectSize46SecC=(mu4SecC-mu6SecC)/(sqrt(((N_4_S-1)*sigma4SqSecC^2+(N_6_S-1)*sigma6SqSecC^2)/(N_4_S+N_6_S-2)))
      effectSize46ThdC=(mu4ThdC-mu6ThdC)/(sqrt(((N_4_S-1)*sigma4SqThdC^2+(N_6_S-1)*sigma6SqThdC^2)/(N_4_S+N_6_S-2)))
      effectSize46FrtC=(mu4FrtC-mu6FrtC)/(sqrt(((N_4_S-1)*sigma4SqFrtC^2+(N_6_S-1)*sigma6SqFrtC^2)/(N_4_S+N_6_S-2)))
      
      effectSize56=(mu5-mu6)/(sqrt(((N_5_S-1)*sigma5Sq^2+(N_6_S-1)*sigma6Sq^2)/(N_5_S+N_6_S-2)))
      effectSize56SecC=(mu5SecC-mu6SecC)/(sqrt(((N_5_S-1)*sigma5SqSecC^2+(N_6_S-1)*sigma6SqSecC^2)/(N_5_S+N_6_S-2)))
      effectSize56ThdC=(mu5ThdC-mu6ThdC)/(sqrt(((N_5_S-1)*sigma5SqThdC^2+(N_6_S-1)*sigma6SqThdC^2)/(N_5_S+N_6_S-2)))
      effectSize56FrtC=(mu5FrtC-mu6FrtC)/(sqrt(((N_5_S-1)*sigma5SqFrtC^2+(N_6_S-1)*sigma6SqFrtC^2)/(N_5_S+N_6_S-2)))
    }
  }
  
  
  ##################################### CONSOLE OUTPUT ############################################
  message("Bayesian ANOVA output:")
  if(length(first)>1 && length(second)>1 && length(third)>1 && length(fourth)>1 && length(fifth)>1 && length(sixth)>1) {
    message("Details: Gaussian-mixture model with six components")
  } else if(length(first)>1 && length(second)>1 && length(third)>1 && length(fourth)>1 && length(fifth)>1) {
    message("Details: Gaussian-mixture model with five components")
  } else if(length(first)>1 && length(second)>1 && length(third)>1 && length(fourth)>1) {
    message("Details: Gaussian-mixture model with four components")
  } else if(length(first)>1 && length(second)>1 && length(third)>1) {
    message("Details: Gaussian-mixture model with three components")
  }
  # Quantiles
  mu1Quants=quantile(mu1,probs=c((1-ci)/2,ci+(1-ci)/2))
  mu2Quants=quantile(mu2,probs=c((1-ci)/2,ci+(1-ci)/2))
  mu3Quants=quantile(mu3,probs=c((1-ci)/2,ci+(1-ci)/2))
  if(length(fourth>1)){
    mu4Quants=quantile(mu4,probs=c((1-ci)/2,ci+(1-ci)/2))
  }
  if(length(fifth>1)){
    mu5Quants=quantile(mu5,probs=c((1-ci)/2,ci+(1-ci)/2))
  }
  if(length(sixth>1)){
    mu6Quants=quantile(mu6,probs=c((1-ci)/2,ci+(1-ci)/2))
  }
  
  sigma1Quants=quantile(sigma1Sq,probs=c((1-ci)/2,ci+(1-ci)/2))
  sigma2Quants=quantile(sigma2Sq,probs=c((1-ci)/2,ci+(1-ci)/2))
  sigma3Quants=quantile(sigma3Sq,probs=c((1-ci)/2,ci+(1-ci)/2))
  if(length(fourth>1)){
    sigma4Quants=quantile(sigma4Sq,probs=c((1-ci)/2,ci+(1-ci)/2))
  }
  if(length(fifth>1)){
    sigma5Quants=quantile(sigma5Sq,probs=c((1-ci)/2,ci+(1-ci)/2))
  }
  if(length(sixth>1)){
    sigma6Quants=quantile(sigma6Sq,probs=c((1-ci)/2,ci+(1-ci)/2))
  }
  
  diffOfMeansMu2MinusMu1Quants=quantile(diffOfMeansMu2MinusMu1,probs=c((1-ci)/2,ci+(1-ci)/2))
  diffOfMeansMu3MinusMu1Quants=quantile(diffOfMeansMu3MinusMu1,probs=c((1-ci)/2,ci+(1-ci)/2))
  diffOfMeansMu3MinusMu2Quants=quantile(diffOfMeansMu3MinusMu2,probs=c((1-ci)/2,ci+(1-ci)/2))
  if(length(fourth>1)){
    diffOfMeansMu4MinusMu1Quants=quantile(diffOfMeansMu4MinusMu1,probs=c((1-ci)/2,ci+(1-ci)/2))
    diffOfMeansMu4MinusMu2Quants=quantile(diffOfMeansMu4MinusMu2,probs=c((1-ci)/2,ci+(1-ci)/2))
    diffOfMeansMu4MinusMu3Quants=quantile(diffOfMeansMu4MinusMu3,probs=c((1-ci)/2,ci+(1-ci)/2))
  }
  if(length(fifth>1)){
    diffOfMeansMu5MinusMu1Quants=quantile(diffOfMeansMu5MinusMu1,probs=c((1-ci)/2,ci+(1-ci)/2))
    diffOfMeansMu5MinusMu2Quants=quantile(diffOfMeansMu5MinusMu2,probs=c((1-ci)/2,ci+(1-ci)/2))
    diffOfMeansMu5MinusMu3Quants=quantile(diffOfMeansMu5MinusMu3,probs=c((1-ci)/2,ci+(1-ci)/2))
    diffOfMeansMu5MinusMu4Quants=quantile(diffOfMeansMu5MinusMu4,probs=c((1-ci)/2,ci+(1-ci)/2))
  }
  if(length(sixth>1)){
    diffOfMeansMu6MinusMu1Quants=quantile(diffOfMeansMu6MinusMu1,probs=c((1-ci)/2,ci+(1-ci)/2))
    diffOfMeansMu6MinusMu2Quants=quantile(diffOfMeansMu6MinusMu2,probs=c((1-ci)/2,ci+(1-ci)/2))
    diffOfMeansMu6MinusMu3Quants=quantile(diffOfMeansMu6MinusMu3,probs=c((1-ci)/2,ci+(1-ci)/2))
    diffOfMeansMu6MinusMu4Quants=quantile(diffOfMeansMu6MinusMu4,probs=c((1-ci)/2,ci+(1-ci)/2))
    diffOfMeansMu6MinusMu5Quants=quantile(diffOfMeansMu6MinusMu5,probs=c((1-ci)/2,ci+(1-ci)/2))
  }
  
  diffOfVariancesS2MinusS1Quants=quantile(diffOfVariancesS2MinusS1,probs=c((1-ci)/2,ci+(1-ci)/2))
  diffOfVariancesS3MinusS1Quants=quantile(diffOfVariancesS3MinusS1,probs=c((1-ci)/2,ci+(1-ci)/2))
  diffOfVariancesS3MinusS2Quants=quantile(diffOfVariancesS3MinusS2,probs=c((1-ci)/2,ci+(1-ci)/2))
  if(length(fourth>1)){
    diffOfVariancesS4MinusS1Quants=quantile(diffOfVariancesS4MinusS1,probs=c((1-ci)/2,ci+(1-ci)/2))
    diffOfVariancesS4MinusS2Quants=quantile(diffOfVariancesS4MinusS2,probs=c((1-ci)/2,ci+(1-ci)/2))
    diffOfVariancesS4MinusS3Quants=quantile(diffOfVariancesS4MinusS3,probs=c((1-ci)/2,ci+(1-ci)/2))
  }
  if(length(fifth>1)){
    diffOfVariancesS5MinusS1Quants=quantile(diffOfVariancesS5MinusS1,probs=c((1-ci)/2,ci+(1-ci)/2))
    diffOfVariancesS5MinusS2Quants=quantile(diffOfVariancesS5MinusS2,probs=c((1-ci)/2,ci+(1-ci)/2))
    diffOfVariancesS5MinusS3Quants=quantile(diffOfVariancesS5MinusS3,probs=c((1-ci)/2,ci+(1-ci)/2))
    diffOfVariancesS5MinusS4Quants=quantile(diffOfVariancesS5MinusS4,probs=c((1-ci)/2,ci+(1-ci)/2))
  }
  if(length(sixth>1)){
    diffOfVariancesS6MinusS1Quants=quantile(diffOfVariancesS6MinusS1,probs=c((1-ci)/2,ci+(1-ci)/2))
    diffOfVariancesS6MinusS2Quants=quantile(diffOfVariancesS6MinusS2,probs=c((1-ci)/2,ci+(1-ci)/2))
    diffOfVariancesS6MinusS3Quants=quantile(diffOfVariancesS6MinusS3,probs=c((1-ci)/2,ci+(1-ci)/2))
    diffOfVariancesS6MinusS4Quants=quantile(diffOfVariancesS6MinusS4,probs=c((1-ci)/2,ci+(1-ci)/2))
    diffOfVariancesS6MinusS5Quants=quantile(diffOfVariancesS6MinusS5,probs=c((1-ci)/2,ci+(1-ci)/2))
  }
  
  effectSize12Quants=quantile(effectSize12,probs=c((1-ci)/2,ci+(1-ci)/2))
  effectSize13Quants=quantile(effectSize13,probs=c((1-ci)/2,ci+(1-ci)/2))
  effectSize23Quants=quantile(effectSize23,probs=c((1-ci)/2,ci+(1-ci)/2))
  if(length(fourth>1)){
    effectSize14Quants=quantile(effectSize14,probs=c((1-ci)/2,ci+(1-ci)/2))
    effectSize24Quants=quantile(effectSize24,probs=c((1-ci)/2,ci+(1-ci)/2))
    effectSize34Quants=quantile(effectSize34,probs=c((1-ci)/2,ci+(1-ci)/2))
  }
  if(length(fifth>1)){
    effectSize15Quants=quantile(effectSize15,probs=c((1-ci)/2,ci+(1-ci)/2))
    effectSize25Quants=quantile(effectSize25,probs=c((1-ci)/2,ci+(1-ci)/2))
    effectSize35Quants=quantile(effectSize35,probs=c((1-ci)/2,ci+(1-ci)/2))
    effectSize45Quants=quantile(effectSize45,probs=c((1-ci)/2,ci+(1-ci)/2))
  }
  if(length(sixth>1)){
    effectSize16Quants=quantile(effectSize16,probs=c((1-ci)/2,ci+(1-ci)/2))
    effectSize26Quants=quantile(effectSize26,probs=c((1-ci)/2,ci+(1-ci)/2))
    effectSize36Quants=quantile(effectSize36,probs=c((1-ci)/2,ci+(1-ci)/2))
    effectSize46Quants=quantile(effectSize46,probs=c((1-ci)/2,ci+(1-ci)/2))
    effectSize56Quants=quantile(effectSize56,probs=c((1-ci)/2,ci+(1-ci)/2))
  }
  
  # Console Output
  if(length(first)>1 && length(second)>1 && length(third)>1 && length(fourth)>1 && length(fifth)>1 && length(sixth)>1) {
    Parameter<-c("mu1",
                 "mu2",
                 "mu3",
                 "mu4",
                 "mu5",
                 "mu6",
                 "sigma1",
                 "sigma2",
                 "sigma3",
                 "sigma4",
                 "sigma5",
                 "sigma6",
                 "mu2-mu1",
                 "mu3-mu1",
                 "mu4-mu1", 
                 "mu5-mu1", 
                 "mu6-mu1", 
                 "mu3-mu2", 
                 "mu4-mu2", 
                 "mu5-mu2", 
                 "mu6-mu2", 
                 "mu4-mu3", 
                 "mu5-mu3",
                 "mu6-mu3", 
                 "mu5-mu4", 
                 "mu6-mu4",
                 "mu6-mu5",
                 "sigma2-sigma1",
                 "sigma3-sigma1",
                 "sigma4-sigma1",
                 "sigma5-sigma1", 
                 "sigma6-sigma1", 
                 "sigma3-sigma2", 
                 "sigma4-sigma2", 
                 "sigma5-sigma2", 
                 "sigma6-sigma2", 
                 "sigma4-sigma3",
                 "sigma5-sigma3", 
                 "sigma6-sigma3", 
                 "sigma5-sigma4", 
                 "sigma6-sigma4", 
                 "sigma6-sigma5",
                 "delta12",
                 "delta13", 
                 "delta14", 
                 "delta15", 
                 "delta16", 
                 "delta23", 
                 "delta24", 
                 "delta25", 
                 "delta26", 
                 "delta34", 
                 "delta35", 
                 "delta36", 
                 "delta45",
                 "delta46", 
                 "delta56")
    LQ<-round(c(as.numeric(mu1Quants[1]),
                as.numeric(mu2Quants[1]),
                as.numeric(mu3Quants[1]), 
                as.numeric(mu4Quants[1]),
                as.numeric(mu5Quants[1]),
                as.numeric(mu6Quants[1]), 
                as.numeric(sigma1Quants[1]),
                as.numeric(sigma2Quants[1]),
                as.numeric(sigma3Quants[1]), 
                as.numeric(sigma4Quants[1]),
                as.numeric(sigma5Quants[1]), 
                as.numeric(sigma6Quants[1]),
                as.numeric(diffOfMeansMu2MinusMu1Quants[1]), 
                as.numeric(diffOfMeansMu3MinusMu1Quants[1]),
                as.numeric(diffOfMeansMu4MinusMu1Quants[1]),
                as.numeric(diffOfMeansMu5MinusMu1Quants[1]),
                as.numeric(diffOfMeansMu6MinusMu1Quants[1]), 
                as.numeric(diffOfMeansMu3MinusMu2Quants[1]), 
                as.numeric(diffOfMeansMu4MinusMu2Quants[1]), 
                as.numeric(diffOfMeansMu5MinusMu2Quants[1]), 
                as.numeric(diffOfMeansMu6MinusMu2Quants[1]), 
                as.numeric(diffOfMeansMu4MinusMu3Quants[1]), 
                as.numeric(diffOfMeansMu5MinusMu3Quants[1]), 
                as.numeric(diffOfMeansMu6MinusMu3Quants[1]),
                as.numeric(diffOfMeansMu5MinusMu4Quants[1]), 
                as.numeric(diffOfMeansMu6MinusMu4Quants[1]), 
                as.numeric(diffOfMeansMu6MinusMu5Quants[1]), 
                as.numeric(diffOfVariancesS2MinusS1Quants[1]),
                as.numeric(diffOfVariancesS3MinusS1Quants[1]),
                as.numeric(diffOfVariancesS4MinusS1Quants[1]),
                as.numeric(diffOfVariancesS5MinusS1Quants[1]), 
                as.numeric(diffOfVariancesS6MinusS1Quants[1]), 
                as.numeric(diffOfVariancesS3MinusS2Quants[1]), 
                as.numeric(diffOfVariancesS4MinusS2Quants[1]), 
                as.numeric(diffOfVariancesS5MinusS2Quants[1]), 
                as.numeric(diffOfVariancesS6MinusS2Quants[1]), 
                as.numeric(diffOfVariancesS4MinusS3Quants[1]),
                as.numeric(diffOfVariancesS5MinusS3Quants[1]), 
                as.numeric(diffOfVariancesS6MinusS3Quants[1]), 
                as.numeric(diffOfVariancesS5MinusS4Quants[1]), 
                as.numeric(diffOfVariancesS6MinusS4Quants[1]), 
                as.numeric(diffOfVariancesS6MinusS5Quants[1]),
                as.numeric(effectSize12Quants[1]), 
                as.numeric(effectSize13Quants[1]), 
                as.numeric(effectSize14Quants[1]), 
                as.numeric(effectSize15Quants[1]), 
                as.numeric(effectSize16Quants[1]), 
                as.numeric(effectSize23Quants[1]), 
                as.numeric(effectSize24Quants[1]), 
                as.numeric(effectSize25Quants[1]), 
                as.numeric(effectSize26Quants[1]), 
                as.numeric(effectSize34Quants[1]), 
                as.numeric(effectSize35Quants[1]), 
                as.numeric(effectSize36Quants[1]), 
                as.numeric(effectSize45Quants[1]), 
                as.numeric(effectSize46Quants[1]), 
                as.numeric(effectSize56Quants[1])),2)
    Mean<-round(c(mean(mu1), 
                  mean(mu2), 
                  mean(mu3), 
                  mean(mu4), 
                  mean(mu5), 
                  mean(mu6),
                  mean(sigma1Sq),
                  mean(sigma2Sq),
                  mean(sigma3Sq), 
                  mean(sigma4Sq), 
                  mean(sigma5Sq),
                  mean(sigma6Sq), 
                  mean(diffOfMeansMu2MinusMu1), 
                  mean(diffOfMeansMu3MinusMu1), 
                  mean(diffOfMeansMu4MinusMu1), 
                  mean(diffOfMeansMu5MinusMu1), 
                  mean(diffOfMeansMu6MinusMu1), 
                  mean(diffOfMeansMu3MinusMu2), 
                  mean(diffOfMeansMu4MinusMu2), 
                  mean(diffOfMeansMu5MinusMu2), 
                  mean(diffOfMeansMu6MinusMu2), 
                  mean(diffOfMeansMu4MinusMu3), 
                  mean(diffOfMeansMu5MinusMu3), 
                  mean(diffOfMeansMu6MinusMu3), 
                  mean(diffOfMeansMu5MinusMu4), 
                  mean(diffOfMeansMu6MinusMu4), 
                  mean(diffOfMeansMu6MinusMu5), 
                  mean(diffOfVariancesS2MinusS1),
                  mean(diffOfVariancesS3MinusS1),
                  mean(diffOfVariancesS4MinusS1),
                  mean(diffOfVariancesS5MinusS1), 
                  mean(diffOfVariancesS6MinusS1), 
                  mean(diffOfVariancesS3MinusS2), 
                  mean(diffOfVariancesS4MinusS2), 
                  mean(diffOfVariancesS5MinusS2), 
                  mean(diffOfVariancesS6MinusS2),
                  mean(diffOfVariancesS4MinusS3), 
                  mean(diffOfVariancesS5MinusS3), 
                  mean(diffOfVariancesS6MinusS3), 
                  mean(diffOfVariancesS5MinusS4), 
                  mean(diffOfVariancesS6MinusS4), 
                  mean(diffOfVariancesS6MinusS5), 
                  mean(effectSize12), 
                  mean(effectSize13), 
                  mean(effectSize14), 
                  mean(effectSize15), 
                  mean(effectSize16), 
                  mean(effectSize23),
                  mean(effectSize24), 
                  mean(effectSize25), 
                  mean(effectSize26), 
                  mean(effectSize34), 
                  mean(effectSize35), 
                  mean(effectSize36), 
                  mean(effectSize45), 
                  mean(effectSize46), 
                  mean(effectSize56)),2)
    UQ<-round(c(as.numeric(mu1Quants[2]), 
                as.numeric(mu2Quants[2]), 
                as.numeric(mu3Quants[2]), 
                as.numeric(mu4Quants[2]), 
                as.numeric(mu5Quants[2]), 
                as.numeric(mu6Quants[2]), 
                as.numeric(sigma1Quants[2]), 
                as.numeric(sigma2Quants[2]), 
                as.numeric(sigma3Quants[2]), 
                as.numeric(sigma4Quants[2]),
                as.numeric(sigma5Quants[2]),
                as.numeric(sigma6Quants[2]),
                as.numeric(diffOfMeansMu2MinusMu1Quants[2]),
                as.numeric(diffOfMeansMu3MinusMu1Quants[2]), 
                as.numeric(diffOfMeansMu4MinusMu1Quants[2]), 
                as.numeric(diffOfMeansMu5MinusMu1Quants[2]), 
                as.numeric(diffOfMeansMu6MinusMu1Quants[2]), 
                as.numeric(diffOfMeansMu3MinusMu2Quants[2]), 
                as.numeric(diffOfMeansMu4MinusMu2Quants[2]), 
                as.numeric(diffOfMeansMu5MinusMu2Quants[2]),
                as.numeric(diffOfMeansMu6MinusMu2Quants[2]), 
                as.numeric(diffOfMeansMu4MinusMu3Quants[2]), 
                as.numeric(diffOfMeansMu5MinusMu3Quants[2]), 
                as.numeric(diffOfMeansMu6MinusMu3Quants[2]), 
                as.numeric(diffOfMeansMu5MinusMu4Quants[2]), 
                as.numeric(diffOfMeansMu6MinusMu4Quants[2]), 
                as.numeric(diffOfMeansMu6MinusMu5Quants[2]), 
                as.numeric(diffOfVariancesS2MinusS1Quants[2]), 
                as.numeric(diffOfVariancesS3MinusS1Quants[2]), 
                as.numeric(diffOfVariancesS4MinusS1Quants[2]), 
                as.numeric(diffOfVariancesS5MinusS1Quants[2]), 
                as.numeric(diffOfVariancesS6MinusS1Quants[2]), 
                as.numeric(diffOfVariancesS3MinusS2Quants[2]), 
                as.numeric(diffOfVariancesS4MinusS2Quants[2]), 
                as.numeric(diffOfVariancesS5MinusS2Quants[2]), 
                as.numeric(diffOfVariancesS6MinusS2Quants[2]), 
                as.numeric(diffOfVariancesS4MinusS3Quants[2]), 
                as.numeric(diffOfVariancesS5MinusS3Quants[2]), 
                as.numeric(diffOfVariancesS6MinusS3Quants[2]), 
                as.numeric(diffOfVariancesS5MinusS4Quants[2]),
                as.numeric(diffOfVariancesS6MinusS4Quants[2]),
                as.numeric(diffOfVariancesS6MinusS5Quants[2]),
                as.numeric(effectSize12Quants[2]),
                as.numeric(effectSize13Quants[2]),
                as.numeric(effectSize14Quants[2]), 
                as.numeric(effectSize15Quants[2]), 
                as.numeric(effectSize16Quants[2]), 
                as.numeric(effectSize23Quants[2]), 
                as.numeric(effectSize24Quants[2]), 
                as.numeric(effectSize25Quants[2]), 
                as.numeric(effectSize26Quants[2]), 
                as.numeric(effectSize34Quants[2]), 
                as.numeric(effectSize35Quants[2]), 
                as.numeric(effectSize36Quants[2]), 
                as.numeric(effectSize45Quants[2]), 
                as.numeric(effectSize46Quants[2]),
                as.numeric(effectSize56Quants[2])),2)
    Std.Err<-round(c(sd(mu1), 
                     sd(mu2), 
                     sd(mu3), 
                     sd(mu4),
                     sd(mu5), 
                     sd(mu6),
                     sd(sigma1Sq), 
                     sd(sigma2Sq),
                     sd(sigma3Sq), 
                     sd(sigma4Sq), 
                     sd(sigma5Sq), 
                     sd(sigma6Sq), 
                     sd(diffOfMeansMu2MinusMu1), 
                     sd(diffOfMeansMu3MinusMu1), 
                     sd(diffOfMeansMu4MinusMu1), 
                     sd(diffOfMeansMu5MinusMu1),
                     sd(diffOfMeansMu6MinusMu1), 
                     sd(diffOfMeansMu3MinusMu2), 
                     sd(diffOfMeansMu4MinusMu2), 
                     sd(diffOfMeansMu5MinusMu2), 
                     sd(diffOfMeansMu6MinusMu2), 
                     sd(diffOfMeansMu4MinusMu3), 
                     sd(diffOfMeansMu5MinusMu3), 
                     sd(diffOfMeansMu6MinusMu3), 
                     sd(diffOfMeansMu5MinusMu4), 
                     sd(diffOfMeansMu6MinusMu4), 
                     sd(diffOfMeansMu6MinusMu5), 
                     sd(diffOfVariancesS2MinusS1), 
                     sd(diffOfVariancesS3MinusS1), 
                     sd(diffOfVariancesS4MinusS1), 
                     sd(diffOfVariancesS5MinusS1),
                     sd(diffOfVariancesS6MinusS1), 
                     sd(diffOfVariancesS3MinusS2), 
                     sd(diffOfVariancesS4MinusS2), 
                     sd(diffOfVariancesS5MinusS2), 
                     sd(diffOfVariancesS6MinusS2), 
                     sd(diffOfVariancesS4MinusS3), 
                     sd(diffOfVariancesS5MinusS3),
                     sd(diffOfVariancesS6MinusS3), 
                     sd(diffOfVariancesS5MinusS4), 
                     sd(diffOfVariancesS6MinusS4), 
                     sd(diffOfVariancesS6MinusS5), 
                     sd(effectSize12), 
                     sd(effectSize13), 
                     sd(effectSize14), 
                     sd(effectSize15), 
                     sd(effectSize16), 
                     sd(effectSize23), 
                     sd(effectSize24), 
                     sd(effectSize25), 
                     sd(effectSize26), 
                     sd(effectSize34), 
                     sd(effectSize35), 
                     sd(effectSize36), 
                     sd(effectSize45), 
                     sd(effectSize46), 
                     sd(effectSize56)),2)
    #dataframe output
    out <- cbind(Parameter,LQ,Mean,UQ,Std.Err)
    
    # posterior draw dataframe
    df=data.frame(mu1, mu1SecC, mu1ThdC, mu1FrtC, 
               mu2, mu2SecC, mu2ThdC, mu2FrtC, 
               mu3, mu3SecC, mu3ThdC, mu3FrtC,
               mu4, mu4SecC, mu4ThdC, mu4FrtC,
               mu5, mu5SecC, mu5ThdC, mu5FrtC,
               mu6, mu6SecC, mu6ThdC, mu6FrtC,
               sigma1Sq, sigma1SqSecC, sigma1SqThdC, sigma1SqFrtC, 
               sigma2Sq, sigma2SqSecC, sigma2SqThdC, sigma2SqFrtC, 
               sigma3Sq, sigma3SqSecC, sigma3SqThdC, sigma3SqFrtC,
               sigma4Sq, sigma4SqSecC, sigma4SqThdC, sigma4SqFrtC,
               sigma5Sq, sigma5SqSecC, sigma5SqThdC, sigma5SqFrtC,
               sigma6Sq, sigma6SqSecC, sigma6SqThdC, sigma6SqFrtC,
               diffOfMeansMu2MinusMu1, diffOfMeansMu2MinusMu1SecC, diffOfMeansMu2MinusMu1ThdC, diffOfMeansMu2MinusMu1FrtC,
               diffOfMeansMu3MinusMu1, diffOfMeansMu3MinusMu1SecC, diffOfMeansMu3MinusMu1ThdC, diffOfMeansMu3MinusMu1FrtC,
               diffOfMeansMu4MinusMu1, diffOfMeansMu4MinusMu1SecC, diffOfMeansMu4MinusMu1ThdC, diffOfMeansMu4MinusMu1FrtC,
               diffOfMeansMu5MinusMu1, diffOfMeansMu5MinusMu1SecC, diffOfMeansMu5MinusMu1ThdC, diffOfMeansMu5MinusMu1FrtC,
               diffOfMeansMu6MinusMu1, diffOfMeansMu6MinusMu1SecC, diffOfMeansMu6MinusMu1ThdC, diffOfMeansMu6MinusMu1FrtC,
               
               diffOfMeansMu3MinusMu2, diffOfMeansMu3MinusMu2SecC, diffOfMeansMu3MinusMu2ThdC, diffOfMeansMu3MinusMu2FrtC,
               diffOfMeansMu4MinusMu2, diffOfMeansMu4MinusMu2SecC, diffOfMeansMu4MinusMu2ThdC, diffOfMeansMu4MinusMu2FrtC,
               diffOfMeansMu5MinusMu2, diffOfMeansMu5MinusMu2SecC, diffOfMeansMu5MinusMu2ThdC, diffOfMeansMu5MinusMu2FrtC,
               diffOfMeansMu6MinusMu2, diffOfMeansMu6MinusMu2SecC, diffOfMeansMu6MinusMu2ThdC, diffOfMeansMu6MinusMu2FrtC,
               
               diffOfMeansMu4MinusMu3, diffOfMeansMu4MinusMu3SecC, diffOfMeansMu4MinusMu3ThdC, diffOfMeansMu4MinusMu3FrtC,
               diffOfMeansMu5MinusMu3, diffOfMeansMu5MinusMu3SecC, diffOfMeansMu5MinusMu3ThdC, diffOfMeansMu5MinusMu3FrtC,
               diffOfMeansMu6MinusMu3, diffOfMeansMu6MinusMu3SecC, diffOfMeansMu6MinusMu3ThdC, diffOfMeansMu6MinusMu3FrtC,
               
               diffOfMeansMu5MinusMu4, diffOfMeansMu5MinusMu4SecC, diffOfMeansMu5MinusMu4ThdC, diffOfMeansMu5MinusMu4FrtC,
               diffOfMeansMu6MinusMu4, diffOfMeansMu6MinusMu4SecC, diffOfMeansMu6MinusMu4ThdC, diffOfMeansMu6MinusMu4FrtC,
               
               diffOfMeansMu6MinusMu5, diffOfMeansMu6MinusMu5SecC, diffOfMeansMu6MinusMu5ThdC, diffOfMeansMu6MinusMu5FrtC,
               
               
               diffOfVariancesS2MinusS1, diffOfVariancesS2MinusS1SecC, diffOfVariancesS2MinusS1ThdC, diffOfVariancesS2MinusS1FrtC,
               diffOfVariancesS3MinusS1, diffOfVariancesS3MinusS1SecC, diffOfVariancesS3MinusS1ThdC, diffOfVariancesS3MinusS1FrtC,
               diffOfVariancesS4MinusS1, diffOfVariancesS4MinusS1SecC, diffOfVariancesS4MinusS1ThdC, diffOfVariancesS4MinusS1FrtC,
               diffOfVariancesS5MinusS1, diffOfVariancesS5MinusS1SecC, diffOfVariancesS5MinusS1ThdC, diffOfVariancesS5MinusS1FrtC,
               diffOfVariancesS6MinusS1, diffOfVariancesS6MinusS1SecC, diffOfVariancesS6MinusS1ThdC, diffOfVariancesS6MinusS1FrtC,
               
               diffOfVariancesS3MinusS2, diffOfVariancesS3MinusS2SecC, diffOfVariancesS3MinusS2ThdC, diffOfVariancesS3MinusS2FrtC,
               diffOfVariancesS4MinusS2, diffOfVariancesS4MinusS2SecC, diffOfVariancesS4MinusS2ThdC, diffOfVariancesS4MinusS2FrtC,
               diffOfVariancesS5MinusS2, diffOfVariancesS5MinusS2SecC, diffOfVariancesS5MinusS2ThdC, diffOfVariancesS5MinusS2FrtC,
               diffOfVariancesS6MinusS2, diffOfVariancesS6MinusS2SecC, diffOfVariancesS6MinusS2ThdC, diffOfVariancesS6MinusS2FrtC,
               
               diffOfVariancesS4MinusS3, diffOfVariancesS4MinusS3SecC, diffOfVariancesS4MinusS3ThdC, diffOfVariancesS4MinusS3FrtC,
               diffOfVariancesS5MinusS3, diffOfVariancesS5MinusS3SecC, diffOfVariancesS5MinusS3ThdC, diffOfVariancesS5MinusS3FrtC,
               diffOfVariancesS6MinusS3, diffOfVariancesS6MinusS3SecC, diffOfVariancesS6MinusS3ThdC, diffOfVariancesS6MinusS3FrtC,
               
               diffOfVariancesS5MinusS4, diffOfVariancesS5MinusS4SecC, diffOfVariancesS5MinusS4ThdC, diffOfVariancesS5MinusS4FrtC,
               diffOfVariancesS6MinusS4, diffOfVariancesS6MinusS4SecC, diffOfVariancesS6MinusS4ThdC, diffOfVariancesS6MinusS4FrtC,
               
               diffOfVariancesS6MinusS5, diffOfVariancesS6MinusS5SecC, diffOfVariancesS6MinusS5ThdC, diffOfVariancesS6MinusS5FrtC,
               
               effectSize12, effectSize12SecC, effectSize12ThdC, effectSize12FrtC,
               effectSize13, effectSize13SecC, effectSize13ThdC, effectSize13FrtC,
               effectSize14, effectSize14SecC, effectSize14ThdC, effectSize14FrtC,
               effectSize15, effectSize15SecC, effectSize15ThdC, effectSize15FrtC,
               effectSize16, effectSize16SecC, effectSize16ThdC, effectSize16FrtC,
               
               effectSize23, effectSize23SecC, effectSize23ThdC, effectSize23FrtC,
               effectSize24, effectSize24SecC, effectSize24ThdC, effectSize24FrtC,
               effectSize25, effectSize25SecC, effectSize25ThdC, effectSize25FrtC,
               effectSize26, effectSize26SecC, effectSize26ThdC, effectSize26FrtC,
               
               effectSize34, effectSize34SecC, effectSize34ThdC, effectSize34FrtC,
               effectSize35, effectSize35SecC, effectSize35ThdC, effectSize35FrtC,
               effectSize36, effectSize36SecC, effectSize36ThdC, effectSize36FrtC,
               
               effectSize45, effectSize45SecC, effectSize45ThdC, effectSize45FrtC,
               effectSize46, effectSize46SecC, effectSize46ThdC, effectSize46FrtC,
               
               effectSize56, effectSize56SecC, effectSize56ThdC, effectSize56FrtC)
    
    #pretty table
    print(knitr::kable(out))
    
    # return dataframe
    return(df)

  } else if(length(first)>1 && length(second)>1 && length(third)>1 && length(fourth)>1 && length(fifth)>1 && is.null(sixth)) {
    Parameter<-c("mu1",
                 "mu2",
                 "mu3",
                 "mu4",
                 "mu5",
                 "sigma1",
                 "sigma2",
                 "sigma3",
                 "sigma4",
                 "sigma5",
                 "mu2-mu1",
                 "mu3-mu1",
                 "mu4-mu1", 
                 "mu5-mu1", 
                 "mu3-mu2", 
                 "mu4-mu2", 
                 "mu5-mu2", 
                 "mu4-mu3", 
                 "mu5-mu3", 
                 "mu5-mu4",
                 "sigma2-sigma1",
                 "sigma3-sigma1",
                 "sigma4-sigma1",
                 "sigma5-sigma1",
                 "sigma3-sigma2",
                 "sigma4-sigma2",
                 "sigma5-sigma2",
                 "sigma4-sigma3", 
                 "sigma5-sigma3", 
                 "sigma5-sigma4",
                 "delta12",
                 "delta13", 
                 "delta14",
                 "delta15",
                 "delta23", 
                 "delta24", 
                 "delta25", 
                 "delta34",
                 "delta35", 
                 "delta45")
    LQ<-round(c(as.numeric(mu1Quants[1]), 
                as.numeric(mu2Quants[1]), 
                as.numeric(mu3Quants[1]), 
                as.numeric(mu4Quants[1]), 
                as.numeric(mu5Quants[1]), 
                as.numeric(sigma1Quants[1]), 
                as.numeric(sigma2Quants[1]), 
                as.numeric(sigma3Quants[1]), 
                as.numeric(sigma4Quants[1]), 
                as.numeric(sigma5Quants[1]),
                as.numeric(diffOfMeansMu2MinusMu1Quants[1]), 
                as.numeric(diffOfMeansMu3MinusMu1Quants[1]), 
                as.numeric(diffOfMeansMu4MinusMu1Quants[1]), 
                as.numeric(diffOfMeansMu5MinusMu1Quants[1]), 
                as.numeric(diffOfMeansMu3MinusMu2Quants[1]), 
                as.numeric(diffOfMeansMu4MinusMu2Quants[1]), 
                as.numeric(diffOfMeansMu5MinusMu2Quants[1]), 
                as.numeric(diffOfMeansMu4MinusMu3Quants[1]), 
                as.numeric(diffOfMeansMu5MinusMu3Quants[1]), 
                as.numeric(diffOfMeansMu5MinusMu4Quants[1]), 
                as.numeric(diffOfVariancesS2MinusS1Quants[1]), 
                as.numeric(diffOfVariancesS3MinusS1Quants[1]), 
                as.numeric(diffOfVariancesS4MinusS1Quants[1]), 
                as.numeric(diffOfVariancesS5MinusS1Quants[1]), 
                as.numeric(diffOfVariancesS3MinusS2Quants[1]), 
                as.numeric(diffOfVariancesS4MinusS2Quants[1]), 
                as.numeric(diffOfVariancesS5MinusS2Quants[1]), 
                as.numeric(diffOfVariancesS4MinusS3Quants[1]), 
                as.numeric(diffOfVariancesS5MinusS3Quants[1]),  
                as.numeric(diffOfVariancesS5MinusS4Quants[1]),
                as.numeric(effectSize12Quants[1]), 
                as.numeric(effectSize13Quants[1]), 
                as.numeric(effectSize14Quants[1]), 
                as.numeric(effectSize15Quants[1]), 
                as.numeric(effectSize23Quants[1]), 
                as.numeric(effectSize24Quants[1]), 
                as.numeric(effectSize25Quants[1]), 
                as.numeric(effectSize34Quants[1]), 
                as.numeric(effectSize35Quants[1]), 
                as.numeric(effectSize45Quants[1])),2)
    Mean<-round(c(mean(mu1), 
                  mean(mu2), 
                  mean(mu3), 
                  mean(mu4), 
                  mean(mu5), 
                  mean(sigma1Sq), 
                  mean(sigma2Sq), 
                  mean(sigma3Sq), 
                  mean(sigma4Sq), 
                  mean(sigma5Sq), 
                  mean(diffOfMeansMu2MinusMu1), 
                  mean(diffOfMeansMu3MinusMu1), 
                  mean(diffOfMeansMu4MinusMu1), 
                  mean(diffOfMeansMu5MinusMu1), 
                  mean(diffOfMeansMu3MinusMu2), 
                  mean(diffOfMeansMu4MinusMu2), 
                  mean(diffOfMeansMu5MinusMu2), 
                  mean(diffOfMeansMu4MinusMu3), 
                  mean(diffOfMeansMu5MinusMu3), 
                  mean(diffOfMeansMu5MinusMu4), 
                  mean(diffOfVariancesS2MinusS1), 
                  mean(diffOfVariancesS3MinusS1), 
                  mean(diffOfVariancesS4MinusS1), 
                  mean(diffOfVariancesS5MinusS1), 
                  mean(diffOfVariancesS3MinusS2), 
                  mean(diffOfVariancesS4MinusS2), 
                  mean(diffOfVariancesS5MinusS2), 
                  mean(diffOfVariancesS4MinusS3), 
                  mean(diffOfVariancesS5MinusS3), 
                  mean(diffOfVariancesS5MinusS4), 
                  mean(effectSize12), 
                  mean(effectSize13), 
                  mean(effectSize14), 
                  mean(effectSize15), 
                  mean(effectSize23), 
                  mean(effectSize24), 
                  mean(effectSize25), 
                  mean(effectSize34), 
                  mean(effectSize35), 
                  mean(effectSize45)),2)
    UQ<-round(c(as.numeric(mu1Quants[2]), 
                as.numeric(mu2Quants[2]), 
                as.numeric(mu3Quants[2]), 
                as.numeric(mu4Quants[2]), 
                as.numeric(mu5Quants[2]),
                as.numeric(sigma1Quants[2]), 
                as.numeric(sigma2Quants[2]), 
                as.numeric(sigma3Quants[2]), 
                as.numeric(sigma4Quants[2]), 
                as.numeric(sigma5Quants[2]),
                as.numeric(diffOfMeansMu2MinusMu1Quants[2]), 
                as.numeric(diffOfMeansMu3MinusMu1Quants[2]), 
                as.numeric(diffOfMeansMu4MinusMu1Quants[2]), 
                as.numeric(diffOfMeansMu5MinusMu1Quants[2]), 
                as.numeric(diffOfMeansMu3MinusMu2Quants[2]), 
                as.numeric(diffOfMeansMu4MinusMu2Quants[2]), 
                as.numeric(diffOfMeansMu5MinusMu2Quants[2]), 
                as.numeric(diffOfMeansMu4MinusMu3Quants[2]), 
                as.numeric(diffOfMeansMu5MinusMu3Quants[2]), 
                as.numeric(diffOfMeansMu5MinusMu4Quants[2]), 
                as.numeric(diffOfVariancesS2MinusS1Quants[2]), 
                as.numeric(diffOfVariancesS3MinusS1Quants[2]), 
                as.numeric(diffOfVariancesS4MinusS1Quants[2]), 
                as.numeric(diffOfVariancesS5MinusS1Quants[2]), 
                as.numeric(diffOfVariancesS3MinusS2Quants[2]), 
                as.numeric(diffOfVariancesS4MinusS2Quants[2]), 
                as.numeric(diffOfVariancesS5MinusS2Quants[2]), 
                as.numeric(diffOfVariancesS4MinusS3Quants[2]), 
                as.numeric(diffOfVariancesS5MinusS3Quants[2]),  
                as.numeric(diffOfVariancesS5MinusS4Quants[2]),
                as.numeric(effectSize12Quants[2]), 
                as.numeric(effectSize13Quants[2]), 
                as.numeric(effectSize14Quants[2]), 
                as.numeric(effectSize15Quants[2]), 
                as.numeric(effectSize23Quants[2]), 
                as.numeric(effectSize24Quants[2]), 
                as.numeric(effectSize25Quants[2]), 
                as.numeric(effectSize34Quants[2]), 
                as.numeric(effectSize35Quants[2]), 
                as.numeric(effectSize45Quants[2])),2)
    Std.Err<-round(c(sd(mu1), 
                     sd(mu2), 
                     sd(mu3), 
                     sd(mu4), 
                     sd(mu5),
                     sd(sigma1Sq), 
                     sd(sigma2Sq), 
                     sd(sigma3Sq), 
                     sd(sigma4Sq), 
                     sd(sigma5Sq),
                     sd(diffOfMeansMu2MinusMu1),
                     sd(diffOfMeansMu3MinusMu1),
                     sd(diffOfMeansMu4MinusMu1), 
                     sd(diffOfMeansMu5MinusMu1), 
                     sd(diffOfMeansMu3MinusMu2), 
                     sd(diffOfMeansMu4MinusMu2), 
                     sd(diffOfMeansMu5MinusMu2), 
                     sd(diffOfMeansMu4MinusMu3), 
                     sd(diffOfMeansMu5MinusMu3), 
                     sd(diffOfMeansMu5MinusMu4), 
                     sd(diffOfVariancesS2MinusS1),
                     sd(diffOfVariancesS3MinusS1), 
                     sd(diffOfVariancesS4MinusS1), 
                     sd(diffOfVariancesS5MinusS1), 
                     sd(diffOfVariancesS3MinusS2), 
                     sd(diffOfVariancesS4MinusS2), 
                     sd(diffOfVariancesS5MinusS2), 
                     sd(diffOfVariancesS4MinusS3), 
                     sd(diffOfVariancesS5MinusS3), 
                     sd(diffOfVariancesS5MinusS4), 
                     sd(effectSize12), 
                     sd(effectSize13), 
                     sd(effectSize14), 
                     sd(effectSize15), 
                     sd(effectSize23), 
                     sd(effectSize24), 
                     sd(effectSize25), 
                     sd(effectSize34), 
                     sd(effectSize35), 
                     sd(effectSize45)),2)
    #dataframe output
    out <- cbind(Parameter,LQ,Mean,UQ,Std.Err)
    
    # posterior draw dataframe
    df=data.frame(mu1, mu1SecC, mu1ThdC, mu1FrtC, 
                  mu2, mu2SecC, mu2ThdC, mu2FrtC, 
                  mu3, mu3SecC, mu3ThdC, mu3FrtC,
                  mu4, mu4SecC, mu4ThdC, mu4FrtC,
                  mu5, mu5SecC, mu5ThdC, mu5FrtC,
                  sigma1Sq, sigma1SqSecC, sigma1SqThdC, sigma1SqFrtC, 
                  sigma2Sq, sigma2SqSecC, sigma2SqThdC, sigma2SqFrtC, 
                  sigma3Sq, sigma3SqSecC, sigma3SqThdC, sigma3SqFrtC,
                  sigma4Sq, sigma4SqSecC, sigma4SqThdC, sigma4SqFrtC,
                  sigma5Sq, sigma5SqSecC, sigma5SqThdC, sigma5SqFrtC,
                  diffOfMeansMu2MinusMu1, diffOfMeansMu2MinusMu1SecC, diffOfMeansMu2MinusMu1ThdC, diffOfMeansMu2MinusMu1FrtC,
                  diffOfMeansMu3MinusMu1, diffOfMeansMu3MinusMu1SecC, diffOfMeansMu3MinusMu1ThdC, diffOfMeansMu3MinusMu1FrtC,
                  diffOfMeansMu4MinusMu1, diffOfMeansMu4MinusMu1SecC, diffOfMeansMu4MinusMu1ThdC, diffOfMeansMu4MinusMu1FrtC,
                  diffOfMeansMu5MinusMu1, diffOfMeansMu5MinusMu1SecC, diffOfMeansMu5MinusMu1ThdC, diffOfMeansMu5MinusMu1FrtC,

                  diffOfMeansMu3MinusMu2, diffOfMeansMu3MinusMu2SecC, diffOfMeansMu3MinusMu2ThdC, diffOfMeansMu3MinusMu2FrtC,
                  diffOfMeansMu4MinusMu2, diffOfMeansMu4MinusMu2SecC, diffOfMeansMu4MinusMu2ThdC, diffOfMeansMu4MinusMu2FrtC,
                  diffOfMeansMu5MinusMu2, diffOfMeansMu5MinusMu2SecC, diffOfMeansMu5MinusMu2ThdC, diffOfMeansMu5MinusMu2FrtC,

                  diffOfMeansMu4MinusMu3, diffOfMeansMu4MinusMu3SecC, diffOfMeansMu4MinusMu3ThdC, diffOfMeansMu4MinusMu3FrtC,
                  diffOfMeansMu5MinusMu3, diffOfMeansMu5MinusMu3SecC, diffOfMeansMu5MinusMu3ThdC, diffOfMeansMu5MinusMu3FrtC,

                  diffOfMeansMu5MinusMu4, diffOfMeansMu5MinusMu4SecC, diffOfMeansMu5MinusMu4ThdC, diffOfMeansMu5MinusMu4FrtC,
                  
                  
                  diffOfVariancesS2MinusS1, diffOfVariancesS2MinusS1SecC, diffOfVariancesS2MinusS1ThdC, diffOfVariancesS2MinusS1FrtC,
                  diffOfVariancesS3MinusS1, diffOfVariancesS3MinusS1SecC, diffOfVariancesS3MinusS1ThdC, diffOfVariancesS3MinusS1FrtC,
                  diffOfVariancesS4MinusS1, diffOfVariancesS4MinusS1SecC, diffOfVariancesS4MinusS1ThdC, diffOfVariancesS4MinusS1FrtC,
                  diffOfVariancesS5MinusS1, diffOfVariancesS5MinusS1SecC, diffOfVariancesS5MinusS1ThdC, diffOfVariancesS5MinusS1FrtC,

                  diffOfVariancesS3MinusS2, diffOfVariancesS3MinusS2SecC, diffOfVariancesS3MinusS2ThdC, diffOfVariancesS3MinusS2FrtC,
                  diffOfVariancesS4MinusS2, diffOfVariancesS4MinusS2SecC, diffOfVariancesS4MinusS2ThdC, diffOfVariancesS4MinusS2FrtC,
                  diffOfVariancesS5MinusS2, diffOfVariancesS5MinusS2SecC, diffOfVariancesS5MinusS2ThdC, diffOfVariancesS5MinusS2FrtC,

                  diffOfVariancesS4MinusS3, diffOfVariancesS4MinusS3SecC, diffOfVariancesS4MinusS3ThdC, diffOfVariancesS4MinusS3FrtC,
                  diffOfVariancesS5MinusS3, diffOfVariancesS5MinusS3SecC, diffOfVariancesS5MinusS3ThdC, diffOfVariancesS5MinusS3FrtC,
                  
                  diffOfVariancesS5MinusS4, diffOfVariancesS5MinusS4SecC, diffOfVariancesS5MinusS4ThdC, diffOfVariancesS5MinusS4FrtC,


                  effectSize12, effectSize12SecC, effectSize12ThdC, effectSize12FrtC,
                  effectSize13, effectSize13SecC, effectSize13ThdC, effectSize13FrtC,
                  effectSize14, effectSize14SecC, effectSize14ThdC, effectSize14FrtC,
                  effectSize15, effectSize15SecC, effectSize15ThdC, effectSize15FrtC,

                  effectSize23, effectSize23SecC, effectSize23ThdC, effectSize23FrtC,
                  effectSize24, effectSize24SecC, effectSize24ThdC, effectSize24FrtC,
                  effectSize25, effectSize25SecC, effectSize25ThdC, effectSize25FrtC,

                  effectSize34, effectSize34SecC, effectSize34ThdC, effectSize34FrtC,
                  effectSize35, effectSize35SecC, effectSize35ThdC, effectSize35FrtC,

                  effectSize45, effectSize45SecC, effectSize45ThdC, effectSize45FrtC)
    
    
    #pretty table
    print(knitr::kable(out))
    
    # return posterior draw dataframe
    return(df)
    
  } else if(length(first)>1 && length(second)>1 && length(third)>1 && length(fourth)>1 && is.null(fifth) && is.null(sixth)) {
    Parameter<-c("mu1",
                 "mu2",
                 "mu3",
                 "mu4",
                 "sigma1",
                 "sigma2",
                 "sigma3",
                 "sigma4",
                 "mu2-mu1",
                 "mu3-mu1",
                 "mu4-mu1", 
                 "mu3-mu2", 
                 "mu4-mu2", 
                 "mu4-mu3",
                 "sigma2-sigma1",
                 "sigma3-sigma1",
                 "sigma4-sigma1",
                 "sigma3-sigma2",
                 "sigma4-sigma2",
                 "sigma4-sigma3",
                 "delta12",
                 "delta13", 
                 "delta14",
                 "delta23", 
                 "delta24", 
                 "delta34")
    LQ<-round(c(as.numeric(mu1Quants[1]), 
                as.numeric(mu2Quants[1]), 
                as.numeric(mu3Quants[1]), 
                as.numeric(mu4Quants[1]), 
                as.numeric(sigma1Quants[1]), 
                as.numeric(sigma2Quants[1]), 
                as.numeric(sigma3Quants[1]), 
                as.numeric(sigma4Quants[1]),
                as.numeric(diffOfMeansMu2MinusMu1Quants[1]), 
                as.numeric(diffOfMeansMu3MinusMu1Quants[1]), 
                as.numeric(diffOfMeansMu4MinusMu1Quants[1]), 
                as.numeric(diffOfMeansMu3MinusMu2Quants[1]), 
                as.numeric(diffOfMeansMu4MinusMu2Quants[1]), 
                as.numeric(diffOfMeansMu4MinusMu3Quants[1]), 
                as.numeric(diffOfVariancesS2MinusS1Quants[1]), 
                as.numeric(diffOfVariancesS3MinusS1Quants[1]), 
                as.numeric(diffOfVariancesS4MinusS1Quants[1]), 
                as.numeric(diffOfVariancesS3MinusS2Quants[1]), 
                as.numeric(diffOfVariancesS4MinusS2Quants[1]), 
                as.numeric(diffOfVariancesS4MinusS3Quants[1]),
                as.numeric(effectSize12Quants[1]), 
                as.numeric(effectSize13Quants[1]), 
                as.numeric(effectSize14Quants[1]), 
                as.numeric(effectSize23Quants[1]), 
                as.numeric(effectSize24Quants[1]), 
                as.numeric(effectSize34Quants[1])),2)
    Mean<-round(c(mean(mu1), 
                  mean(mu2), 
                  mean(mu3), 
                  mean(mu4), 
                  mean(sigma1Sq), 
                  mean(sigma2Sq), 
                  mean(sigma3Sq), 
                  mean(sigma4Sq), 
                  mean(diffOfMeansMu2MinusMu1), 
                  mean(diffOfMeansMu3MinusMu1), 
                  mean(diffOfMeansMu4MinusMu1), 
                  mean(diffOfMeansMu3MinusMu2), 
                  mean(diffOfMeansMu4MinusMu2), 
                  mean(diffOfMeansMu4MinusMu3), 
                  mean(diffOfVariancesS2MinusS1), 
                  mean(diffOfVariancesS3MinusS1), 
                  mean(diffOfVariancesS4MinusS1), 
                  mean(diffOfVariancesS3MinusS2), 
                  mean(diffOfVariancesS4MinusS2), 
                  mean(diffOfVariancesS4MinusS3), 
                  mean(effectSize12), 
                  mean(effectSize13), 
                  mean(effectSize14), 
                  mean(effectSize23), 
                  mean(effectSize24), 
                  mean(effectSize34)),2)
    UQ<-round(c(as.numeric(mu1Quants[2]), 
                as.numeric(mu2Quants[2]), 
                as.numeric(mu3Quants[2]), 
                as.numeric(mu4Quants[2]),
                as.numeric(sigma1Quants[2]), 
                as.numeric(sigma2Quants[2]), 
                as.numeric(sigma3Quants[2]), 
                as.numeric(sigma4Quants[2]),
                as.numeric(diffOfMeansMu2MinusMu1Quants[2]), 
                as.numeric(diffOfMeansMu3MinusMu1Quants[2]), 
                as.numeric(diffOfMeansMu4MinusMu1Quants[2]), 
                as.numeric(diffOfMeansMu3MinusMu2Quants[2]), 
                as.numeric(diffOfMeansMu4MinusMu2Quants[2]), 
                as.numeric(diffOfMeansMu4MinusMu3Quants[2]), 
                as.numeric(diffOfVariancesS2MinusS1Quants[2]), 
                as.numeric(diffOfVariancesS3MinusS1Quants[2]), 
                as.numeric(diffOfVariancesS4MinusS1Quants[2]), 
                as.numeric(diffOfVariancesS3MinusS2Quants[2]), 
                as.numeric(diffOfVariancesS4MinusS2Quants[2]), 
                as.numeric(diffOfVariancesS4MinusS3Quants[2]),
                as.numeric(effectSize12Quants[2]), 
                as.numeric(effectSize13Quants[2]), 
                as.numeric(effectSize14Quants[2]), 
                as.numeric(effectSize23Quants[2]), 
                as.numeric(effectSize24Quants[2]), 
                as.numeric(effectSize34Quants[2])),2)
    Std.Err<-round(c(sd(mu1), 
                     sd(mu2), 
                     sd(mu3), 
                     sd(mu4),
                     sd(sigma1Sq), 
                     sd(sigma2Sq), 
                     sd(sigma3Sq), 
                     sd(sigma4Sq),
                     sd(diffOfMeansMu2MinusMu1),
                     sd(diffOfMeansMu3MinusMu1),
                     sd(diffOfMeansMu4MinusMu1), 
                     sd(diffOfMeansMu3MinusMu2), 
                     sd(diffOfMeansMu4MinusMu2), 
                     sd(diffOfMeansMu4MinusMu3), 
                     sd(diffOfVariancesS2MinusS1),
                     sd(diffOfVariancesS3MinusS1), 
                     sd(diffOfVariancesS4MinusS1), 
                     sd(diffOfVariancesS3MinusS2), 
                     sd(diffOfVariancesS4MinusS2), 
                     sd(diffOfVariancesS4MinusS3), 
                     sd(effectSize12), 
                     sd(effectSize13), 
                     sd(effectSize14), 
                     sd(effectSize23), 
                     sd(effectSize24), 
                     sd(effectSize34)),2)
    #dataframe output
    out <- cbind(Parameter,LQ,Mean,UQ,Std.Err)
    
    #pretty table
    print(knitr::kable(out))
    
    df=data.frame(mu1, mu1SecC, mu1ThdC, mu1FrtC, 
                  mu2, mu2SecC, mu2ThdC, mu2FrtC, 
                  mu3, mu3SecC, mu3ThdC, mu3FrtC,
                  mu4, mu4SecC, mu4ThdC, mu4FrtC,
                  sigma1Sq, sigma1SqSecC, sigma1SqThdC, sigma1SqFrtC, 
                  sigma2Sq, sigma2SqSecC, sigma2SqThdC, sigma2SqFrtC, 
                  sigma3Sq, sigma3SqSecC, sigma3SqThdC, sigma3SqFrtC,
                  sigma4Sq, sigma4SqSecC, sigma4SqThdC, sigma4SqFrtC,
                  diffOfMeansMu2MinusMu1, diffOfMeansMu2MinusMu1SecC, diffOfMeansMu2MinusMu1ThdC, diffOfMeansMu2MinusMu1FrtC,
                  diffOfMeansMu3MinusMu1, diffOfMeansMu3MinusMu1SecC, diffOfMeansMu3MinusMu1ThdC, diffOfMeansMu3MinusMu1FrtC,
                  diffOfMeansMu4MinusMu1, diffOfMeansMu4MinusMu1SecC, diffOfMeansMu4MinusMu1ThdC, diffOfMeansMu4MinusMu1FrtC,

                  diffOfMeansMu3MinusMu2, diffOfMeansMu3MinusMu2SecC, diffOfMeansMu3MinusMu2ThdC, diffOfMeansMu3MinusMu2FrtC,
                  diffOfMeansMu4MinusMu2, diffOfMeansMu4MinusMu2SecC, diffOfMeansMu4MinusMu2ThdC, diffOfMeansMu4MinusMu2FrtC,

                  diffOfMeansMu4MinusMu3, diffOfMeansMu4MinusMu3SecC, diffOfMeansMu4MinusMu3ThdC, diffOfMeansMu4MinusMu3FrtC,

                  
                  diffOfVariancesS2MinusS1, diffOfVariancesS2MinusS1SecC, diffOfVariancesS2MinusS1ThdC, diffOfVariancesS2MinusS1FrtC,
                  diffOfVariancesS3MinusS1, diffOfVariancesS3MinusS1SecC, diffOfVariancesS3MinusS1ThdC, diffOfVariancesS3MinusS1FrtC,
                  diffOfVariancesS4MinusS1, diffOfVariancesS4MinusS1SecC, diffOfVariancesS4MinusS1ThdC, diffOfVariancesS4MinusS1FrtC,

                  diffOfVariancesS3MinusS2, diffOfVariancesS3MinusS2SecC, diffOfVariancesS3MinusS2ThdC, diffOfVariancesS3MinusS2FrtC,
                  diffOfVariancesS4MinusS2, diffOfVariancesS4MinusS2SecC, diffOfVariancesS4MinusS2ThdC, diffOfVariancesS4MinusS2FrtC,

                  diffOfVariancesS4MinusS3, diffOfVariancesS4MinusS3SecC, diffOfVariancesS4MinusS3ThdC, diffOfVariancesS4MinusS3FrtC,

                  
                  effectSize12, effectSize12SecC, effectSize12ThdC, effectSize12FrtC,
                  effectSize13, effectSize13SecC, effectSize13ThdC, effectSize13FrtC,
                  effectSize14, effectSize14SecC, effectSize14ThdC, effectSize14FrtC,

                  effectSize23, effectSize23SecC, effectSize23ThdC, effectSize23FrtC,
                  effectSize24, effectSize24SecC, effectSize24ThdC, effectSize24FrtC,

                  effectSize34, effectSize34SecC, effectSize34ThdC, effectSize34FrtC)
    
    return(df)

  } else if(length(first)>1 && length(second)>1 && length(third)>1 && is.null(fourth) && is.null(fifth) && is.null(sixth)) {
    Parameter<-c("mu1",
                 "mu2",
                 "mu3",
                 "sigma1",
                 "sigma2",
                 "sigma3",
                 "mu2-mu1",
                 "mu3-mu1", 
                 "mu3-mu2",
                 "sigma2-sigma1",
                 "sigma3-sigma1",
                 "sigma3-sigma2",
                 "delta12",
                 "delta13",
                 "delta23")
    LQ<-round(c(as.numeric(mu1Quants[1]), 
                as.numeric(mu2Quants[1]), 
                as.numeric(mu3Quants[1]), 
                as.numeric(sigma1Quants[1]), 
                as.numeric(sigma2Quants[1]), 
                as.numeric(sigma3Quants[1]),
                as.numeric(diffOfMeansMu2MinusMu1Quants[1]), 
                as.numeric(diffOfMeansMu3MinusMu1Quants[1]), 
                as.numeric(diffOfMeansMu3MinusMu2Quants[1]), 
                as.numeric(diffOfVariancesS2MinusS1Quants[1]), 
                as.numeric(diffOfVariancesS3MinusS1Quants[1]), 
                as.numeric(diffOfVariancesS3MinusS2Quants[1]),
                as.numeric(effectSize12Quants[1]), 
                as.numeric(effectSize13Quants[1]), 
                as.numeric(effectSize23Quants[1])),2)
    Mean<-round(c(mean(mu1), 
                  mean(mu2), 
                  mean(mu3), 
                  mean(sigma1Sq), 
                  mean(sigma2Sq), 
                  mean(sigma3Sq), 
                  mean(diffOfMeansMu2MinusMu1), 
                  mean(diffOfMeansMu3MinusMu1), 
                  mean(diffOfMeansMu3MinusMu2), 
                  mean(diffOfVariancesS2MinusS1), 
                  mean(diffOfVariancesS3MinusS1), 
                  mean(diffOfVariancesS3MinusS2), 
                  mean(effectSize12), 
                  mean(effectSize13), 
                  mean(effectSize23)),2)
    UQ<-round(c(as.numeric(mu1Quants[2]), 
                as.numeric(mu2Quants[2]), 
                as.numeric(mu3Quants[2]),
                as.numeric(sigma1Quants[2]), 
                as.numeric(sigma2Quants[2]), 
                as.numeric(sigma3Quants[2]),
                as.numeric(diffOfMeansMu2MinusMu1Quants[2]), 
                as.numeric(diffOfMeansMu3MinusMu1Quants[2]), 
                as.numeric(diffOfMeansMu3MinusMu2Quants[2]), 
                as.numeric(diffOfVariancesS2MinusS1Quants[2]), 
                as.numeric(diffOfVariancesS3MinusS1Quants[2]), 
                as.numeric(diffOfVariancesS3MinusS2Quants[2]),
                as.numeric(effectSize12Quants[2]), 
                as.numeric(effectSize13Quants[2]), 
                as.numeric(effectSize23Quants[2])),2)
    Std.Err<-round(c(sd(mu1), 
                     sd(mu2), 
                     sd(mu3),
                     sd(sigma1Sq), 
                     sd(sigma2Sq), 
                     sd(sigma3Sq),
                     sd(diffOfMeansMu2MinusMu1),
                     sd(diffOfMeansMu3MinusMu1), 
                     sd(diffOfMeansMu3MinusMu2), 
                     sd(diffOfVariancesS2MinusS1),
                     sd(diffOfVariancesS3MinusS1), 
                     sd(diffOfVariancesS3MinusS2), 
                     sd(effectSize12), 
                     sd(effectSize13), 
                     sd(effectSize23)),2)
    #dataframe output
    out <- cbind(Parameter,LQ,Mean,UQ,Std.Err)
    
    #pretty table
    print(knitr::kable(out))
    
    df=data.frame(mu1, mu1SecC, mu1ThdC, mu1FrtC, 
               mu2, mu2SecC, mu2ThdC, mu2FrtC, 
               mu3, mu3SecC, mu3ThdC, mu3FrtC, 
               sigma1Sq, sigma1SqSecC, sigma1SqThdC, sigma1SqFrtC, 
               sigma2Sq, sigma2SqSecC, sigma2SqThdC, sigma2SqFrtC, 
               sigma3Sq, sigma3SqSecC, sigma3SqThdC, sigma3SqFrtC, 
               diffOfMeansMu2MinusMu1, diffOfMeansMu2MinusMu1SecC, diffOfMeansMu2MinusMu1ThdC, diffOfMeansMu2MinusMu1FrtC,
               diffOfMeansMu3MinusMu1, diffOfMeansMu3MinusMu1SecC, diffOfMeansMu3MinusMu1ThdC, diffOfMeansMu3MinusMu1FrtC,
               diffOfMeansMu3MinusMu2, diffOfMeansMu3MinusMu2SecC, diffOfMeansMu3MinusMu2ThdC, diffOfMeansMu3MinusMu2FrtC,
               diffOfVariancesS2MinusS1, diffOfVariancesS2MinusS1SecC, diffOfVariancesS2MinusS1ThdC, diffOfVariancesS2MinusS1FrtC,
               diffOfVariancesS3MinusS1, diffOfVariancesS3MinusS1SecC, diffOfVariancesS3MinusS1ThdC, diffOfVariancesS3MinusS1FrtC,
               diffOfVariancesS3MinusS2, diffOfVariancesS3MinusS2SecC, diffOfVariancesS3MinusS2ThdC, diffOfVariancesS3MinusS2FrtC,
               effectSize12, effectSize12SecC, effectSize12ThdC, effectSize12FrtC,
               effectSize13, effectSize13SecC, effectSize13ThdC, effectSize13FrtC,
               effectSize23, effectSize23SecC, effectSize23ThdC, effectSize23FrtC)
    return(df)
  }
}



anovaplot = function(dataframe, type="rope", sd="sd", ci=0.95){
  # load parameter chains from dataframe
  mu1=dataframe$mu1; mu1SecC=dataframe$mu1SecC; mu1ThdC=dataframe$mu1ThdC; mu1FrtC=dataframe$mu1FrtC;
  mu2=dataframe$mu2; mu2SecC=dataframe$mu2SecC; mu2ThdC=dataframe$mu2ThdC; mu2FrtC=dataframe$mu2FrtC;
  mu3=dataframe$mu3; mu3SecC=dataframe$mu3SecC; mu3ThdC=dataframe$mu3ThdC; mu3FrtC=dataframe$mu3FrtC;
  if(!is.null(dataframe$mu4)){
    mu4=dataframe$mu4; mu4SecC=dataframe$mu4SecC; mu4ThdC=dataframe$mu4ThdC; mu4FrtC=dataframe$mu4FrtC;
  }
  if(!is.null(dataframe$mu5)){
    mu5=dataframe$mu5; mu5SecC=dataframe$mu5SecC; mu5ThdC=dataframe$mu5ThdC; mu5FrtC=dataframe$mu5FrtC;
  }
  if(!is.null(dataframe$mu6)){
    mu6=dataframe$mu6; mu6SecC=dataframe$mu6SecC; mu6ThdC=dataframe$mu6ThdC; mu6FrtC=dataframe$mu6FrtC;
  }
  sigma1Sq=dataframe$sigma1Sq; sigma1SqSecC=dataframe$sigma1SqSecC; sigma1SqThdC=dataframe$sigma1SqThdC; sigma1SqFrtC=dataframe$sigma1SqFrtC;
  sigma2Sq=dataframe$sigma2Sq; sigma2SqSecC=dataframe$sigma2SqSecC; sigma2SqThdC=dataframe$sigma2SqThdC; sigma2SqFrtC=dataframe$sigma2SqFrtC;
  sigma3Sq=dataframe$sigma3Sq; sigma3SqSecC=dataframe$sigma3SqSecC; sigma3SqThdC=dataframe$sigma3SqThdC; sigma3SqFrtC=dataframe$sigma3SqFrtC;
  if(!is.null(dataframe$sigma4Sq)){
    sigma4Sq=dataframe$sigma4Sq; sigma4SqSecC=dataframe$sigma4SqSecC; sigma4SqThdC=dataframe$sigma4SqThdC; sigma4SqFrtC=dataframe$sigma4SqFrtC;
  }
  if(!is.null(dataframe$sigma5Sq)){
    sigma5Sq=dataframe$sigma5Sq; sigma5SqSecC=dataframe$sigma5SqSecC; sigma5SqThdC=dataframe$sigma5SqThdC; sigma5SqFrtC=dataframe$sigma5SqFrtC;
  }
  if(!is.null(dataframe$sigma6Sq)){
    sigma6Sq=dataframe$sigma6Sq; sigma6SqSecC=dataframe$sigma6SqSecC; sigma6SqThdC=dataframe$sigma6SqThdC; sigma6SqFrtC=dataframe$sigma6SqFrtC;
  }
  
  diffOfMeansMu2MinusMu1=dataframe$diffOfMeansMu2MinusMu1; diffOfMeansMu2MinusMu1SecC=dataframe$diffOfMeansMu2MinusMu1SecC; diffOfMeansMu2MinusMu1ThdC=dataframe$diffOfMeansMu2MinusMu1ThdC; diffOfMeansMu2MinusMu1FrtC=dataframe$diffOfMeansMu2MinusMu1FrtC;
  diffOfMeansMu3MinusMu1=dataframe$diffOfMeansMu3MinusMu1; diffOfMeansMu3MinusMu1SecC=dataframe$diffOfMeansMu3MinusMu1SecC; diffOfMeansMu3MinusMu1ThdC=dataframe$diffOfMeansMu3MinusMu1ThdC; diffOfMeansMu3MinusMu1FrtC=dataframe$diffOfMeansMu3MinusMu1FrtC;
  diffOfMeansMu3MinusMu2=dataframe$diffOfMeansMu3MinusMu2;diffOfMeansMu3MinusMu2SecC=dataframe$diffOfMeansMu3MinusMu2SecC; diffOfMeansMu3MinusMu2ThdC=dataframe$diffOfMeansMu3MinusMu2ThdC; diffOfMeansMu3MinusMu2FrtC=dataframe$diffOfMeansMu3MinusMu2FrtC;
  if(!is.null(dataframe$sigma4Sq)){
    diffOfMeansMu4MinusMu1=dataframe$diffOfMeansMu4MinusMu1; diffOfMeansMu4MinusMu1SecC=dataframe$diffOfMeansMu4MinusMu1SecC; diffOfMeansMu4MinusMu1ThdC=dataframe$diffOfMeansMu4MinusMu1ThdC; diffOfMeansMu4MinusMu1FrtC=dataframe$diffOfMeansMu4MinusMu1FrtC;
    diffOfMeansMu4MinusMu2=dataframe$diffOfMeansMu4MinusMu2; diffOfMeansMu4MinusMu2SecC=dataframe$diffOfMeansMu4MinusMu2SecC; diffOfMeansMu4MinusMu2ThdC=dataframe$diffOfMeansMu4MinusMu2ThdC; diffOfMeansMu4MinusMu2FrtC=dataframe$diffOfMeansMu4MinusMu2FrtC;
    diffOfMeansMu4MinusMu3=dataframe$diffOfMeansMu4MinusMu3; diffOfMeansMu4MinusMu3SecC=dataframe$diffOfMeansMu4MinusMu3SecC; diffOfMeansMu4MinusMu3ThdC=dataframe$diffOfMeansMu4MinusMu3ThdC; diffOfMeansMu4MinusMu3FrtC=dataframe$diffOfMeansMu4MinusMu3FrtC;
  }
  if(!is.null(dataframe$sigma5Sq)){
    diffOfMeansMu5MinusMu1=dataframe$diffOfMeansMu5MinusMu1; diffOfMeansMu5MinusMu1SecC=dataframe$diffOfMeansMu5MinusMu1SecC; diffOfMeansMu5MinusMu1ThdC=dataframe$diffOfMeansMu5MinusMu1ThdC; diffOfMeansMu5MinusMu1FrtC=dataframe$diffOfMeansMu5MinusMu1FrtC;
    diffOfMeansMu5MinusMu2=dataframe$diffOfMeansMu5MinusMu2; diffOfMeansMu5MinusMu2SecC=dataframe$diffOfMeansMu5MinusMu2SecC; diffOfMeansMu5MinusMu2ThdC=dataframe$diffOfMeansMu5MinusMu2ThdC; diffOfMeansMu5MinusMu2FrtC=dataframe$diffOfMeansMu5MinusMu2FrtC;
    diffOfMeansMu5MinusMu3=dataframe$diffOfMeansMu5MinusMu3; diffOfMeansMu5MinusMu3SecC=dataframe$diffOfMeansMu5MinusMu3SecC; diffOfMeansMu5MinusMu3ThdC=dataframe$diffOfMeansMu5MinusMu3ThdC; diffOfMeansMu5MinusMu3FrtC=dataframe$diffOfMeansMu5MinusMu3FrtC;
    diffOfMeansMu5MinusMu4=dataframe$diffOfMeansMu5MinusMu4; diffOfMeansMu5MinusMu4SecC=dataframe$diffOfMeansMu5MinusMu4SecC; diffOfMeansMu5MinusMu4ThdC=dataframe$diffOfMeansMu5MinusMu4ThdC; diffOfMeansMu5MinusMu4FrtC=dataframe$diffOfMeansMu5MinusMu4FrtC;
  }
  if(!is.null(dataframe$sigma6Sq)){
    diffOfMeansMu6MinusMu1=dataframe$diffOfMeansMu6MinusMu1; diffOfMeansMu6MinusMu1SecC=dataframe$diffOfMeansMu6MinusMu1SecC; diffOfMeansMu6MinusMu1ThdC=dataframe$diffOfMeansMu6MinusMu1ThdC; diffOfMeansMu6MinusMu1FrtC=dataframe$diffOfMeansMu6MinusMu1FrtC;
    diffOfMeansMu6MinusMu2=dataframe$diffOfMeansMu6MinusMu2; diffOfMeansMu6MinusMu2SecC=dataframe$diffOfMeansMu6MinusMu2SecC; diffOfMeansMu6MinusMu2ThdC=dataframe$diffOfMeansMu6MinusMu2ThdC; diffOfMeansMu6MinusMu2FrtC=dataframe$diffOfMeansMu6MinusMu2FrtC;
    diffOfMeansMu6MinusMu3=dataframe$diffOfMeansMu6MinusMu3; diffOfMeansMu6MinusMu3SecC=dataframe$diffOfMeansMu6MinusMu3SecC; diffOfMeansMu6MinusMu3ThdC=dataframe$diffOfMeansMu6MinusMu3ThdC; diffOfMeansMu6MinusMu3FrtC=dataframe$diffOfMeansMu6MinusMu3FrtC;
    diffOfMeansMu6MinusMu4=dataframe$diffOfMeansMu6MinusMu4; diffOfMeansMu6MinusMu4SecC=dataframe$diffOfMeansMu6MinusMu4SecC; diffOfMeansMu6MinusMu4ThdC=dataframe$diffOfMeansMu6MinusMu4ThdC; diffOfMeansMu6MinusMu4FrtC=dataframe$diffOfMeansMu6MinusMu4FrtC;
    diffOfMeansMu6MinusMu5=dataframe$diffOfMeansMu6MinusMu5; diffOfMeansMu6MinusMu5SecC=dataframe$diffOfMeansMu6MinusMu5SecC; diffOfMeansMu6MinusMu5ThdC=dataframe$diffOfMeansMu6MinusMu5ThdC; diffOfMeansMu6MinusMu5FrtC=dataframe$diffOfMeansMu6MinusMu5FrtC;
  }
  
  diffOfVariancesS2MinusS1=dataframe$diffOfVariancesS2MinusS1; diffOfVariancesS2MinusS1SecC=dataframe$diffOfVariancesS2MinusS1SecC; diffOfVariancesS2MinusS1ThdC=dataframe$diffOfVariancesS2MinusS1ThdC; diffOfVariancesS2MinusS1FrtC=dataframe$diffOfVariancesS2MinusS1FrtC;
  diffOfVariancesS3MinusS1=dataframe$diffOfVariancesS3MinusS1; diffOfVariancesS3MinusS1SecC=dataframe$diffOfVariancesS3MinusS1SecC; diffOfVariancesS3MinusS1ThdC=dataframe$diffOfVariancesS3MinusS1ThdC; diffOfVariancesS3MinusS1FrtC=dataframe$diffOfVariancesS3MinusS1FrtC;
  diffOfVariancesS3MinusS2=dataframe$diffOfVariancesS3MinusS2;diffOfVariancesS3MinusS2SecC=dataframe$diffOfVariancesS3MinusS2SecC; diffOfVariancesS3MinusS2ThdC=dataframe$diffOfVariancesS3MinusS2ThdC; diffOfVariancesS3MinusS2FrtC=dataframe$diffOfVariancesS3MinusS2FrtC;
  if(!is.null(dataframe$sigma4Sq)){
    diffOfVariancesS4MinusS1=dataframe$diffOfVariancesS4MinusS1; diffOfVariancesS4MinusS1SecC=dataframe$diffOfVariancesS4MinusS1SecC; diffOfVariancesS4MinusS1ThdC=dataframe$diffOfVariancesS4MinusS1ThdC; diffOfVariancesS4MinusS1FrtC=dataframe$diffOfVariancesS4MinusS1FrtC;
    diffOfVariancesS4MinusS2=dataframe$diffOfVariancesS4MinusS2; diffOfVariancesS4MinusS2SecC=dataframe$diffOfVariancesS4MinusS2SecC; diffOfVariancesS4MinusS2ThdC=dataframe$diffOfVariancesS4MinusS2ThdC; diffOfVariancesS4MinusS2FrtC=dataframe$diffOfVariancesS4MinusS2FrtC;
    diffOfVariancesS4MinusS3=dataframe$diffOfVariancesS4MinusS3; diffOfVariancesS4MinusS3SecC=dataframe$diffOfVariancesS4MinusS3SecC; diffOfVariancesS4MinusS3ThdC=dataframe$diffOfVariancesS4MinusS3ThdC; diffOfVariancesS4MinusS3FrtC=dataframe$diffOfVariancesS4MinusS3FrtC;
  }
  if(!is.null(dataframe$sigma5Sq)){
    diffOfVariancesS5MinusS1=dataframe$diffOfVariancesS5MinusS1; diffOfVariancesS5MinusS1SecC=dataframe$diffOfVariancesS5MinusS1SecC; diffOfVariancesS5MinusS1ThdC=dataframe$diffOfVariancesS5MinusS1ThdC; diffOfVariancesS5MinusS1FrtC=dataframe$diffOfVariancesS5MinusS1FrtC;
    diffOfVariancesS5MinusS2=dataframe$diffOfVariancesS5MinusS2; diffOfVariancesS5MinusS2SecC=dataframe$diffOfVariancesS5MinusS2SecC; diffOfVariancesS5MinusS2ThdC=dataframe$diffOfVariancesS5MinusS2ThdC; diffOfVariancesS5MinusS2FrtC=dataframe$diffOfVariancesS5MinusS2FrtC;
    diffOfVariancesS5MinusS3=dataframe$diffOfVariancesS5MinusS3; diffOfVariancesS5MinusS3SecC=dataframe$diffOfVariancesS5MinusS3SecC; diffOfVariancesS5MinusS3ThdC=dataframe$diffOfVariancesS5MinusS3ThdC; diffOfVariancesS5MinusS3FrtC=dataframe$diffOfVariancesS5MinusS3FrtC;
    diffOfVariancesS5MinusS4=dataframe$diffOfVariancesS5MinusS4; diffOfVariancesS5MinusS4SecC=dataframe$diffOfVariancesS5MinusS4SecC; diffOfVariancesS5MinusS4ThdC=dataframe$diffOfVariancesS5MinusS4ThdC; diffOfVariancesS5MinusS4FrtC=dataframe$diffOfVariancesS5MinusS4FrtC;
  }
  if(!is.null(dataframe$sigma6Sq)){
    diffOfVariancesS6MinusS1=dataframe$diffOfVariancesS6MinusS1; diffOfVariancesS6MinusS1SecC=dataframe$diffOfVariancesS6MinusS1SecC; diffOfVariancesS6MinusS1ThdC=dataframe$diffOfVariancesS6MinusS1ThdC; diffOfVariancesS6MinusS1FrtC=dataframe$diffOfVariancesS6MinusS1FrtC;
    diffOfVariancesS6MinusS2=dataframe$diffOfVariancesS6MinusS2; diffOfVariancesS6MinusS2SecC=dataframe$diffOfVariancesS6MinusS2SecC; diffOfVariancesS6MinusS2ThdC=dataframe$diffOfVariancesS6MinusS2ThdC; diffOfVariancesS6MinusS2FrtC=dataframe$diffOfVariancesS6MinusS2FrtC;
    diffOfVariancesS6MinusS3=dataframe$diffOfVariancesS6MinusS3; diffOfVariancesS6MinusS3SecC=dataframe$diffOfVariancesS6MinusS3SecC; diffOfVariancesS6MinusS3ThdC=dataframe$diffOfVariancesS6MinusS3ThdC; diffOfVariancesS6MinusS3FrtC=dataframe$diffOfVariancesS6MinusS3FrtC;
    diffOfVariancesS6MinusS4=dataframe$diffOfVariancesS6MinusS4; diffOfVariancesS6MinusS4SecC=dataframe$diffOfVariancesS6MinusS4SecC; diffOfVariancesS6MinusS4ThdC=dataframe$diffOfVariancesS6MinusS4ThdC; diffOfVariancesS6MinusS4FrtC=dataframe$diffOfVariancesS6MinusS4FrtC;
    diffOfVariancesS6MinusS5=dataframe$diffOfVariancesS6MinusS5; diffOfVariancesS6MinusS5SecC=dataframe$diffOfVariancesS6MinusS5SecC; diffOfVariancesS6MinusS5ThdC=dataframe$diffOfVariancesS6MinusS5ThdC; diffOfVariancesS6MinusS5FrtC=dataframe$diffOfVariancesS6MinusS5FrtC;
  }
  
  effectSize12=dataframe$effectSize12; effectSize12SecC=dataframe$effectSize12SecC; effectSize12ThdC=dataframe$effectSize12ThdC; effectSize12FrtC=dataframe$effectSize12FrtC;
  effectSize13=dataframe$effectSize13; effectSize13SecC=dataframe$effectSize13SecC; effectSize13ThdC=dataframe$effectSize13ThdC; effectSize13FrtC=dataframe$effectSize13FrtC;
  effectSize23=dataframe$effectSize23; effectSize23SecC=dataframe$effectSize23SecC; effectSize23ThdC=dataframe$effectSize23ThdC; effectSize23FrtC=dataframe$effectSize23FrtC;
  if(!is.null(dataframe$sigma4Sq)){
    effectSize14=dataframe$effectSize14; effectSize14SecC=dataframe$effectSize14SecC; effectSize14ThdC=dataframe$effectSize14ThdC; effectSize14FrtC=dataframe$effectSize14FrtC;
    effectSize24=dataframe$effectSize24; effectSize24SecC=dataframe$effectSize24SecC; effectSize24ThdC=dataframe$effectSize24ThdC; effectSize24FrtC=dataframe$effectSize24FrtC;
    effectSize34=dataframe$effectSize34; effectSize34SecC=dataframe$effectSize34SecC; effectSize34ThdC=dataframe$effectSize34ThdC; effectSize34FrtC=dataframe$effectSize34FrtC;
  }
  if(!is.null(dataframe$sigma5Sq)){
    effectSize15=dataframe$effectSize15; effectSize15SecC=dataframe$effectSize15SecC; effectSize15ThdC=dataframe$effectSize15ThdC; effectSize15FrtC=dataframe$effectSize15FrtC;
    effectSize25=dataframe$effectSize25; effectSize25SecC=dataframe$effectSize25SecC; effectSize25ThdC=dataframe$effectSize25ThdC; effectSize25FrtC=dataframe$effectSize25FrtC;
    effectSize35=dataframe$effectSize35; effectSize35SecC=dataframe$effectSize35SecC; effectSize35ThdC=dataframe$effectSize35ThdC; effectSize35FrtC=dataframe$effectSize35FrtC;
    effectSize45=dataframe$effectSize45; effectSize45SecC=dataframe$effectSize45SecC; effectSize45ThdC=dataframe$effectSize45ThdC; effectSize45FrtC=dataframe$effectSize45FrtC;
  }
  if(!is.null(dataframe$sigma6Sq)){
    effectSize16=dataframe$effectSize16; effectSize16SecC=dataframe$effectSize16SecC; effectSize16ThdC=dataframe$effectSize16ThdC; effectSize16FrtC=dataframe$effectSize16FrtC;
    effectSize26=dataframe$effectSize26; effectSize26SecC=dataframe$effectSize26SecC; effectSize26ThdC=dataframe$effectSize26ThdC; effectSize26FrtC=dataframe$effectSize26FrtC;
    effectSize36=dataframe$effectSize36; effectSize36SecC=dataframe$effectSize36SecC; effectSize36ThdC=dataframe$effectSize36ThdC; effectSize36FrtC=dataframe$effectSize36FrtC;
    effectSize46=dataframe$effectSize46; effectSize46SecC=dataframe$effectSize46SecC; effectSize46ThdC=dataframe$effectSize46ThdC; effectSize46FrtC=dataframe$effectSize46FrtC;
    effectSize56=dataframe$effectSize56; effectSize56SecC=dataframe$effectSize56SecC; effectSize56ThdC=dataframe$effectSize56ThdC; effectSize56FrtC=dataframe$effectSize56FrtC;
  }
  

  
  # posterior plots for each parameter
  if(type=="pars"){
    # MU1
    dev.new();
    plotpar<-par(mfrow=c(2,2),mai = c(1,0.5,0.5,0.5), oma=c(1,1,1,1))
    # Reset par only after first change of par
    on.exit(par(plotpar))
    
    hist(mu1,freq=FALSE,main="",col="cornflowerblue",border="white",xlab=expression(mu[1]))
    plot(mu1,ty="l",col="cornflowerblue",xlab=expression(mu[1]),ylab="")
    mcmcObj1=coda::mcmc(data=mu1)
    mcmcObj2=coda::mcmc(data=mu1SecC)
    mcmcObj3=coda::mcmc(data=mu1ThdC)
    mcmcObj4=coda::mcmc(data=mu1FrtC)
    mcmcListObj=coda::mcmc.list(mcmcObj1,mcmcObj2,mcmcObj3,mcmcObj4)
    coda::gelman.plot(mcmcListObj, auto.layout = FALSE)
    stats::acf(mu1)
    
    
    # MU2
    dev.new();
    plotpar<-par(mfrow=c(2,2),mai = c(1,0.5,0.5,0.5), oma=c(1,1,1,1))
    # Reset par only after first change of par
    on.exit(par(plotpar))
    
    hist(mu2,freq=FALSE,main="",col="cornflowerblue",border="white",xlab=expression(mu[2]))
    plot(mu2,ty="l",col="cornflowerblue",xlab=expression(mu[2]),ylab="")
    mcmcObj1=coda::mcmc(data=mu2)
    mcmcObj2=coda::mcmc(data=mu2SecC)
    mcmcObj3=coda::mcmc(data=mu2ThdC)
    mcmcObj4=coda::mcmc(data=mu2FrtC)
    mcmcListObj=coda::mcmc.list(mcmcObj1,mcmcObj2,mcmcObj3,mcmcObj4)
    coda::gelman.plot(mcmcListObj,auto.layout = FALSE)
    stats::acf(mu2)
    
    # MU3
    dev.new();
    plotpar<-par(mfrow=c(2,2),mai = c(1,0.5,0.5,0.5), oma=c(1,1,1,1))
    # Reset par only after first change of par
    on.exit(par(plotpar))
    
    hist(mu3,freq=FALSE,main="",col="cornflowerblue",border="white",xlab=expression(mu[3]))
    plot(mu3,ty="l",col="cornflowerblue",xlab=expression(mu[3]),ylab="")
    mcmcObj1=coda::mcmc(data=mu3)
    mcmcObj2=coda::mcmc(data=mu3SecC)
    mcmcObj3=coda::mcmc(data=mu3ThdC)
    mcmcObj4=coda::mcmc(data=mu3FrtC)
    mcmcListObj=coda::mcmc.list(mcmcObj1,mcmcObj2,mcmcObj3,mcmcObj4)
    coda::gelman.plot(mcmcListObj,auto.layout = FALSE)
    stats::acf(mu3)
    
    
    # MU4
    if(!is.null(dataframe$mu4)){
      dev.new();
      plotpar<-par(mfrow=c(2,2),mai = c(1,0.5,0.5,0.5), oma=c(1,1,1,1))
      # Reset par only after first change of par
      on.exit(par(plotpar))
      
      hist(mu4,freq=FALSE,main="",col="cornflowerblue",border="white",xlab=expression(mu[4]))
      plot(mu4,ty="l",col="cornflowerblue",xlab=expression(mu[4]),ylab="")
      mcmcObj1=coda::mcmc(data=mu4)
      mcmcObj2=coda::mcmc(data=mu4SecC)
      mcmcObj3=coda::mcmc(data=mu4ThdC)
      mcmcObj4=coda::mcmc(data=mu4FrtC)
      mcmcListObj=coda::mcmc.list(mcmcObj1,mcmcObj2,mcmcObj3,mcmcObj4)
      coda::gelman.plot(mcmcListObj, auto.layout = FALSE)
      stats::acf(mu4)
    }
    
    
    # MU5
    if(!is.null(dataframe$mu5)){
      dev.new();
      plotpar<-par(mfrow=c(2,2),mai = c(1,0.5,0.5,0.5), oma=c(1,1,1,1))
      # Reset par only after first change of par
      on.exit(par(plotpar))
      
      hist(mu5,freq=FALSE,main="",col="cornflowerblue",border="white",xlab=expression(mu[5]))
      plot(mu5,ty="l",col="cornflowerblue",xlab=expression(mu[5]),ylab="")
      mcmcObj1=coda::mcmc(data=mu5)
      mcmcObj2=coda::mcmc(data=mu5SecC)
      mcmcObj3=coda::mcmc(data=mu5ThdC)
      mcmcObj4=coda::mcmc(data=mu5FrtC)
      mcmcListObj=coda::mcmc.list(mcmcObj1,mcmcObj2,mcmcObj3,mcmcObj4)
      coda::gelman.plot(mcmcListObj, auto.layout = FALSE)
      stats::acf(mu5)
    }
    
    
    # MU6
    if(!is.null(dataframe$mu6)){
      dev.new();
      plotpar<-par(mfrow=c(2,2),mai = c(1,0.5,0.5,0.5), oma=c(1,1,1,1))
      # Reset par only after first change of par
      on.exit(par(plotpar))
      
      hist(mu6,freq=FALSE,main="",col="cornflowerblue",border="white",xlab=expression(mu[6]))
      plot(mu6,ty="l",col="cornflowerblue",xlab=expression(mu[6]),ylab="")
      mcmcObj1=coda::mcmc(data=mu6)
      mcmcObj2=coda::mcmc(data=mu6SecC)
      mcmcObj3=coda::mcmc(data=mu6ThdC)
      mcmcObj4=coda::mcmc(data=mu6FrtC)
      mcmcListObj=coda::mcmc.list(mcmcObj1,mcmcObj2,mcmcObj3,mcmcObj4)
      coda::gelman.plot(mcmcListObj, auto.layout = FALSE)
      stats::acf(mu6)
    }
    
    
    # SD1
    dev.new();
    plotpar<-par(mfrow=c(2,2),mai = c(1,0.5,0.5,0.5), oma=c(1,1,1,1))
    # Reset par only after first change of par
    on.exit(par(plotpar))
    
    if(sd=="var"){
      hist(sigma1Sq,freq=FALSE,main="",col="cornflowerblue",border="white",xlab=expression(sigma[1]^2))
      plot(sigma1Sq,ty="l",col="cornflowerblue",xlab=expression(sigma[1]^2),ylab="")
    }
    if(sd=="sd"){
      hist(sigma1Sq,freq=FALSE,main="",col="cornflowerblue",border="white",xlab=expression(sigma[1]))
      plot(sigma1Sq,ty="l",col="cornflowerblue",xlab=expression(sigma[1]),ylab="")
    }
    mcmcObj1=coda::mcmc(data=sigma1Sq)
    mcmcObj2=coda::mcmc(data=sigma1SqSecC)
    mcmcObj3=coda::mcmc(data=sigma1SqThdC)
    mcmcObj4=coda::mcmc(data=sigma1SqFrtC)
    mcmcListObj=coda::mcmc.list(mcmcObj1,mcmcObj2,mcmcObj3,mcmcObj4)
    coda::gelman.plot(mcmcListObj, auto.layout = FALSE)
    stats::acf(sigma1Sq)
    
    
    # SD2
    dev.new()
    plotpar<-par(mfrow=c(2,2),mai = c(1,0.5,0.5,0.5), oma=c(1,1,1,1))
    # Reset par only after first change of par
    on.exit(par(plotpar))
    
    if(sd=="var"){
      hist(sigma2Sq,freq=FALSE,main="",col="cornflowerblue",border="white",xlab=expression(sigma[2]^2))
      plot(sigma2Sq,ty="l",col="cornflowerblue",xlab=expression(sigma[2]^2),ylab="")
    }
    if(sd=="sd"){
      hist(sigma2Sq,freq=FALSE,main="",col="cornflowerblue",border="white",xlab=expression(sigma[2]))
      plot(sigma2Sq,ty="l",col="cornflowerblue",xlab=expression(sigma[2]),ylab="")
    }
    mcmcObj1=coda::mcmc(data=sigma2Sq)
    mcmcObj2=coda::mcmc(data=sigma2SqSecC)
    mcmcObj3=coda::mcmc(data=sigma2SqThdC)
    mcmcObj4=coda::mcmc(data=sigma2SqFrtC)
    mcmcListObj=coda::mcmc.list(mcmcObj1,mcmcObj2,mcmcObj3,mcmcObj4)
    coda::gelman.plot(mcmcListObj, auto.layout = FALSE)
    stats::acf(sigma2Sq)
    
    
    # SD3
    dev.new()
    plotpar<-par(mfrow=c(2,2),mai = c(1,0.5,0.5,0.5), oma=c(1,1,1,1))
    # Reset par only after first change of par
    on.exit(par(plotpar))
    
    if(sd=="var"){
      hist(sigma3Sq,freq=FALSE,main="",col="cornflowerblue",border="white",xlab=expression(sigma[3]^2))
      plot(sigma3Sq,ty="l",col="cornflowerblue",xlab=expression(sigma[3]^2),ylab="")
    }
    if(sd=="sd"){
      hist(sigma3Sq,freq=FALSE,main="",col="cornflowerblue",border="white",xlab=expression(sigma[3]))
      plot(sigma3Sq,ty="l",col="cornflowerblue",xlab=expression(sigma[3]),ylab="")
    }
    mcmcObj1=coda::mcmc(data=sigma3Sq)
    mcmcObj2=coda::mcmc(data=sigma3SqSecC)
    mcmcObj3=coda::mcmc(data=sigma3SqThdC)
    mcmcObj4=coda::mcmc(data=sigma3SqFrtC)
    mcmcListObj=coda::mcmc.list(mcmcObj1,mcmcObj2,mcmcObj3,mcmcObj4)
    coda::gelman.plot(mcmcListObj, auto.layout = FALSE)
    stats::acf(sigma3Sq)
  
  
    # SD4
    if(!is.null(dataframe$mu4)){
      dev.new();
      plotpar<-par(mfrow=c(2,2),mai = c(1,0.5,0.5,0.5), oma=c(1,1,1,1))
      # Reset par only after first change of par
      on.exit(par(plotpar))
      
      if(sd=="var"){
        hist(sigma4Sq,freq=FALSE,main="",col="cornflowerblue",border="white",xlab=expression(sigma[4]^2))
        plot(sigma4Sq,ty="l",col="cornflowerblue",xlab=expression(sigma[4]^2),ylab="")
      }
      if(sd=="sd"){
        hist(sigma4Sq,freq=FALSE,main="",col="cornflowerblue",border="white",xlab=expression(sigma[4]))
        plot(sigma4Sq,ty="l",col="cornflowerblue",xlab=expression(sigma[4]),ylab="")
      }
      mcmcObj1=coda::mcmc(data=sigma4Sq)
      mcmcObj2=coda::mcmc(data=sigma4SqSecC)
      mcmcObj3=coda::mcmc(data=sigma4SqThdC)
      mcmcObj4=coda::mcmc(data=sigma4SqFrtC)
      mcmcListObj=coda::mcmc.list(mcmcObj1,mcmcObj2,mcmcObj3,mcmcObj4)
      coda::gelman.plot(mcmcListObj, auto.layout = FALSE)
      stats::acf(sigma4Sq)
    }
    
    
    # SD5
    if(!is.null(dataframe$mu5)){
      dev.new();
      plotpar<-par(mfrow=c(2,2),mai = c(1,0.5,0.5,0.5), oma=c(1,1,1,1))
      # Reset par only after first change of par
      on.exit(par(plotpar))
      
      if(sd=="var"){
        hist(sigma5Sq,freq=FALSE,main="",col="cornflowerblue",border="white",xlab=expression(sigma[5]^2))
        plot(sigma5Sq,ty="l",col="cornflowerblue",xlab=expression(sigma[5]^2),ylab="")
      }
      if(sd=="sd"){
        hist(sigma5Sq,freq=FALSE,main="",col="cornflowerblue",border="white",xlab=expression(sigma[5]))
        plot(sigma5Sq,ty="l",col="cornflowerblue",xlab=expression(sigma[5]),ylab="")
      }
      mcmcObj1=coda::mcmc(data=sigma5Sq)
      mcmcObj2=coda::mcmc(data=sigma5SqSecC)
      mcmcObj3=coda::mcmc(data=sigma5SqThdC)
      mcmcObj4=coda::mcmc(data=sigma5SqFrtC)
      mcmcListObj=coda::mcmc.list(mcmcObj1,mcmcObj2,mcmcObj3,mcmcObj4)
      coda::gelman.plot(mcmcListObj, auto.layout = FALSE)
      stats::acf(sigma5Sq)
    }
    
    
    # SD6
    if(!is.null(dataframe$mu6)){
      dev.new();
      plotpar<-par(mfrow=c(2,2),mai = c(1,0.5,0.5,0.5), oma=c(1,1,1,1))
      # Reset par only after first change of par
      on.exit(par(plotpar))
      
      if(sd=="var"){
        hist(sigma6Sq,freq=FALSE,main="",col="cornflowerblue",border="white",xlab=expression(sigma[6]^2))
        plot(sigma6Sq,ty="l",col="cornflowerblue",xlab=expression(sigma[6]^2),ylab="")
      }
      if(sd=="sd"){
        hist(sigma6Sq,freq=FALSE,main="",col="cornflowerblue",border="white",xlab=expression(sigma[6]))
        plot(sigma6Sq,ty="l",col="cornflowerblue",xlab=expression(sigma[6]),ylab="")
      }
      mcmcObj1=coda::mcmc(data=sigma6Sq)
      mcmcObj2=coda::mcmc(data=sigma6SqSecC)
      mcmcObj3=coda::mcmc(data=sigma6SqThdC)
      mcmcObj4=coda::mcmc(data=sigma6SqFrtC)
      mcmcListObj=coda::mcmc.list(mcmcObj1,mcmcObj2,mcmcObj3,mcmcObj4)
      coda::gelman.plot(mcmcListObj, auto.layout = FALSE)
      stats::acf(sigma6Sq)
    }
  }
  
  
  # posterior plot of difference of means & variances
  if(type=="diff"){
    dev.new()
    plotpar<-par(mfrow=c(2,2),mai = c(1,0.5,0.5,0.5), oma=c(1,1,1,1))
    # Reset par only after first change of par
    on.exit(par(plotpar))
    hist(diffOfMeansMu2MinusMu1,freq=FALSE,main="Difference of means Group2-Group1",col="cornflowerblue",border="white",xlab=expression(mu[2]-mu[1]))
    plot(diffOfMeansMu2MinusMu1,ty="l",col="cornflowerblue",xlab=expression(mu[2]-mu[1]),ylab="")
    mcmcObj1=coda::mcmc(data=diffOfMeansMu2MinusMu1)
    mcmcObj2=coda::mcmc(data=diffOfMeansMu2MinusMu1SecC)
    mcmcObj3=coda::mcmc(data=diffOfMeansMu2MinusMu1ThdC)
    mcmcObj4=coda::mcmc(data=diffOfMeansMu2MinusMu1FrtC)
    mcmcListObj=coda::mcmc.list(mcmcObj1,mcmcObj2,mcmcObj3,mcmcObj4)
    coda::gelman.plot(mcmcListObj, auto.layout = FALSE)
    stats::acf(diffOfMeansMu2MinusMu1)
    
    
    dev.new()
    plotpar<-par(mfrow=c(2,2),mai = c(1,0.5,0.5,0.5), oma=c(1,1,1,1))
    # Reset par only after first change of par
    on.exit(par(plotpar))
    hist(diffOfMeansMu3MinusMu2,freq=FALSE,main="Difference of means Group3-Group2",col="cornflowerblue",border="white",xlab=expression(mu[3]-mu[2]))
    plot(diffOfMeansMu3MinusMu2,ty="l",col="cornflowerblue",xlab=expression(mu[3]-mu[2]),ylab="")
    mcmcObj1=coda::mcmc(data=diffOfMeansMu3MinusMu2)
    mcmcObj2=coda::mcmc(data=diffOfMeansMu3MinusMu2SecC)
    mcmcObj3=coda::mcmc(data=diffOfMeansMu3MinusMu2ThdC)
    mcmcObj4=coda::mcmc(data=diffOfMeansMu3MinusMu2FrtC)
    mcmcListObj=coda::mcmc.list(mcmcObj1,mcmcObj2,mcmcObj3,mcmcObj4)
    coda::gelman.plot(mcmcListObj, auto.layout = FALSE)
    stats::acf(diffOfMeansMu3MinusMu2)
    
    
    dev.new()
    plotpar<-par(mfrow=c(2,2),mai = c(1,0.5,0.5,0.5), oma=c(1,1,1,1))
    # Reset par only after first change of par
    on.exit(par(plotpar))
    hist(diffOfMeansMu3MinusMu1,freq=FALSE,main="Difference of means Group3-Group1",col="cornflowerblue",border="white",xlab=expression(mu[3]-mu[1]))
    plot(diffOfMeansMu3MinusMu1,ty="l",col="cornflowerblue",xlab=expression(mu[3]-mu[1]),ylab="")
    mcmcObj1=coda::mcmc(data=diffOfMeansMu3MinusMu1)
    mcmcObj2=coda::mcmc(data=diffOfMeansMu3MinusMu1SecC)
    mcmcObj3=coda::mcmc(data=diffOfMeansMu3MinusMu1ThdC)
    mcmcObj4=coda::mcmc(data=diffOfMeansMu3MinusMu1FrtC)
    mcmcListObj=coda::mcmc.list(mcmcObj1,mcmcObj2,mcmcObj3,mcmcObj4)
    coda::gelman.plot(mcmcListObj, auto.layout = FALSE)
    stats::acf(diffOfMeansMu3MinusMu1)
    
    
    if(!is.null(dataframe$mu4)){
      dev.new()
      plotpar<-par(mfrow=c(2,2),mai = c(1,0.5,0.5,0.5), oma=c(1,1,1,1))
      # Reset par only after first change of par
      on.exit(par(plotpar))
      hist(diffOfMeansMu4MinusMu1,freq=FALSE,main="Difference of means Group4-Group1",col="cornflowerblue",border="white",xlab=expression(mu[4]-mu[1]))
      plot(diffOfMeansMu4MinusMu1,ty="l",col="cornflowerblue",xlab=expression(mu[4]-mu[1]),ylab="")
      mcmcObj1=coda::mcmc(data=diffOfMeansMu4MinusMu1)
      mcmcObj2=coda::mcmc(data=diffOfMeansMu4MinusMu1SecC)
      mcmcObj3=coda::mcmc(data=diffOfMeansMu4MinusMu1ThdC)
      mcmcObj4=coda::mcmc(data=diffOfMeansMu4MinusMu1FrtC)
      mcmcListObj=coda::mcmc.list(mcmcObj1,mcmcObj2,mcmcObj3,mcmcObj4)
      coda::gelman.plot(mcmcListObj, auto.layout = FALSE)
      stats::acf(diffOfMeansMu4MinusMu1)
      
      
      dev.new()
      plotpar<-par(mfrow=c(2,2),mai = c(1,0.5,0.5,0.5), oma=c(1,1,1,1))
      # Reset par only after first change of par
      on.exit(par(plotpar))
      hist(diffOfMeansMu4MinusMu2,freq=FALSE,main="Difference of means Group2-Group2",col="cornflowerblue",border="white",xlab=expression(mu[4]-mu[2]))
      plot(diffOfMeansMu4MinusMu2,ty="l",col="cornflowerblue",xlab=expression(mu[4]-mu[2]),ylab="")
      mcmcObj1=coda::mcmc(data=diffOfMeansMu4MinusMu2)
      mcmcObj2=coda::mcmc(data=diffOfMeansMu4MinusMu2SecC)
      mcmcObj3=coda::mcmc(data=diffOfMeansMu4MinusMu2ThdC)
      mcmcObj4=coda::mcmc(data=diffOfMeansMu4MinusMu2FrtC)
      mcmcListObj=coda::mcmc.list(mcmcObj1,mcmcObj2,mcmcObj3,mcmcObj4)
      coda::gelman.plot(mcmcListObj, auto.layout = FALSE)
      stats::acf(diffOfMeansMu4MinusMu2)
      
      
      dev.new()
      plotpar<-par(mfrow=c(2,2),mai = c(1,0.5,0.5,0.5), oma=c(1,1,1,1))
      # Reset par only after first change of par
      on.exit(par(plotpar))
      hist(diffOfMeansMu4MinusMu3,freq=FALSE,main="Difference of means Group4-Group3",col="cornflowerblue",border="white",xlab=expression(mu[4]-mu[3]))
      plot(diffOfMeansMu4MinusMu3,ty="l",col="cornflowerblue",xlab=expression(mu[4]-mu[3]),ylab="")
      mcmcObj1=coda::mcmc(data=diffOfMeansMu4MinusMu3)
      mcmcObj2=coda::mcmc(data=diffOfMeansMu4MinusMu3SecC)
      mcmcObj3=coda::mcmc(data=diffOfMeansMu4MinusMu3ThdC)
      mcmcObj4=coda::mcmc(data=diffOfMeansMu4MinusMu3FrtC)
      mcmcListObj=coda::mcmc.list(mcmcObj1,mcmcObj2,mcmcObj3,mcmcObj4)
      coda::gelman.plot(mcmcListObj, auto.layout = FALSE)
      stats::acf(diffOfMeansMu4MinusMu3)
    }
    
    if(!is.null(dataframe$mu5)){
      dev.new()
      plotpar<-par(mfrow=c(2,2),mai = c(1,0.5,0.5,0.5), oma=c(1,1,1,1))
      # Reset par only after first change of par
      on.exit(par(plotpar))
      hist(diffOfMeansMu5MinusMu1,freq=FALSE,main="Difference of means Group5-Group1",col="cornflowerblue",border="white",xlab=expression(mu[5]-mu[1]))
      plot(diffOfMeansMu5MinusMu1,ty="l",col="cornflowerblue",xlab=expression(mu[5]-mu[1]),ylab="")
      mcmcObj1=coda::mcmc(data=diffOfMeansMu5MinusMu1)
      mcmcObj2=coda::mcmc(data=diffOfMeansMu5MinusMu1SecC)
      mcmcObj3=coda::mcmc(data=diffOfMeansMu5MinusMu1ThdC)
      mcmcObj4=coda::mcmc(data=diffOfMeansMu5MinusMu1FrtC)
      mcmcListObj=coda::mcmc.list(mcmcObj1,mcmcObj2,mcmcObj3,mcmcObj4)
      coda::gelman.plot(mcmcListObj, auto.layout = FALSE)
      stats::acf(diffOfMeansMu5MinusMu1)
      
      
      dev.new()
      plotpar<-par(mfrow=c(2,2),mai = c(1,0.5,0.5,0.5), oma=c(1,1,1,1))
      # Reset par only after first change of par
      on.exit(par(plotpar))
      hist(diffOfMeansMu5MinusMu2,freq=FALSE,main="Difference of means Group5-Group2",col="cornflowerblue",border="white",xlab=expression(mu[5]-mu[2]))
      plot(diffOfMeansMu5MinusMu2,ty="l",col="cornflowerblue",xlab=expression(mu[5]-mu[2]),ylab="")
      mcmcObj1=coda::mcmc(data=diffOfMeansMu5MinusMu2)
      mcmcObj2=coda::mcmc(data=diffOfMeansMu5MinusMu2SecC)
      mcmcObj3=coda::mcmc(data=diffOfMeansMu5MinusMu2ThdC)
      mcmcObj4=coda::mcmc(data=diffOfMeansMu5MinusMu2FrtC)
      mcmcListObj=coda::mcmc.list(mcmcObj1,mcmcObj2,mcmcObj3,mcmcObj4)
      coda::gelman.plot(mcmcListObj, auto.layout = FALSE)
      stats::acf(diffOfMeansMu5MinusMu2)
      
      
      dev.new()
      plotpar<-par(mfrow=c(2,2),mai = c(1,0.5,0.5,0.5), oma=c(1,1,1,1))
      # Reset par only after first change of par
      on.exit(par(plotpar))
      hist(diffOfMeansMu5MinusMu3,freq=FALSE,main="Difference of means Group5-Group3",col="cornflowerblue",border="white",xlab=expression(mu[5]-mu[3]))
      plot(diffOfMeansMu5MinusMu3,ty="l",col="cornflowerblue",xlab=expression(mu[5]-mu[3]),ylab="")
      mcmcObj1=coda::mcmc(data=diffOfMeansMu5MinusMu3)
      mcmcObj2=coda::mcmc(data=diffOfMeansMu5MinusMu3SecC)
      mcmcObj3=coda::mcmc(data=diffOfMeansMu5MinusMu3ThdC)
      mcmcObj4=coda::mcmc(data=diffOfMeansMu5MinusMu3FrtC)
      mcmcListObj=coda::mcmc.list(mcmcObj1,mcmcObj2,mcmcObj3,mcmcObj4)
      coda::gelman.plot(mcmcListObj, auto.layout = FALSE)
      stats::acf(diffOfMeansMu5MinusMu3)
      
      
      dev.new()
      plotpar<-par(mfrow=c(2,2),mai = c(1,0.5,0.5,0.5), oma=c(1,1,1,1))
      # Reset par only after first change of par
      on.exit(par(plotpar))
      hist(diffOfMeansMu5MinusMu4,freq=FALSE,main="Difference of means Group5-Group4",col="cornflowerblue",border="white",xlab=expression(mu[5]-mu[4]))
      plot(diffOfMeansMu5MinusMu4,ty="l",col="cornflowerblue",xlab=expression(mu[5]-mu[4]),ylab="")
      mcmcObj1=coda::mcmc(data=diffOfMeansMu5MinusMu4)
      mcmcObj2=coda::mcmc(data=diffOfMeansMu5MinusMu4SecC)
      mcmcObj3=coda::mcmc(data=diffOfMeansMu5MinusMu4ThdC)
      mcmcObj4=coda::mcmc(data=diffOfMeansMu5MinusMu4FrtC)
      mcmcListObj=coda::mcmc.list(mcmcObj1,mcmcObj2,mcmcObj3,mcmcObj4)
      coda::gelman.plot(mcmcListObj, auto.layout = FALSE)
      stats::acf(diffOfMeansMu5MinusMu4)
    }
    
    if(!is.null(dataframe$mu6)){
      dev.new()
      plotpar<-par(mfrow=c(2,2),mai = c(1,0.5,0.5,0.5), oma=c(1,1,1,1))
      # Reset par only after first change of par
      on.exit(par(plotpar))
      hist(diffOfMeansMu6MinusMu1,freq=FALSE,main="Difference of means Group6-Group1",col="cornflowerblue",border="white",xlab=expression(mu[6]-mu[1]))
      plot(diffOfMeansMu6MinusMu1,ty="l",col="cornflowerblue",xlab=expression(mu[6]-mu[1]),ylab="")
      mcmcObj1=coda::mcmc(data=diffOfMeansMu6MinusMu1)
      mcmcObj2=coda::mcmc(data=diffOfMeansMu6MinusMu1SecC)
      mcmcObj3=coda::mcmc(data=diffOfMeansMu6MinusMu1ThdC)
      mcmcObj4=coda::mcmc(data=diffOfMeansMu6MinusMu1FrtC)
      mcmcListObj=coda::mcmc.list(mcmcObj1,mcmcObj2,mcmcObj3,mcmcObj4)
      coda::gelman.plot(mcmcListObj, auto.layout = FALSE)
      stats::acf(diffOfMeansMu6MinusMu1)
      
      
      dev.new()
      plotpar<-par(mfrow=c(2,2),mai = c(1,0.5,0.5,0.5), oma=c(1,1,1,1))
      # Reset par only after first change of par
      on.exit(par(plotpar))
      hist(diffOfMeansMu6MinusMu2,freq=FALSE,main="Difference of means Group6-Group2",col="cornflowerblue",border="white",xlab=expression(mu[6]-mu[2]))
      plot(diffOfMeansMu6MinusMu2,ty="l",col="cornflowerblue",xlab=expression(mu[6]-mu[2]),ylab="")
      mcmcObj1=coda::mcmc(data=diffOfMeansMu6MinusMu2)
      mcmcObj2=coda::mcmc(data=diffOfMeansMu6MinusMu2SecC)
      mcmcObj3=coda::mcmc(data=diffOfMeansMu6MinusMu2ThdC)
      mcmcObj4=coda::mcmc(data=diffOfMeansMu6MinusMu2FrtC)
      mcmcListObj=coda::mcmc.list(mcmcObj1,mcmcObj2,mcmcObj3,mcmcObj4)
      coda::gelman.plot(mcmcListObj, auto.layout = FALSE)
      stats::acf(diffOfMeansMu6MinusMu2)
      
      
      dev.new()
      plotpar<-par(mfrow=c(2,2),mai = c(1,0.5,0.5,0.5), oma=c(1,1,1,1))
      # Reset par only after first change of par
      on.exit(par(plotpar))
      hist(diffOfMeansMu6MinusMu3,freq=FALSE,main="Difference of means Group6-Group3",col="cornflowerblue",border="white",xlab=expression(mu[6]-mu[3]))
      plot(diffOfMeansMu6MinusMu3,ty="l",col="cornflowerblue",xlab=expression(mu[6]-mu[3]),ylab="")
      mcmcObj1=coda::mcmc(data=diffOfMeansMu6MinusMu3)
      mcmcObj2=coda::mcmc(data=diffOfMeansMu6MinusMu3SecC)
      mcmcObj3=coda::mcmc(data=diffOfMeansMu6MinusMu3ThdC)
      mcmcObj4=coda::mcmc(data=diffOfMeansMu6MinusMu3FrtC)
      mcmcListObj=coda::mcmc.list(mcmcObj1,mcmcObj2,mcmcObj3,mcmcObj4)
      coda::gelman.plot(mcmcListObj, auto.layout = FALSE)
      stats::acf(diffOfMeansMu6MinusMu3)
      
      
      dev.new()
      plotpar<-par(mfrow=c(2,2),mai = c(1,0.5,0.5,0.5), oma=c(1,1,1,1))
      # Reset par only after first change of par
      on.exit(par(plotpar))
      hist(diffOfMeansMu6MinusMu4,freq=FALSE,main="Difference of means Group6-Group4",col="cornflowerblue",border="white",xlab=expression(mu[6]-mu[4]))
      plot(diffOfMeansMu6MinusMu4,ty="l",col="cornflowerblue",xlab=expression(mu[6]-mu[4]),ylab="")
      mcmcObj1=coda::mcmc(data=diffOfMeansMu6MinusMu4)
      mcmcObj2=coda::mcmc(data=diffOfMeansMu6MinusMu4SecC)
      mcmcObj3=coda::mcmc(data=diffOfMeansMu6MinusMu4ThdC)
      mcmcObj4=coda::mcmc(data=diffOfMeansMu6MinusMu4FrtC)
      mcmcListObj=coda::mcmc.list(mcmcObj1,mcmcObj2,mcmcObj3,mcmcObj4)
      coda::gelman.plot(mcmcListObj, auto.layout = FALSE)
      stats::acf(diffOfMeansMu6MinusMu4)
      
      
      dev.new()
      plotpar<-par(mfrow=c(2,2),mai = c(1,0.5,0.5,0.5), oma=c(1,1,1,1))
      # Reset par only after first change of par
      on.exit(par(plotpar))
      hist(diffOfMeansMu6MinusMu5,freq=FALSE,main="Difference of means Group6-Group5",col="cornflowerblue",border="white",xlab=expression(mu[6]-mu[5]))
      plot(diffOfMeansMu6MinusMu5,ty="l",col="cornflowerblue",xlab=expression(mu[6]-mu[5]),ylab="")
      mcmcObj1=coda::mcmc(data=diffOfMeansMu6MinusMu5)
      mcmcObj2=coda::mcmc(data=diffOfMeansMu6MinusMu5SecC)
      mcmcObj3=coda::mcmc(data=diffOfMeansMu6MinusMu5ThdC)
      mcmcObj4=coda::mcmc(data=diffOfMeansMu6MinusMu5FrtC)
      mcmcListObj=coda::mcmc.list(mcmcObj1,mcmcObj2,mcmcObj3,mcmcObj4)
      coda::gelman.plot(mcmcListObj, auto.layout = FALSE)
      stats::acf(diffOfMeansMu6MinusMu5)
    }
    
    
    
    # DIFFERENCES OF VARIANCES / STANDARD DEVIATIONS
    dev.new()
    plotpar<-par(mfrow=c(2,2),mai = c(1,0.5,0.5,0.5), oma=c(1,1,1,1))
    # Reset par only after first change of par
    on.exit(par(plotpar))
    if(sd=="var"){
      hist(diffOfVariancesS2MinusS1 ,freq=FALSE,main="Difference of variances Group2-Group1",col="cornflowerblue",border="white",
           xlab=expression(sigma[2]^2-sigma[1]^2))
      plot(diffOfVariancesS2MinusS1,ty="l",col="cornflowerblue",xlab=expression(sigma[2]^2-sigma[1]^2),ylab="")
    }
    if(sd=="sd"){
      hist(diffOfVariancesS2MinusS1,freq=FALSE,main="Difference of standard deviations Group2-Group1",col="cornflowerblue",border="white",
           xlab=expression(sigma[2]-sigma[1]))
      plot(diffOfVariancesS2MinusS1,ty="l",col="cornflowerblue",xlab=expression(sigma[2]-sigma[1]),ylab="")
    }
    mcmcObj1=coda::mcmc(data=diffOfVariancesS2MinusS1)
    mcmcObj2=coda::mcmc(data=diffOfVariancesS2MinusS1SecC)
    mcmcObj3=coda::mcmc(data=diffOfVariancesS2MinusS1ThdC)
    mcmcObj4=coda::mcmc(data=diffOfVariancesS2MinusS1FrtC)
    mcmcListObj=coda::mcmc.list(mcmcObj1,mcmcObj2,mcmcObj3,mcmcObj4)
    coda::gelman.plot(mcmcListObj, auto.layout = FALSE)
    stats::acf(diffOfVariancesS2MinusS1)
    
    
    dev.new()
    plotpar<-par(mfrow=c(2,2),mai = c(1,0.5,0.5,0.5), oma=c(1,1,1,1))
    # Reset par only after first change of par
    on.exit(par(plotpar))
    if(sd=="var"){
      hist(diffOfVariancesS3MinusS2,freq=FALSE,main="Difference of variances Group3-Group2",col="cornflowerblue",border="white",
           xlab=expression(sigma[3]^2-sigma[2]^2))
      plot(diffOfVariancesS3MinusS2,ty="l",col="cornflowerblue",xlab=expression(sigma[3]^2-sigma[2]^2),ylab="")
    }
    if(sd=="sd"){
      hist(diffOfVariancesS3MinusS2,freq=FALSE,main="Difference of standard deviations Group3-Group2",col="cornflowerblue",border="white",
           xlab=expression(sigma[3]-sigma[2]))
      plot(diffOfVariancesS3MinusS2,ty="l",col="cornflowerblue",xlab=expression(sigma[3]-sigma[2]),ylab="")
    }
    mcmcObj1=coda::mcmc(data=diffOfVariancesS3MinusS2)
    mcmcObj2=coda::mcmc(data=diffOfVariancesS3MinusS2SecC)
    mcmcObj3=coda::mcmc(data=diffOfVariancesS3MinusS2ThdC)
    mcmcObj4=coda::mcmc(data=diffOfVariancesS3MinusS2FrtC)
    mcmcListObj=coda::mcmc.list(mcmcObj1,mcmcObj2,mcmcObj3,mcmcObj4)
    coda::gelman.plot(mcmcListObj, auto.layout = FALSE)
    stats::acf(diffOfVariancesS3MinusS2)
    
    
    dev.new()
    plotpar<-par(mfrow=c(2,2),mai = c(1,0.5,0.5,0.5), oma=c(1,1,1,1))
    # Reset par only after first change of par
    on.exit(par(plotpar))
    if(sd=="var"){
      hist(diffOfVariancesS3MinusS1,freq=FALSE,main="Difference of variances Group3-Group1",col="cornflowerblue",border="white",
           xlab=expression(sigma[3]^2-sigma[1]^2))
      plot(diffOfVariancesS3MinusS1,ty="l",col="cornflowerblue",xlab=expression(sigma[3]^2-sigma[1]^2),ylab="")
    }
    if(sd=="sd"){
      hist(diffOfVariancesS3MinusS1,freq=FALSE,main="Difference of standard deviations Group3-Group1",col="cornflowerblue",border="white",
           xlab=expression(sigma[3]-sigma[1]))
      plot(diffOfVariancesS3MinusS1,ty="l",col="cornflowerblue",xlab=expression(sigma[3]-sigma[1]),ylab="")
    }
    mcmcObj1=coda::mcmc(data=diffOfVariancesS3MinusS1)
    mcmcObj2=coda::mcmc(data=diffOfVariancesS3MinusS1SecC)
    mcmcObj3=coda::mcmc(data=diffOfVariancesS3MinusS1ThdC)
    mcmcObj4=coda::mcmc(data=diffOfVariancesS3MinusS1FrtC)
    mcmcListObj=coda::mcmc.list(mcmcObj1,mcmcObj2,mcmcObj3,mcmcObj4)
    coda::gelman.plot(mcmcListObj, auto.layout = FALSE)
    stats::acf(diffOfVariancesS3MinusS1)
    
    
    if(!is.null(dataframe$mu4)){
      dev.new()
      plotpar<-par(mfrow=c(2,2),mai = c(1,0.5,0.5,0.5), oma=c(1,1,1,1))
      # Reset par only after first change of par
      on.exit(par(plotpar))
      if(sd=="var"){
        hist(diffOfVariancesS4MinusS1,freq=FALSE,main="Difference of variances Group4-Group1",col="cornflowerblue",border="white",
             xlab=expression(sigma[4]^2-sigma[1]^2))
        plot(diffOfVariancesS4MinusS1,ty="l",col="cornflowerblue",xlab=expression(sigma[4]^2-sigma[1]^2),ylab="")
      }
      if(sd=="sd"){
        hist(diffOfVariancesS4MinusS1,freq=FALSE,main="Difference of standard deviations Group4-Group1",col="cornflowerblue",border="white",
             xlab=expression(sigma[4]-sigma[1]))
        plot(diffOfVariancesS4MinusS1,ty="l",col="cornflowerblue",xlab=expression(sigma[4]-sigma[1]),ylab="")
      }
      mcmcObj1=coda::mcmc(data=diffOfVariancesS4MinusS1)
      mcmcObj2=coda::mcmc(data=diffOfVariancesS4MinusS1SecC)
      mcmcObj3=coda::mcmc(data=diffOfVariancesS4MinusS1ThdC)
      mcmcObj4=coda::mcmc(data=diffOfVariancesS4MinusS1FrtC)
      mcmcListObj=coda::mcmc.list(mcmcObj1,mcmcObj2,mcmcObj3,mcmcObj4)
      coda::gelman.plot(mcmcListObj, auto.layout = FALSE)
      stats::acf(diffOfVariancesS4MinusS1)
      
      
      dev.new()
      plotpar<-par(mfrow=c(2,2),mai = c(1,0.5,0.5,0.5), oma=c(1,1,1,1))
      # Reset par only after first change of par
      on.exit(par(plotpar))
      if(sd=="var"){
        hist(diffOfVariancesS4MinusS2,freq=FALSE,main="Difference of variances Group4-Group2",col="cornflowerblue",border="white",
             xlab=expression(sigma[4]^2-sigma[2]^2))
        plot(diffOfVariancesS4MinusS2,ty="l",col="cornflowerblue",xlab=expression(sigma[4]^2-sigma[2]^2),ylab="")
      }
      if(sd=="sd"){
        hist(diffOfVariancesS4MinusS2,freq=FALSE,main="Difference of standard deviations Group4-Group2",col="cornflowerblue",border="white",
             xlab=expression(sigma[4]-sigma[2]))
        plot(diffOfVariancesS4MinusS2,ty="l",col="cornflowerblue",xlab=expression(sigma[4]-sigma[2]),ylab="")
      }
      mcmcObj1=coda::mcmc(data=diffOfVariancesS4MinusS2)
      mcmcObj2=coda::mcmc(data=diffOfVariancesS4MinusS2SecC)
      mcmcObj3=coda::mcmc(data=diffOfVariancesS4MinusS2ThdC)
      mcmcObj4=coda::mcmc(data=diffOfVariancesS4MinusS2FrtC)
      mcmcListObj=coda::mcmc.list(mcmcObj1,mcmcObj2,mcmcObj3,mcmcObj4)
      coda::gelman.plot(mcmcListObj, auto.layout = FALSE)
      stats::acf(diffOfVariancesS4MinusS2)
      
      
      dev.new()
      plotpar<-par(mfrow=c(2,2),mai = c(1,0.5,0.5,0.5), oma=c(1,1,1,1))
      # Reset par only after first change of par
      on.exit(par(plotpar))
      if(sd=="var"){
        hist(diffOfVariancesS4MinusS3,freq=FALSE,main="Difference of variances Group4-Group3",col="cornflowerblue",border="white",
             xlab=expression(sigma[4]^2-sigma[3]^2))
        plot(diffOfVariancesS4MinusS3,ty="l",col="cornflowerblue",xlab=expression(sigma[4]^2-sigma[3]^2),ylab="")
      }
      if(sd=="sd"){
        hist(diffOfVariancesS4MinusS3,freq=FALSE,main="Difference of standard deviations Group4-Group3",col="cornflowerblue",border="white",
             xlab=expression(sigma[4]-sigma[3]))
        plot(diffOfVariancesS4MinusS3,ty="l",col="cornflowerblue",xlab=expression(sigma[4]-sigma[3]),ylab="")
      }
      mcmcObj1=coda::mcmc(data=diffOfVariancesS4MinusS3)
      mcmcObj2=coda::mcmc(data=diffOfVariancesS4MinusS3SecC)
      mcmcObj3=coda::mcmc(data=diffOfVariancesS4MinusS3ThdC)
      mcmcObj4=coda::mcmc(data=diffOfVariancesS4MinusS3FrtC)
      mcmcListObj=coda::mcmc.list(mcmcObj1,mcmcObj2,mcmcObj3,mcmcObj4)
      coda::gelman.plot(mcmcListObj, auto.layout = FALSE)
      stats::acf(diffOfVariancesS4MinusS3)
    }
    
    if(!is.null(dataframe$mu5)){
      dev.new()
      plotpar<-par(mfrow=c(2,2),mai = c(1,0.5,0.5,0.5), oma=c(1,1,1,1))
      # Reset par only after first change of par
      on.exit(par(plotpar))
      if(sd=="var"){
        hist(diffOfVariancesS5MinusS1,freq=FALSE,main="Difference of variances Group5-Group1",col="cornflowerblue",border="white",
             xlab=expression(sigma[5]^2-sigma[1]^2))
        plot(diffOfVariancesS5MinusS1,ty="l",col="cornflowerblue",xlab=expression(sigma[5]^2-sigma[1]^2),ylab="")
      }
      if(sd=="sd"){
        hist(diffOfVariancesS5MinusS1,freq=FALSE,main="Difference of standard deviations Group5-Group1",col="cornflowerblue",border="white",
             xlab=expression(sigma[5]-sigma[1]))
        plot(diffOfVariancesS5MinusS1,ty="l",col="cornflowerblue",xlab=expression(sigma[5]-sigma[1]),ylab="")
      }
      mcmcObj1=coda::mcmc(data=diffOfVariancesS5MinusS1)
      mcmcObj2=coda::mcmc(data=diffOfVariancesS5MinusS1SecC)
      mcmcObj3=coda::mcmc(data=diffOfVariancesS5MinusS1ThdC)
      mcmcObj4=coda::mcmc(data=diffOfVariancesS5MinusS1FrtC)
      mcmcListObj=coda::mcmc.list(mcmcObj1,mcmcObj2,mcmcObj3,mcmcObj4)
      coda::gelman.plot(mcmcListObj, auto.layout = FALSE)
      stats::acf(diffOfVariancesS5MinusS1)
      
      
      
      dev.new()
      plotpar<-par(mfrow=c(2,2),mai = c(1,0.5,0.5,0.5), oma=c(1,1,1,1))
      # Reset par only after first change of par
      on.exit(par(plotpar))
      if(sd=="var"){
        hist(diffOfVariancesS5MinusS2,freq=FALSE,main="Difference of variances Group5-Group2",col="cornflowerblue",border="white",
             xlab=expression(sigma[5]^2-sigma[2]^2))
        plot(diffOfVariancesS5MinusS2,ty="l",col="cornflowerblue",xlab=expression(sigma[5]^2-sigma[2]^2),ylab="")
      }
      if(sd=="sd"){
        hist(diffOfVariancesS5MinusS2,freq=FALSE,main="Difference of standard deviations Group5-Group2",col="cornflowerblue",border="white",
             xlab=expression(sigma[5]-sigma[2]))
        plot(diffOfVariancesS5MinusS2,ty="l",col="cornflowerblue",xlab=expression(sigma[5]-sigma[2]),ylab="")
      }
      mcmcObj1=coda::mcmc(data=diffOfVariancesS5MinusS2)
      mcmcObj2=coda::mcmc(data=diffOfVariancesS5MinusS2SecC)
      mcmcObj3=coda::mcmc(data=diffOfVariancesS5MinusS2ThdC)
      mcmcObj4=coda::mcmc(data=diffOfVariancesS5MinusS2FrtC)
      mcmcListObj=coda::mcmc.list(mcmcObj1,mcmcObj2,mcmcObj3,mcmcObj4)
      coda::gelman.plot(mcmcListObj, auto.layout = FALSE)
      stats::acf(diffOfVariancesS5MinusS2)
      
      
      
      dev.new()
      plotpar<-par(mfrow=c(2,2),mai = c(1,0.5,0.5,0.5), oma=c(1,1,1,1))
      # Reset par only after first change of par
      on.exit(par(plotpar))
      if(sd=="var"){
        hist(diffOfVariancesS5MinusS3,freq=FALSE,main="Difference of variances Group5-Group3",col="cornflowerblue",border="white",
             xlab=expression(sigma[5]^2-sigma[3]^2))
        plot(diffOfVariancesS5MinusS3,ty="l",col="cornflowerblue",xlab=expression(sigma[5]^2-sigma[3]^2),ylab="")
      }
      if(sd=="sd"){
        hist(diffOfVariancesS5MinusS3,freq=FALSE,main="Difference of standard deviations Group5-Group3",col="cornflowerblue",border="white",
             xlab=expression(sigma[5]-sigma[3]))
        plot(diffOfVariancesS5MinusS3,ty="l",col="cornflowerblue",xlab=expression(sigma[5]-sigma[3]),ylab="")
      }
      mcmcObj1=coda::mcmc(data=diffOfVariancesS5MinusS3)
      mcmcObj2=coda::mcmc(data=diffOfVariancesS5MinusS3SecC)
      mcmcObj3=coda::mcmc(data=diffOfVariancesS5MinusS3ThdC)
      mcmcObj4=coda::mcmc(data=diffOfVariancesS5MinusS3FrtC)
      mcmcListObj=coda::mcmc.list(mcmcObj1,mcmcObj2,mcmcObj3,mcmcObj4)
      coda::gelman.plot(mcmcListObj, auto.layout = FALSE)
      stats::acf(diffOfVariancesS5MinusS3)
      
      
      
      dev.new()
      plotpar<-par(mfrow=c(2,2),mai = c(1,0.5,0.5,0.5), oma=c(1,1,1,1))
      # Reset par only after first change of par
      on.exit(par(plotpar))
      if(sd=="var"){
        hist(diffOfVariancesS5MinusS4,freq=FALSE,main="Difference of variances Group5-Group4",col="cornflowerblue",border="white",
             xlab=expression(sigma[5]^2-sigma[4]^2))
        plot(diffOfVariancesS5MinusS4,ty="l",col="cornflowerblue",xlab=expression(sigma[5]^2-sigma[4]^2),ylab="")
      }
      if(sd=="sd"){
        hist(diffOfVariancesS5MinusS4,freq=FALSE,main="Difference of standard deviations Group5-Group4",col="cornflowerblue",border="white",
             xlab=expression(sigma[5]-sigma[4]))
        plot(diffOfVariancesS5MinusS4,ty="l",col="cornflowerblue",xlab=expression(sigma[5]-sigma[4]),ylab="")
      }
      mcmcObj1=coda::mcmc(data=diffOfVariancesS5MinusS4)
      mcmcObj2=coda::mcmc(data=diffOfVariancesS5MinusS4SecC)
      mcmcObj3=coda::mcmc(data=diffOfVariancesS5MinusS4ThdC)
      mcmcObj4=coda::mcmc(data=diffOfVariancesS5MinusS4FrtC)
      mcmcListObj=coda::mcmc.list(mcmcObj1,mcmcObj2,mcmcObj3,mcmcObj4)
      coda::gelman.plot(mcmcListObj, auto.layout = FALSE)
      stats::acf(diffOfVariancesS5MinusS4)
    }
    
    if(!is.null(dataframe$mu6)){
      dev.new()
      plotpar<-par(mfrow=c(2,2),mai = c(1,0.5,0.5,0.5), oma=c(1,1,1,1))
      # Reset par only after first change of par
      on.exit(par(plotpar))
      if(sd=="var"){
        hist(diffOfVariancesS6MinusS1,freq=FALSE,main="Difference of variances Group6-Group1",col="cornflowerblue",border="white",
             xlab=expression(sigma[5]^2-sigma[1]^2))
        plot(diffOfVariancesS6MinusS1,ty="l",col="cornflowerblue",xlab=expression(sigma[6]^2-sigma[1]^2),ylab="")
      }
      if(sd=="sd"){
        hist(diffOfVariancesS6MinusS1,freq=FALSE,main="Difference of standard deviations Group6-Group1",col="cornflowerblue",border="white",
             xlab=expression(sigma[6]-sigma[1]))
        plot(diffOfVariancesS6MinusS1,ty="l",col="cornflowerblue",xlab=expression(sigma[6]-sigma[1]),ylab="")
      }
      mcmcObj1=coda::mcmc(data=diffOfVariancesS6MinusS1)
      mcmcObj2=coda::mcmc(data=diffOfVariancesS6MinusS1SecC)
      mcmcObj3=coda::mcmc(data=diffOfVariancesS6MinusS1ThdC)
      mcmcObj4=coda::mcmc(data=diffOfVariancesS6MinusS1FrtC)
      mcmcListObj=coda::mcmc.list(mcmcObj1,mcmcObj2,mcmcObj3,mcmcObj4)
      coda::gelman.plot(mcmcListObj, auto.layout = FALSE)
      stats::acf(diffOfVariancesS6MinusS1)
      
      
      dev.new()
      plotpar<-par(mfrow=c(2,2),mai = c(1,0.5,0.5,0.5), oma=c(1,1,1,1))
      # Reset par only after first change of par
      on.exit(par(plotpar))
      if(sd=="var"){
        hist(diffOfVariancesS6MinusS2,freq=FALSE,main="Difference of variances Group6-Group2",col="cornflowerblue",border="white",
             xlab=expression(sigma[5]^2-sigma[2]^2))
        plot(diffOfVariancesS6MinusS2,ty="l",col="cornflowerblue",xlab=expression(sigma[6]^2-sigma[2]^2),ylab="")
      }
      if(sd=="sd"){
        hist(diffOfVariancesS6MinusS2,freq=FALSE,main="Difference of standard deviations Group6-Group2",col="cornflowerblue",border="white",
             xlab=expression(sigma[6]-sigma[2]))
        plot(diffOfVariancesS6MinusS2,ty="l",col="cornflowerblue",xlab=expression(sigma[6]-sigma[2]),ylab="")
      }
      mcmcObj1=coda::mcmc(data=diffOfVariancesS6MinusS2)
      mcmcObj2=coda::mcmc(data=diffOfVariancesS6MinusS2SecC)
      mcmcObj3=coda::mcmc(data=diffOfVariancesS6MinusS2ThdC)
      mcmcObj4=coda::mcmc(data=diffOfVariancesS6MinusS2FrtC)
      mcmcListObj=coda::mcmc.list(mcmcObj1,mcmcObj2,mcmcObj3,mcmcObj4)
      coda::gelman.plot(mcmcListObj, auto.layout = FALSE)
      stats::acf(diffOfVariancesS6MinusS2)
      
      
      dev.new()
      plotpar<-par(mfrow=c(2,2),mai = c(1,0.5,0.5,0.5), oma=c(1,1,1,1))
      # Reset par only after first change of par
      on.exit(par(plotpar))
      if(sd=="var"){
        hist(diffOfVariancesS6MinusS3,freq=FALSE,main="Difference of variances Group6-Group3",col="cornflowerblue",border="white",
             xlab=expression(sigma[5]^2-sigma[3]^2))
        plot(diffOfVariancesS6MinusS3,ty="l",col="cornflowerblue",xlab=expression(sigma[6]^2-sigma[3]^2),ylab="")
      }
      if(sd=="sd"){
        hist(diffOfVariancesS6MinusS3,freq=FALSE,main="Difference of standard deviations Group6-Group3",col="cornflowerblue",border="white",
             xlab=expression(sigma[6]-sigma[3]))
        plot(diffOfVariancesS6MinusS3,ty="l",col="cornflowerblue",xlab=expression(sigma[6]-sigma[3]),ylab="")
      }
      mcmcObj1=coda::mcmc(data=diffOfVariancesS6MinusS3)
      mcmcObj2=coda::mcmc(data=diffOfVariancesS6MinusS3SecC)
      mcmcObj3=coda::mcmc(data=diffOfVariancesS6MinusS3ThdC)
      mcmcObj4=coda::mcmc(data=diffOfVariancesS6MinusS3FrtC)
      mcmcListObj=coda::mcmc.list(mcmcObj1,mcmcObj2,mcmcObj3,mcmcObj4)
      coda::gelman.plot(mcmcListObj, auto.layout = FALSE)
      stats::acf(diffOfVariancesS6MinusS3)
      
      
      dev.new()
      plotpar<-par(mfrow=c(2,2),mai = c(1,0.5,0.5,0.5), oma=c(1,1,1,1))
      # Reset par only after first change of par
      on.exit(par(plotpar))
      if(sd=="var"){
        hist(diffOfVariancesS6MinusS4,freq=FALSE,main="Difference of variances Group6-Group4",col="cornflowerblue",border="white",
             xlab=expression(sigma[5]^2-sigma[4]^2))
        plot(diffOfVariancesS6MinusS4,ty="l",col="cornflowerblue",xlab=expression(sigma[6]^2-sigma[4]^2),ylab="")
      }
      if(sd=="sd"){
        hist(diffOfVariancesS6MinusS4,freq=FALSE,main="Difference of standard deviations Group6-Group4",col="cornflowerblue",border="white",
             xlab=expression(sigma[6]-sigma[4]))
        plot(diffOfVariancesS6MinusS4,ty="l",col="cornflowerblue",xlab=expression(sigma[6]-sigma[4]),ylab="")
      }
      mcmcObj1=coda::mcmc(data=diffOfVariancesS6MinusS4)
      mcmcObj2=coda::mcmc(data=diffOfVariancesS6MinusS4SecC)
      mcmcObj3=coda::mcmc(data=diffOfVariancesS6MinusS4ThdC)
      mcmcObj4=coda::mcmc(data=diffOfVariancesS6MinusS4FrtC)
      mcmcListObj=coda::mcmc.list(mcmcObj1,mcmcObj2,mcmcObj3,mcmcObj4)
      coda::gelman.plot(mcmcListObj, auto.layout = FALSE)
      stats::acf(diffOfVariancesS6MinusS4)
      
      
      
      dev.new()
      plotpar<-par(mfrow=c(2,2),mai = c(1,0.5,0.5,0.5), oma=c(1,1,1,1))
      # Reset par only after first change of par
      on.exit(par(plotpar))
      if(sd=="var"){
        hist(diffOfVariancesS6MinusS5,freq=FALSE,main="Difference of variances Group6-Group5",col="cornflowerblue",border="white",
             xlab=expression(sigma[5]^2-sigma[5]^2))
        plot(diffOfVariancesS6MinusS5,ty="l",col="cornflowerblue",xlab=expression(sigma[6]^2-sigma[5]^2),ylab="")
      }
      if(sd=="sd"){
        hist(diffOfVariancesS6MinusS5,freq=FALSE,main="Difference of standard deviations Group6-Group5",col="cornflowerblue",border="white",
             xlab=expression(sigma[6]-sigma[5]))
        plot(diffOfVariancesS6MinusS5,ty="l",col="cornflowerblue",xlab=expression(sigma[6]-sigma[5]),ylab="")
      }
      mcmcObj1=coda::mcmc(data=diffOfVariancesS6MinusS5)
      mcmcObj2=coda::mcmc(data=diffOfVariancesS6MinusS5SecC)
      mcmcObj3=coda::mcmc(data=diffOfVariancesS6MinusS5ThdC)
      mcmcObj4=coda::mcmc(data=diffOfVariancesS6MinusS5FrtC)
      mcmcListObj=coda::mcmc.list(mcmcObj1,mcmcObj2,mcmcObj3,mcmcObj4)
      coda::gelman.plot(mcmcListObj, auto.layout = FALSE)
      stats::acf(diffOfVariancesS6MinusS5)
    }
  }

  
  # posterior plot of effect size
  if(type=="effect"){
    dev.new()
    plotpar<-par(mfrow=c(2,2),mai = c(1,0.5,0.5,0.5), oma=c(1,1,1,1))
    # Reset par only after first change of par
    on.exit(par(plotpar))
    hist(effectSize12,freq=FALSE,main="Effect size Group1-Group2",col="cornflowerblue",border="white",xlab=expression(delta[12]))
    plot(effectSize12,ty="l",col="cornflowerblue",xlab=expression(delta[12]),ylab="")
    mcmcObj1=coda::mcmc(data=effectSize12)
    mcmcObj2=coda::mcmc(data=effectSize12SecC)
    mcmcObj3=coda::mcmc(data=effectSize12ThdC)
    mcmcObj4=coda::mcmc(data=effectSize12FrtC)
    mcmcListObj=coda::mcmc.list(mcmcObj1,mcmcObj2,mcmcObj3,mcmcObj4)
    coda::gelman.plot(mcmcListObj, auto.layout = FALSE)
    stats::acf(effectSize12)
    
    
    dev.new()
    plotpar<-par(mfrow=c(2,2),mai = c(1,0.5,0.5,0.5), oma=c(1,1,1,1))
    # Reset par only after first change of par
    on.exit(par(plotpar))
    hist(effectSize23,freq=FALSE,main="Effect size Group2-Group3",col="cornflowerblue",border="white",xlab=expression(delta[23]))
    plot(effectSize23,ty="l",col="cornflowerblue",xlab=expression(delta[23]),ylab="")
    mcmcObj1=coda::mcmc(data=effectSize23)
    mcmcObj2=coda::mcmc(data=effectSize23SecC)
    mcmcObj3=coda::mcmc(data=effectSize23ThdC)
    mcmcObj4=coda::mcmc(data=effectSize23FrtC)
    mcmcListObj=coda::mcmc.list(mcmcObj1,mcmcObj2,mcmcObj3,mcmcObj4)
    coda::gelman.plot(mcmcListObj, auto.layout = FALSE)
    stats::acf(effectSize23)
    
    
    dev.new()
    plotpar<-par(mfrow=c(2,2),mai = c(1,0.5,0.5,0.5), oma=c(1,1,1,1))
    # Reset par only after first change of par
    on.exit(par(plotpar))
    hist(effectSize13,freq=FALSE,main="Effect size Group1-Group3",col="cornflowerblue",border="white",xlab=expression(delta[13]))
    plot(effectSize13,ty="l",col="cornflowerblue",xlab=expression(delta[13]),ylab="")
    mcmcObj1=coda::mcmc(data=effectSize13)
    mcmcObj2=coda::mcmc(data=effectSize13SecC)
    mcmcObj3=coda::mcmc(data=effectSize13ThdC)
    mcmcObj4=coda::mcmc(data=effectSize13FrtC)
    mcmcListObj=coda::mcmc.list(mcmcObj1,mcmcObj2,mcmcObj3,mcmcObj4)
    coda::gelman.plot(mcmcListObj, auto.layout = FALSE)
    stats::acf(effectSize13)
    
    
    if(!is.null(dataframe$mu4)){
      dev.new()
      plotpar<-par(mfrow=c(2,2),mai = c(1,0.5,0.5,0.5), oma=c(1,1,1,1))
      # Reset par only after first change of par
      on.exit(par(plotpar))
      hist(effectSize14,freq=FALSE,main="Effect size Group1-Group4",col="cornflowerblue",border="white",xlab=expression(delta[14]))
      plot(effectSize14,ty="l",col="cornflowerblue",xlab=expression(delta[14]),ylab="")
      mcmcObj1=coda::mcmc(data=effectSize14)
      mcmcObj2=coda::mcmc(data=effectSize14SecC)
      mcmcObj3=coda::mcmc(data=effectSize14ThdC)
      mcmcObj4=coda::mcmc(data=effectSize14FrtC)
      mcmcListObj=coda::mcmc.list(mcmcObj1,mcmcObj2,mcmcObj3,mcmcObj4)
      coda::gelman.plot(mcmcListObj, auto.layout = FALSE)
      stats::acf(effectSize14)
      
      
      dev.new()
      plotpar<-par(mfrow=c(2,2),mai = c(1,0.5,0.5,0.5), oma=c(1,1,1,1))
      # Reset par only after first change of par
      on.exit(par(plotpar))
      hist(effectSize24,freq=FALSE,main="Effect size Group2-Group4",col="cornflowerblue",border="white",xlab=expression(delta[24]))
      plot(effectSize24,ty="l",col="cornflowerblue",xlab=expression(delta[24]),ylab="")
      mcmcObj1=coda::mcmc(data=effectSize24)
      mcmcObj2=coda::mcmc(data=effectSize24SecC)
      mcmcObj3=coda::mcmc(data=effectSize24ThdC)
      mcmcObj4=coda::mcmc(data=effectSize24FrtC)
      mcmcListObj=coda::mcmc.list(mcmcObj1,mcmcObj2,mcmcObj3,mcmcObj4)
      coda::gelman.plot(mcmcListObj, auto.layout = FALSE)
      stats::acf(effectSize24)
      
      
      dev.new()
      plotpar<-par(mfrow=c(2,2),mai = c(1,0.5,0.5,0.5), oma=c(1,1,1,1))
      # Reset par only after first change of par
      on.exit(par(plotpar))
      hist(effectSize34,freq=FALSE,main="Effect size Group3-Group4",col="cornflowerblue",border="white",xlab=expression(delta[34]))
      plot(effectSize34,ty="l",col="cornflowerblue",xlab=expression(delta[34]),ylab="")
      mcmcObj1=coda::mcmc(data=effectSize34)
      mcmcObj2=coda::mcmc(data=effectSize34SecC)
      mcmcObj3=coda::mcmc(data=effectSize34ThdC)
      mcmcObj4=coda::mcmc(data=effectSize34FrtC)
      mcmcListObj=coda::mcmc.list(mcmcObj1,mcmcObj2,mcmcObj3,mcmcObj4)
      coda::gelman.plot(mcmcListObj, auto.layout = FALSE)
      stats::acf(effectSize34)
    }
    
    if(!is.null(dataframe$mu5)){
      dev.new()
      plotpar<-par(mfrow=c(2,2),mai = c(1,0.5,0.5,0.5), oma=c(1,1,1,1))
      # Reset par only after first change of par
      on.exit(par(plotpar))
      hist(effectSize15,freq=FALSE,main="Effect size Group1-Group5",col="cornflowerblue",border="white",xlab=expression(delta[15]))
      plot(effectSize15,ty="l",col="cornflowerblue",xlab=expression(delta[15]),ylab="")
      mcmcObj1=coda::mcmc(data=effectSize15)
      mcmcObj2=coda::mcmc(data=effectSize15SecC)
      mcmcObj3=coda::mcmc(data=effectSize15ThdC)
      mcmcObj4=coda::mcmc(data=effectSize15FrtC)
      mcmcListObj=coda::mcmc.list(mcmcObj1,mcmcObj2,mcmcObj3,mcmcObj4)
      coda::gelman.plot(mcmcListObj, auto.layout = FALSE)
      stats::acf(effectSize15)
      
      
      dev.new()
      plotpar<-par(mfrow=c(2,2),mai = c(1,0.5,0.5,0.5), oma=c(1,1,1,1))
      # Reset par only after first change of par
      on.exit(par(plotpar))
      hist(effectSize25,freq=FALSE,main="Effect size Group2-Group5",col="cornflowerblue",border="white",xlab=expression(delta[25]))
      plot(effectSize25,ty="l",col="cornflowerblue",xlab=expression(delta[25]),ylab="")
      mcmcObj1=coda::mcmc(data=effectSize25)
      mcmcObj2=coda::mcmc(data=effectSize25SecC)
      mcmcObj3=coda::mcmc(data=effectSize25ThdC)
      mcmcObj4=coda::mcmc(data=effectSize25FrtC)
      mcmcListObj=coda::mcmc.list(mcmcObj1,mcmcObj2,mcmcObj3,mcmcObj4)
      coda::gelman.plot(mcmcListObj, auto.layout = FALSE)
      stats::acf(effectSize25)
      
      
      dev.new()
      plotpar<-par(mfrow=c(2,2),mai = c(1,0.5,0.5,0.5), oma=c(1,1,1,1))
      # Reset par only after first change of par
      on.exit(par(plotpar))
      hist(effectSize35,freq=FALSE,main="Effect size Group3-Group5",col="cornflowerblue",border="white",xlab=expression(delta[35]))
      plot(effectSize35,ty="l",col="cornflowerblue",xlab=expression(delta[35]),ylab="")
      mcmcObj1=coda::mcmc(data=effectSize35)
      mcmcObj2=coda::mcmc(data=effectSize35SecC)
      mcmcObj3=coda::mcmc(data=effectSize35ThdC)
      mcmcObj4=coda::mcmc(data=effectSize35FrtC)
      mcmcListObj=coda::mcmc.list(mcmcObj1,mcmcObj2,mcmcObj3,mcmcObj4)
      coda::gelman.plot(mcmcListObj, auto.layout = FALSE)
      stats::acf(effectSize35)
      
      
      dev.new()
      plotpar<-par(mfrow=c(2,2),mai = c(1,0.5,0.5,0.5), oma=c(1,1,1,1))
      # Reset par only after first change of par
      on.exit(par(plotpar))
      hist(effectSize45,freq=FALSE,main="Effect size Group4-Group5",col="cornflowerblue",border="white",xlab=expression(delta[45]))
      plot(effectSize45,ty="l",col="cornflowerblue",xlab=expression(delta[45]),ylab="")
      mcmcObj1=coda::mcmc(data=effectSize45)
      mcmcObj2=coda::mcmc(data=effectSize45SecC)
      mcmcObj3=coda::mcmc(data=effectSize45ThdC)
      mcmcObj4=coda::mcmc(data=effectSize45FrtC)
      mcmcListObj=coda::mcmc.list(mcmcObj1,mcmcObj2,mcmcObj3,mcmcObj4)
      coda::gelman.plot(mcmcListObj, auto.layout = FALSE)
      stats::acf(effectSize45)
    }
    
    if(!is.null(dataframe$mu6)){
      dev.new()
      plotpar<-par(mfrow=c(2,2),mai = c(1,0.5,0.5,0.5), oma=c(1,1,1,1))
      # Reset par only after first change of par
      on.exit(par(plotpar))
      hist(effectSize16,freq=FALSE,main="Effect size Group1-Group6",col="cornflowerblue",border="white",xlab=expression(delta[16]))
      plot(effectSize16,ty="l",col="cornflowerblue",xlab=expression(delta[16]),ylab="")
      mcmcObj1=coda::mcmc(data=effectSize16)
      mcmcObj2=coda::mcmc(data=effectSize16SecC)
      mcmcObj3=coda::mcmc(data=effectSize16ThdC)
      mcmcObj4=coda::mcmc(data=effectSize16FrtC)
      mcmcListObj=coda::mcmc.list(mcmcObj1,mcmcObj2,mcmcObj3,mcmcObj4)
      coda::gelman.plot(mcmcListObj, auto.layout = FALSE)
      stats::acf(effectSize16)
      
      
      dev.new()
      plotpar<-par(mfrow=c(2,2),mai = c(1,0.5,0.5,0.5), oma=c(1,1,1,1))
      # Reset par only after first change of par
      on.exit(par(plotpar))
      hist(effectSize26,freq=FALSE,main="Effect size Group2-Group6",col="cornflowerblue",border="white",xlab=expression(delta[26]))
      plot(effectSize26,ty="l",col="cornflowerblue",xlab=expression(delta[26]),ylab="")
      mcmcObj1=coda::mcmc(data=effectSize26)
      mcmcObj2=coda::mcmc(data=effectSize26SecC)
      mcmcObj3=coda::mcmc(data=effectSize26ThdC)
      mcmcObj4=coda::mcmc(data=effectSize26FrtC)
      mcmcListObj=coda::mcmc.list(mcmcObj1,mcmcObj2,mcmcObj3,mcmcObj4)
      coda::gelman.plot(mcmcListObj, auto.layout = FALSE)
      stats::acf(effectSize26)
      
      
      dev.new()
      plotpar<-par(mfrow=c(2,2),mai = c(1,0.5,0.5,0.5), oma=c(1,1,1,1))
      # Reset par only after first change of par
      on.exit(par(plotpar))
      hist(effectSize36,freq=FALSE,main="Effect size Group3-Group6",col="cornflowerblue",border="white",xlab=expression(delta[36]))
      plot(effectSize36,ty="l",col="cornflowerblue",xlab=expression(delta[36]),ylab="")
      mcmcObj1=coda::mcmc(data=effectSize36)
      mcmcObj2=coda::mcmc(data=effectSize36SecC)
      mcmcObj3=coda::mcmc(data=effectSize36ThdC)
      mcmcObj4=coda::mcmc(data=effectSize36FrtC)
      mcmcListObj=coda::mcmc.list(mcmcObj1,mcmcObj2,mcmcObj3,mcmcObj4)
      coda::gelman.plot(mcmcListObj, auto.layout = FALSE)
      stats::acf(effectSize36)
      
      
      dev.new()
      plotpar<-par(mfrow=c(2,2),mai = c(1,0.5,0.5,0.5), oma=c(1,1,1,1))
      # Reset par only after first change of par
      on.exit(par(plotpar))
      hist(effectSize46,freq=FALSE,main="Effect size Group4-Group6",col="cornflowerblue",border="white",xlab=expression(delta[46]))
      plot(effectSize46,ty="l",col="cornflowerblue",xlab=expression(delta[46]),ylab="")
      mcmcObj1=coda::mcmc(data=effectSize46)
      mcmcObj2=coda::mcmc(data=effectSize46SecC)
      mcmcObj3=coda::mcmc(data=effectSize46ThdC)
      mcmcObj4=coda::mcmc(data=effectSize46FrtC)
      mcmcListObj=coda::mcmc.list(mcmcObj1,mcmcObj2,mcmcObj3,mcmcObj4)
      coda::gelman.plot(mcmcListObj, auto.layout = FALSE)
      stats::acf(effectSize46)
      
      

      dev.new()
      plotpar<-par(mfrow=c(2,2),mai = c(1,0.5,0.5,0.5), oma=c(1,1,1,1))
      # Reset par only after first change of par
      on.exit(par(plotpar))
      hist(effectSize56,freq=FALSE,main="Effect size Group5-Group6",col="cornflowerblue",border="white",xlab=expression(delta[56]))
      plot(effectSize56,ty="l",col="cornflowerblue",xlab=expression(delta[56]),ylab="")
      mcmcObj1=coda::mcmc(data=effectSize56)
      mcmcObj2=coda::mcmc(data=effectSize56SecC)
      mcmcObj3=coda::mcmc(data=effectSize56ThdC)
      mcmcObj4=coda::mcmc(data=effectSize56FrtC)
      mcmcListObj=coda::mcmc.list(mcmcObj1,mcmcObj2,mcmcObj3,mcmcObj4)
      coda::gelman.plot(mcmcListObj, auto.layout = FALSE)
      stats::acf(effectSize56)
    }
  }
  
  # ROPE posterior analysis plot
  if(type=="rope"){
    
    # utility function for computing the mode
    getmode <- function(v) {
      uniqv <- unique(v)
      uniqv[which.max(tabulate(match(v, uniqv)))]
    }

    
    ################# POSTERIOR ROPE ANALYSIS FOR DELTA_12 ######################
    # Set 2x3 par
    dev.new()
    newpar <- par(mfrow=c(2,3),mai = c(1,0.5,0.5,0.5), oma=c(1,1,1,1))
    # Reset par only after first change of par
    on.exit(par(newpar))
    # Effect size posterior histogram
    # Make a vector of values to draw ticks at:
    hist(effectSize12,freq=FALSE,col="cornflowerblue",border="white",xlab=expression(delta[12]),main=mtext("Posterior distribution of Effect size "~delta[12], side=3,line=2.2,col="black",cex=0.8,font=2),xaxt="n")
    ticks <- seq(from=-10, to=10, by=0.05)
    # And draw the axis:
    axis(1, at=ticks)
    
    # Posterior mode
    #abline(v=getmode(effectSize),col="blue",lwd=2,lty="solid")
    lines(x=c(getmode(effectSize12),getmode(effectSize12)),y=c(-0.025,0.025),type="l",lwd=2.5,col="blue")
    
    # Posterior CI
    postCI=quantile(effectSize12,probs=c((1-ci)/2,ci+(1-ci)/2))
    lines(x=c(postCI[1],postCI[2]),y=c(0,0),type="l",lwd=2,col="black")
    lines(x=c(postCI[1],postCI[1]),y=c(-0.02,0.02),type="l",lwd=2,col="black")
    lines(x=c(postCI[2],postCI[2]),y=c(-0.02,0.02),type="l",lwd=2,col="black")
    
    # Visualization of typical ROPEs
    # ROPE for no effect
    lines(x=c(-0.2,0.2),y=c(0.04,0.04),type="l",lwd=2,col="red")
    # ROPE for small effect
    lines(x=c(-0.5,-0.2),y=c(0.04,0.04),type="l",lwd=2,col="orange")
    lines(x=c(0.2,0.5),y=c(0.04,0.04),type="l",lwd=2,col="orange")
    
    # ROPE for medium effect
    lines(x=c(-0.8,-0.5),y=c(0.04,0.04),type="l",lwd=2,col="green")
    lines(x=c(0.5,0.8),y=c(0.04,0.04),type="l",lwd=2,col="green")
    # ROPE for large effect
    lines(x=c(-2,-0.8),y=c(0.04,0.04),type="l",lwd=2,col="purple")
    lines(x=c(0.8,2),y=c(0.04,0.04),type="l",lwd=2,col="purple")
    
    abline(v=0.8,col="purple",lwd=0.5,lty="dashed")
    abline(v=-0.8,col="purple",lwd=0.5,lty="dashed")
    abline(v=0.5,col="green",lwd=0.5,lty="dashed")
    abline(v=-0.5,col="green",lwd=0.5,lty="dashed")
    abline(v=0.2,col="orange",lwd=0.5,lty="dashed")
    abline(v=-0.2,col="orange",lwd=0.5,lty="dashed")
    
    # Positions of text above the plot
    #posLeft=quantile(effectSize,probs=c(0.1,0.5,0.95))[1]
    #posMid=quantile(effectSize,probs=c(0.1,0.5,0.95))[2]
    #posRight=quantile(effectSize,probs=c(0.1,0.5,0.95))[3]
    
    # Text above the plot
    mtext(paste0(ci*100,"% ","CI: [", round(postCI[1], 3),",",round(postCI[2], 3),"]"), side=3,line=0.15,at=min(effectSize12),col="black",cex=0.8,font=2)
    mtext(paste0("Mean: ", round(mean(effectSize12),3)),side=3,line=1.15,at=min(effectSize12),col="black",cex=0.8,font=2)
    mtext(paste0("Mode: ", round(getmode(effectSize12),3)),side=3,line=1.15,at=max(effectSize12),col="black",cex=0.8,font=2)
    
    
    
    ########################################### HISTOGRAM DELTA_13 ###########################################
    # Effect size posterior histogram
    # Make a vector of values to draw ticks at:
    hist(effectSize13,freq=FALSE,col="cornflowerblue",border="white",xlab=expression(delta[13]),main=mtext("Posterior distribution of Effect size "~delta[13], side=3,line=2.2,col="black",cex=0.8,font=2),xaxt="n")
    ticks <- seq(from=-10, to=10, by=0.05)
    # And draw the axis:
    axis(1, at=ticks)
    
    # Posterior mode
    #abline(v=getmode(effectSize),col="blue",lwd=2,lty="solid")
    lines(x=c(getmode(effectSize13),getmode(effectSize13)),y=c(-0.025,0.025),type="l",lwd=2.5,col="blue")
    
    # Posterior CI
    postCI=quantile(effectSize13,probs=c((1-ci)/2,ci+(1-ci)/2))
    lines(x=c(postCI[1],postCI[2]),y=c(0,0),type="l",lwd=2,col="black")
    lines(x=c(postCI[1],postCI[1]),y=c(-0.02,0.02),type="l",lwd=2,col="black")
    lines(x=c(postCI[2],postCI[2]),y=c(-0.02,0.02),type="l",lwd=2,col="black")
    
    # Visualization of typical ROPEs
    # ROPE for no effect
    lines(x=c(-0.2,0.2),y=c(0.04,0.04),type="l",lwd=2,col="red")
    # ROPE for small effect
    lines(x=c(-0.5,-0.2),y=c(0.04,0.04),type="l",lwd=2,col="orange")
    lines(x=c(0.2,0.5),y=c(0.04,0.04),type="l",lwd=2,col="orange")
    
    # ROPE for medium effect
    lines(x=c(-0.8,-0.5),y=c(0.04,0.04),type="l",lwd=2,col="green")
    lines(x=c(0.5,0.8),y=c(0.04,0.04),type="l",lwd=2,col="green")
    # ROPE for large effect
    lines(x=c(-2,-0.8),y=c(0.04,0.04),type="l",lwd=2,col="purple")
    lines(x=c(0.8,2),y=c(0.04,0.04),type="l",lwd=2,col="purple")
    
    abline(v=0.8,col="purple",lwd=0.5,lty="dashed")
    abline(v=-0.8,col="purple",lwd=0.5,lty="dashed")
    abline(v=0.5,col="green",lwd=0.5,lty="dashed")
    abline(v=-0.5,col="green",lwd=0.5,lty="dashed")
    abline(v=0.2,col="orange",lwd=0.5,lty="dashed")
    abline(v=-0.2,col="orange",lwd=0.5,lty="dashed")
    
    # Positions of text above the plot
    #posLeft=quantile(effectSize,probs=c(0.1,0.5,0.95))[1]
    #posMid=quantile(effectSize,probs=c(0.1,0.5,0.95))[2]
    #posRight=quantile(effectSize,probs=c(0.1,0.5,0.95))[3]
    
    # Text above the plot
    mtext(paste0(ci*100,"% ","CI: [", round(postCI[1], 3),",",round(postCI[2], 3),"]"), side=3,line=0.15,at=min(effectSize13),col="black",cex=0.8,font=2)
    mtext(paste0("Mean: ", round(mean(effectSize13),3)),side=3,line=1.15,at=min(effectSize13),col="black",cex=0.8,font=2)
    mtext(paste0("Mode: ", round(getmode(effectSize13),3)),side=3,line=1.15,at=max(effectSize13),col="black",cex=0.8,font=2)
    
    
    ####################################### HISTOGRAM DELTA_23 ###############################################
    # Effect size posterior histogram
    # Make a vector of values to draw ticks at:
    hist(effectSize23,freq=FALSE,col="cornflowerblue",border="white",xlab=expression(delta[23]),main=mtext("Posterior distribution of Effect size "~delta[23], side=3,line=2.2,col="black",cex=0.8,font=2),xaxt="n")
    ticks <- seq(from=-10, to=10, by=0.05)
    # And draw the axis:
    axis(1, at=ticks)
    
    # Posterior mode
    #abline(v=getmode(effectSize),col="blue",lwd=2,lty="solid")
    lines(x=c(getmode(effectSize23),getmode(effectSize23)),y=c(-0.025,0.025),type="l",lwd=2.5,col="blue")
    
    # Posterior CI
    postCI=quantile(effectSize23,probs=c((1-ci)/2,ci+(1-ci)/2))
    lines(x=c(postCI[1],postCI[2]),y=c(0,0),type="l",lwd=2,col="black")
    lines(x=c(postCI[1],postCI[1]),y=c(-0.02,0.02),type="l",lwd=2,col="black")
    lines(x=c(postCI[2],postCI[2]),y=c(-0.02,0.02),type="l",lwd=2,col="black")
    
    # Visualization of typical ROPEs
    # ROPE for no effect
    lines(x=c(-0.2,0.2),y=c(0.04,0.04),type="l",lwd=2,col="red")
    # ROPE for small effect
    lines(x=c(-0.5,-0.2),y=c(0.04,0.04),type="l",lwd=2,col="orange")
    lines(x=c(0.2,0.5),y=c(0.04,0.04),type="l",lwd=2,col="orange")
    
    # ROPE for medium effect
    lines(x=c(-0.8,-0.5),y=c(0.04,0.04),type="l",lwd=2,col="green")
    lines(x=c(0.5,0.8),y=c(0.04,0.04),type="l",lwd=2,col="green")
    # ROPE for large effect
    lines(x=c(-2,-0.8),y=c(0.04,0.04),type="l",lwd=2,col="purple")
    lines(x=c(0.8,2),y=c(0.04,0.04),type="l",lwd=2,col="purple")
    
    abline(v=0.8,col="purple",lwd=0.5,lty="dashed")
    abline(v=-0.8,col="purple",lwd=0.5,lty="dashed")
    abline(v=0.5,col="green",lwd=0.5,lty="dashed")
    abline(v=-0.5,col="green",lwd=0.5,lty="dashed")
    abline(v=0.2,col="orange",lwd=0.5,lty="dashed")
    abline(v=-0.2,col="orange",lwd=0.5,lty="dashed")
    
    # Positions of text above the plot
    #posLeft=quantile(effectSize,probs=c(0.1,0.5,0.95))[1]
    #posMid=quantile(effectSize,probs=c(0.1,0.5,0.95))[2]
    #posRight=quantile(effectSize,probs=c(0.1,0.5,0.95))[3]
    
    # Text above the plot
    mtext(paste0(ci*100,"% ","CI: [", round(postCI[1], 3),",",round(postCI[2], 3),"]"), side=3,line=0.15,at=min(effectSize23),col="black",cex=0.8,font=2)
    mtext(paste0("Mean: ", round(mean(effectSize23),3)),side=3,line=1.15,at=min(effectSize23),col="black",cex=0.8,font=2)
    mtext(paste0("Mode: ", round(getmode(effectSize23),3)),side=3,line=1.15,at=max(effectSize23),col="black",cex=0.8,font=2)
    
    
    ######################## POSTERIOR ROPE ANALYSIS DELTA_12 #################################################
    # Bar plot
    #effSizeCIValues = effectSize[{q<-rank(effectSize)/length(effectSize);q>=((1-ci)/2) & q <= (ci+(1-ci)/2)}]
    effSizeCI=quantile(effectSize12,probs=c(((1-ci)/2),(ci+(1-ci)/2)))
    effectSizeCIValues = effectSize12[effectSize12 >effSizeCI[1] & effectSize12 < effSizeCI[2]]
    largePosEffectIterations = length(which(effectSizeCIValues >= 0.8))/length(effectSizeCIValues)
    largeNegEffectIterations = length(which(effectSizeCIValues <= -0.8))/length(effectSizeCIValues)
    mediumPosEffectIterations = length(which(effectSizeCIValues >= 0.5 & effectSizeCIValues < 0.8))/length(effectSizeCIValues)
    mediumNegEffectIterations = length(which(effectSizeCIValues <= -0.5 & effectSizeCIValues > -0.8))/length(effectSizeCIValues)
    smallPosEffectIterations = length(which(effectSizeCIValues >= 0.2 & effectSizeCIValues < 0.5))/length(effectSizeCIValues)
    smallNegEffectIterations = length(which(effectSizeCIValues <= -0.2 & effectSizeCIValues > -0.5))/length(effectSizeCIValues)
    noEffectIterations = length(which(effectSizeCIValues < 0.2 & effectSizeCIValues > -0.2))/length(effectSizeCIValues)
    
    lbls=c("large positive:", "medium positive:", "small positive:", "no effect:", "small negative:", "medium negative:","large negative:")
    pct=c(largeNegEffectIterations,mediumNegEffectIterations,smallNegEffectIterations,noEffectIterations,smallPosEffectIterations,mediumPosEffectIterations,largePosEffectIterations)
    pct=round(pct,digits=4)
    pct=pct*100
    lbls=paste(lbls,pct)
    lbls <- paste(pct,"%",sep="") # ad % to labels
    #pie(pct,labels = lbls, col=rainbow(length(lbls)),
    #main="Posterior probabilities of varying effect sizes")
    
    # Barplot
    b<-barplot(pct,
               main = mtext("Posterior ROPE-analysis for"~delta[12]),
               xlab = "ROPEs of standard effect sizes",
               ylab = "Posterior-Percentage included inside ROPE",
               names.arg = c("L-", "M-", "S-", "No Eff.", "S+", "M+", "L+"),
               col = "cornflowerblue",
               border="white",
               horiz = FALSE,
               ylim=c(0,107.5))
    b
    text(x=b,y=pct+4,labels=lbls)
    
    
    ################################## POSTERIOR ROPE ANALYSIS FOR DELTA_[13] ############################
    # Bar plot
    #effSizeCIValues = effectSize[{q<-rank(effectSize)/length(effectSize);q>=((1-ci)/2) & q <= (ci+(1-ci)/2)}]
    effSizeCI=quantile(effectSize13,probs=c(((1-ci)/2),(ci+(1-ci)/2)))
    effectSizeCIValues = effectSize13[effectSize13 >effSizeCI[1] & effectSize13 < effSizeCI[2]]
    largePosEffectIterations = length(which(effectSizeCIValues >= 0.8))/length(effectSizeCIValues)
    largeNegEffectIterations = length(which(effectSizeCIValues <= -0.8))/length(effectSizeCIValues)
    mediumPosEffectIterations = length(which(effectSizeCIValues >= 0.5 & effectSizeCIValues < 0.8))/length(effectSizeCIValues)
    mediumNegEffectIterations = length(which(effectSizeCIValues <= -0.5 & effectSizeCIValues > -0.8))/length(effectSizeCIValues)
    smallPosEffectIterations = length(which(effectSizeCIValues >= 0.2 & effectSizeCIValues < 0.5))/length(effectSizeCIValues)
    smallNegEffectIterations = length(which(effectSizeCIValues <= -0.2 & effectSizeCIValues > -0.5))/length(effectSizeCIValues)
    noEffectIterations = length(which(effectSizeCIValues < 0.2 & effectSizeCIValues > -0.2))/length(effectSizeCIValues)
    
    lbls=c("large positive:", "medium positive:", "small positive:", "no effect:", "small negative:", "medium negative:","large negative:")
    pct=c(largeNegEffectIterations,mediumNegEffectIterations,smallNegEffectIterations,noEffectIterations,smallPosEffectIterations,mediumPosEffectIterations,largePosEffectIterations)
    pct=round(pct,digits=4)
    pct=pct*100
    lbls=paste(lbls,pct)
    lbls <- paste(pct,"%",sep="") # ad % to labels
    #pie(pct,labels = lbls, col=rainbow(length(lbls)),
    #main="Posterior probabilities of varying effect sizes")
    
    # Barplot
    b<-barplot(pct,
               main = mtext("Posterior ROPE-analysis for"~delta[13]),
               xlab = "ROPEs of standard effect sizes",
               ylab = "Posterior-Percentage included inside ROPE",
               names.arg = c("L-", "M-", "S-", "No Eff.", "S+", "M+", "L+"),
               col = "cornflowerblue",
               border="white",
               horiz = FALSE,
               ylim=c(0,107.5))
    b
    text(x=b,y=pct+4,labels=lbls)
    
    
    
    
    ####################### POSTERIOR ROPE ANALYSIS FOR DELTA_23 #######################################
    # Bar plot
    #effSizeCIValues = effectSize[{q<-rank(effectSize)/length(effectSize);q>=((1-ci)/2) & q <= (ci+(1-ci)/2)}]
    effSizeCI=quantile(effectSize23,probs=c(((1-ci)/2),(ci+(1-ci)/2)))
    effectSizeCIValues = effectSize23[effectSize23 >effSizeCI[1] & effectSize23 < effSizeCI[2]]
    largePosEffectIterations = length(which(effectSizeCIValues >= 0.8))/length(effectSizeCIValues)
    largeNegEffectIterations = length(which(effectSizeCIValues <= -0.8))/length(effectSizeCIValues)
    mediumPosEffectIterations = length(which(effectSizeCIValues >= 0.5 & effectSizeCIValues < 0.8))/length(effectSizeCIValues)
    mediumNegEffectIterations = length(which(effectSizeCIValues <= -0.5 & effectSizeCIValues > -0.8))/length(effectSizeCIValues)
    smallPosEffectIterations = length(which(effectSizeCIValues >= 0.2 & effectSizeCIValues < 0.5))/length(effectSizeCIValues)
    smallNegEffectIterations = length(which(effectSizeCIValues <= -0.2 & effectSizeCIValues > -0.5))/length(effectSizeCIValues)
    noEffectIterations = length(which(effectSizeCIValues < 0.2 & effectSizeCIValues > -0.2))/length(effectSizeCIValues)
    
    lbls=c("large positive:", "medium positive:", "small positive:", "no effect:", "small negative:", "medium negative:","large negative:")
    pct=c(largeNegEffectIterations,mediumNegEffectIterations,smallNegEffectIterations,noEffectIterations,smallPosEffectIterations,mediumPosEffectIterations,largePosEffectIterations)
    pct=round(pct,digits=4)
    pct=pct*100
    lbls=paste(lbls,pct)
    lbls <- paste(pct,"%",sep="") # ad % to labels
    #pie(pct,labels = lbls, col=rainbow(length(lbls)),
    #main="Posterior probabilities of varying effect sizes")
    
    # Barplot
    b<-barplot(pct,
               main = mtext("Posterior ROPE-analysis for"~delta[23]),
               xlab = "ROPEs of standard effect sizes",
               ylab = "Posterior-Percentage included inside ROPE",
               names.arg = c("L-", "M-", "S-", "No Eff.", "S+", "M+", "L+"),
               col = "cornflowerblue",
               border="white",
               horiz = FALSE,
               ylim=c(0,107.5))
    b
    text(x=b,y=pct+4,labels=lbls)
    
    ###################################################### ADDITIONAL EFFECT SIZES FOR FOURTH GROUP #########################################################
    if(!is.null(dataframe$mu4)){
      # Set 2x3 par
      dev.new()
      newpar <- par(mfrow=c(2,3),mai = c(1,0.5,0.5,0.5), oma=c(1,1,1,1))
      # Reset par only after first change of par
      on.exit(par(newpar))
      
      ########################################### HISTOGRAM DELTA_14 ###########################################
      # Effect size posterior histogram
      # Make a vector of values to draw ticks at:
      hist(effectSize14,freq=FALSE,col="cornflowerblue",border="white",xlab=expression(delta[14]),main=mtext("Posterior distribution of Effect size "~delta[14], side=3,line=2.2,col="black",cex=0.8,font=2),xaxt="n")
      ticks <- seq(from=-10, to=10, by=0.05)
      # And draw the axis:
      axis(1, at=ticks)
      
      # Posterior mode
      #abline(v=getmode(effectSize),col="blue",lwd=2,lty="solid")
      lines(x=c(getmode(effectSize14),getmode(effectSize14)),y=c(-0.025,0.025),type="l",lwd=2.5,col="blue")
      
      # Posterior CI
      postCI=quantile(effectSize14,probs=c((1-ci)/2,ci+(1-ci)/2))
      lines(x=c(postCI[1],postCI[2]),y=c(0,0),type="l",lwd=2,col="black")
      lines(x=c(postCI[1],postCI[1]),y=c(-0.02,0.02),type="l",lwd=2,col="black")
      lines(x=c(postCI[2],postCI[2]),y=c(-0.02,0.02),type="l",lwd=2,col="black")
      
      # Visualization of typical ROPEs
      # ROPE for no effect
      lines(x=c(-0.2,0.2),y=c(0.04,0.04),type="l",lwd=2,col="red")
      # ROPE for small effect
      lines(x=c(-0.5,-0.2),y=c(0.04,0.04),type="l",lwd=2,col="orange")
      lines(x=c(0.2,0.5),y=c(0.04,0.04),type="l",lwd=2,col="orange")
      
      # ROPE for medium effect
      lines(x=c(-0.8,-0.5),y=c(0.04,0.04),type="l",lwd=2,col="green")
      lines(x=c(0.5,0.8),y=c(0.04,0.04),type="l",lwd=2,col="green")
      # ROPE for large effect
      lines(x=c(-2,-0.8),y=c(0.04,0.04),type="l",lwd=2,col="purple")
      lines(x=c(0.8,2),y=c(0.04,0.04),type="l",lwd=2,col="purple")
      
      abline(v=0.8,col="purple",lwd=0.5,lty="dashed")
      abline(v=-0.8,col="purple",lwd=0.5,lty="dashed")
      abline(v=0.5,col="green",lwd=0.5,lty="dashed")
      abline(v=-0.5,col="green",lwd=0.5,lty="dashed")
      abline(v=0.2,col="orange",lwd=0.5,lty="dashed")
      abline(v=-0.2,col="orange",lwd=0.5,lty="dashed")
      
      # Positions of text above the plot
      #posLeft=quantile(effectSize,probs=c(0.1,0.5,0.95))[1]
      #posMid=quantile(effectSize,probs=c(0.1,0.5,0.95))[2]
      #posRight=quantile(effectSize,probs=c(0.1,0.5,0.95))[3]
      
      # Text above the plot
      mtext(paste0(ci*100,"% ","CI: [", round(postCI[1], 3),",",round(postCI[2], 3),"]"), side=3,line=0.15,at=min(effectSize14),col="black",cex=0.8,font=2)
      mtext(paste0("Mean: ", round(mean(effectSize14),3)),side=3,line=1.15,at=min(effectSize14),col="black",cex=0.8,font=2)
      mtext(paste0("Mode: ", round(getmode(effectSize14),3)),side=3,line=1.15,at=max(effectSize14),col="black",cex=0.8,font=2)
      
      
      
      
      ########################################### HISTOGRAM DELTA_24 ###########################################
      # Effect size posterior histogram
      # Make a vector of values to draw ticks at:
      hist(effectSize24,freq=FALSE,col="cornflowerblue",border="white",xlab=expression(delta[24]),main=mtext("Posterior distribution of Effect size "~delta[24], side=3,line=2.2,col="black",cex=0.8,font=2),xaxt="n")
      ticks <- seq(from=-10, to=10, by=0.05)
      # And draw the axis:
      axis(1, at=ticks)
      
      # Posterior mode
      #abline(v=getmode(effectSize),col="blue",lwd=2,lty="solid")
      lines(x=c(getmode(effectSize24),getmode(effectSize24)),y=c(-0.025,0.025),type="l",lwd=2.5,col="blue")
      
      # Posterior CI
      postCI=quantile(effectSize24,probs=c((1-ci)/2,ci+(1-ci)/2))
      lines(x=c(postCI[1],postCI[2]),y=c(0,0),type="l",lwd=2,col="black")
      lines(x=c(postCI[1],postCI[1]),y=c(-0.02,0.02),type="l",lwd=2,col="black")
      lines(x=c(postCI[2],postCI[2]),y=c(-0.02,0.02),type="l",lwd=2,col="black")
      
      # Visualization of typical ROPEs
      # ROPE for no effect
      lines(x=c(-0.2,0.2),y=c(0.04,0.04),type="l",lwd=2,col="red")
      # ROPE for small effect
      lines(x=c(-0.5,-0.2),y=c(0.04,0.04),type="l",lwd=2,col="orange")
      lines(x=c(0.2,0.5),y=c(0.04,0.04),type="l",lwd=2,col="orange")
      
      # ROPE for medium effect
      lines(x=c(-0.8,-0.5),y=c(0.04,0.04),type="l",lwd=2,col="green")
      lines(x=c(0.5,0.8),y=c(0.04,0.04),type="l",lwd=2,col="green")
      # ROPE for large effect
      lines(x=c(-2,-0.8),y=c(0.04,0.04),type="l",lwd=2,col="purple")
      lines(x=c(0.8,2),y=c(0.04,0.04),type="l",lwd=2,col="purple")
      
      abline(v=0.8,col="purple",lwd=0.5,lty="dashed")
      abline(v=-0.8,col="purple",lwd=0.5,lty="dashed")
      abline(v=0.5,col="green",lwd=0.5,lty="dashed")
      abline(v=-0.5,col="green",lwd=0.5,lty="dashed")
      abline(v=0.2,col="orange",lwd=0.5,lty="dashed")
      abline(v=-0.2,col="orange",lwd=0.5,lty="dashed")
      
      # Positions of text above the plot
      #posLeft=quantile(effectSize,probs=c(0.1,0.5,0.95))[1]
      #posMid=quantile(effectSize,probs=c(0.1,0.5,0.95))[2]
      #posRight=quantile(effectSize,probs=c(0.1,0.5,0.95))[3]
      
      # Text above the plot
      mtext(paste0(ci*100,"% ","CI: [", round(postCI[1], 3),",",round(postCI[2], 3),"]"), side=3,line=0.15,at=min(effectSize24),col="black",cex=0.8,font=2)
      mtext(paste0("Mean: ", round(mean(effectSize24),3)),side=3,line=1.15,at=min(effectSize24),col="black",cex=0.8,font=2)
      mtext(paste0("Mode: ", round(getmode(effectSize24),3)),side=3,line=1.15,at=max(effectSize24),col="black",cex=0.8,font=2)
      
      
      
      ########################################### HISTOGRAM DELTA_34 ###########################################
      # Effect size posterior histogram
      # Make a vector of values to draw ticks at:
      hist(effectSize34,freq=FALSE,col="cornflowerblue",border="white",xlab=expression(delta[34]),main=mtext("Posterior distribution of Effect size "~delta[34], side=3,line=2.2,col="black",cex=0.8,font=2),xaxt="n")
      ticks <- seq(from=-10, to=10, by=0.05)
      # And draw the axis:
      axis(1, at=ticks)
      
      # Posterior mode
      #abline(v=getmode(effectSize),col="blue",lwd=2,lty="solid")
      lines(x=c(getmode(effectSize34),getmode(effectSize34)),y=c(-0.025,0.025),type="l",lwd=2.5,col="blue")
      
      # Posterior CI
      postCI=quantile(effectSize34,probs=c((1-ci)/2,ci+(1-ci)/2))
      lines(x=c(postCI[1],postCI[2]),y=c(0,0),type="l",lwd=2,col="black")
      lines(x=c(postCI[1],postCI[1]),y=c(-0.02,0.02),type="l",lwd=2,col="black")
      lines(x=c(postCI[2],postCI[2]),y=c(-0.02,0.02),type="l",lwd=2,col="black")
      
      # Visualization of typical ROPEs
      # ROPE for no effect
      lines(x=c(-0.2,0.2),y=c(0.04,0.04),type="l",lwd=2,col="red")
      # ROPE for small effect
      lines(x=c(-0.5,-0.2),y=c(0.04,0.04),type="l",lwd=2,col="orange")
      lines(x=c(0.2,0.5),y=c(0.04,0.04),type="l",lwd=2,col="orange")
      
      # ROPE for medium effect
      lines(x=c(-0.8,-0.5),y=c(0.04,0.04),type="l",lwd=2,col="green")
      lines(x=c(0.5,0.8),y=c(0.04,0.04),type="l",lwd=2,col="green")
      # ROPE for large effect
      lines(x=c(-2,-0.8),y=c(0.04,0.04),type="l",lwd=2,col="purple")
      lines(x=c(0.8,2),y=c(0.04,0.04),type="l",lwd=2,col="purple")
      
      abline(v=0.8,col="purple",lwd=0.5,lty="dashed")
      abline(v=-0.8,col="purple",lwd=0.5,lty="dashed")
      abline(v=0.5,col="green",lwd=0.5,lty="dashed")
      abline(v=-0.5,col="green",lwd=0.5,lty="dashed")
      abline(v=0.2,col="orange",lwd=0.5,lty="dashed")
      abline(v=-0.2,col="orange",lwd=0.5,lty="dashed")
      
      # Positions of text above the plot
      #posLeft=quantile(effectSize,probs=c(0.1,0.5,0.95))[1]
      #posMid=quantile(effectSize,probs=c(0.1,0.5,0.95))[2]
      #posRight=quantile(effectSize,probs=c(0.1,0.5,0.95))[3]
      
      # Text above the plot
      mtext(paste0(ci*100,"% ","CI: [", round(postCI[1], 3),",",round(postCI[2], 3),"]"), side=3,line=0.15,at=min(effectSize34),col="black",cex=0.8,font=2)
      mtext(paste0("Mean: ", round(mean(effectSize34),3)),side=3,line=1.15,at=min(effectSize34),col="black",cex=0.8,font=2)
      mtext(paste0("Mode: ", round(getmode(effectSize34),3)),side=3,line=1.15,at=max(effectSize34),col="black",cex=0.8,font=2)
      
      
      ####################### POSTERIOR ROPE ANALYSIS FOR DELTA_14 #######################################
      # Bar plot
      #effSizeCIValues = effectSize[{q<-rank(effectSize)/length(effectSize);q>=((1-ci)/2) & q <= (ci+(1-ci)/2)}]
      effSizeCI=quantile(effectSize14,probs=c(((1-ci)/2),(ci+(1-ci)/2)))
      effectSizeCIValues = effectSize14[effectSize14 >effSizeCI[1] & effectSize14 < effSizeCI[2]]
      largePosEffectIterations = length(which(effectSizeCIValues >= 0.8))/length(effectSizeCIValues)
      largeNegEffectIterations = length(which(effectSizeCIValues <= -0.8))/length(effectSizeCIValues)
      mediumPosEffectIterations = length(which(effectSizeCIValues >= 0.5 & effectSizeCIValues < 0.8))/length(effectSizeCIValues)
      mediumNegEffectIterations = length(which(effectSizeCIValues <= -0.5 & effectSizeCIValues > -0.8))/length(effectSizeCIValues)
      smallPosEffectIterations = length(which(effectSizeCIValues >= 0.2 & effectSizeCIValues < 0.5))/length(effectSizeCIValues)
      smallNegEffectIterations = length(which(effectSizeCIValues <= -0.2 & effectSizeCIValues > -0.5))/length(effectSizeCIValues)
      noEffectIterations = length(which(effectSizeCIValues < 0.2 & effectSizeCIValues > -0.2))/length(effectSizeCIValues)
      
      lbls=c("large positive:", "medium positive:", "small positive:", "no effect:", "small negative:", "medium negative:","large negative:")
      pct=c(largeNegEffectIterations,mediumNegEffectIterations,smallNegEffectIterations,noEffectIterations,smallPosEffectIterations,mediumPosEffectIterations,largePosEffectIterations)
      pct=round(pct,digits=4)
      pct=pct*100
      lbls=paste(lbls,pct)
      lbls <- paste(pct,"%",sep="") # ad % to labels
      #pie(pct,labels = lbls, col=rainbow(length(lbls)),
      #main="Posterior probabilities of varying effect sizes")
      
      # Barplot
      b<-barplot(pct,
                 main = mtext("Posterior ROPE-analysis for"~delta[14]),
                 xlab = "ROPEs of standard effect sizes",
                 ylab = "Posterior-Percentage included inside ROPE",
                 names.arg = c("L-", "M-", "S-", "No Eff.", "S+", "M+", "L+"),
                 col = "cornflowerblue",
                 border="white",
                 horiz = FALSE,
                 ylim=c(0,107.5))
      b
      text(x=b,y=pct+4,labels=lbls)
      
      
      
      ####################### POSTERIOR ROPE ANALYSIS FOR DELTA_24 #######################################
      # Bar plot
      #effSizeCIValues = effectSize[{q<-rank(effectSize)/length(effectSize);q>=((1-ci)/2) & q <= (ci+(1-ci)/2)}]
      effSizeCI=quantile(effectSize24,probs=c(((1-ci)/2),(ci+(1-ci)/2)))
      effectSizeCIValues = effectSize24[effectSize24 >effSizeCI[1] & effectSize24 < effSizeCI[2]]
      largePosEffectIterations = length(which(effectSizeCIValues >= 0.8))/length(effectSizeCIValues)
      largeNegEffectIterations = length(which(effectSizeCIValues <= -0.8))/length(effectSizeCIValues)
      mediumPosEffectIterations = length(which(effectSizeCIValues >= 0.5 & effectSizeCIValues < 0.8))/length(effectSizeCIValues)
      mediumNegEffectIterations = length(which(effectSizeCIValues <= -0.5 & effectSizeCIValues > -0.8))/length(effectSizeCIValues)
      smallPosEffectIterations = length(which(effectSizeCIValues >= 0.2 & effectSizeCIValues < 0.5))/length(effectSizeCIValues)
      smallNegEffectIterations = length(which(effectSizeCIValues <= -0.2 & effectSizeCIValues > -0.5))/length(effectSizeCIValues)
      noEffectIterations = length(which(effectSizeCIValues < 0.2 & effectSizeCIValues > -0.2))/length(effectSizeCIValues)
      
      lbls=c("large positive:", "medium positive:", "small positive:", "no effect:", "small negative:", "medium negative:","large negative:")
      pct=c(largeNegEffectIterations,mediumNegEffectIterations,smallNegEffectIterations,noEffectIterations,smallPosEffectIterations,mediumPosEffectIterations,largePosEffectIterations)
      pct=round(pct,digits=4)
      pct=pct*100
      lbls=paste(lbls,pct)
      lbls <- paste(pct,"%",sep="") # ad % to labels
      #pie(pct,labels = lbls, col=rainbow(length(lbls)),
      #main="Posterior probabilities of varying effect sizes")
      
      # Barplot
      b<-barplot(pct,
                 main = mtext("Posterior ROPE-analysis for"~delta[24]),
                              xlab = "ROPEs of standard effect sizes",
                              ylab = "Posterior-Percentage included inside ROPE",
                              names.arg = c("L-", "M-", "S-", "No Eff.", "S+", "M+", "L+"),
                              col = "cornflowerblue",
                              border="white",
                              horiz = FALSE,
                              ylim=c(0,107.5))
      b
      text(x=b,y=pct+4,labels=lbls)
      
      
      
      ####################### POSTERIOR ROPE ANALYSIS FOR DELTA_34 #######################################
      # Bar plot
      #effSizeCIValues = effectSize[{q<-rank(effectSize)/length(effectSize);q>=((1-ci)/2) & q <= (ci+(1-ci)/2)}]
      effSizeCI=quantile(effectSize34,probs=c(((1-ci)/2),(ci+(1-ci)/2)))
      effectSizeCIValues = effectSize34[effectSize34 >effSizeCI[1] & effectSize34 < effSizeCI[2]]
      largePosEffectIterations = length(which(effectSizeCIValues >= 0.8))/length(effectSizeCIValues)
      largeNegEffectIterations = length(which(effectSizeCIValues <= -0.8))/length(effectSizeCIValues)
      mediumPosEffectIterations = length(which(effectSizeCIValues >= 0.5 & effectSizeCIValues < 0.8))/length(effectSizeCIValues)
      mediumNegEffectIterations = length(which(effectSizeCIValues <= -0.5 & effectSizeCIValues > -0.8))/length(effectSizeCIValues)
      smallPosEffectIterations = length(which(effectSizeCIValues >= 0.2 & effectSizeCIValues < 0.5))/length(effectSizeCIValues)
      smallNegEffectIterations = length(which(effectSizeCIValues <= -0.2 & effectSizeCIValues > -0.5))/length(effectSizeCIValues)
      noEffectIterations = length(which(effectSizeCIValues < 0.2 & effectSizeCIValues > -0.2))/length(effectSizeCIValues)
      
      lbls=c("large positive:", "medium positive:", "small positive:", "no effect:", "small negative:", "medium negative:","large negative:")
      pct=c(largeNegEffectIterations,mediumNegEffectIterations,smallNegEffectIterations,noEffectIterations,smallPosEffectIterations,mediumPosEffectIterations,largePosEffectIterations)
      pct=round(pct,digits=4)
      pct=pct*100
      lbls=paste(lbls,pct)
      lbls <- paste(pct,"%",sep="") # ad % to labels
      #pie(pct,labels = lbls, col=rainbow(length(lbls)),
      #main="Posterior probabilities of varying effect sizes")
      b<-barplot(pct,
                 main = mtext("Posterior ROPE-analysis for"~delta[34]),
                 xlab = "ROPEs of standard effect sizes",
                 ylab = "Posterior-Percentage included inside ROPE",
                 names.arg = c("L-", "M-", "S-", "No Eff.", "S+", "M+", "L+"),
                 col = "cornflowerblue",
                 border="white",
                 horiz = FALSE,
                 ylim=c(0,107.5))
      b
      text(x=b,y=pct+4,labels=lbls)  
    }
    
    ###################################################### ADDITIONAL EFFECT SIZES FOR FIFTH GROUP #########################################################
    if(!is.null(dataframe$mu5)){
      # Set 2x4 par
      dev.new()
      newpar <- par(mfrow=c(2,4),mai = c(1,0.5,0.5,0.5), oma=c(1,1,1,1))
      # Reset par only after first change of par
      on.exit(par(newpar))
      
      ########################################### HISTOGRAM DELTA_15 ###########################################
      # Effect size posterior histogram
      # Make a vector of values to draw ticks at:
      hist(effectSize15,freq=FALSE,col="cornflowerblue",border="white",xlab=expression(delta[15]),main=mtext("Posterior distribution of Effect size "~delta[15], side=3,line=2.2,col="black",cex=0.8,font=2),xaxt="n")
      ticks <- seq(from=-10, to=10, by=0.05)
      # And draw the axis:
      axis(1, at=ticks)
      
      # Posterior mode
      #abline(v=getmode(effectSize),col="blue",lwd=2,lty="solid")
      lines(x=c(getmode(effectSize15),getmode(effectSize15)),y=c(-0.025,0.025),type="l",lwd=2.5,col="blue")
      
      # Posterior CI
      postCI=quantile(effectSize15,probs=c((1-ci)/2,ci+(1-ci)/2))
      lines(x=c(postCI[1],postCI[2]),y=c(0,0),type="l",lwd=2,col="black")
      lines(x=c(postCI[1],postCI[1]),y=c(-0.02,0.02),type="l",lwd=2,col="black")
      lines(x=c(postCI[2],postCI[2]),y=c(-0.02,0.02),type="l",lwd=2,col="black")
      
      # Visualization of typical ROPEs
      # ROPE for no effect
      lines(x=c(-0.2,0.2),y=c(0.04,0.04),type="l",lwd=2,col="red")
      # ROPE for small effect
      lines(x=c(-0.5,-0.2),y=c(0.04,0.04),type="l",lwd=2,col="orange")
      lines(x=c(0.2,0.5),y=c(0.04,0.04),type="l",lwd=2,col="orange")
      
      # ROPE for medium effect
      lines(x=c(-0.8,-0.5),y=c(0.04,0.04),type="l",lwd=2,col="green")
      lines(x=c(0.5,0.8),y=c(0.04,0.04),type="l",lwd=2,col="green")
      # ROPE for large effect
      lines(x=c(-2,-0.8),y=c(0.04,0.04),type="l",lwd=2,col="purple")
      lines(x=c(0.8,2),y=c(0.04,0.04),type="l",lwd=2,col="purple")
      
      abline(v=0.8,col="purple",lwd=0.5,lty="dashed")
      abline(v=-0.8,col="purple",lwd=0.5,lty="dashed")
      abline(v=0.5,col="green",lwd=0.5,lty="dashed")
      abline(v=-0.5,col="green",lwd=0.5,lty="dashed")
      abline(v=0.2,col="orange",lwd=0.5,lty="dashed")
      abline(v=-0.2,col="orange",lwd=0.5,lty="dashed")
      
      # Positions of text above the plot
      #posLeft=quantile(effectSize,probs=c(0.1,0.5,0.95))[1]
      #posMid=quantile(effectSize,probs=c(0.1,0.5,0.95))[2]
      #posRight=quantile(effectSize,probs=c(0.1,0.5,0.95))[3]
      
      # Text above the plot
      mtext(paste0(ci*100,"% ","CI: [", round(postCI[1], 3),",",round(postCI[2], 3),"]"), side=3,line=0.15,at=min(effectSize15),col="black",cex=0.8,font=2)
      mtext(paste0("Mean: ", round(mean(effectSize15),3)),side=3,line=1.15,at=min(effectSize15),col="black",cex=0.8,font=2)
      mtext(paste0("Mode: ", round(getmode(effectSize15),3)),side=3,line=1.15,at=max(effectSize15),col="black",cex=0.8,font=2)
      
      
      
      
      ########################################### HISTOGRAM DELTA_25 ###########################################
      # Effect size posterior histogram
      # Make a vector of values to draw ticks at:
      hist(effectSize25,freq=FALSE,col="cornflowerblue",border="white",xlab=expression(delta[25]),main=mtext("Posterior distribution of Effect size "~delta[25], side=3,line=2.2,col="black",cex=0.8,font=2),xaxt="n")
      ticks <- seq(from=-10, to=10, by=0.05)
      # And draw the axis:
      axis(1, at=ticks)
      
      # Posterior mode
      #abline(v=getmode(effectSize),col="blue",lwd=2,lty="solid")
      lines(x=c(getmode(effectSize25),getmode(effectSize25)),y=c(-0.025,0.025),type="l",lwd=2.5,col="blue")
      
      # Posterior CI
      postCI=quantile(effectSize25,probs=c((1-ci)/2,ci+(1-ci)/2))
      lines(x=c(postCI[1],postCI[2]),y=c(0,0),type="l",lwd=2,col="black")
      lines(x=c(postCI[1],postCI[1]),y=c(-0.02,0.02),type="l",lwd=2,col="black")
      lines(x=c(postCI[2],postCI[2]),y=c(-0.02,0.02),type="l",lwd=2,col="black")
      
      # Visualization of typical ROPEs
      # ROPE for no effect
      lines(x=c(-0.2,0.2),y=c(0.04,0.04),type="l",lwd=2,col="red")
      # ROPE for small effect
      lines(x=c(-0.5,-0.2),y=c(0.04,0.04),type="l",lwd=2,col="orange")
      lines(x=c(0.2,0.5),y=c(0.04,0.04),type="l",lwd=2,col="orange")
      
      # ROPE for medium effect
      lines(x=c(-0.8,-0.5),y=c(0.04,0.04),type="l",lwd=2,col="green")
      lines(x=c(0.5,0.8),y=c(0.04,0.04),type="l",lwd=2,col="green")
      # ROPE for large effect
      lines(x=c(-2,-0.8),y=c(0.04,0.04),type="l",lwd=2,col="purple")
      lines(x=c(0.8,2),y=c(0.04,0.04),type="l",lwd=2,col="purple")
      
      abline(v=0.8,col="purple",lwd=0.5,lty="dashed")
      abline(v=-0.8,col="purple",lwd=0.5,lty="dashed")
      abline(v=0.5,col="green",lwd=0.5,lty="dashed")
      abline(v=-0.5,col="green",lwd=0.5,lty="dashed")
      abline(v=0.2,col="orange",lwd=0.5,lty="dashed")
      abline(v=-0.2,col="orange",lwd=0.5,lty="dashed")
      
      # Positions of text above the plot
      #posLeft=quantile(effectSize,probs=c(0.1,0.5,0.95))[1]
      #posMid=quantile(effectSize,probs=c(0.1,0.5,0.95))[2]
      #posRight=quantile(effectSize,probs=c(0.1,0.5,0.95))[3]
      
      # Text above the plot
      mtext(paste0(ci*100,"% ","CI: [", round(postCI[1], 3),",",round(postCI[2], 3),"]"), side=3,line=0.15,at=min(effectSize25),col="black",cex=0.8,font=2)
      mtext(paste0("Mean: ", round(mean(effectSize25),3)),side=3,line=1.15,at=min(effectSize25),col="black",cex=0.8,font=2)
      mtext(paste0("Mode: ", round(getmode(effectSize25),3)),side=3,line=1.15,at=max(effectSize25),col="black",cex=0.8,font=2)
      
      
      
      ########################################### HISTOGRAM DELTA_35 ###########################################
      # Effect size posterior histogram
      # Make a vector of values to draw ticks at:
      hist(effectSize35,freq=FALSE,col="cornflowerblue",border="white",xlab=expression(delta[35]),main=mtext("Posterior distribution of Effect size "~delta[35], side=3,line=2.2,col="black",cex=0.8,font=2),xaxt="n")
      ticks <- seq(from=-10, to=10, by=0.05)
      # And draw the axis:
      axis(1, at=ticks)
      
      # Posterior mode
      #abline(v=getmode(effectSize),col="blue",lwd=2,lty="solid")
      lines(x=c(getmode(effectSize35),getmode(effectSize35)),y=c(-0.025,0.025),type="l",lwd=2.5,col="blue")
      
      # Posterior CI
      postCI=quantile(effectSize35,probs=c((1-ci)/2,ci+(1-ci)/2))
      lines(x=c(postCI[1],postCI[2]),y=c(0,0),type="l",lwd=2,col="black")
      lines(x=c(postCI[1],postCI[1]),y=c(-0.02,0.02),type="l",lwd=2,col="black")
      lines(x=c(postCI[2],postCI[2]),y=c(-0.02,0.02),type="l",lwd=2,col="black")
      
      # Visualization of typical ROPEs
      # ROPE for no effect
      lines(x=c(-0.2,0.2),y=c(0.04,0.04),type="l",lwd=2,col="red")
      # ROPE for small effect
      lines(x=c(-0.5,-0.2),y=c(0.04,0.04),type="l",lwd=2,col="orange")
      lines(x=c(0.2,0.5),y=c(0.04,0.04),type="l",lwd=2,col="orange")
      
      # ROPE for medium effect
      lines(x=c(-0.8,-0.5),y=c(0.04,0.04),type="l",lwd=2,col="green")
      lines(x=c(0.5,0.8),y=c(0.04,0.04),type="l",lwd=2,col="green")
      # ROPE for large effect
      lines(x=c(-2,-0.8),y=c(0.04,0.04),type="l",lwd=2,col="purple")
      lines(x=c(0.8,2),y=c(0.04,0.04),type="l",lwd=2,col="purple")
      
      abline(v=0.8,col="purple",lwd=0.5,lty="dashed")
      abline(v=-0.8,col="purple",lwd=0.5,lty="dashed")
      abline(v=0.5,col="green",lwd=0.5,lty="dashed")
      abline(v=-0.5,col="green",lwd=0.5,lty="dashed")
      abline(v=0.2,col="orange",lwd=0.5,lty="dashed")
      abline(v=-0.2,col="orange",lwd=0.5,lty="dashed")
      
      # Positions of text above the plot
      #posLeft=quantile(effectSize,probs=c(0.1,0.5,0.95))[1]
      #posMid=quantile(effectSize,probs=c(0.1,0.5,0.95))[2]
      #posRight=quantile(effectSize,probs=c(0.1,0.5,0.95))[3]
      
      # Text above the plot
      mtext(paste0(ci*100,"% ","CI: [", round(postCI[1], 3),",",round(postCI[2], 3),"]"), side=3,line=0.15,at=min(effectSize35),col="black",cex=0.8,font=2)
      mtext(paste0("Mean: ", round(mean(effectSize35),3)),side=3,line=1.15,at=min(effectSize35),col="black",cex=0.8,font=2)
      mtext(paste0("Mode: ", round(getmode(effectSize35),3)),side=3,line=1.15,at=max(effectSize35),col="black",cex=0.8,font=2)
      
      
      
      ########################################### HISTOGRAM DELTA_45 ###########################################
      # Effect size posterior histogram
      # Make a vector of values to draw ticks at:
      hist(effectSize45,freq=FALSE,col="cornflowerblue",border="white",xlab=expression(delta[45]),main=mtext("Posterior distribution of Effect size "~delta[45], side=3,line=2.2,col="black",cex=0.8,font=2),xaxt="n")
      ticks <- seq(from=-10, to=10, by=0.05)
      # And draw the axis:
      axis(1, at=ticks)
      
      # Posterior mode
      #abline(v=getmode(effectSize),col="blue",lwd=2,lty="solid")
      lines(x=c(getmode(effectSize45),getmode(effectSize45)),y=c(-0.025,0.025),type="l",lwd=2.5,col="blue")
      
      # Posterior CI
      postCI=quantile(effectSize45,probs=c((1-ci)/2,ci+(1-ci)/2))
      lines(x=c(postCI[1],postCI[2]),y=c(0,0),type="l",lwd=2,col="black")
      lines(x=c(postCI[1],postCI[1]),y=c(-0.02,0.02),type="l",lwd=2,col="black")
      lines(x=c(postCI[2],postCI[2]),y=c(-0.02,0.02),type="l",lwd=2,col="black")
      
      # Visualization of typical ROPEs
      # ROPE for no effect
      lines(x=c(-0.2,0.2),y=c(0.04,0.04),type="l",lwd=2,col="red")
      # ROPE for small effect
      lines(x=c(-0.5,-0.2),y=c(0.04,0.04),type="l",lwd=2,col="orange")
      lines(x=c(0.2,0.5),y=c(0.04,0.04),type="l",lwd=2,col="orange")
      
      # ROPE for medium effect
      lines(x=c(-0.8,-0.5),y=c(0.04,0.04),type="l",lwd=2,col="green")
      lines(x=c(0.5,0.8),y=c(0.04,0.04),type="l",lwd=2,col="green")
      # ROPE for large effect
      lines(x=c(-2,-0.8),y=c(0.04,0.04),type="l",lwd=2,col="purple")
      lines(x=c(0.8,2),y=c(0.04,0.04),type="l",lwd=2,col="purple")
      
      abline(v=0.8,col="purple",lwd=0.5,lty="dashed")
      abline(v=-0.8,col="purple",lwd=0.5,lty="dashed")
      abline(v=0.5,col="green",lwd=0.5,lty="dashed")
      abline(v=-0.5,col="green",lwd=0.5,lty="dashed")
      abline(v=0.2,col="orange",lwd=0.5,lty="dashed")
      abline(v=-0.2,col="orange",lwd=0.5,lty="dashed")
      
      # Positions of text above the plot
      #posLeft=quantile(effectSize,probs=c(0.1,0.5,0.95))[1]
      #posMid=quantile(effectSize,probs=c(0.1,0.5,0.95))[2]
      #posRight=quantile(effectSize,probs=c(0.1,0.5,0.95))[3]
      
      # Text above the plot
      mtext(paste0(ci*100,"% ","CI: [", round(postCI[1], 3),",",round(postCI[2], 3),"]"), side=3,line=0.15,at=min(effectSize45),col="black",cex=0.8,font=2)
      mtext(paste0("Mean: ", round(mean(effectSize45),3)),side=3,line=1.15,at=min(effectSize45),col="black",cex=0.8,font=2)
      mtext(paste0("Mode: ", round(getmode(effectSize45),3)),side=3,line=1.15,at=max(effectSize45),col="black",cex=0.8,font=2)
      
      
      ####################### POSTERIOR ROPE ANALYSIS FOR DELTA_15 #######################################
      # Bar plot
      #effSizeCIValues = effectSize[{q<-rank(effectSize)/length(effectSize);q>=((1-ci)/2) & q <= (ci+(1-ci)/2)}]
      effSizeCI=quantile(effectSize15,probs=c(((1-ci)/2),(ci+(1-ci)/2)))
      effectSizeCIValues = effectSize15[effectSize15 >effSizeCI[1] & effectSize15 < effSizeCI[2]]
      largePosEffectIterations = length(which(effectSizeCIValues >= 0.8))/length(effectSizeCIValues)
      largeNegEffectIterations = length(which(effectSizeCIValues <= -0.8))/length(effectSizeCIValues)
      mediumPosEffectIterations = length(which(effectSizeCIValues >= 0.5 & effectSizeCIValues < 0.8))/length(effectSizeCIValues)
      mediumNegEffectIterations = length(which(effectSizeCIValues <= -0.5 & effectSizeCIValues > -0.8))/length(effectSizeCIValues)
      smallPosEffectIterations = length(which(effectSizeCIValues >= 0.2 & effectSizeCIValues < 0.5))/length(effectSizeCIValues)
      smallNegEffectIterations = length(which(effectSizeCIValues <= -0.2 & effectSizeCIValues > -0.5))/length(effectSizeCIValues)
      noEffectIterations = length(which(effectSizeCIValues < 0.2 & effectSizeCIValues > -0.2))/length(effectSizeCIValues)
      
      lbls=c("large positive:", "medium positive:", "small positive:", "no effect:", "small negative:", "medium negative:","large negative:")
      pct=c(largeNegEffectIterations,mediumNegEffectIterations,smallNegEffectIterations,noEffectIterations,smallPosEffectIterations,mediumPosEffectIterations,largePosEffectIterations)
      pct=round(pct,digits=4)
      pct=pct*100
      lbls=paste(lbls,pct)
      lbls <- paste(pct,"%",sep="") # ad % to labels
      #pie(pct,labels = lbls, col=rainbow(length(lbls)),
      #main="Posterior probabilities of varying effect sizes")
      
      # Barplot
      b<-barplot(pct,
                 main = mtext("Posterior ROPE-analysis for"~delta[15]),
                 xlab = "ROPEs of standard effect sizes",
                 ylab = "Posterior-Percentage included inside ROPE",
                 names.arg = c("L-", "M-", "S-", "No Eff.", "S+", "M+", "L+"),
                 col = "cornflowerblue",
                 border="white",
                 horiz = FALSE,
                 ylim=c(0,107.5))
      b
      text(x=b,y=pct+4,labels=lbls)
      
      
      
      ####################### POSTERIOR ROPE ANALYSIS FOR DELTA_25 #######################################
      # Bar plot
      #effSizeCIValues = effectSize[{q<-rank(effectSize)/length(effectSize);q>=((1-ci)/2) & q <= (ci+(1-ci)/2)}]
      effSizeCI=quantile(effectSize25,probs=c(((1-ci)/2),(ci+(1-ci)/2)))
      effectSizeCIValues = effectSize25[effectSize25 >effSizeCI[1] & effectSize25 < effSizeCI[2]]
      largePosEffectIterations = length(which(effectSizeCIValues >= 0.8))/length(effectSizeCIValues)
      largeNegEffectIterations = length(which(effectSizeCIValues <= -0.8))/length(effectSizeCIValues)
      mediumPosEffectIterations = length(which(effectSizeCIValues >= 0.5 & effectSizeCIValues < 0.8))/length(effectSizeCIValues)
      mediumNegEffectIterations = length(which(effectSizeCIValues <= -0.5 & effectSizeCIValues > -0.8))/length(effectSizeCIValues)
      smallPosEffectIterations = length(which(effectSizeCIValues >= 0.2 & effectSizeCIValues < 0.5))/length(effectSizeCIValues)
      smallNegEffectIterations = length(which(effectSizeCIValues <= -0.2 & effectSizeCIValues > -0.5))/length(effectSizeCIValues)
      noEffectIterations = length(which(effectSizeCIValues < 0.2 & effectSizeCIValues > -0.2))/length(effectSizeCIValues)
      
      lbls=c("large positive:", "medium positive:", "small positive:", "no effect:", "small negative:", "medium negative:","large negative:")
      pct=c(largeNegEffectIterations,mediumNegEffectIterations,smallNegEffectIterations,noEffectIterations,smallPosEffectIterations,mediumPosEffectIterations,largePosEffectIterations)
      pct=round(pct,digits=4)
      pct=pct*100
      lbls=paste(lbls,pct)
      lbls <- paste(pct,"%",sep="") # ad % to labels
      #pie(pct,labels = lbls, col=rainbow(length(lbls)),
      #main="Posterior probabilities of varying effect sizes")
      
      # Barplot
      b<-barplot(pct,
                 main = mtext("Posterior ROPE-analysis for"~delta[25]),
                 xlab = "ROPEs of standard effect sizes",
                 ylab = "Posterior-Percentage included inside ROPE",
                 names.arg = c("L-", "M-", "S-", "No Eff.", "S+", "M+", "L+"),
                 col = "cornflowerblue",
                 border="white",
                 horiz = FALSE,
                 ylim=c(0,107.5))
      b
      text(x=b,y=pct+4,labels=lbls)
      
      
      
      ####################### POSTERIOR ROPE ANALYSIS FOR DELTA_35 #######################################
      # Bar plot
      #effSizeCIValues = effectSize[{q<-rank(effectSize)/length(effectSize);q>=((1-ci)/2) & q <= (ci+(1-ci)/2)}]
      effSizeCI=quantile(effectSize35,probs=c(((1-ci)/2),(ci+(1-ci)/2)))
      effectSizeCIValues = effectSize35[effectSize35 >effSizeCI[1] & effectSize35 < effSizeCI[2]]
      largePosEffectIterations = length(which(effectSizeCIValues >= 0.8))/length(effectSizeCIValues)
      largeNegEffectIterations = length(which(effectSizeCIValues <= -0.8))/length(effectSizeCIValues)
      mediumPosEffectIterations = length(which(effectSizeCIValues >= 0.5 & effectSizeCIValues < 0.8))/length(effectSizeCIValues)
      mediumNegEffectIterations = length(which(effectSizeCIValues <= -0.5 & effectSizeCIValues > -0.8))/length(effectSizeCIValues)
      smallPosEffectIterations = length(which(effectSizeCIValues >= 0.2 & effectSizeCIValues < 0.5))/length(effectSizeCIValues)
      smallNegEffectIterations = length(which(effectSizeCIValues <= -0.2 & effectSizeCIValues > -0.5))/length(effectSizeCIValues)
      noEffectIterations = length(which(effectSizeCIValues < 0.2 & effectSizeCIValues > -0.2))/length(effectSizeCIValues)
      
      lbls=c("large positive:", "medium positive:", "small positive:", "no effect:", "small negative:", "medium negative:","large negative:")
      pct=c(largeNegEffectIterations,mediumNegEffectIterations,smallNegEffectIterations,noEffectIterations,smallPosEffectIterations,mediumPosEffectIterations,largePosEffectIterations)
      pct=round(pct,digits=4)
      pct=pct*100
      lbls=paste(lbls,pct)
      lbls <- paste(pct,"%",sep="") # ad % to labels
      #pie(pct,labels = lbls, col=rainbow(length(lbls)),
      #main="Posterior probabilities of varying effect sizes")
      b<-barplot(pct,
                 main = mtext("Posterior ROPE-analysis for"~delta[35]),
                 xlab = "ROPEs of standard effect sizes",
                 ylab = "Posterior-Percentage included inside ROPE",
                 names.arg = c("L-", "M-", "S-", "No Eff.", "S+", "M+", "L+"),
                 col = "cornflowerblue",
                 border="white",
                 horiz = FALSE,
                 ylim=c(0,107.5))
      b
      text(x=b,y=pct+4,labels=lbls)  
      
      
      
      ####################### POSTERIOR ROPE ANALYSIS FOR DELTA_45 #######################################
      # Bar plot
      #effSizeCIValues = effectSize[{q<-rank(effectSize)/length(effectSize);q>=((1-ci)/2) & q <= (ci+(1-ci)/2)}]
      effSizeCI=quantile(effectSize45,probs=c(((1-ci)/2),(ci+(1-ci)/2)))
      effectSizeCIValues = effectSize45[effectSize45 >effSizeCI[1] & effectSize45 < effSizeCI[2]]
      largePosEffectIterations = length(which(effectSizeCIValues >= 0.8))/length(effectSizeCIValues)
      largeNegEffectIterations = length(which(effectSizeCIValues <= -0.8))/length(effectSizeCIValues)
      mediumPosEffectIterations = length(which(effectSizeCIValues >= 0.5 & effectSizeCIValues < 0.8))/length(effectSizeCIValues)
      mediumNegEffectIterations = length(which(effectSizeCIValues <= -0.5 & effectSizeCIValues > -0.8))/length(effectSizeCIValues)
      smallPosEffectIterations = length(which(effectSizeCIValues >= 0.2 & effectSizeCIValues < 0.5))/length(effectSizeCIValues)
      smallNegEffectIterations = length(which(effectSizeCIValues <= -0.2 & effectSizeCIValues > -0.5))/length(effectSizeCIValues)
      noEffectIterations = length(which(effectSizeCIValues < 0.2 & effectSizeCIValues > -0.2))/length(effectSizeCIValues)
      
      lbls=c("large positive:", "medium positive:", "small positive:", "no effect:", "small negative:", "medium negative:","large negative:")
      pct=c(largeNegEffectIterations,mediumNegEffectIterations,smallNegEffectIterations,noEffectIterations,smallPosEffectIterations,mediumPosEffectIterations,largePosEffectIterations)
      pct=round(pct,digits=4)
      pct=pct*100
      lbls=paste(lbls,pct)
      lbls <- paste(pct,"%",sep="") # ad % to labels
      #pie(pct,labels = lbls, col=rainbow(length(lbls)),
      #main="Posterior probabilities of varying effect sizes")
      b<-barplot(pct,
                 main = mtext("Posterior ROPE-analysis for"~delta[45]),
                 xlab = "ROPEs of standard effect sizes",
                 ylab = "Posterior-Percentage included inside ROPE",
                 names.arg = c("L-", "M-", "S-", "No Eff.", "S+", "M+", "L+"),
                 col = "cornflowerblue",
                 border="white",
                 horiz = FALSE,
                 ylim=c(0,107.5))
      b
      text(x=b,y=pct+4,labels=lbls) 
    }
    
    ###################################################### ADDITIONAL EFFECT SIZES FOR SIXTH GROUP #########################################################
    if(!is.null(dataframe$mu6)){
      # Set 2x4 par
      dev.new()
      newpar <- par(mfrow=c(2,5),mai = c(1,0.5,0.5,0.5), oma=c(1,1,1,1))
      # Reset par only after first change of par
      on.exit(par(newpar))
      
      ########################################### HISTOGRAM DELTA_16 ###########################################
      # Effect size posterior histogram
      # Make a vector of values to draw ticks at:
      hist(effectSize16,freq=FALSE,col="cornflowerblue",border="white",xlab=expression(delta[16]),main=mtext("Posterior distribution of Effect size "~delta[16], side=3,line=2.2,col="black",cex=0.8,font=2),xaxt="n")
      ticks <- seq(from=-10, to=10, by=0.05)
      # And draw the axis:
      axis(1, at=ticks)
      
      # Posterior mode
      #abline(v=getmode(effectSize),col="blue",lwd=2,lty="solid")
      lines(x=c(getmode(effectSize16),getmode(effectSize16)),y=c(-0.025,0.025),type="l",lwd=2.5,col="blue")
      
      # Posterior CI
      postCI=quantile(effectSize16,probs=c((1-ci)/2,ci+(1-ci)/2))
      lines(x=c(postCI[1],postCI[2]),y=c(0,0),type="l",lwd=2,col="black")
      lines(x=c(postCI[1],postCI[1]),y=c(-0.02,0.02),type="l",lwd=2,col="black")
      lines(x=c(postCI[2],postCI[2]),y=c(-0.02,0.02),type="l",lwd=2,col="black")
      
      # Visualization of typical ROPEs
      # ROPE for no effect
      lines(x=c(-0.2,0.2),y=c(0.04,0.04),type="l",lwd=2,col="red")
      # ROPE for small effect
      lines(x=c(-0.5,-0.2),y=c(0.04,0.04),type="l",lwd=2,col="orange")
      lines(x=c(0.2,0.5),y=c(0.04,0.04),type="l",lwd=2,col="orange")
      
      # ROPE for medium effect
      lines(x=c(-0.8,-0.5),y=c(0.04,0.04),type="l",lwd=2,col="green")
      lines(x=c(0.5,0.8),y=c(0.04,0.04),type="l",lwd=2,col="green")
      # ROPE for large effect
      lines(x=c(-2,-0.8),y=c(0.04,0.04),type="l",lwd=2,col="purple")
      lines(x=c(0.8,2),y=c(0.04,0.04),type="l",lwd=2,col="purple")
      
      abline(v=0.8,col="purple",lwd=0.5,lty="dashed")
      abline(v=-0.8,col="purple",lwd=0.5,lty="dashed")
      abline(v=0.5,col="green",lwd=0.5,lty="dashed")
      abline(v=-0.5,col="green",lwd=0.5,lty="dashed")
      abline(v=0.2,col="orange",lwd=0.5,lty="dashed")
      abline(v=-0.2,col="orange",lwd=0.5,lty="dashed")
      
      # Positions of text above the plot
      #posLeft=quantile(effectSize,probs=c(0.1,0.5,0.95))[1]
      #posMid=quantile(effectSize,probs=c(0.1,0.5,0.95))[2]
      #posRight=quantile(effectSize,probs=c(0.1,0.5,0.95))[3]
      
      # Text above the plot
      mtext(paste0(ci*100,"% ","CI: [", round(postCI[1], 3),",",round(postCI[2], 3),"]"), side=3,line=0.15,at=min(effectSize16),col="black",cex=0.8,font=2)
      mtext(paste0("Mean: ", round(mean(effectSize16),3)),side=3,line=1.15,at=min(effectSize16),col="black",cex=0.8,font=2)
      mtext(paste0("Mode: ", round(getmode(effectSize16),3)),side=3,line=1.15,at=max(effectSize16),col="black",cex=0.8,font=2)
      
      
      
      
      ########################################### HISTOGRAM DELTA_26 ###########################################
      # Effect size posterior histogram
      # Make a vector of values to draw ticks at:
      hist(effectSize26,freq=FALSE,col="cornflowerblue",border="white",xlab=expression(delta[26]),main=mtext("Posterior distribution of Effect size "~delta[26], side=3,line=2.2,col="black",cex=0.8,font=2),xaxt="n")
      ticks <- seq(from=-10, to=10, by=0.05)
      # And draw the axis:
      axis(1, at=ticks)
      
      # Posterior mode
      #abline(v=getmode(effectSize),col="blue",lwd=2,lty="solid")
      lines(x=c(getmode(effectSize26),getmode(effectSize26)),y=c(-0.025,0.025),type="l",lwd=2.5,col="blue")
      
      # Posterior CI
      postCI=quantile(effectSize26,probs=c((1-ci)/2,ci+(1-ci)/2))
      lines(x=c(postCI[1],postCI[2]),y=c(0,0),type="l",lwd=2,col="black")
      lines(x=c(postCI[1],postCI[1]),y=c(-0.02,0.02),type="l",lwd=2,col="black")
      lines(x=c(postCI[2],postCI[2]),y=c(-0.02,0.02),type="l",lwd=2,col="black")
      
      # Visualization of typical ROPEs
      # ROPE for no effect
      lines(x=c(-0.2,0.2),y=c(0.04,0.04),type="l",lwd=2,col="red")
      # ROPE for small effect
      lines(x=c(-0.5,-0.2),y=c(0.04,0.04),type="l",lwd=2,col="orange")
      lines(x=c(0.2,0.5),y=c(0.04,0.04),type="l",lwd=2,col="orange")
      
      # ROPE for medium effect
      lines(x=c(-0.8,-0.5),y=c(0.04,0.04),type="l",lwd=2,col="green")
      lines(x=c(0.5,0.8),y=c(0.04,0.04),type="l",lwd=2,col="green")
      # ROPE for large effect
      lines(x=c(-2,-0.8),y=c(0.04,0.04),type="l",lwd=2,col="purple")
      lines(x=c(0.8,2),y=c(0.04,0.04),type="l",lwd=2,col="purple")
      
      abline(v=0.8,col="purple",lwd=0.5,lty="dashed")
      abline(v=-0.8,col="purple",lwd=0.5,lty="dashed")
      abline(v=0.5,col="green",lwd=0.5,lty="dashed")
      abline(v=-0.5,col="green",lwd=0.5,lty="dashed")
      abline(v=0.2,col="orange",lwd=0.5,lty="dashed")
      abline(v=-0.2,col="orange",lwd=0.5,lty="dashed")
      
      # Positions of text above the plot
      #posLeft=quantile(effectSize,probs=c(0.1,0.5,0.95))[1]
      #posMid=quantile(effectSize,probs=c(0.1,0.5,0.95))[2]
      #posRight=quantile(effectSize,probs=c(0.1,0.5,0.95))[3]
      
      # Text above the plot
      mtext(paste0(ci*100,"% ","CI: [", round(postCI[1], 3),",",round(postCI[2], 3),"]"), side=3,line=0.15,at=min(effectSize26),col="black",cex=0.8,font=2)
      mtext(paste0("Mean: ", round(mean(effectSize26),3)),side=3,line=1.15,at=min(effectSize26),col="black",cex=0.8,font=2)
      mtext(paste0("Mode: ", round(getmode(effectSize26),3)),side=3,line=1.15,at=max(effectSize26),col="black",cex=0.8,font=2)
      
      
      
      ########################################### HISTOGRAM DELTA_36 ###########################################
      # Effect size posterior histogram
      # Make a vector of values to draw ticks at:
      hist(effectSize36,freq=FALSE,col="cornflowerblue",border="white",xlab=expression(delta[36]),main=mtext("Posterior distribution of Effect size "~delta[36], side=3,line=2.2,col="black",cex=0.8,font=2),xaxt="n")
      ticks <- seq(from=-10, to=10, by=0.05)
      # And draw the axis:
      axis(1, at=ticks)
      
      # Posterior mode
      #abline(v=getmode(effectSize),col="blue",lwd=2,lty="solid")
      lines(x=c(getmode(effectSize36),getmode(effectSize36)),y=c(-0.025,0.025),type="l",lwd=2.5,col="blue")
      
      # Posterior CI
      postCI=quantile(effectSize36,probs=c((1-ci)/2,ci+(1-ci)/2))
      lines(x=c(postCI[1],postCI[2]),y=c(0,0),type="l",lwd=2,col="black")
      lines(x=c(postCI[1],postCI[1]),y=c(-0.02,0.02),type="l",lwd=2,col="black")
      lines(x=c(postCI[2],postCI[2]),y=c(-0.02,0.02),type="l",lwd=2,col="black")
      
      # Visualization of typical ROPEs
      # ROPE for no effect
      lines(x=c(-0.2,0.2),y=c(0.04,0.04),type="l",lwd=2,col="red")
      # ROPE for small effect
      lines(x=c(-0.5,-0.2),y=c(0.04,0.04),type="l",lwd=2,col="orange")
      lines(x=c(0.2,0.5),y=c(0.04,0.04),type="l",lwd=2,col="orange")
      
      # ROPE for medium effect
      lines(x=c(-0.8,-0.5),y=c(0.04,0.04),type="l",lwd=2,col="green")
      lines(x=c(0.5,0.8),y=c(0.04,0.04),type="l",lwd=2,col="green")
      # ROPE for large effect
      lines(x=c(-2,-0.8),y=c(0.04,0.04),type="l",lwd=2,col="purple")
      lines(x=c(0.8,2),y=c(0.04,0.04),type="l",lwd=2,col="purple")
      
      abline(v=0.8,col="purple",lwd=0.5,lty="dashed")
      abline(v=-0.8,col="purple",lwd=0.5,lty="dashed")
      abline(v=0.5,col="green",lwd=0.5,lty="dashed")
      abline(v=-0.5,col="green",lwd=0.5,lty="dashed")
      abline(v=0.2,col="orange",lwd=0.5,lty="dashed")
      abline(v=-0.2,col="orange",lwd=0.5,lty="dashed")
      
      # Positions of text above the plot
      #posLeft=quantile(effectSize,probs=c(0.1,0.5,0.95))[1]
      #posMid=quantile(effectSize,probs=c(0.1,0.5,0.95))[2]
      #posRight=quantile(effectSize,probs=c(0.1,0.5,0.95))[3]
      
      # Text above the plot
      mtext(paste0(ci*100,"% ","CI: [", round(postCI[1], 3),",",round(postCI[2], 3),"]"), side=3,line=0.15,at=min(effectSize36),col="black",cex=0.8,font=2)
      mtext(paste0("Mean: ", round(mean(effectSize36),3)),side=3,line=1.15,at=min(effectSize36),col="black",cex=0.8,font=2)
      mtext(paste0("Mode: ", round(getmode(effectSize36),3)),side=3,line=1.15,at=max(effectSize36),col="black",cex=0.8,font=2)
      
      
      
      ########################################### HISTOGRAM DELTA_46 ###########################################
      # Effect size posterior histogram
      # Make a vector of values to draw ticks at:
      hist(effectSize46,freq=FALSE,col="cornflowerblue",border="white",xlab=expression(delta[46]),main=mtext("Posterior distribution of Effect size "~delta[46], side=3,line=2.2,col="black",cex=0.8,font=2),xaxt="n")
      ticks <- seq(from=-10, to=10, by=0.05)
      # And draw the axis:
      axis(1, at=ticks)
      
      # Posterior mode
      #abline(v=getmode(effectSize),col="blue",lwd=2,lty="solid")
      lines(x=c(getmode(effectSize46),getmode(effectSize46)),y=c(-0.025,0.025),type="l",lwd=2.5,col="blue")
      
      # Posterior CI
      postCI=quantile(effectSize46,probs=c((1-ci)/2,ci+(1-ci)/2))
      lines(x=c(postCI[1],postCI[2]),y=c(0,0),type="l",lwd=2,col="black")
      lines(x=c(postCI[1],postCI[1]),y=c(-0.02,0.02),type="l",lwd=2,col="black")
      lines(x=c(postCI[2],postCI[2]),y=c(-0.02,0.02),type="l",lwd=2,col="black")
      
      # Visualization of typical ROPEs
      # ROPE for no effect
      lines(x=c(-0.2,0.2),y=c(0.04,0.04),type="l",lwd=2,col="red")
      # ROPE for small effect
      lines(x=c(-0.5,-0.2),y=c(0.04,0.04),type="l",lwd=2,col="orange")
      lines(x=c(0.2,0.5),y=c(0.04,0.04),type="l",lwd=2,col="orange")
      
      # ROPE for medium effect
      lines(x=c(-0.8,-0.5),y=c(0.04,0.04),type="l",lwd=2,col="green")
      lines(x=c(0.5,0.8),y=c(0.04,0.04),type="l",lwd=2,col="green")
      # ROPE for large effect
      lines(x=c(-2,-0.8),y=c(0.04,0.04),type="l",lwd=2,col="purple")
      lines(x=c(0.8,2),y=c(0.04,0.04),type="l",lwd=2,col="purple")
      
      abline(v=0.8,col="purple",lwd=0.5,lty="dashed")
      abline(v=-0.8,col="purple",lwd=0.5,lty="dashed")
      abline(v=0.5,col="green",lwd=0.5,lty="dashed")
      abline(v=-0.5,col="green",lwd=0.5,lty="dashed")
      abline(v=0.2,col="orange",lwd=0.5,lty="dashed")
      abline(v=-0.2,col="orange",lwd=0.5,lty="dashed")
      
      # Positions of text above the plot
      #posLeft=quantile(effectSize,probs=c(0.1,0.5,0.95))[1]
      #posMid=quantile(effectSize,probs=c(0.1,0.5,0.95))[2]
      #posRight=quantile(effectSize,probs=c(0.1,0.5,0.95))[3]
      
      # Text above the plot
      mtext(paste0(ci*100,"% ","CI: [", round(postCI[1], 3),",",round(postCI[2], 3),"]"), side=3,line=0.15,at=min(effectSize46),col="black",cex=0.8,font=2)
      mtext(paste0("Mean: ", round(mean(effectSize46),3)),side=3,line=1.15,at=min(effectSize46),col="black",cex=0.8,font=2)
      mtext(paste0("Mode: ", round(getmode(effectSize46),3)),side=3,line=1.15,at=max(effectSize46),col="black",cex=0.8,font=2)
      
      
      
      ########################################### HISTOGRAM DELTA_56 ###########################################
      # Effect size posterior histogram
      # Make a vector of values to draw ticks at:
      hist(effectSize56,freq=FALSE,col="cornflowerblue",border="white",xlab=expression(delta[56]),main=mtext("Posterior distribution of Effect size "~delta[56], side=3,line=2.2,col="black",cex=0.8,font=2),xaxt="n")
      ticks <- seq(from=-10, to=10, by=0.05)
      # And draw the axis:
      axis(1, at=ticks)
      
      # Posterior mode
      #abline(v=getmode(effectSize),col="blue",lwd=2,lty="solid")
      lines(x=c(getmode(effectSize56),getmode(effectSize56)),y=c(-0.025,0.025),type="l",lwd=2.5,col="blue")
      
      # Posterior CI
      postCI=quantile(effectSize56,probs=c((1-ci)/2,ci+(1-ci)/2))
      lines(x=c(postCI[1],postCI[2]),y=c(0,0),type="l",lwd=2,col="black")
      lines(x=c(postCI[1],postCI[1]),y=c(-0.02,0.02),type="l",lwd=2,col="black")
      lines(x=c(postCI[2],postCI[2]),y=c(-0.02,0.02),type="l",lwd=2,col="black")
      
      # Visualization of typical ROPEs
      # ROPE for no effect
      lines(x=c(-0.2,0.2),y=c(0.04,0.04),type="l",lwd=2,col="red")
      # ROPE for small effect
      lines(x=c(-0.5,-0.2),y=c(0.04,0.04),type="l",lwd=2,col="orange")
      lines(x=c(0.2,0.5),y=c(0.04,0.04),type="l",lwd=2,col="orange")
      
      # ROPE for medium effect
      lines(x=c(-0.8,-0.5),y=c(0.04,0.04),type="l",lwd=2,col="green")
      lines(x=c(0.5,0.8),y=c(0.04,0.04),type="l",lwd=2,col="green")
      # ROPE for large effect
      lines(x=c(-2,-0.8),y=c(0.04,0.04),type="l",lwd=2,col="purple")
      lines(x=c(0.8,2),y=c(0.04,0.04),type="l",lwd=2,col="purple")
      
      abline(v=0.8,col="purple",lwd=0.5,lty="dashed")
      abline(v=-0.8,col="purple",lwd=0.5,lty="dashed")
      abline(v=0.5,col="green",lwd=0.5,lty="dashed")
      abline(v=-0.5,col="green",lwd=0.5,lty="dashed")
      abline(v=0.2,col="orange",lwd=0.5,lty="dashed")
      abline(v=-0.2,col="orange",lwd=0.5,lty="dashed")
      
      # Positions of text above the plot
      #posLeft=quantile(effectSize,probs=c(0.1,0.5,0.95))[1]
      #posMid=quantile(effectSize,probs=c(0.1,0.5,0.95))[2]
      #posRight=quantile(effectSize,probs=c(0.1,0.5,0.95))[3]
      
      # Text above the plot
      mtext(paste0(ci*100,"% ","CI: [", round(postCI[1], 3),",",round(postCI[2], 3),"]"), side=3,line=0.15,at=min(effectSize56),col="black",cex=0.8,font=2)
      mtext(paste0("Mean: ", round(mean(effectSize56),3)),side=3,line=1.15,at=min(effectSize56),col="black",cex=0.8,font=2)
      mtext(paste0("Mode: ", round(getmode(effectSize56),3)),side=3,line=1.15,at=max(effectSize56),col="black",cex=0.8,font=2)
      
      
      
      ####################### POSTERIOR ROPE ANALYSIS FOR DELTA_15 #######################################
      # Bar plot
      #effSizeCIValues = effectSize[{q<-rank(effectSize)/length(effectSize);q>=((1-ci)/2) & q <= (ci+(1-ci)/2)}]
      effSizeCI=quantile(effectSize16,probs=c(((1-ci)/2),(ci+(1-ci)/2)))
      effectSizeCIValues = effectSize16[effectSize16 >effSizeCI[1] & effectSize16 < effSizeCI[2]]
      largePosEffectIterations = length(which(effectSizeCIValues >= 0.8))/length(effectSizeCIValues)
      largeNegEffectIterations = length(which(effectSizeCIValues <= -0.8))/length(effectSizeCIValues)
      mediumPosEffectIterations = length(which(effectSizeCIValues >= 0.5 & effectSizeCIValues < 0.8))/length(effectSizeCIValues)
      mediumNegEffectIterations = length(which(effectSizeCIValues <= -0.5 & effectSizeCIValues > -0.8))/length(effectSizeCIValues)
      smallPosEffectIterations = length(which(effectSizeCIValues >= 0.2 & effectSizeCIValues < 0.5))/length(effectSizeCIValues)
      smallNegEffectIterations = length(which(effectSizeCIValues <= -0.2 & effectSizeCIValues > -0.5))/length(effectSizeCIValues)
      noEffectIterations = length(which(effectSizeCIValues < 0.2 & effectSizeCIValues > -0.2))/length(effectSizeCIValues)
      
      lbls=c("large positive:", "medium positive:", "small positive:", "no effect:", "small negative:", "medium negative:","large negative:")
      pct=c(largeNegEffectIterations,mediumNegEffectIterations,smallNegEffectIterations,noEffectIterations,smallPosEffectIterations,mediumPosEffectIterations,largePosEffectIterations)
      pct=round(pct,digits=4)
      pct=pct*100
      lbls=paste(lbls,pct)
      lbls <- paste(pct,"%",sep="") # ad % to labels
      #pie(pct,labels = lbls, col=rainbow(length(lbls)),
      #main="Posterior probabilities of varying effect sizes")
      
      # Barplot
      b<-barplot(pct,
                 main = mtext("Posterior ROPE-analysis for"~delta[16]),
                 xlab = "ROPEs of standard effect sizes",
                 ylab = "Posterior-Percentage included inside ROPE",
                 names.arg = c("L-", "M-", "S-", "No Eff.", "S+", "M+", "L+"),
                 col = "cornflowerblue",
                 border="white",
                 horiz = FALSE,
                 ylim=c(0,107.5))
      b
      text(x=b,y=pct+4,labels=lbls)
      
      
      
      ####################### POSTERIOR ROPE ANALYSIS FOR DELTA_26 #######################################
      # Bar plot
      #effSizeCIValues = effectSize[{q<-rank(effectSize)/length(effectSize);q>=((1-ci)/2) & q <= (ci+(1-ci)/2)}]
      effSizeCI=quantile(effectSize26,probs=c(((1-ci)/2),(ci+(1-ci)/2)))
      effectSizeCIValues = effectSize26[effectSize26 >effSizeCI[1] & effectSize26 < effSizeCI[2]]
      largePosEffectIterations = length(which(effectSizeCIValues >= 0.8))/length(effectSizeCIValues)
      largeNegEffectIterations = length(which(effectSizeCIValues <= -0.8))/length(effectSizeCIValues)
      mediumPosEffectIterations = length(which(effectSizeCIValues >= 0.5 & effectSizeCIValues < 0.8))/length(effectSizeCIValues)
      mediumNegEffectIterations = length(which(effectSizeCIValues <= -0.5 & effectSizeCIValues > -0.8))/length(effectSizeCIValues)
      smallPosEffectIterations = length(which(effectSizeCIValues >= 0.2 & effectSizeCIValues < 0.5))/length(effectSizeCIValues)
      smallNegEffectIterations = length(which(effectSizeCIValues <= -0.2 & effectSizeCIValues > -0.5))/length(effectSizeCIValues)
      noEffectIterations = length(which(effectSizeCIValues < 0.2 & effectSizeCIValues > -0.2))/length(effectSizeCIValues)
      
      lbls=c("large positive:", "medium positive:", "small positive:", "no effect:", "small negative:", "medium negative:","large negative:")
      pct=c(largeNegEffectIterations,mediumNegEffectIterations,smallNegEffectIterations,noEffectIterations,smallPosEffectIterations,mediumPosEffectIterations,largePosEffectIterations)
      pct=round(pct,digits=4)
      pct=pct*100
      lbls=paste(lbls,pct)
      lbls <- paste(pct,"%",sep="") # ad % to labels
      #pie(pct,labels = lbls, col=rainbow(length(lbls)),
      #main="Posterior probabilities of varying effect sizes")
      
      # Barplot
      b<-barplot(pct,
                 main = mtext("Posterior ROPE-analysis for"~delta[26]),
                 xlab = "ROPEs of standard effect sizes",
                 ylab = "Posterior-Percentage included inside ROPE",
                 names.arg = c("L-", "M-", "S-", "No Eff.", "S+", "M+", "L+"),
                 col = "cornflowerblue",
                 border="white",
                 horiz = FALSE,
                 ylim=c(0,107.5))
      b
      text(x=b,y=pct+4,labels=lbls)
      
      
      
      ####################### POSTERIOR ROPE ANALYSIS FOR DELTA_36 #######################################
      # Bar plot
      #effSizeCIValues = effectSize[{q<-rank(effectSize)/length(effectSize);q>=((1-ci)/2) & q <= (ci+(1-ci)/2)}]
      effSizeCI=quantile(effectSize36,probs=c(((1-ci)/2),(ci+(1-ci)/2)))
      effectSizeCIValues = effectSize36[effectSize36 >effSizeCI[1] & effectSize36 < effSizeCI[2]]
      largePosEffectIterations = length(which(effectSizeCIValues >= 0.8))/length(effectSizeCIValues)
      largeNegEffectIterations = length(which(effectSizeCIValues <= -0.8))/length(effectSizeCIValues)
      mediumPosEffectIterations = length(which(effectSizeCIValues >= 0.5 & effectSizeCIValues < 0.8))/length(effectSizeCIValues)
      mediumNegEffectIterations = length(which(effectSizeCIValues <= -0.5 & effectSizeCIValues > -0.8))/length(effectSizeCIValues)
      smallPosEffectIterations = length(which(effectSizeCIValues >= 0.2 & effectSizeCIValues < 0.5))/length(effectSizeCIValues)
      smallNegEffectIterations = length(which(effectSizeCIValues <= -0.2 & effectSizeCIValues > -0.5))/length(effectSizeCIValues)
      noEffectIterations = length(which(effectSizeCIValues < 0.2 & effectSizeCIValues > -0.2))/length(effectSizeCIValues)
      
      lbls=c("large positive:", "medium positive:", "small positive:", "no effect:", "small negative:", "medium negative:","large negative:")
      pct=c(largeNegEffectIterations,mediumNegEffectIterations,smallNegEffectIterations,noEffectIterations,smallPosEffectIterations,mediumPosEffectIterations,largePosEffectIterations)
      pct=round(pct,digits=4)
      pct=pct*100
      lbls=paste(lbls,pct)
      lbls <- paste(pct,"%",sep="") # ad % to labels
      #pie(pct,labels = lbls, col=rainbow(length(lbls)),
      #main="Posterior probabilities of varying effect sizes")
      b<-barplot(pct,
                 main = mtext("Posterior ROPE-analysis for"~delta[36]),
                 xlab = "ROPEs of standard effect sizes",
                 ylab = "Posterior-Percentage included inside ROPE",
                 names.arg = c("L-", "M-", "S-", "No Eff.", "S+", "M+", "L+"),
                 col = "cornflowerblue",
                 border="white",
                 horiz = FALSE,
                 ylim=c(0,107.5))
      b
      text(x=b,y=pct+4,labels=lbls)  
      
      
      
      ####################### POSTERIOR ROPE ANALYSIS FOR DELTA_46 #######################################
      # Bar plot
      #effSizeCIValues = effectSize[{q<-rank(effectSize)/length(effectSize);q>=((1-ci)/2) & q <= (ci+(1-ci)/2)}]
      effSizeCI=quantile(effectSize46,probs=c(((1-ci)/2),(ci+(1-ci)/2)))
      effectSizeCIValues = effectSize46[effectSize46 >effSizeCI[1] & effectSize46 < effSizeCI[2]]
      largePosEffectIterations = length(which(effectSizeCIValues >= 0.8))/length(effectSizeCIValues)
      largeNegEffectIterations = length(which(effectSizeCIValues <= -0.8))/length(effectSizeCIValues)
      mediumPosEffectIterations = length(which(effectSizeCIValues >= 0.5 & effectSizeCIValues < 0.8))/length(effectSizeCIValues)
      mediumNegEffectIterations = length(which(effectSizeCIValues <= -0.5 & effectSizeCIValues > -0.8))/length(effectSizeCIValues)
      smallPosEffectIterations = length(which(effectSizeCIValues >= 0.2 & effectSizeCIValues < 0.5))/length(effectSizeCIValues)
      smallNegEffectIterations = length(which(effectSizeCIValues <= -0.2 & effectSizeCIValues > -0.5))/length(effectSizeCIValues)
      noEffectIterations = length(which(effectSizeCIValues < 0.2 & effectSizeCIValues > -0.2))/length(effectSizeCIValues)
      
      lbls=c("large positive:", "medium positive:", "small positive:", "no effect:", "small negative:", "medium negative:","large negative:")
      pct=c(largeNegEffectIterations,mediumNegEffectIterations,smallNegEffectIterations,noEffectIterations,smallPosEffectIterations,mediumPosEffectIterations,largePosEffectIterations)
      pct=round(pct,digits=4)
      pct=pct*100
      lbls=paste(lbls,pct)
      lbls <- paste(pct,"%",sep="") # ad % to labels
      #pie(pct,labels = lbls, col=rainbow(length(lbls)),
      #main="Posterior probabilities of varying effect sizes")
      b<-barplot(pct,
                 main = mtext("Posterior ROPE-analysis for"~delta[46]),
                 xlab = "ROPEs of standard effect sizes",
                 ylab = "Posterior-Percentage included inside ROPE",
                 names.arg = c("L-", "M-", "S-", "No Eff.", "S+", "M+", "L+"),
                 col = "cornflowerblue",
                 border="white",
                 horiz = FALSE,
                 ylim=c(0,107.5))
      b
      text(x=b,y=pct+4,labels=lbls) 
      
      
      ####################### POSTERIOR ROPE ANALYSIS FOR DELTA_56 #######################################
      # Bar plot
      #effSizeCIValues = effectSize[{q<-rank(effectSize)/length(effectSize);q>=((1-ci)/2) & q <= (ci+(1-ci)/2)}]
      effSizeCI=quantile(effectSize56,probs=c(((1-ci)/2),(ci+(1-ci)/2)))
      effectSizeCIValues = effectSize56[effectSize56 >effSizeCI[1] & effectSize56 < effSizeCI[2]]
      largePosEffectIterations = length(which(effectSizeCIValues >= 0.8))/length(effectSizeCIValues)
      largeNegEffectIterations = length(which(effectSizeCIValues <= -0.8))/length(effectSizeCIValues)
      mediumPosEffectIterations = length(which(effectSizeCIValues >= 0.5 & effectSizeCIValues < 0.8))/length(effectSizeCIValues)
      mediumNegEffectIterations = length(which(effectSizeCIValues <= -0.5 & effectSizeCIValues > -0.8))/length(effectSizeCIValues)
      smallPosEffectIterations = length(which(effectSizeCIValues >= 0.2 & effectSizeCIValues < 0.5))/length(effectSizeCIValues)
      smallNegEffectIterations = length(which(effectSizeCIValues <= -0.2 & effectSizeCIValues > -0.5))/length(effectSizeCIValues)
      noEffectIterations = length(which(effectSizeCIValues < 0.2 & effectSizeCIValues > -0.2))/length(effectSizeCIValues)
      
      lbls=c("large positive:", "medium positive:", "small positive:", "no effect:", "small negative:", "medium negative:","large negative:")
      pct=c(largeNegEffectIterations,mediumNegEffectIterations,smallNegEffectIterations,noEffectIterations,smallPosEffectIterations,mediumPosEffectIterations,largePosEffectIterations)
      pct=round(pct,digits=4)
      pct=pct*100
      lbls=paste(lbls,pct)
      lbls <- paste(pct,"%",sep="") # ad % to labels
      #pie(pct,labels = lbls, col=rainbow(length(lbls)),
      #main="Posterior probabilities of varying effect sizes")
      b<-barplot(pct,
                 main = mtext("Posterior ROPE-analysis for"~delta[56]),
                 xlab = "ROPEs of standard effect sizes",
                 ylab = "Posterior-Percentage included inside ROPE",
                 names.arg = c("L-", "M-", "S-", "No Eff.", "S+", "M+", "L+"),
                 col = "cornflowerblue",
                 border="white",
                 horiz = FALSE,
                 ylim=c(0,107.5))
      b
      text(x=b,y=pct+4,labels=lbls) 
    }
  }
}


assumption.check = function(x1,x2,x3,x4=NULL,x5=NULL,x6=NULL,conf.level=0.95){
  check=TRUE;
  if(shapiro.test(x1)$p.value < 1-conf.level | shapiro.test(x2)$p.value < 1-conf.level | shapiro.test(x3)$p.value < 1-conf.level){
    check=FALSE;
  }
  
  if(!is.null(x4)){
    if(shapiro.test(x4)$p.value < 1-conf.level){
      check=FALSE;
    }
  }
  
  if(!is.null(x5)){
    if(shapiro.test(x5)$p.value < 1-conf.level){
      check=FALSE; 
    }
  }
  
  if(!is.null(x6)){
    if(shapiro.test(x6)$p.value < 1-conf.level){
      check=FALSE;
 }
  }
  
  if(check==TRUE){
    message("Model assumptions checked. No significant deviations from normality detected. Bayesian ANOVA can be run safely.")
  } else {
    warning("Model assumption of normally distributed data in each group is violated.\n All results of the Bayesian ANOVA based on a three-component Gaussian mixture could therefore be unreliable and not trustworthy.")
    warning("Run further diagnostics (like Quantile-Quantile-plots) to check if the Bayesian ANOVA can be expected to be robust to the violations of normality\n")
  }
  # Set 2x3 par
  dev.new()
  newpar <- par(mfrow=c(2,3))
  # Effect size posterior histogram
  # Make a vector of values to draw ticks at:
  hist(x1,freq=FALSE,col="cornflowerblue",border="white",xlab="First group",main="Histogram of first group")
  lines(density(x1))
  hist(x2,freq=FALSE,col="cornflowerblue",border="white",xlab="Second group",main="Histogram of second group")
  lines(density(x2))
  hist(x3,freq=FALSE,col="cornflowerblue",border="white",xlab="Third group",main="Histogram of third group")
  lines(density(x3))
  qqnorm(x1, col="cornflowerblue", main="Q-Q Plot for first group"); qqline(x1, col = "black")
  qqnorm(x2, col="cornflowerblue", main="Q-Q Plot for second group"); qqline(x2, col = "black")
  qqnorm(x3, col="cornflowerblue", main="Q-Q Plot for third group"); qqline(x3, col = "black")
  
  
  if(length(x4)>1 && length(x5)>1 && length(x6)>1){
    dev.new()
    newpar <- par(mfrow=c(2,3))
    hist(x4,freq=FALSE,col="cornflowerblue",border="white",xlab="Fourth group",main="Histogram of fourth group")
    lines(density(x4))
    hist(x5,freq=FALSE,col="cornflowerblue",border="white",xlab="Fifth group",main="Histogram of fifth group")
    lines(density(x5))
    hist(x6,freq=FALSE,col="cornflowerblue",border="white",xlab="Sixth group",main="Histogram of sixth group")
    lines(density(x6))
    qqnorm(x4, col="cornflowerblue", main="Q-Q Plot for fourth group"); qqline(x4, col = "black")
    qqnorm(x5, col="cornflowerblue", main="Q-Q Plot for fifth group"); qqline(x5, col = "black")
    qqnorm(x6, col="cornflowerblue", main="Q-Q Plot for sixth group"); qqline(x6, col = "black")
  } else if(length(x4)>1 && length(x5)>1 && is.null(x6)){
    dev.new()
    newpar <- par(mfrow=c(2,2))
    hist(x4,freq=FALSE,col="cornflowerblue",border="white",xlab="Fourth group",main="Histogram of fourth group")
    lines(density(x4))
    hist(x5,freq=FALSE,col="cornflowerblue",border="white",xlab="Fifth group",main="Histogram of fifth group")
    lines(density(x5))
    qqnorm(x4, col="cornflowerblue", main="Q-Q Plot for fourth group"); qqline(x4, col = "black")
    qqnorm(x5, col="cornflowerblue", main="Q-Q Plot for fifth group"); qqline(x5, col = "black")
  } else if(length(x4)>1 && is.null(x5) && is.null(x6)){
    dev.new()
    newpar <- par(mfrow=c(2,1))
    hist(x4,freq=FALSE,col="cornflowerblue",border="white",xlab="Fourth group",main="Histogram of fourth group")
    lines(density(x4))
    qqnorm(x4, col="cornflowerblue", main="Q-Q Plot for fourth group"); qqline(x4, col = "black")
  }
}

post.fit = function(dataframe,x1,x2,x3,x4=NULL,x5=NULL,x6=NULL){
  
  x1Seq=seq(min(x1),max(x1),length.out = 100)
  x2Seq=seq(min(x2),max(x2),length.out = 100)
  x3Seq=seq(min(x3),max(x3),length.out = 100)
  
  if(length(x4)>1){
    x4Seq=seq(min(x4),max(x4),length.out = 100)
  }
  if(length(x5)){
    x5Seq=seq(min(x5),max(x5),length.out = 100)
  }
  if(length(x6)){
    x6Seq=seq(min(x6),max(x6),length.out = 100)
  }
  
  if(length(x4)>1 && length(x5)>1 && length(x6)>1){
    hist(c(x1,x2,x3,x4,x5,x6),main="Posterior fit",xlab="Original data",ylab="Posterior density",freq=FALSE,ylim=c(0,1))
    lines(x1Seq,dnorm(x1Seq,mean=mean(dataframe$mu1),sd=mean(dataframe$sigma1Sq))/6,col="blue")
    lines(x2Seq,dnorm(x2Seq,mean=mean(dataframe$mu2),sd=mean(dataframe$sigma2Sq))/6,col="red")
    lines(x3Seq,dnorm(x3Seq,mean=mean(dataframe$mu3),sd=mean(dataframe$sigma3Sq))/6,col="purple")
    lines(x4Seq,dnorm(x4Seq,mean=mean(dataframe$mu4),sd=mean(dataframe$sigma4Sq))/6,col="orange")
    lines(x5Seq,dnorm(x5Seq,mean=mean(dataframe$mu5),sd=mean(dataframe$sigma5Sq))/6,col="green")
    lines(x6Seq,dnorm(x6Seq,mean=mean(dataframe$mu6),sd=mean(dataframe$sigma6Sq))/6,col="cornflowerblue")
  }
  if(length(x4)>1 && length(x5)>1 && is.null(x6)){
    hist(c(x1,x2,x3,x4,x5),main="Posterior fit",xlab="Original data",ylab="Posterior density",freq=FALSE,ylim=c(0,1))
    lines(x1Seq,dnorm(x1Seq,mean=mean(dataframe$mu1),sd=mean(dataframe$sigma1Sq))/5,col="blue")
    lines(x2Seq,dnorm(x2Seq,mean=mean(dataframe$mu2),sd=mean(dataframe$sigma2Sq))/5,col="red")
    lines(x3Seq,dnorm(x3Seq,mean=mean(dataframe$mu3),sd=mean(dataframe$sigma3Sq))/5,col="purple")
    lines(x4Seq,dnorm(x4Seq,mean=mean(dataframe$mu4),sd=mean(dataframe$sigma4Sq))/5,col="orange")
    lines(x5Seq,dnorm(x5Seq,mean=mean(dataframe$mu5),sd=mean(dataframe$sigma5Sq))/5,col="green")
  }
  if(length(x4)>1 && is.null(x5) && is.null(x6)){
    hist(c(x1,x2,x3,x4),main="Posterior fit",xlab="Original data",ylab="Posterior density",freq=FALSE,ylim=c(0,1))
    lines(x1Seq,dnorm(x1Seq,mean=mean(dataframe$mu1),sd=mean(dataframe$sigma1Sq))/4,col="blue")
    lines(x2Seq,dnorm(x2Seq,mean=mean(dataframe$mu2),sd=mean(dataframe$sigma2Sq))/4,col="red")
    lines(x3Seq,dnorm(x3Seq,mean=mean(dataframe$mu3),sd=mean(dataframe$sigma3Sq))/4,col="purple")
    lines(x4Seq,dnorm(x4Seq,mean=mean(dataframe$mu4),sd=mean(dataframe$sigma4Sq))/4,col="orange")
  }
  if(is.null(x4) && is.null(x5) && is.null(x6)){
    hist(c(x1,x2,x3),main="Posterior fit",xlab="Original data",ylab="Posterior density",freq=FALSE,ylim=c(0,1))
    lines(x1Seq,dnorm(x1Seq,mean=mean(dataframe$mu1),sd=mean(dataframe$sigma1Sq))/3,col="blue")
    lines(x2Seq,dnorm(x2Seq,mean=mean(dataframe$mu2),sd=mean(dataframe$sigma2Sq))/3,col="red")
    lines(x3Seq,dnorm(x3Seq,mean=mean(dataframe$mu3),sd=mean(dataframe$sigma3Sq))/3,col="purple")
  }
}

