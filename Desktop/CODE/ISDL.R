#
# This master R file includes all the R functions for implementing 
#Innovated Scalable Dynamic Learning for Graphical Models (ISDL)
#
# The ISDLcode was rewritten by Liwan Li(emmaly@mail.ustc.edu.cn), April 10, 2020
#on the basis of ISEE Written by Yingying Fan (fanyingy@marshall.usc.edu) and Jinchi Lv 
# (jinchilv@marshall.usc.edu) 
# 
#
# List of R functions for ISDL:
# Part I: main functions
#   1) isdl: main function calculating hat matrix for Omega
#   2) isdl.X: function calculating hat matrix for X tilde
#   3) slasso: function implementing tuning-free sparse regression with scaled Lasso
# Part II: additional functions
#   1) beta.block: function calculating beta hat matrix for each small block
#   2) isdl.cv: cv function selecting threshold tau
# Part III: functions summarizing results and generating data matrix
#   1) rate.func: function summarizing selection results
#   2) make.data1: function generating data matrix from band Omega
#	3) make.data2: function generating data matrix from block diagonal Omega
#   4) band.mat: function generating band Omega
#   5) blk.mat: function generating block diagonal Omega
#

# Part I: main functions
# Part I 1) isdl: main function calculating hat matrix for Omega^{t}
isdl<- function(X, regfactor = "log", npermu = 1,t, sis.use = 1, bia.cor = 0){
  n = nrow(X)
  p = ncol(X)
  
  if (regfactor == "log") {
    reg.fac = log(n)
  }
  if (regfactor == "one") {
    reg.fac = 1
  }
  if (regfactor == "sqrt") {
    reg.fac = sqrt(n)
  }
  
  Omega.isdl = matrix(0, p, p)
  # bias corrected isdl estimator
  Omega.isdl.c = matrix(0, p, p)
  temp = matrix(0, p, p)
  temp1 = matrix(0, p, p)
  for (i in 1:npermu) {
    if (i == 1) {
      permu = c(1:p)
    } else {
      permu = sample.int(p,size=p)
    }
    
    obj = isdl.X(X, permu, reg.fac,t, sis.use, bia.cor)
    X.tilde = obj$X.tilde
    
    obj1 = isdl.cv(X.tilde)
    Omega.isdl = Omega.isdl + obj1$Omega.isdl
    temp = temp + (obj1$Omega.isdl != 0)
    
    if (bia.cor == 1){
      # bias corrected isdl estimator with correction on its support
      Omega.isdl.c0 = obj1$Omega.isdl + (obj1$Omega.isdl != 0)*(abs(cov2cor(obj$biascor.mat)) > obj1$tau)*(obj$biascor.mat - obj1$Omega.isdl)/2
      Omega.isdl.c = Omega.isdl.c + Omega.isdl.c0
      temp1 = temp1 + (Omega.isdl.c0 != 0)
    }
  }
  
  	Omega.isdl = Omega.isdl/(temp + matrix(1e-8, p, p))
  	if (bia.cor == 1){
    Omega.isdl.c = (Omega.isdl != 0)*Omega.isdl.c/(temp1 + matrix(1e-8, p, p))
  	}
  
  	if (bia.cor == 0){
    obj=list(Omega.isdl=Omega.isdl)
  	}
  	if (bia.cor == 1){
    obj=list(Omega.isdl=Omega.isdl, Omega.isdl.c=Omega.isdl.c)
 	 }
  return(obj)
}
###########################################################################

# Part I 2) isdl.X: function calculating hat matrix for X tilde
isdl.X <- function(X, permu, reg.fac, t,sis.use = 0, bia.cor = 0){
  # detect the OS
  if (Sys.info()['sysname'] == 'Windows') {
    windows.use = 1
  }else{
    windows.use = 0}
  
  X = X[,permu] 
  # K is the blocksize
  K <-2 #or 3,4,5
  p <- ncol(X)
  n <- nrow(X)
  # bias corrected matrix
  biascor.mat <- matrix(0, p, p)
  # beta coefficient matrix
  beta.mat <- matrix(0, p, p)
  
  gata2 <- sqrt(n)/(log(p)*2*p)
  B2 <- qt(1-gata2, n-1)
  lam.coef<-1
  lambda <- lam.coef*B2/sqrt(n-1+B2^2)
  #####reweighted matrix
    h=5.848/(n^(1/3)) #h<-geth() from a cv
    w=c()
    wt=c()
    for (i in 1:n){
    u=(i-t)/n*h
    w[i]=exp(-0.5*u^2)
   }
    wt=w/sum(w)
    Xs=matrix(0,n,p)
  for (ii in 1:n){
     Xs[ii,]=X[ii,]*sqrt(wt[ii])/sqrt(wt[t])
   }
  #####

  stan_X = scale(Xs, center=F, scale=T)/sqrt(n-1)
  y.norm = attr(stan_X,"scale")*sqrt(n-1) 
  i.index <- seq(1,p-1,by=K)
  
  X.tilde = matrix(0, nrow=n, ncol=p)
  for (k in 1:length(i.index)){ 
    i = i.index[k]
    # for K=2
    if(i==p-2){
      j=c(p-1,p)
    }else{
      j=i+1
    }
    
    block.ind = c(i,j)
    temp= beta.block(stan_X, block.ind, lambda, reg.fac, windows.use, sis.use)
    beta.coef = temp$beta.j
    
    temp2 = scale(stan_X[,-block.ind] %*% beta.coef, center=F, scale=1/y.norm[block.ind])
    epsi <- Xs[, block.ind] - temp2
    omega <- solve(t(epsi)%*%epsi/n)
    X.tilde[,block.ind] =  as.matrix(epsi%*%omega)
    
    if (bia.cor == 1){
      biascor.mat[block.ind,block.ind] = as.matrix(omega)
      beta.mat[block.ind,-block.ind] = t(scale(beta.coef, center=F, scale=1/y.norm[block.ind]))
    }
  }

  if (bia.cor == 1){
    beta.mat = scale(beta.mat, center=F, scale=y.norm)
    Omega.ini = t(X.tilde)%*%X.tilde/n
    #for K=2
    for (ii in 1:(length(i.index)-1)){
      for (jj in (ii+1):length(i.index)){
        i = i.index[ii]
        indset1 = i:(i+K-1)
        j = i.index[jj]
        indset2 = j:(j+K-1)
        
        biascor.mat[indset1,indset2] = -(Omega.ini[indset1,indset2] + t(beta.mat[indset1,indset2])%*%Omega.ini[indset1,indset1] + t(beta.mat[indset2,indset1])%*%Omega.ini[indset2,indset2])
        biascor.mat[indset2,indset1] = t(biascor.mat[indset1,indset2])
      }
    }
    biascor.mat = biascor.mat[order(permu),order(permu)]
  }
  
  X.tilde = X.tilde[,order(permu)]
  if (bia.cor == 0){
    obj=list(X.tilde=X.tilde)
  }
  if (bia.cor == 1){
    obj=list(X.tilde=X.tilde, biascor.mat=biascor.mat)
  }
  return(obj)
}
###########################################################################

# Part I 3) slasso: function implementing tuning-free sparse regression with scaled Lasso

# Windows version; Mac or Linux version
# Lasso implemented with coordinate descent
# each column of n x p design matrix X is rescaled to have unit L2-norm

slasso <- function(X, y, lam, betap, windows.use = 1, maxite = 50, tol = 1e-2){
  # detect the OS
  #if (Sys.info()['sysname'] == 'Windows') {
  #  windows.use = 1
  #}else{
  #  windows.use = 0}
  
  nr = nrow(X)
  nc = ncol(X)
  if (windows.use) {
    dyn.load("slasso.dll")
  } else {
    dyn.load("slasso.so")
  }
  # for output variables, use outvar = some argument
  out = .C("slasso", as.vector(as.numeric(X)), as.integer(nr), as.integer(nc), as.double(y), as.double(lam), as.integer(maxite), as.double(tol), betap = as.double(betap))
  # take out betap component of the returned list
  if (is.nan(out$betap[1])) stop("NaN values")
  return(out$betap)
}
###########################################################################


###########################################################################
# Part II: additional functions
# Part II 1) beta.block: function calculating beta hat matrix for each small block

beta.block <- function(X, block.ind, lambda, reg.fac, windows.use, sis.use){
  p <- ncol(X)
  n <- nrow(X)
  tol = 1e-3
  block.size = length(block.ind)
  precision <- matrix(0, ncol = block.size , nrow = block.size)
  
  beta.j<- matrix(0, nrow = p-block.size, ncol = block.size)
  
  for (j in 1:block.size){
    if (sis.use) {
      org.set = c(1:p)[ -block.ind]
      set0<-order(abs(cor(X[, -block.ind], X[, block.ind[j]])),decreasing=T)
      sis.size = min(floor(n/log(n)),p-block.size)
      sis.set0<-sort(set0[1:sis.size])
      sis.set= org.set[sis.set0]
    }else {
      sis.set = c(1:p)[ -block.ind]
      sis.set0 = c(1:length(sis.set)) 
    }
    
    sigma.old <- 1
    beta.old <- c(rep(0, length(sis.set)))
    last_beta_diff <- 0
    inner <- 0
    while (inner < 5) {
      reg <- sigma.old*lambda*reg.fac
      beta <- slasso(as.matrix(X[, sis.set]), as.matrix(X[, block.ind[j]]),reg, c(rep(0, length(sis.set))), windows.use)
      if (reg.fac == 1) {
        beta = beta*(abs(beta) > 1*lambda)
      }else{
        beta = beta*(abs(beta) > 0.5*lambda)
      }
      
      sel.mod = which(beta != 0)
      beta[sel.mod] = slasso(as.matrix(X[, sis.set[sel.mod]]), as.matrix(X[, block.ind[j]]), reg, c(rep(0, length(sel.mod))), windows.use)
      beta_diff <- max(abs(beta - beta.old))
      A <- X[, block.ind[j]] - X[, sis.set] %*% beta
      sigma <- sqrt(sum(A^2))/sqrt(max(n-sum(beta!=0),0)+1e-20)
      sigma_diff <- abs(sigma - sigma.old)
      
      if (sigma_diff < tol) 
        break
      else if (sigma_diff < tol * 0.1& 
               abs(beta_diff - last_beta_diff) < tol) 
        break
      else {
        sigma.old <- sigma
        beta.old <- beta
        last_beta_diff <- beta_diff
      }
      inner <- inner + 1
    }
    beta.j[sis.set0,j] <- as.vector(beta)
  }
  obj=list(beta.j=beta.j)
  return(obj)
}
###########################################################################

# Part II 2) isdl.cv: cv function selecting threshold tau

isdl.cv <- function(X.tilde, ntau = 20, split.ratio = 0.9, n.split=5, criterion="Frob"){
  n=nrow(X.tilde)
  p=ncol(X.tilde)
  n.train = ceiling(n*split.ratio)
  n.test = n-n.train
  
  tau.min=0.5 #tau.min<-gettau.min()
  tau.max=0.8#tau.max<-gettau.max() 
  tau.path = seq(tau.min,tau.max,length.out=ntau)*sqrt(log(p)/n)
  
  # loss function used in calculating prediction error
  loss.func <- function(Sig, Sig.tau, criterion){
    if(criterion=="Frob")err= sum((Sig-Sig.tau)^2)
    if(criterion=="Infi")err= max(abs(Sig-Sig.tau))
    if(criterion=="L1")err = sum(abs(Sig-Sig.tau))
    if(criterion=="L1L1")err = sum(abs(Sig-Sig.tau))+0.5*sum(abs(Sig.tau))
    return(err)
  }
  
  err = matrix(0, n.split, length(tau.path))
  for(j in 1:n.split){
    train.id = sample.int(n, size = n.train, replace = FALSE)
    test.id = setdiff(c(1:n), train.id)
    Sig.train = cor(X.tilde[train.id,])
    Sig.test = cor(X.tilde[test.id,])
    
    for(i in 1:length(tau.path)){
      Sig.train.tau = Sig.train*(abs(Sig.train)>tau.path[i])
      err[j,i] = loss.func(Sig.test,Sig.train.tau, criterion)
    }
  }
  
  pe.vec = apply(err, 2, mean)
  tau.id = which.min(pe.vec)
  
  Omega.ini = (t(X.tilde)%*%X.tilde/n)
  Omega.isdl = Omega.ini*(abs(cov2cor(Omega.ini)) > tau.path[tau.id])
  
  obj=list(Omega.isdl=Omega.isdl, tau=tau.path[tau.id])
  return(obj)
}


###########################################################################
# Part III: functions summarizing results and generating data matrix
# Part III 1) rate.func: function summarizing selection results
rate.func <- function(Omega.hat, Omega){
  TP <- sum((Omega.hat!=0)*(Omega!=0))
  FP <- sum((Omega.hat!=0)*(Omega==0))
  TN <- sum((Omega.hat == 0)*(Omega==0))
  FN <- sum((Omega.hat == 0)*(Omega!=0))
  
  TPR <- TP/(TP+FN)
  FPR <- FP/(FP+TN) 
  PRE <- TP/(TP+FP)
  F1<-2*TPR*PRE/(PRE+TPR)
  obj <-list(TPR=TPR, FPR=FPR,PRE=PRE,F1=F1)
  return(obj)
}

# Part III 3) make.data1: function generating data matrix from band Omega
# acoording to the setting of simulation example 1
make.data1 <- function(n,p,initial,nedge,oedge,max=0.3,min=0.1,re=0.25, seed){
  	set.seed(seed)  
    initial=p;nedge=p/10;oedge=p/10
    mu=c();
    for(j in 1:p)
    {
    	mu[j]=0;
    }
    X=c();
    omega.final=array(0,c(p,p,n))
    sigma.final=array(0,c(p,p,n))
    omega.final=band.mat(n,p,initial,nedge,oedge,max=0.3,min=0.1,re=0.25)
    for(i in 1:n){
    omega=omega.final[,,i]
    sigma=solve(omega)
    sigma.final[,,i]=sigma
    X=rbind(X,mvrnorm(n=1, mu, sigma))
	}
	obj<-list(X=X,omega.final=omega.final,sigma.final=sigma.final)
	return(obj)
}

# Part III 3) make.data1: function generating data matrix from block diagonal Omega
# or acoording to the the setting of simulation example 2
#n = 200；p = 2000；obj=make.data2(n,p,seed=1000) 
make.data2 <- function(n,p,dim=20,initial=20,oedge=2,nedge=2,min=0.1,max=0.3,seed){ 
	set.seed(seed) 
  	ok=n*p/20
 	num=p/20
	omega.final=array(0,c(p,p,n))
	sigma.final=array(0,c(p,p,n))
	medomega=array(0,c(20,20,ok))
	dim=20
	for(xy in 1:num){
	  ll=(xy-1)*n+1
	  ww=xy*n
 	  medomega[,,ll:ww]=band.mat(n=200,p=dim,initial=dim,nedge=dim/10,oedge=dim/10,max=0.3,min=0.1,re=0.25)
	}
	dif=c(1:p)
	permu=sample(dif,length(dif))
 	for(t in 1:n){
    obj=blk.mat(t,n,num,permu,medomega)
    yxy=obj$Omega
    sigma.final=obj$Sigma.half
    omega.final[,,t]=as.numeric(unlist(yxy))
    X = rbind(X,matrix(rnorm(p),1,p)%*%(obj$Sigma.half))	
    }
  obj<-list(X=X,omega.final=omega.final,sigma.final=sigma.final)
  return(obj)
}



###########################################################################
# Part III 4) band.mat: function generating band Omega
band.mat<-function(n=200,p=50,initial=50,nedge=5,oedge=5,max=0.3,min=0.1,re=0.25){
  omega0=diag(rep(re, p))
  ####pick up 50 edges  
  edge=matrix(0,initial,2)
  e=c()
  edge[1,]=sample(1:p,2)
  e[1]=sum(exp(edge[1,]))
  i<-2
  while(i<=initial){
    edge[i,]=sample(1:p,2)
    e[i]=sum(exp(edge[i,]))
    if(e[i] %in% e[1:(i-1)]){
      next
    }else{i=i+1}
  }
  a=c()
  a=runif(initial,min,max)
  for (i in 1:initial){
    omega0[edge[i,1],edge[i,2]]=-a[i]
    omega0[edge[i,2],edge[i,1]]=-a[i]
    omega0[edge[i,1],edge[i,1]]=omega0[edge[i,1],edge[i,1]]+a[i]
    omega0[edge[i,2],edge[i,2]]=omega0[edge[i,2],edge[i,2]]+a[i]
  }
  
  d1=sample(1:initial,oedge)
  edge1=matrix(0,oedge,2)
  edge2=matrix(0,nedge,2)
  for(i in 1:oedge){
    edge1[i,]=edge[d1[i],]#####select 5 old edges
  }
  ####select 5 new edges
  ii<-1
  while(ii<=nedge){
    edge2[ii,]=sample(1:p,2)
    e[initial+ii]=sum(exp(edge2[ii,]))
    if(e[initial+ii] %in% e[1:initial+ii-1]){
      next
    }else{ii=ii+1}
  }
  a2=c()
  a2=runif(nedge,min,max)
  omega.final=array(0,c(p,p,n))
  #####rule 
  for (j in 1:n){
    for(i in 1:nedge){
      omega0[edge2[i,1],edge2[i,2]]=omega0[edge2[i,1],edge2[i,2]]-a2[i]/n 
      omega0[edge2[i,2],edge2[i,1]]=omega0[edge2[i,2],edge2[i,1]]-a2[i]/n
      omega0[edge1[i,1],edge1[i,2]]=omega0[edge1[i,1],edge1[i,2]]+a[d1[i]]/n 
      omega0[edge1[i,2],edge1[i,1]]=omega0[edge1[i,2],edge1[i,1]]+a[d1[i]]/n
      omega0[edge2[i,1],edge2[i,1]]=omega0[edge2[i,1],edge2[i,1]]+a2[i]/n
      omega0[edge2[i,2],edge2[i,2]]=omega0[edge2[i,2],edge2[i,2]]+a2[i]/n
      omega0[edge1[i,1],edge1[i,1]]=omega0[edge1[i,1],edge1[i,1]]-a[d1[i]]/n
      omega0[edge1[i,2],edge1[i,2]]=omega0[edge1[i,2],edge1[i,2]]-a[d1[i]]/n
    }
    omega.final[,,j]=omega0
  }
  return(omega.final)
}

###########################################################################
# Part III 5) blk.mat: function generating block diagonal Omega
blk.mat <- function(t,n,num,permu, medomega){
	Omega.blk =matrix(0,20,20)
	Omega.blk = medomega[,,t]
    Sigma.blk = solve(Omega.blk)
    Sigma.half.blk =  chol(Sigma.blk)
    Omega = Omega.blk
    Sigma = Sigma.blk
    Sigma.half = Sigma.half.blk
    for(dada in 1:(num-1)){
    	yy=t+dada*n
    	Omega.blk = medomega[,,yy]
    	Sigma.blk = solve(Omega.blk)
    	Sigma.half.blk =  chol(Sigma.blk)
   
   	    Omega=bdiag(Omega, Omega.blk)
   	    Sigma = bdiag(Sigma, Sigma.blk)
        Sigma.half = bdiag(Sigma.half, Sigma.half.blk)
  	}
    Sigma = Sigma[permu,permu]
    Omega = Omega[permu,permu]
    Sigma.half = Sigma.half[permu,permu]
    obj = list(Sigma = Sigma, Omega=Omega, Sigma.half = Sigma.half)
    return(obj)
}

################################## end ######################################
