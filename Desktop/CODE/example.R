#
# This is an example of applying ISDL.
#
# This master R file includes all the R functions for implementing 
#Innovated Scalable Dynamic Learning for Graphical Models (ISDL)
#
# The ISDLcode was rewritten by Liwan Li(emmaly@mail.ustc.edu.cn), April 10, 2020
# on the basis of ISEE Written by Yingying Fan (fanyingy@marshall.usc.edu) and Jinchi Lv 
# (jinchilv@marshall.usc.edu) 
#
# Reference:
# Zheng, Z., Li, L. ,Zhou, J.,Kong, Y.(2020) Innovated Scalable Dynamic Learning for
#Graphical Models. Preprint submitted to Statistics and Probability Letters.
# Fan, Y. and Lv, J. (2015). Innovated scalable efficient estimation in ultra-large 
# Gaussian graphical models. The Annals of Statistics.
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

## slasso.so or slasso.dll
if (windows.use) {
  dyn.load("slasso.dll")
} else {
  dyn.load("slasso.so")
}

# load packages
library(MASS) 
library(Matrix)

# simulation parameters
n = 200
p = 100
Omega.final=array(0,c(p,p,n))
Sigma.final=array(0,c(p,p,n))
obj= make.data1(n,p,initial,nedge,oedge,max,min,re,seed=1000)
Omega.final=obj$omega.final
Sigma.final= obj$sigma.final
X.mat= obj$X

# initialization parameters
regfactor = "log"  # or "one", "sqrt"
npermu = 1         # or >= 2
sis.use = 1        # or 0, whether to use SIS for screening
bia.cor = 1        # or 0, whether to apply bias correction for ISDL

# ISDL 
Omega.hat.c=array(0,c(p,p,n))
t = proc.time()
for(t in 1:n){
Omega.hat0= isdl(X.mat, regfactor, npermu,t,sis.use, bia.cor)$Omega.isdl.c
Omega.hat0=as.numeric(unlist(Omega.hat0))
Omega.hat.c[,,t]=matrix(Omega.hat0,ncol=p)
}
t.c = proc.time()-t
T=t.c/n

# squared Frobenius norm
# ISDL estimator with bias correction
sum((Omega.hat.c[,,t] - Omega.final[,,t])^2)/sum(Omega.final[,,t]^2)
##t in (1,...,n)
