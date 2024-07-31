####### Lasso implementation for a capture-recapture model
####### using white stork data as a case study
####### check out this paper http://bit.ly/2ik437S

# load detection/non-detection capture-recapture data and covariates
load('whitestork.RData')

# initial run with all covariates to compute weights for the adaptive lasso
library(RMark)
storks=convert.inp("cigognes")
storks.processed=process.data(storks)
storks.ddl=make.design.data(storks.processed)
storks.ddl$Phi$x1=0
storks.ddl$Phi$x2=0
storks.ddl$Phi$x3=0
storks.ddl$Phi$x4=0
storks.ddl$Phi$x5=0
storks.ddl$Phi$x6=0
storks.ddl$Phi$x7=0
storks.ddl$Phi$x8=0
storks.ddl$Phi$x9=0
storks.ddl$Phi$x10=0
for (i in 1:nrow(cov)){
storks.ddl$Phi$x1[storks.ddl$Phi$time==i]=cov[i,1]
storks.ddl$Phi$x2[storks.ddl$Phi$time==i]=cov[i,2]
storks.ddl$Phi$x3[storks.ddl$Phi$time==i]=cov[i,3]
storks.ddl$Phi$x4[storks.ddl$Phi$time==i]=cov[i,4]
storks.ddl$Phi$x5[storks.ddl$Phi$time==i]=cov[i,5]
storks.ddl$Phi$x6[storks.ddl$Phi$time==i]=cov[i,6]
storks.ddl$Phi$x7[storks.ddl$Phi$time==i]=cov[i,7]
storks.ddl$Phi$x8[storks.ddl$Phi$time==i]=cov[i,8]
storks.ddl$Phi$x9[storks.ddl$Phi$time==i]=cov[i,9]
storks.ddl$Phi$x10[storks.ddl$Phi$time==i]=cov[i,10]
}
PhiCov=list(formula=~x1+x2+x3+x4+x5+x6+x7+x8+x9+x10)
pdot=list(formula=~1)
phixp = mark(storks.processed,storks.ddl,model.parameters=list(Phi=PhiCov,p=pdot))
phixp$results$beta

# compute weights for adaptive lasso
a.weights = c(0,1/abs(phixp$results$beta$estimate[2:11]),0) # intercept and logit(p) have 0 weights

# specify various quantities
nkp1 <- dim(data)[2]
nh <- dim(data)[1] # nb of individuals
k <- nkp1-1 # nb of capture occasions
km1 <- k-1 
eff <- data[,nkp1] # counts
# select ecounter histories only - remove counts
data <- data[,1:k]
# compute the date of first capture fc, and state at initial capture init.state
fc <- NULL
init.state <- NULL
for (i in 1:nh){
temp <- 1:k
fc <- c(fc,min(temp[data[i,]==1]))
init.state <- c(init.state,data[i,fc[i]])
}

# form design matrix for survival (detection is constant)
X = cbind(1,cov)

apply(X,2,sum)
apply(X,2,sd)
cor(cov(X))

##################
### Lasso path ###
##################

maxlamb = 2
penalty.min = 1.e-3
n.fits = 50

penalty = exp(seq(log(penalty.min), log(maxlamb), length.out = n.fits))
lassopath = matrix(NA, ncol(X)+1, n.fits)
lassolike = rep(NA, n.fits)
par.init  = rep(0, ncol(X)+1)
data <- t(data) # transpose data
source('wtlasso_devCJSphicovp_intercept.R') # deviance capture-recapture model with lasso regularisation

no.shrink = c(1,ncol(X)+1)

for (i in 1:n.fits)
{
	lasso.vec = rep(penalty[i], length(par.init))*a.weights
	lasso.vec[no.shrink] = 0
	binit <- rep(0,ncol(X)+1) # init values
	fit.lasso = optim(binit,wtlasso_devCJSphicovp_intercept,NULL,hessian=FALSE,data,X,eff,fc,init.state,
	nh,km1,lasso=lasso.vec,no.shrink=no.shrink,method="BFGS",control=list(trace=1, REPORT=1,reltol=1e-12))
	lassopath[,i] = fit.lasso$par
	lassolike[i]  = fit.lasso$value
	cat(paste(i, "\n"))
	flush.console()
}

save(lassopath,lassolike,file='ws_lasso.RData')
load('ws_lasso.RData')
N = nh # to be checked
num.v = apply(abs(lassopath) > 1.e-3, 2, sum) - length(no.shrink)
BIC = 2*lassolike + log(N)*num.v
#BIC = 2*lassolike + 2*num.v

setEPS()
postscript("bic.eps")
plot(penalty, BIC, log = "x", type = "o", xlab = "Lasso Penalty", ylab = "BIC")
dev.off()

pen = penalty[which.min(BIC)] # get 'best' penalisation
lasso.vec = rep(pen, length(par.init))
lasso.vec[no.shrink] = 0

source('All Functions.R')
# Regularization Path for CR model and white stork data
setEPS()
postscript("path.eps")
vis.path(lassopath[-no.shrink,], penalty, bic = pen, main = "")
dev.off()

# get estimates from model with selected lasso penalisation 
lasso.occ = lassopath[,which.min(BIC)]

# display results
x <- lasso.occ
cov.names = names(X)
cbind(x[1:ncol(X)],cov.names)

                              cov.names   
 [1,] "0.677197390019177"     "1"         
 [2,] "7.47114934885098e-05"  "Diourbel"  
 [3,] "-2.99586996679983e-05" "Gao"       
 [4,] "0.000128869648004468"  "Kayes"     
 [5,] "0.241027971579547"     "Kita"      
 [6,] "-0.00012945774462468"  "Maradi"    
 [7,] "0.000351596996537262"  "Mopti"     
 [8,] "-5.91219509372506e-05" "Ouahigouya"
 [9,] "1.6626341903615e-05"   "Segou"     
[10,] "0.000116112255058697"  "Tahoua"    
[11,] "-0.000267682591336625" "Tombouctou"

matrix(x[(ncol(X)+1)])
         [,1]
[1,] 2.371982

# survival probs
phi <- exp(as.matrix(X) %*% x[1:ncol(X)])/(1+exp(as.matrix(X) %*% x[1:ncol(X)]))
# detection probs
p <- exp(x[(ncol(X)+1)])/(1+exp(x[(ncol(X)+1)]))

phi
           [,1]
 [1,] 0.7043902
 [2,] 0.7627816
 [3,] 0.7166949
 [4,] 0.6552528
 [5,] 0.6625752
 [6,] 0.6365405
 [7,] 0.7484468
 [8,] 0.6167897
 [9,] 0.6621625
[10,] 0.6707916
[11,] 0.6238433
[12,] 0.6892665
[13,] 0.6288057
[14,] 0.6141608
[15,] 0.6075211
[16,] 0.5771152

p
[1] 0.9146657


