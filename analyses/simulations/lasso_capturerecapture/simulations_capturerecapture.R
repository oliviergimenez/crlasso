####### simulations to compare lasso to AIC

### step 1 ### Build function to run simulations

simul_lasso_cmr <- function(maxlamb, nboot){

#---- load home-made functions and R packages
source('All Functions.R')
library(RMark)

res_mdot = list()
res_m1 = list()
res_m2 = list()
res_m12 = list()
res_lassolike = list()
res_lassopath = list()
nboot = nboot

for (iter in 1:nboot){

#---- generate capture-recapture data and covariates

# define various quantities
n.occasions = 15 # nb recapture occasionsmarked = rep(15, n.occasions-1) # newly marked indmean.phi = 0.8p = rep(0.9, n.occasions-1)beta1 = -0.6 beta2 = 0 

# simulate covariates
X1 = rnorm(n.occasions-1, 0, 1^0.5)
X1 = scale(X1)X2 = rnorm(n.occasions-1, 0, 1^0.5)
X2 = scale(X2)
logit.phi = qlogis(mean.phi) + beta1*X1 + beta2*X2phi = plogis(logit.phi)
PHI = matrix(phi, ncol = n.occasions-1, nrow = sum(marked), byrow = TRUE)P = matrix(p, ncol = n.occasions-1, nrow = sum(marked))

# simulate capture-recapture histories
data = simul.cjs(PHI, P, marked)

#---- standard analyses with AIC

# read in data
storks = data.frame(ch=pasty(data))
storks.processed = process.data(storks)
storks.ddl = make.design.data(storks.processed)
storks.ddl$Phi$x1 = 0
storks.ddl$Phi$x2 = 0
for (i in 1:length(X1)){
	storks.ddl$Phi$x1[storks.ddl$Phi$time==i] = X1[i]
	storks.ddl$Phi$x2[storks.ddl$Phi$time==i] = X2[i]
}

# define model
phidot=list(formula=~1)
phi1=list(formula=~x1)
phi2=list(formula=~x2)
phi12=list(formula=~x1+x2)
pdot=list(formula=~1)

# fit all possible models
mdot = mark(storks.processed,storks.ddl,model.parameters=list(Phi=phidot,p=pdot),output = FALSE,delete=T)
m1 = mark(storks.processed,storks.ddl,model.parameters=list(Phi=phi1,p=pdot),output = FALSE,delete=T)
m2 = mark(storks.processed,storks.ddl,model.parameters=list(Phi=phi2,p=pdot),output = FALSE,delete=T)
m12 = mark(storks.processed,storks.ddl,model.parameters=list(Phi=phi12,p=pdot),output = FALSE,delete=T)

#---- Adaptive lasso analyses

# specify various quantities
nh = dim(data)[1] # nb of individuals
k = dim(data)[2] # nb of capture occasions
km1 = k-1 
eff = rep(1,nh) # counts
fc = NULL
init.state = NULL
for (i in 1:nh){
	temp = 1:k
	fc = c(fc,min(temp[data[i,]==1]))
	init.state = c(init.state,data[i,fc[i]])
}

# form design matrix for survival (detection is constant)
X = cbind(1,X1,X2)

# compute weights for adaptive lasso
a.weights = c(0,1/abs(m12$results$beta$estimate[2:3]),0) # intercept and logit(p) have 0 weights

# lasso path
maxlamb = maxlamb
penalty.min = 1.e-3
n.fits = 50
penalty = exp(seq(log(penalty.min), log(maxlamb), length.out = n.fits))
lassopath = matrix(NA, ncol(X)+1, n.fits)
lassolike = rep(NA, n.fits)

# prepare fitting
par.init  = rep(0, ncol(X)+1) # inits
data = t(data) # transpose data
source('wtlasso_devCJSphicovp_intercept.R') # deviance capture-recapture model with lasso regularisation

# parameters not subject to lasso (intercept and detection prob)
no.shrink = c(1,ncol(X)+1)

for (i in 1:n.fits){
	lasso.vec = rep(penalty[i], length(par.init))*a.weights
	lasso.vec[no.shrink] = 0 # not needed here, but in general yes
	if (i == 1){
		binit = rep(0,ncol(X)+1) # inits
	} else {
		binit = fit.lasso$par # initial values		
	}
	fit.lasso = optim(binit,wtlasso_devCJSphicovp_intercept,NULL,hessian=FALSE,data,X,eff,fc,init.state,nh,km1,lasso=lasso.vec,no.shrink=no.shrink,method="BFGS",control=list(trace=1, REPORT=1,reltol=1e-12))
	lassopath[,i] = fit.lasso$par
	lassolike[i]  = fit.lasso$value
	cat(paste(i, "\n"))
	flush.console()
}

res_mdot[[iter]] = mdot
res_m1[[iter]] = m1
res_m2[[iter]] = m2
res_m12[[iter]] = m12
res_lassolike[[iter]] = lassolike
res_lassopath[[iter]] = lassopath

}

save(res_mdot,res_m1,res_m2,res_m12,res_lassolike,res_lassopath,file='sim_lasso_cmr.RData')

ind_AIC = rep(NA,nboot)
for (i in 1:nboot){
	temp = c(res_mdot[[i]]$results$AICc,res_m1[[i]]$results$AICc,res_m2[[i]]$results$AICc,res_m12[[i]]$results$AICc)
	ind_AIC[i] = which.min(temp)
}

res1 = mean(ind_AIC==2)
# 0.85

lasso.occ = NULL
for (i in 1:nboot){
N = nh # to be checked
num.v = apply(abs(res_lassopath[[i]]) > 1.e-3, 2, sum) - length(no.shrink)
BIC = 2*res_lassolike[[i]] + log(N)*num.v
AIClasso = 2*res_lassolike[[i]] + 2*num.v
lasso.occ = rbind(lasso.occ,c(res_lassopath[[i]][,which.min(BIC)],res_lassopath[[i]][,which.min(AIClasso)]))
}
res_BIC = (abs(lasso.occ[,2]) > 1e-3) & (abs(lasso.occ[,3]) < 1e-3)
res2 = mean(res_BIC)
res_AIC = (abs(lasso.occ[,6]) > 1e-3) & (abs(lasso.occ[,7]) < 1e-3)
res3 = mean(res_AIC)
restot=data.frame(AIC=res1,lassoBIC=res2,lassoAIC=res3)
return(restot)
}

### step 2 ### decide which value to use for lambdamax

res = list()
res[[1]]=simul_lasso_cmr(0.5,20)
res[[2]]=simul_lasso_cmr(1,20)
res[[3]]=simul_lasso_cmr(1.5,20)
res[[4]]=simul_lasso_cmr(2,20)
res[[5]]=simul_lasso_cmr(2.5,20)
res[[6]]=simul_lasso_cmr(3,20)
res[[7]]=simul_lasso_cmr(3.5,20)
res[[8]]=simul_lasso_cmr(4,20)
res[[9]]=simul_lasso_cmr(4.5,20)
res[[10]]=simul_lasso_cmr(5,20)
res

> res
[[1]]
   AIC lassoBIC lassoAIC
1 0.85     0.85     0.75

[[2]]
  AIC lassoBIC lassoAIC
1 0.9      0.9     0.45

[[3]]
   AIC lassoBIC lassoAIC
1 0.85      0.8      0.7

[[4]]
   AIC lassoBIC lassoAIC
1 0.85      0.8      0.6

[[5]]
  AIC lassoBIC lassoAIC
1 0.8     0.85     0.65

[[6]]
  AIC lassoBIC lassoAIC
1   1        1      0.9

[[7]]
   AIC lassoBIC lassoAIC
1 0.85     0.85      0.7

[[8]]
   AIC lassoBIC lassoAIC
1 0.75      0.8     0.55

[[9]]
  AIC lassoBIC lassoAIC
1 0.8      0.8      0.7

[[10]]
   AIC lassoBIC lassoAIC
1 0.95     0.95      0.9

### step 3 ### do simulations with lambdamax = 3

res = simul_lasso_cmr(3,100)
res

#   AIC       lassoBIC     lassoAIC
#   0.8      0.8          0.65
