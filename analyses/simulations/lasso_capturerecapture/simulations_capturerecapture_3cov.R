####### Lasso implementation for a capture-recapture model
####### using simulated data

source('All Functions.R')
library(RMark)

# generate detection/non-detection capture-recapture data and covariates

n.occasions <- 15 # Number of recapture occasionsmarked <- rep(15, n.occasions-1) # Annual number of newly marked individualsmean.phi <- 0.8p <- rep(0.9, n.occasions-1)beta1 <- -0.6 # Slope of survival-winter relationshipbeta2 <- 0 # Slope of survival-winter relationship
beta3 <- 0 # Slope of survival-winter relationship
# Draw annual survival probabilitiesX1 <- rnorm(n.occasions-1, 0, 1^0.5)
X1 = scale(X1)X2 <- rnorm(n.occasions-1, 0, 1^0.5)
X2 = scale(X2)
X3 <- rnorm(n.occasions-1, 0, 1^0.5)
X3 = scale(X3)
logit.phi <- qlogis(mean.phi) + beta1*X1 + beta2*X2 + beta3*X3phi <- plogis(logit.phi)

clever_inits <- c(logit.phi,log(.9/(1-.9)))

# Define matrices with survival and recapture probabilitiesPHI <- matrix(phi, ncol = n.occasions-1, nrow = sum(marked), byrow = TRUE)P <- matrix(p, ncol = n.occasions-1, nrow = sum(marked))# Simulate capture historiesdata <- simul.cjs(PHI, P, marked)

#------- MARK analysis

storks<-data.frame(ch=pasty(data))
storks.processed=process.data(storks)
storks.ddl=make.design.data(storks.processed)
storks.ddl$Phi$x1=0
storks.ddl$Phi$x2=0
storks.ddl$Phi$x3=0
for (i in 1:length(X1))
{
storks.ddl$Phi$x1[storks.ddl$Phi$time==i]=X1[i]
storks.ddl$Phi$x2[storks.ddl$Phi$time==i]=X2[i]
storks.ddl$Phi$x3[storks.ddl$Phi$time==i]=X3[i]
}

#  Define model
PhiCov=list(formula=~x1+x2+x3)
pdot=list(formula=~1)

# phix,p
phixp = mark(storks.processed,storks.ddl,model.parameters=list(Phi=PhiCov,p=pdot),output = FALSE,delete=T)

phixp$results$beta
phixp$results$real

#---------------------

# specify various quantities
nh <- dim(data)[1] # nb of individuals
k <- dim(data)[2] # nb of capture occasions
km1 <- k-1 
eff <- rep(1,nh) # counts

# compute the date of first capture fc, and state at initial capture init.state
fc <- NULL
init.state <- NULL
for (i in 1:nh){
temp <- 1:k
fc <- c(fc,min(temp[data[i,]==1]))
init.state <- c(init.state,data[i,fc[i]])
}

# form design matrix for survival (detection is constant)
X = cbind(1,X1,X2,X3)
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
	lasso.vec = rep(penalty[i], length(par.init))
	lasso.vec[no.shrink] = 0
	if (i == 1){
		binit <- rep(0,ncol(X)+1) # inits
	} else {
		binit <- fit.lasso$par # initial values		
	}
	fit.lasso = optim(binit,wtlasso_devCJSphicovp_intercept,NULL,hessian=FALSE,data,X,eff,fc,init.state,
	nh,km1,lasso=lasso.vec,no.shrink=no.shrink,method="BFGS",control=list(trace=1, REPORT=1,reltol=1e-12))
	lassopath[,i] = fit.lasso$par
	lassolike[i]  = fit.lasso$value
	cat(paste(i, "\n"))
	flush.console()
}

#save(lassopath,lassolike,file='ws_lasso.RData')

N = nh # to be checked
num.v = apply(abs(lassopath) > 1.e-3, 2, sum) - length(no.shrink)
BIC = 2*lassolike + log(N)*num.v
#BIC = 2*lassolike + 2*num.v
plot(penalty, BIC, log = "x", type = "o", xlab = "Lasso Penalty", ylab = "BIC")

pen = penalty[which.min(BIC)] # get 'best' penalisation
lasso.vec = rep(pen, length(par.init))
lasso.vec[no.shrink] = 0

vis.path(lassopath[-no.shrink,], penalty, bic = pen, main = "Regularization Path")

# get estimates from model with selected lasso penalisation 
lasso.occ = lassopath[,which.min(BIC)]

# display results
x <- lasso.occ
data.frame(true=c(qlogis(mean.phi),-0.5,0,0,log(0.9/(1-.9))),lasso=c(x[1],x[2:(ncol(X)+1)]),mark=phixp$results$beta[,'estimate'])

# survival probs
phi_estim <- exp(as.matrix(X) %*% x[1:ncol(X)])/(1+exp(as.matrix(X) %*% x[1:ncol(X)]))
# detection probs
p_estim <- exp(x[(ncol(X)+1)])/(1+exp(x[(ncol(X)+1)]))
data.frame(true=phi,lasso=phi_estim,mark=phixp$results$real[-nrow(phixp$results$real),'estimate'])
data.frame(true=p,lasso=p_estim,mark=phixp$results$real[nrow(phixp$results$real),'estimate'])

