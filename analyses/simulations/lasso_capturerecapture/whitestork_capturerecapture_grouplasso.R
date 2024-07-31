####### Lasso implementation for a capture-recapture model
####### using white stork data as a case study
####### check out this paper http://bit.ly/2ik437S

# load detection/non-detection capture-recapture data and covariates
load('whitestork.RData')

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
lassopath = matrix(NA, ncol(X)+km1, n.fits)
lassolike = rep(NA, n.fits)
par.init  = rep(0, ncol(X)+km1)
data <- t(data) # transpose data
source('wtlasso_devCJSphicovpt_intercept.R') # deviance capture-recapture model with lasso regularisation

no.shrink = c(1,ncol(X)+1)

for (i in 1:n.fits)
{
	lasso.vec = rep(penalty[i], length(par.init))
	lasso.vec[no.shrink] = 0
	binit <- rep(0,ncol(X)+km1) # init values
	fit.lasso = optim(binit,wtlasso_devCJSphicovpt_intercept,NULL,hessian=FALSE,data,X,eff,fc,init.state,
	nh,km1,lasso=lasso.vec,method="BFGS",control=list(trace=1, REPORT=1))
	lassopath[,i] = fit.lasso$par
	lassolike[i]  = fit.lasso$value
	cat(paste(i, "\n"))
	flush.console()
}

save(lassopath,lassolike,file='ws_lasso.RData')

N = nh # to be checked
num.v = apply(abs(lassopath) > 1.e-7, 2, sum) - length(no.shrink)
BIC = 2*lassolike + log(N)*num.v
plot(penalty, BIC, log = "x", type = "o", xlab = "Lasso Penalty", ylab = "BIC")

pen = penalty[which.min(BIC)] # get 'best' penalisation
lasso.vec = rep(pen, length(par.init))
lasso.vec[no.shrink] = 0

source('All Functions.R')
vis.path(lassopath[-no.shrink,], penalty, bic = pen, main = "Regularization Path for CR model and white stork data")

# get estimates from model with selected lasso penalisation 
lasso.occ = lassopath[,which.min(BIC)]

# display results
x <- lasso.occ
cov.names = names(X)
cbind(x[1:ncol(X)],cov.names)

                            cov.names   
 [1,] "0.703373327472566"   "1"         
 [2,] "0.152287661833539"   "Diourbel"  
 [3,] "-0.193040616382062"  "Gao"       
 [4,] "0.0568281921800258"  "Kayes"     
 [5,] "0.400845288103173"   "Kita"      
 [6,] "-0.270756459622895"  "Maradi"    
 [7,] "0.179752787623818"   "Mopti"     
 [8,] "0.100335189019884"   "Ouahigouya"
 [9,] "-0.0939706236706635" "Segou"     
[10,] "0.196626368055611"   "Tahoua"    
[11,] "0.0622333427522267"  "Tombouctou"

matrix(x[(ncol(X)+1):(ncol(X)+km1)])

 [1,]  2.2479555068
 [2,] -0.1805916512
 [3,] -0.3412256657
 [4,]  0.7675258962
 [5,]  0.6276187320
 [6,] -0.2251504108
 [7,]  0.0006379526
 [8,] -0.0670712394
 [9,]  0.6290252680
[10,]  0.5011853194
[11,]  0.2152453918
[12,]  0.0004945467
[13,] -0.5779829395
[14,]  2.4258155323
[15,] -0.0803505306
[16,]  1.3201443256


# survival probs
phi <- exp(as.matrix(X) %*% x[1:ncol(X)])/(1+exp(as.matrix(X) %*% x[1:ncol(X)]))
# detection probs
p <- exp(model.matrix(~ as.factor(1:km1)) %*% x[(ncol(X)+1):(ncol(X)+km1)])/(1+exp(model.matrix(~ as.factor(1:km1)) %*% x[(ncol(X)+1):(ncol(X)+km1)]))
phi

           [,1]
 [1,] 0.7438028
 [2,] 0.8029669
 [3,] 0.7638501
 [4,] 0.7043928
 [5,] 0.6328182
 [6,] 0.5435576
 [7,] 0.7343664
 [8,] 0.5347651
 [9,] 0.6541267
[10,] 0.7656508
[11,] 0.6173899
[12,] 0.7574634
[13,] 0.6326602
[14,] 0.5911539
[15,] 0.5404448
[16,] 0.5873095

p

        [,1]
1  0.9044740
2  0.8876904
3  0.8706513
4  0.9532686
5  0.9466257
6  0.8831708
7  0.9045291
8  0.8985197
9  0.9466967
10 0.9398648
11 0.9215215
12 0.9045168
13 0.8415722
14 0.9907494
15 0.8973025
16 0.9725645

