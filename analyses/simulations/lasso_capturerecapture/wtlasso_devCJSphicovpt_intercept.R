wtlasso_devCJSphicovpt_intercept <- function(b,data,X,eff,e,garb,nh,km1,lasso=0,no.shrink=NULL){
# data encounter histories, eff counts
# e vector of dates of first captures
# garb vector of initial states 
# km1 nb of recapture occasions (nb of capture occ - 1)
# nh nb ind

# logit link
bphi = b[1:ncol(X)]
bp = b[(ncol(X)+1):(ncol(X)+km1)]
phi <- 1/(1+exp(- as.matrix(X) %*% bphi))
p <- 1/(1+exp(- model.matrix(~ as.factor(1:km1)) %*% bp))

# prob of obs (rows) cond on states (col)
B <- array(0,c(2,2,km1)) # allows for time-dep in the parameters
for (t in 1:km1){
	B[,,t] <- matrix(c(1-p[t],1,p[t],0),nrow=2,ncol=2,byrow=T)
}

# first encounter
BE <- matrix(c(0,1,1,0),nrow=2,ncol=2,byrow=T) 
   
# prob of states at t+1 given states at t
A <- array(0,c(2,2,km1))
	for (t in 1:km1){
		A[,,t] <- matrix(c(phi[t], 1 - phi[t],0,1),nrow=2,ncol=2,byrow=T)
}

# init states
PI <- array(0,c(1,2,km1+1))
for (t in 1:(km1+1)){
		PI[,,t] <- c(1,0)
}

# likelihood
   l <- 0
   for (i in 1:nh) # loop on ind
   {
      ei <- e[i] # date of first det
      oe <- garb[i] + 1 # init obs
      evennt <- data[,i] + 1 # non-det become 1s, det 2s
      ALPHA <- PI[,,ei]*BE[oe,]
      for (j in (ei+1):(km1+1)) # cond on first capture
      {
      	if ((ei+1)>(km1+1)) {break} # sous MATLAB la commande >> 8:7 rend >> null, alors que sous R, ça rend le vecteur c(8,7)!
        ALPHA <- (ALPHA %*% A[,,j-1])*B[evennt[j],,j-1]
      }
      l <- l + log(sum(ALPHA))*eff[i] - sum(as.vector(lasso[-no.shrink])*abs(b[-no.shrink]))
   }
    l <- -l
    l
}

