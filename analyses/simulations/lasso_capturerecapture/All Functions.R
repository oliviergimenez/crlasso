pasty<-function(x){
	k<-ncol(x)
	n<-nrow(x)
	out<-array(dim=n)
	for (i in 1:n){
		y<-(x[i,]>0)*1
		out[i]<-paste(y,collapse="")
	}
return(out)
}


# Define function to simulate a capture-history (CH) matrixsimul.cjs <- function(PHI, P, marked){n.occasions <- dim(PHI)[2] + 1CH <- matrix(0, ncol = n.occasions, nrow = sum(marked))# Define a vector with the occasion of markingmark.occ <- rep(1:length(marked), marked[1:length(marked)])# Fill the CH matrixfor (i in 1:sum(marked)){CH[i, mark.occ[i]] <- 1 # Write an 1 at the release occasionif (mark.occ[i]==n.occasions) nextfor (t in (mark.occ[i]+1):n.occasions){# Bernoulli trial: has individual survived occasion?sur <- rbinom(1, 1, PHI[i,t-1])if (sur==0) break # If dead, move to next individual# Bernoulli trial: has individual been recaptured?rp <- rbinom(1, 1, P[i,t-1])if (rp==1) CH[i,t] <- 1} #t} #ireturn(CH)}

# Create vector with occasion of markingget.first <- function(x) min(which(x!=0))



logit = function(pp) { log(pp) - log(1-pp) }
expit = function(eta) {1/(1+exp(-eta))}

negLL.po = function(alpha = NULL, beta, W.back = NULL, W.po = NULL, X.back, X.po, coord = c("x", "y")) {
  
  coord.col = match(coord, names(X.back))
  back.xy   = X.back[, coord.col]
  po.xy     = X.po[, coord.col]
  X.back    = X.back[, -coord.col]
  X.po      = X.po[, -coord.col]
  
  X.po   = add.ones(X.po)
  X.back = add.ones(X.back)
  if (length(beta) < dim(X.back)[2])
  {
  	beta = c(0, beta)
  }

  if (is.null(alpha) == FALSE)
  {
	  W.back = add.ones(W.back)
 	 W.po   = add.ones(W.po)
 	 if (length(alpha) < dim(W.back)[2])
 	 {
		alpha = c(0, alpha)
 	 }
  }
  
  ux = sort(unique(back.xy[,1]))
  uy = sort(unique(back.xy[,2]))
  cell.width = ux[2] - ux[1]
  cell.height = uy[2] - uy[1]
  area.back = cell.width * cell.height
  lambda = exp(as.matrix(X.back) %*% beta)
  mu = lambda * area.back
  if (is.null(alpha) == FALSE)
  {
    p = expit(W.back %*% alpha)
    logL.po = sum(as.matrix(X.po) %*% beta) + sum(W.po %*% alpha) - sum(log(1 + exp(W.po %*% alpha))) - sum(mu*p)
  }
  else
  {
    logL.po = sum(as.matrix(X.po) %*% beta) - sum(mu)
  }
  
  (-1)*sum(logL.po)
}

makemap = function(fit, Xcov, Xpar = NULL, x, y, main = NA, add = FALSE, minplot = NULL, type = "raw", col.scheme = heat.colors(1024)[900:1], zmax = NULL, cex = 1.2)
{
	if (is.null(Xpar))
	{
		Xpar = 1:length(fit$par)
	}
	beta = fit$par[Xpar]
	lambda = exp(as.matrix(Xcov) %*% beta)
	if (is.null(zmax))
	{
		zmax = max(lambda)
	}
	if (is.null(minplot))
	{
		minplot = min(lambda)
	}
	if (type == "raw")
	{
		plot.out = levelplot(lambda ~ x + y, asp = "iso", col.regions = col.scheme, cuts = length(col.scheme), 
		ylab = "", xlab = "", scales = list(y = list(draw = FALSE), x = list(draw = FALSE)), cex = cex,
		colorkey = list(labels = list(cex = cex)),
		at = seq(minplot, zmax, length.out = length(col.scheme)), main = list(label = main, cex = cex), add = add)
	}
	if (type == "log")
	{
		plot.out = levelplot(log(lambda) ~ x + y, asp = "iso", col.regions = col.scheme, cuts = length(col.scheme), 
		ylab = "", xlab = "", scales = list(y = list(draw = FALSE), x = list(draw = FALSE)), cex = cex,
		colorkey = list(labels = list(cex = cex)),	
		at = seq(min(log(lambda)), zmax, length.out = length(col.scheme)), main = main, add = add)
	}
	print(plot.out)
}

negoccupancyLL <-function(alpha, beta, W, X){
  if (W[,1] != rep(1, dim(W)[1]))
  {
    W = cbind(1, W)
    alpha = c(1, alpha)
  }
  if (X[,1] != rep(1, dim(X)[1]))

  {
    X = cbind(1, X)
    beta = c(1, beta)
  }
  p<-1/(1+exp(-W %*% alpha))
  #psi<-1/(1+exp(-(b0*ones + b1*elev)))
  # cloglog(psi) = log(lambda); 
  # replace the logit link by the cloglog link function, then use the exponential to get lambda
  lambda<-exp(X %*% beta)
  psi = 1-exp(-lambda)
  -1*sum(log(dbinom(y,J,p)*psi + ifelse(y==0,1,0)*(1-psi)))
}

se.param = function(fit)
{
	sqrt(diag(solve(fit$hessian)))
}

newwtnegLL.po = function(par, W.back = NULL, W.po = NULL, X.back, X.po, X.des, ob.wt, pres, y, coord = c("x", "y")) {
  coord.col = match(coord, names(X.back))
  back.xy   = X.back[, coord.col]
  po.xy     = X.po[, coord.col]
  X.back    = X.back[, -coord.col]
  X.po      = X.po[, -coord.col]
  #ob.wt     = weights(po.xy, back.xy, coord = coord)
  #pres      = c(rep(1, dim(X.po)[1]), rep(0, dim(X.back)[1]))
  #y         = pres/ob.wt

  if (is.null(W.back) == FALSE)
  {
	 beta    = par[1:dim(X.back)[2]]
         alpha   = par[(dim(X.back)[2] + 1):length(par)]
#	 X.des   = cbind(as.matrix(rbind(X.po, X.back)), as.matrix(rbind(W.po, W.back)))
	 all.par = c(beta, alpha)
  }
  else
  {
	beta = par
	all.par = beta
#	X.des   = as.matrix(rbind(X.po, X.back))
  }

  mu = exp(as.matrix(X.des) %*% all.par)
  
  #logL.po = sum(ob.wt*(y*log(mu) - mu)) - sum(log(1:sum(y > 0))) - sum(alpha * as.vector(lambda)*abs(beta))
  logL.po = sum(ob.wt*(y*log(mu) - mu)) - sum(log(1:sum(y > 0)))
  (-1)*sum(logL.po)
}


lassonegLL.po = function(par, W.back = NULL, W.po = NULL, X.back, X.po, coord = c("x", "y"), lasso = 0) {
  if (is.null(W.back) == FALSE)
  {
	 alpha  = par[1:dim(W.back)[2]]
	 beta   = par[(dim(W.back)[2] + 1):length(par)]
	 #W.back = add.ones(W.back)
 	 #W.po   = add.ones(W.po)
	 #alpha  = aligndim(alpha, W.back)
  }
  else
  {
	beta = par
  }

  coord.col = match(coord, names(X.back))
  back.xy   = X.back[, coord.col]
  po.xy     = X.po[, coord.col]
  X.back    = X.back[, -coord.col]
  X.po      = X.po[, -coord.col]
  
  #X.po   = add.ones(X.po)
  #X.back = add.ones(X.back)
  #beta   = aligndim(beta, X.back)
  
  ux = sort(unique(back.xy[,1]))
  uy = sort(unique(back.xy[,2]))
  cell.width = ux[2] - ux[1]
  cell.height = uy[2] - uy[1]
  area.back = cell.width * cell.height
  lambda = exp(as.matrix(X.back) %*% beta)
  mu = lambda * area.back
  if (is.null(W.back) == FALSE)
  {
    p = expit(as.matrix(W.back) %*% alpha)
    logL.po = sum(as.matrix(X.po) %*% beta) + sum(as.matrix(W.po) %*% alpha) - sum(log(1 + exp(as.matrix(W.po) %*% alpha))) - sum(mu*p) - sum(as.vector(lasso)*abs(beta))
  }
  else
  {
    logL.po = sum(as.matrix(X.po) %*% beta) - sum(mu) - sum(as.vector(lasso)*abs(beta))
  }
  
  (-1)*sum(logL.po)
}

wtlassonegLL.po = function(par, W.back = NULL, W.po = NULL, X.back, X.po, X.des, ob.wt, pres, y, coord = c("x", "y"), lasso = 0) {
  coord.col = match(coord, names(X.back))
  back.xy   = X.back[, coord.col]
  po.xy     = X.po[, coord.col]
  X.back    = X.back[, -coord.col]
  X.po      = X.po[, -coord.col]
  #ob.wt     = weights(po.xy, back.xy, coord = coord)
  #pres      = c(rep(1, dim(X.po)[1]), rep(0, dim(X.back)[1]))
  #y         = pres/ob.wt

  if (is.null(W.back) == FALSE)
  {
	 beta    = par[1:dim(X.back)[2]]
         alpha   = par[(dim(X.back)[2] + 1):length(par)]
#	 X.des   = cbind(as.matrix(rbind(X.po, X.back)), as.matrix(rbind(W.po, W.back)))
	 all.par = c(beta, alpha)
  }
  else
  {
	beta = par
	all.par = beta
#	X.des   = as.matrix(rbind(X.po, X.back))
  }

  mu = exp(as.matrix(X.des) %*% all.par)
  
  #logL.po = sum(ob.wt*(y*log(mu) - mu)) - sum(log(1:sum(y > 0))) - sum(alpha * as.vector(lambda)*abs(beta))
  logL.po = sum(ob.wt*(y*log(mu) - mu)) - sum(log(1:sum(y > 0))) - sum(as.vector(lasso)*abs(all.par))
  (-1)*sum(logL.po)
}

wtlassonegLL.occ = function(par, W, X, y, site.area, lasso = 0){
	beta = par[1:dim(X)[2]]
	alpha = par[(dim(X)[2]+1):length(par)]
	p = 1/(1 + exp(-W %*% alpha))
	M = nrow(y)
	J = ncol(y) - apply(is.na(y), 1, sum)
	y.bin = data.frame(matrix(as.numeric(as.matrix(y)), M, ncol(y)))
	y.bin = apply(y.bin, 1, sum, na.rm = TRUE)

  #psi<-1/(1+exp(-X %*% beta))
  # cloglog(psi) = log(lambda); 
  # replace the logit link by the cloglog link function, then use the exponential to get lambda
  
  lambda = exp(X %*% beta)
  psi = 1 - exp(-lambda*site.area)
  logL.occ = sum(log(dbinom(y.bin, J, p)*psi + ifelse(y.bin == 0, 1, 0)*(1 - psi))) - sum(as.vector(lasso)*abs(c(beta, alpha)))
  (-1)*sum(logL.occ)
}

wtlassonegLL.log <-function(par, X, y.pa, lasso = 0){
  X = add.ones(X)
  if (length(par) < dim(X)[2])
  {
  	par = c(0, par)
  }
  #psi<-1/(1+exp(-(b0*ones + b1*elev)))
  # cloglog(psi) = log(lambda); 
  # replace the logit link by the cloglog link function, then use the exponential to get lambda
  lambda <- exp(X %*% par)
  psi = 1-exp(-lambda)
  logL.log = sum(log(dbinom(y.pa, 1, psi))) - sum(as.vector(lasso)*abs(par))
  (-1)*logL.log
}

neglogisticLL <-function(beta, X, y.pa){
  X = add.ones(X)
  if (length(beta) < dim(X)[2])
  {
  	beta = c(0, beta)
  }
  #psi<-1/(1+exp(-(b0*ones + b1*elev)))
  # cloglog(psi) = log(lambda); 
  # replace the logit link by the cloglog link function, then use the exponential to get lambda
  lambda <- exp(X %*% beta)
  psi = 1-exp(-lambda)
  -1*sum(log(dbinom(y.pa, 1, psi)))
}

negLLcombined <- function(alpha = NULL, beta, X, W = NULL, X.po, X.back, W.po = NULL, W.back = NULL, y.pa, coord = c("x", "y"))
{
  #X = add.ones(X)
  #X.back = add.ones(X.back)
  #if (length(beta) < dim(X)[2])
  #{
  #	beta = c(0, beta)
  #}
  if (is.null(alpha))
  {
    out = neglogisticLL(beta, X, y.pa) + negLL.po(alpha, beta, W.back, W.po, X.back, X.po, coord)
  }
  else
  {
    out = negoccupancyLL(alpha, beta, X, W) + negLL.po(alpha, beta, W.back, W.po, X.back, X.po, coord)
  }
out
}

add.ones = function(X)
{
	if (sum(X[,1] != rep(1, dim(X)[1])) != 0)
  {
    X = cbind(1, X)
  }
  X
 } 

weights = function(sp.xy, quad.xy, coord = c("X", "Y"))
{
	sp.col   = c(which(names(sp.xy) == coord[1]), which(names(sp.xy) == coord[2]))
	quad.col = c(which(names(quad.xy) == coord[1]), which(names(quad.xy) == coord[2]))
	
	X.inc   = sort(unique(quad.xy[,quad.col[1]]))[2] - sort(unique(quad.xy[,quad.col[1]]))[1]
	Y.inc   = sort(unique(quad.xy[,quad.col[2]]))[2] - sort(unique(quad.xy[,quad.col[2]]))[1]
	quad.0X = min(quad.xy[,quad.col[1]]) - floor(min(quad.xy[,quad.col[1]])/X.inc)*X.inc
	quad.0Y = min(quad.xy[,quad.col[2]]) - floor(min(quad.xy[,quad.col[2]])/Y.inc)*Y.inc

	X = c(sp.xy[,quad.col[1]], quad.xy[,quad.col[1]])
	Y = c(sp.xy[,quad.col[2]], quad.xy[,quad.col[2]])

	round.X     = round((X - quad.0X)/X.inc)*X.inc
	round.Y     = round((Y - quad.0Y)/Y.inc)*Y.inc
	round.id    = paste(round.X, round.Y)
	round.table = table(round.id)
	wt          = X.inc*Y.inc/as.numeric(round.table[match(round.id, names(round.table))])

	wt
}

newnegLL.po = function(par, W.back = NULL, W.po = NULL, X.back, X.po, coord = c("x", "y")) {
  if (is.null(W.back) == FALSE)
  {
	 alpha  = par[1:dim(W.back)[2]]
	 beta   = par[(dim(W.back)[2] + 1):length(par)]
	 #W.back = add.ones(W.back)
 	 #W.po   = add.ones(W.po)
	 #alpha  = aligndim(alpha, W.back)
  }
  else
  {
	beta = par
  }

  coord.col = match(coord, names(X.back))
  back.xy   = X.back[, coord.col]
  po.xy     = X.po[, coord.col]
  X.back    = X.back[, -coord.col]
  X.po      = X.po[, -coord.col]
  
  #X.po   = add.ones(X.po)
  #X.back = add.ones(X.back)
  #beta   = aligndim(beta, X.back)
  
  ux = sort(unique(back.xy[,1]))
  uy = sort(unique(back.xy[,2]))
  cell.width = ux[2] - ux[1]
  cell.height = uy[2] - uy[1]
  area.back = cell.width * cell.height
  lambda = exp(as.matrix(X.back) %*% beta)
  mu = lambda * area.back
  if (is.null(W.back) == FALSE)
  {
    p = expit(as.matrix(W.back) %*% alpha)
    logL.po = sum(as.matrix(X.po) %*% beta) + sum(as.matrix(W.po) %*% alpha) - sum(log(1 + exp(as.matrix(W.po) %*% alpha))) - sum(mu*p)
  }
  else
  {
    logL.po = sum(as.matrix(X.po) %*% beta) - sum(mu)
  }
  
  (-1)*sum(logL.po)
}


get.polys = function(X, X.train = NULL, degree)
{
	X.poly       = c()
	degs         = c()
	poly.coefs   = c()
	names.X.poly = c()
	int.degs     = c()
	n = dim(X)[1]
        for(i in 1:dim(X)[2])
        {
		poly.i          = poly(X[,i], degree = degree, coefs = X.train$coefs[[i]])
		X.poly          = cbind(X.poly, poly.i)
		degs            = c(degs, attr(poly.i, "degree"))
		int.degs[[i]]   = 0:(degree - 1)
		poly.coefs[[i]] = attr(poly.i, "coefs")
		names.X.poly    = c(names.X.poly, dimnames(X)[[2]][i], paste(dimnames(X)[[2]][i], "^", 2:degree, sep = ""))
	}
        dimnames(X.poly)[[2]] = names.X.poly
	
	int.grid = expand.grid(int.degs)
	int.grid = int.grid[rowSums(int.grid) <= degree & rowSums(int.grid > 0) >= 2,]
	exc      = 0:(dim(X)[2] - 1)*degree
	cand.int = int.grid + matrix(rep(exc, each = dim(int.grid)[1]), dim(int.grid)[1], dim(X)[2])
	for (int.i in 1:dim(int.grid)[1])
	{
		X.prod = X.poly[,unlist(setdiff(cand.int[int.i,], exc))]
		X.poly = cbind(X.poly, apply(X.prod, 1, prod))
		dimnames(X.poly)[[2]][dim(X.poly)[2]] = paste(dimnames(X.poly)[[2]][unlist(setdiff(cand.int[int.i,], exc))], collapse = "*")
	}
    return(list(X = X.poly, coefs = poly.coefs))
}

standardise.X = function(mat)
{
	X = scale(as.matrix(mat))
	dat.means = apply(as.matrix(mat), 2, mean, na.rm = TRUE)
	dat.sds   = apply(as.matrix(mat), 2, sd, na.rm = TRUE)
	return(list(X = X, dat.means = dat.means, dat.sds = dat.sds))
}

newinterp = function(sp.xy, sp.scale, f, back.xy, coord = c("X","Y"))
{
	options(scipen = 999)
	x.dat   = sp.xy[,which(names(sp.xy) == coord[1])]
	y.dat   = sp.xy[,which(names(sp.xy) == coord[2])]
	x.back  = back.xy[,which(names(back.xy) == coord[1])]
	y.back  = back.xy[,which(names(back.xy) == coord[2])]

	grid    = data.table(x.back, y.back, f, key = c("x.back", "y.back"))

	ux    = sort(unique(x.back))
	uy    = sort(unique(y.back))

	x.col   = which(names(back.xy) == coord[1])
	y.col   = which(names(back.xy) == coord[2])

	x.step  = ux[2] - ux[1]
	y.step  = uy[2] - uy[1]

	x.o = min(back.xy[,x.col]) - floor(min(back.xy[,x.col])/x.step)*x.step
	y.o = min(back.xy[,y.col]) - floor(min(back.xy[,y.col])/y.step)*y.step

	x.1   = floor((x.dat - x.o)/sp.scale)*sp.scale + x.o
	y.1   = floor((y.dat - y.o)/sp.scale)*sp.scale + y.o
	x.2   = pmin(x.1 + sp.scale, max(ux))
	y.2   = pmin(y.1 + sp.scale, max(uy))

	w11   = (x.2 - x.dat)*(y.2 - y.dat)/((x.2 - x.1)*(y.2 - y.1))
	w12   = (x.2 - x.dat)*(y.dat - y.1)/((x.2 - x.1)*(y.2 - y.1))
	w21   = (x.dat - x.1)*(y.2 - y.dat)/((x.2 - x.1)*(y.2 - y.1))
	w22   = (x.dat - x.1)*(y.dat - y.1)/((x.2 - x.1)*(y.2 - y.1))

	f11   = grid[list(x.1, y.1)]$f
	f12   = grid[list(x.1, y.2)]$f
	f21   = grid[list(x.2, y.1)]$f
	f22   = grid[list(x.2, y.2)]$f

	c11 = 1 - is.na(f11)
	c12 = 1 - is.na(f12)
	c21 = 1 - is.na(f21)
	c22 = 1 - is.na(f22)

	env.wt.mat = cbind(f11*w11*c11, f12*w12*c12, f21*w21*c21, f22*w22*c22) 
	
	f.interp = apply(env.wt.mat, 1, sum, na.rm = TRUE)/(w11*c11 + w12*c12 + w21*c21 + w22*c22)
	f.interp
}


newenv.var = function(sp.xy, env.grid, env.scale, coord = c("X","Y"), file.name = NA)
{
	convert = FALSE
	if (any(lapply(env.grid, class) == "factor"))
	{
		convert  = TRUE
		out.grid = CatConvert(env.grid)
		env.grid = out.grid$X
	}
	x.dat   = sp.xy[,which(names(sp.xy) == coord[1])]
	y.dat   = sp.xy[,which(names(sp.xy) == coord[2])]
	x.back  = env.grid[,which(names(env.grid) == coord[1])]
	y.back  = env.grid[,which(names(env.grid) == coord[2])]
	x.col   = which(names(env.grid) == coord[1])
	y.col   = which(names(env.grid) == coord[2])
	var.col = setdiff(1:dim(env.grid)[2], c(x.col, y.col))
	
	sp.dat        = as.data.frame(matrix(NA, length(x.dat), length(var.col)))
	names(sp.dat) = names(env.grid[var.col])

	for (var in 1:length(var.col))
	{
		loop.scale = env.scale
		loc        = which(is.na(sp.dat[,var]))
		while(sum(is.na(sp.dat[,var])) > 0)
		{
			loc = which(is.na(sp.dat[,var]))
			sp.dat[loc, var] = newinterp(sp.xy[loc,], loop.scale, env.grid[,var.col[var]], env.grid, coord = c("X","Y"))
			loop.scale = loop.scale*2
		}
		cat(paste("Calculating species environmental data for variable:", names(sp.dat)[var], "\n"))
		flush.console()
	}
	
	sp.dat = data.frame(x.dat, y.dat, sp.dat)
	names(sp.dat)[1:2] = c("X", "Y")
	if (is.na(file.name) == FALSE)
	{
		save.name = paste(file.name, ".RData", sep = "")
		save(sp.dat, file = save.name)
		print(paste("Output saved in the file", save.name))
	}
	if (convert == TRUE)
	{
		sp.dat = list(X = sp.dat, cat.names = out.grid$cat.names)
	}
	sp.dat
}

negoccupancyLL <-function(par, W, X, y){
	alpha = par[1:dim(W)[2]]
	beta = par[(dim(W)[2]+1):length(par)]
  if (sum(W[,1] != rep(1, dim(W)[1]))!=0)
  {
    W = cbind(1, W)
    alpha = c(1, alpha)
  }
  if (sum(X[,1] != rep(1, dim(X)[1]))!=0)
  {
    X = cbind(1, X)
    beta = c(1, beta)
  }
  p<-1/(1+exp(-W %*% alpha))
  #psi<-1/(1+exp(-(b0*ones + b1*elev)))
  # cloglog(psi) = log(lambda); 
  # replace the logit link by the cloglog link function, then use the exponential to get lambda
  lambda<-exp(X %*% beta)
  psi = 1-exp(-lambda)
  -1*sum(log(dbinom(y,J,p)*psi + ifelse(y==0,1,0)*(1-psi)))
}

lassomaxpo = function(X.po, X.quad, W.po = NULL, W.quad = NULL, ob.wt = rep(1, dim(X.po)[1] + dim(X.quad)[1]), family = poisson())
{
	y = c(rep(1, dim(X.po)[1]), rep(0, dim(X.quad)[1]))
	#y    = data$Pres/data$wt
	if (is.null(W.po) == FALSE)
	{
		X.des = cbind(1, rbind(X.po, X.quad), rbind(as.matrix(W.po), as.matrix(W.quad)))
	}
	else
	{
		X.des = cbind(1, rbind(X.po, X.quad))
	}
	int.mod = glm(y ~ 1, family = family, weights = ob.wt)
	mu     = int.mod$fitted
	vari   = family$variance(mu)
	deriv  = 1/vari
	wt.mat = ob.wt*1/(deriv^2*vari)
	Xw.s   = t(as.vector(wt.mat) * t(t(X.des)))
	score  = t(as.vector(deriv) * t(Xw.s)) %*% (y - mu)
	scoremax = max(abs(score))
	scoremax
}

lassomaxocc = function(W = NULL, X, y, site.area, family = binomial()) #Not working
{
	M = nrow(y)
	J = ncol(y) - apply(is.na(y), 1, sum)
	y.bin = data.frame(matrix(as.numeric(as.matrix(y)), M, ncol(y)))
	y.bin = apply(y.bin, 1, sum, na.rm = TRUE)
	ob.wt = rep(1, length(y.bin))
	if (is.null(W) == FALSE)
	{
		X.des = cbind(X, W)
		intercept.v = c(1, dim(X)[2] + 1)
	}
	else
	{
		X.des = X
		intercept.v = 1
	}
	int.mod = glm(y.bin/J ~ 1, family = family, weights = ob.wt)
	mu     = int.mod$fitted
	vari   = family$variance(mu)
	deriv  = 1/vari
	wt.mat = ob.wt*1/(deriv^2*vari)
	Xw.s   = t(as.vector(wt.mat) * t(t(X.des)))
	score  = t(as.vector(deriv) * t(Xw.s)) %*% (y.bin - mu)
	scoremax = max(abs(score[-intercept.v]))
	scoremax
}


wtlassomaxpo = function(X.po, X.quad, W.po = NULL, W.quad = NULL, ob.wt = rep(1, dim(X.po)[1] + dim(X.quad)[1]), family = poisson(), coord = c("X", "Y"))
{
	coord.col = match(coord, names(X.quad))
 	back.xy   = X.quad[, coord.col]
  	po.xy     = X.po[, coord.col]
	X.po      = X.po[, -coord.col]
	X.quad    = X.quad[, -coord.col]
	ob.wt     = weights(po.xy, back.xy, coord = coord)
  	pres      = c(rep(1, dim(X.po)[1]), rep(0, dim(X.quad)[1]))
  	y         = pres/ob.wt
	#y = c(rep(1, dim(X.po)[1]), rep(0, dim(X.quad)[1]))
	#y    = data$Pres/data$wt
	if (is.null(W.po) == FALSE)
	{
		X.des = cbind(1, rbind(X.po, X.quad), rbind(as.matrix(W.po), as.matrix(W.quad)))
	}
	else
	{
		X.des = cbind(1, rbind(X.po, X.quad))
	}
	int.mod = glm(y ~ 1, family = family, weights = ob.wt)
	mu     = int.mod$fitted
	vari   = family$variance(mu)
	deriv  = 1/vari
	wt.mat = ob.wt*1/(deriv^2*vari)
	Xw.s   = t(as.vector(wt.mat) * t(t(X.des)))
	score  = t(as.vector(deriv) * t(Xw.s)) %*% (y - mu)
	scoremax = max(abs(score[-1]))
	scoremax
}

score.int = function(y, X.des, ob.wt = rep(1, length(y)), area.int = FALSE, int = NA, family = "poisson")
{
	if (area.int == FALSE)
	{
		int.mod = glm(y ~ 1, family = family, weights = ob.wt)
	}
	if (area.int != FALSE)
	{
		int.mod = glm(y ~ 1 + int, family = family, weights = ob.wt)
	}

	mu     = int.mod$fitted
	vari   = family$variance(mu)
	deriv  = 1/vari
	wt.mat = ob.wt*1/(deriv^2*vari)
	Xw.s   = t(as.vector(wt.mat) * t(t(X.des)))
	score  = t(as.vector(deriv) * t(Xw.s)) %*% (y - mu)
	score
}

#y    = data$Pres/data$wt
#score.0    = score.int(y, X, ob.wt = ob.wt, area.int = area.int, int = interactions, family)

#new.score  = abs(score.0/adapt.weights)
#		sub.score  = new.score[-1]
#		max.lambda = max(sub.score[is.infinite(sub.score) == FALSE])


DecimalCount = function(x, max.dec = max(10, max(nchar(x))))
{
	x.mult  = x*10^(0:max.dec)
	num.dec = min(which(x.mult %% 1 == 0)) - 1
	num.dec
}

sample.quad = function(env.grid, sp.scale, coord = c("X", "Y"), file = "Quad")
{
	convert = FALSE
	#if (any(lapply(env.grid, class) == "factor"))
	#{
	#	is.cat    = which(lapply(env.grid, class) == "factor")
	#	out.grid  = CatConvert(env.grid)
	#	cat.names = out.grid$cat.names
	#	env.grid  = out.grid$X
	#	convert   = TRUE
	#}
	f.name = list()
	for(i in 1:length(sp.scale))
	{
		i.scale = sp.scale[i]
		x.col   = which(names(env.grid) == coord[1])
		y.col   = which(names(env.grid) == coord[2])
		x.uq    = sort(unique(env.grid[, x.col]))
		y.uq    = sort(unique(env.grid[, y.col]))
		
		scale.mult = max(c(unlist(lapply(x.uq, DecimalCount)), unlist(lapply(y.uq, DecimalCount))))
		
		i.scale = i.scale*10^scale.mult
		env.grid[, x.col] = env.grid[, x.col]*10^scale.mult
		env.grid[, y.col] = env.grid[, y.col]*10^scale.mult

		#x.step  = diff((x.uq*10^scale.mult)[1:2])
		x.step	= min(diff(x.uq*10^scale.mult))
		#y.step  = diff((x.uq*10^scale.mult)[1:2])
		y.step	= min(diff(y.uq*10^scale.mult))

		x.o = min(env.grid[,x.col]) - floor(min(env.grid[,x.col])/x.step)*x.step
		y.o = min(env.grid[,y.col]) - floor(min(env.grid[,y.col])/y.step)*y.step
		if (x.o/x.step > 0.5)
		{
			x.o = x.o - x.step
		}	
		if (y.o/y.step > 0.5)
		{
			y.o = y.o - y.step
		}

		is.on.scale   = abs((env.grid[,x.col]/i.scale) - round(env.grid[,x.col]/i.scale) - x.o/i.scale) + abs((env.grid[,y.col]/i.scale) - round(env.grid[,y.col]/i.scale) - y.o/i.scale) < 1.e-8
	  	dat.quad      = env.grid[is.on.scale,]
		dat.quad[, x.col] = dat.quad[, x.col]/(10^scale.mult)
		dat.quad[, y.col] = dat.quad[, y.col]/(10^scale.mult)
		if (is.na(file) == FALSE)
		{
			f.name[[i]] = paste(file, sp.scale[i], ".RData", sep = "")
			save(dat.quad, file = f.name[[i]])
			print(paste("Output saved in the file", f.name[[i]]))
		}
	}
	if (convert == TRUE)
	{
		dat.quad = list(X = dat.quad, cat.names = cat.names)
	}
	if (length(sp.scale) == 1)
		dat.quad
	else
		f.name
}

vis.path = function(path, lasso, main = main, bic = NULL)
{
	par(cex = 1.4)
	max.y = max(path)
	min.y = min(path)
	plot(lasso, path[1,], type = "l", log = "x", xlab = "Lasso penalty", ylab = "Coefficient values", ylim = c(min.y, max.y), main = main)
	for (v in 2:dim(path)[1])
	{
		points(lasso, path[v,], type = "l")
	}
	if (is.null(bic) == FALSE)
	{
		abline(v = bic, lwd = 3, col = "red")
		legend("bottomright", legend = "BIC", lwd = 3, col = "red")
	}
}

lasso.po = function(X.po, X.back, W.po = NULL, W.back = NULL, coord = c("x", "y"), method = 'BFGS', hessian = FALSE, ob.wt = rep(1, dim(X.po)[1] + dim(X.back)[1]), n.fits = 50, lasso.min = 1.e-3)
{
	#X.po = as.matrix(X.po)
	#X.back = as.matrix(X.back)
	coord.id = match(coord, dimnames(X.po)[[2]])
	if (is.null(W.po) == TRUE)
	{
		init.par = rep(0, dim(X.po)[2] - 1)
		fit.init = optim(par = init.par, fn = newnegLL.po, W.back = NULL, W.po = NULL, X.back = add.ones(X.back), X.po = add.ones(X.po),  coord = coord, method = method, hessian = hessian)
		w.col = 0
	}
	if (is.null(W.po) == FALSE)
	{
		W.po = as.matrix(W.po)
		W.back = as.matrix(W.back)
		init.par = rep(0, dim(X.po)[2] + dim(W.po)[2] - 1)
		fit.init = optim(par = init.par, fn = newnegLL.po, W.back = as.matrix(W.back), W.po = as.matrix(W.po), X.back = add.ones(X.back), X.po = add.ones(X.po),  coord = coord, method = method, hessian = hessian)
		w.col = dim(W.po)[2]
	}
	lasso.max = lassomaxpo(X.po[,-coord.id], X.back[,-coord.id], W.po, W.back, ob.wt)
	path.penalty = exp(seq(log(lasso.min), log(lasso.max), length.out = (n.fits + 1)))[-1]
	lassopar = matrix(NA, length(fit.init$par), n.fits)
	lassolikes = rep(NA, n.fits)
	for (path in 1:n.fits)
	{
		lasso.vec = c(rep(path.penalty[path], w.col), 0, rep(path.penalty[path], dim(X.po[,-coord.id])[2]))
		if (path == 1)
		{
			fit.lasso = optim(par = fit.init$par, fn = lassonegLL.po, W.back = W.back, W.po = W.po, X.back = add.ones(X.back), X.po = add.ones(X.po),  coord = c("x", "y"), lasso = lasso.vec, method=method, control = list(maxit = 20))
		}
		if (path > 1)
		{
			fit.lasso = optim(par = lassopar[, (path - 1)], fn = lassonegLL.po, W.back = W.back, W.po = W.po, X.back = add.ones(X.back), X.po = add.ones(X.po),  coord = c("x", "y"), lasso = lasso.vec, method=method, control = list(maxit = 20))
		}
		lassopar[, path] = fit.lasso$par
		lassolikes[path]  = fit.lasso$value
		cat(paste(path, "of", n.fits, "\n"))
		flush.console()
	}
	return(list(par = lassopar, like = lassolikes, lasso = path.penalty))
}

wtlasso.po = function(X.po, X.back, W.po = NULL, W.back = NULL, X.des, ob.wt, pres, y, coord = c("x", "y"), method = 'BFGS', hessian = FALSE, n.fits = 50, lasso.min = 1.e-3)
{
	#X.po = as.matrix(X.po)
	#X.back = as.matrix(X.back)
	coord.id = match(coord, dimnames(X.po)[[2]])
	if (is.null(W.po) == TRUE)
	{
		init.par = rep(0, dim(X.po)[2] - 1)
		#fit.init = optim(par = init.par, fn = newnegLL.po, W.back = NULL, W.po = NULL, X.back = add.ones(X.back), X.po = add.ones(X.po),  coord = coord, method = method, hessian = hessian)
		fit.init = optim(par = init.par, fn = newwtnegLL.po, W.back = NULL, W.po = NULL, X.back = add.ones(X.back), X.po = add.ones(X.po), X.des = X.des, ob.wt = ob.wt, pres = pres, y = y, coord = coord, method = method, hessian = hessian)
		w.col = 0
	}
	if (is.null(W.po) == FALSE)
	{
		W.po = as.matrix(W.po)
		W.back = as.matrix(W.back)
		init.par = rep(0, dim(X.po)[2] + dim(W.po)[2] - 1)
		#fit.init = optim(par = init.par, fn = newnegLL.po, W.back = as.matrix(W.back), W.po = as.matrix(W.po), X.back = add.ones(X.back), X.po = add.ones(X.po),  coord = coord, method = method, hessian = hessian)
		fit.init = optim(par = init.par, fn = newwtnegLL.po, W.back = as.matrix(W.back), W.po = as.matrix(W.po), X.back = add.ones(X.back), X.po = add.ones(X.po), X.des = X.des, ob.wt = ob.wt, pres = pres, y = y, coord = coord, method = method, hessian = hessian)
		w.col = dim(W.po)[2]
	}
	cat(paste("Initial fit done", "\n"))
	flush.console()
	lasso.max = wtlassomaxpo(X.po, X.back, W.po, W.back, ob.wt, coord = coord)
	path.penalty = exp(seq(log(lasso.min), log(lasso.max*1.001), length.out = (n.fits + 1)))[-1]
	lassopar = matrix(NA, length(fit.init$par), n.fits)
	lassolikes = rep(NA, n.fits)
	for (path in 1:n.fits)
	{
		lasso.vec = c(0, rep(path.penalty[path], dim(X.des)[2]))
		if (path == 1)
		{
			#fit.lasso = optim(par = fit.init$par, fn = lassonegLL.po, W.back = W.back, W.po = W.po, X.back = add.ones(X.back), X.po = add.ones(X.po),  coord = c("x", "y"), lasso = lasso.vec, method=method, control = list(maxit = 20))
			fit.lasso = optim(par = fit.init$par, fn = wtlassonegLL.po, W.back = W.back, W.po = W.po, X.back = add.ones(X.back), X.po = add.ones(X.po), X.des = X.des, ob.wt = ob.wt, pres = pres, y = y,  coord = coord, lasso = lasso.vec, method=method, control = list(maxit = 100))
		}
		if (path > 1)
		{
			#fit.lasso = optim(par = lassopar[, (path - 1)], fn = lassonegLL.po, W.back = W.back, W.po = W.po, X.back = add.ones(X.back), X.po = add.ones(X.po),  coord = c("x", "y"), lasso = lasso.vec, method=method, control = list(maxit = 20))
			fit.lasso = optim(par = lassopar[, (path - 1)], fn = wtlassonegLL.po, W.back = W.back, W.po = W.po, X.back = add.ones(X.back), X.po = add.ones(X.po), X.des = X.des, ob.wt = ob.wt, pres = pres, y = y, coord = coord, lasso = lasso.vec, method=method, control = list(maxit = 100))
		}
		lassopar[, path] = fit.lasso$par
		lassolikes[path]  = fit.lasso$value
		cat(paste(path, "of", n.fits, "\n"))
		flush.console()
	}
	return(list(par = lassopar, like = lassolikes, lasso = path.penalty))
}

altwtlasso.po = function(X.po, X.back, W.po = NULL, W.back = NULL, X.des, ob.wt, pres, y, coord = c("x", "y"), method = 'BFGS', hessian = FALSE, n.fits = 50, lasso.min = 1.e-3)
{
	#X.po = as.matrix(X.po)
	#X.back = as.matrix(X.back)
	coord.id = match(coord, dimnames(X.po)[[2]])
	if (is.null(W.po) == TRUE)
	{
		init.par = rep(0, dim(X.po)[2] - 1)
		#fit.init = optim(par = init.par, fn = newnegLL.po, W.back = NULL, W.po = NULL, X.back = add.ones(X.back), X.po = add.ones(X.po),  coord = coord, method = method, hessian = hessian)
		fit.init = optim(par = init.par, fn = newwtnegLL.po, W.back = NULL, W.po = NULL, X.back = add.ones(X.back), X.po = add.ones(X.po), X.des = X.des, ob.wt = ob.wt, pres = pres, y = y, coord = coord, method = method, hessian = hessian)
		w.col = 0
	}
	if (is.null(W.po) == FALSE)
	{
		W.po = as.matrix(W.po)
		W.back = as.matrix(W.back)
		init.par = rep(0, dim(X.po)[2] + dim(W.po)[2] - 1)
		#fit.init = optim(par = init.par, fn = newnegLL.po, W.back = as.matrix(W.back), W.po = as.matrix(W.po), X.back = add.ones(X.back), X.po = add.ones(X.po),  coord = coord, method = method, hessian = hessian)
		fit.init = optim(par = init.par, fn = newwtnegLL.po, W.back = as.matrix(W.back), W.po = as.matrix(W.po), X.back = add.ones(X.back), X.po = add.ones(X.po), X.des = X.des, ob.wt = ob.wt, pres = pres, y = y, coord = coord, method = method, hessian = hessian)
		w.col = dim(W.po)[2]
	}
	cat(paste("Initial fit done", "\n"))
	flush.console()
	lasso.max = wtlassomaxpo(X.po, X.back, W.po, W.back, ob.wt, coord = coord)
	path.penalty = exp(seq(log(lasso.min), log(lasso.max*1.001), length.out = (n.fits + 1)))[-1]
	path.penalty = sort(path.penalty, decreasing = TRUE)
	lassopar = matrix(NA, length(fit.init$par), n.fits)
	lassolikes = rep(NA, n.fits)
	for (path in 1:n.fits)
	{
		lasso.vec = c(0, rep(path.penalty[path], dim(X.des)[2]))
		if (path == 1)
		{
			#fit.lasso = optim(par = fit.init$par, fn = lassonegLL.po, W.back = W.back, W.po = W.po, X.back = add.ones(X.back), X.po = add.ones(X.po),  coord = c("x", "y"), lasso = lasso.vec, method=method, control = list(maxit = 20))
			fit.lasso = optim(par = rep(0, length(fit.init$par)), fn = wtlassonegLL.po, W.back = W.back, W.po = W.po, X.back = add.ones(X.back), X.po = add.ones(X.po), X.des = X.des, ob.wt = ob.wt, pres = pres, y = y,  coord = coord, lasso = lasso.vec, method=method, control = list(maxit = 1000))
		}
		if (path > 1)
		{
			#fit.lasso = optim(par = lassopar[, (path - 1)], fn = lassonegLL.po, W.back = W.back, W.po = W.po, X.back = add.ones(X.back), X.po = add.ones(X.po),  coord = c("x", "y"), lasso = lasso.vec, method=method, control = list(maxit = 20))
			fit.lasso = optim(par = lassopar[, (path - 1)], fn = wtlassonegLL.po, W.back = W.back, W.po = W.po, X.back = add.ones(X.back), X.po = add.ones(X.po), X.des = X.des, ob.wt = ob.wt, pres = pres, y = y, coord = coord, lasso = lasso.vec, method=method, control = list(maxit = 1000))
		}
		lassopar[, path] = fit.lasso$par
		lassolikes[path]  = fit.lasso$value
		cat(paste(path, "of", n.fits, "\n"))
		flush.console()
	}
	return(list(par = lassopar, like = lassolikes, lasso = path.penalty))
}

