#################################################################################
##                                                                             ##
##             When using this code and the underlying ideas,                  ##
##                                                                             ##
##                       >>>>>  PLEASE CITE  <<<<<<                            ##
##-----------------------------------------------------------------------------##
##        Fast Bayesian inference using Laplace approximations                 ##
##       in a flexible promotion time cure model based on P-splines            ##
##                                                                             ##
##               by Oswaldo GRESSANI & Philippe LAMBERT                        ##
##          (Computational Statistics and Data Analysis, 2018)                 ##
##-----------------------------------------------------------------------------##
##                                                                             ##
## The following code implements the Laplace-P-spline model to analyse the     ##
## malignant melanoma survival data (Andersen et al. 1993)                     ##
## presented in Section 4.                                                     ##
##                                                                             ##
## The algorithm reproduces the results given in Table 5 in the manuscript.    ##
## A more general package will be released in the coming months                ##
##                                                                             ##
#################################################################################


library("MASS")
library("survival")
library("fda")

tic = proc.time()

###################### MELANOMA DATA (ANDERSEN ET AL. 1993) ####################
Melanomdata = Melanoma
n = dim(Melanomdata)[1] #Sample size

## Transformation of the survival time in years
Melanomdata$time = Melanomdata$time/365.25

##Covariates influencing the probability to be cured
X = cbind(1,Melanomdata$thickness,Melanomdata$ulcer)
nbeta = dim(X)[2]

## Covariates influencing the time for a cell to yield a detectable tumor
W = cbind(Melanomdata$thickness,Melanomdata$ulcer)
ngamma = dim(W)[2]

## Observed survival time in years
tobs = Melanomdata$time

## Event indicator: 1 died from melanoma, 0 censored.
delta = rep(1,n)
censored = which(Melanomdata$status==2 | Melanomdata$status==3)
delta[censored] = 0

## Last B-spline knot is set at largest observed survival time
xtrunk = max(tobs)

## Kaplan-Meier estimate to check if the follow-up is sufficiently long
time = Surv(tobs,delta)
KapMeier = survfit(time~1)
plot(KapMeier,mark.time=TRUE,mark=4,xlab="Time (in years)")


################################################################################
#### B-splines basis and Penalty matrix
##
## Lower and upper bounds to define the B-spline basis, knots and degree
xl = 0; xr = xtrunk; nknots = 47; degree = 3
dx = (xr-xl)/nknots         #Width of the intervals for B-spline basis knots
ddx = degree*dx
t0 = xl-ddx; tm = xr+ddx    #Upper and lower bounds of B-spline basis
B.knots = seq(t0,tm,by=dx)  #Knots sequence
K = nknots+degree           #Number of B-splines in the basis
H = K+nbeta+ngamma          #Length of the latent field
constr = 10   #Constraint on last B-spline coefficient

## B-spline basis evaluated on tobs
Bbasis1 = create.bspline.basis(c(t0,tm),norder=4,breaks=B.knots)
Bobs = getbasismatrix(tobs,Bbasis1,nderiv=0)[,4:53]

## Plot of the B-spline basis
plot(Bbasis1,xlim=c(xl,xr))

## Fine partition for B-spline basis
binsf = 1500                                  #Total number of bins
partitionf = seq(0,xtrunk,length=binsf+1)     #Partition (has binsf+1 elements)
widthf = partitionf[2]-partitionf[1]          #Width of the bin
middleBinsf = partitionf[1:binsf]+(widthf/2)  #Middle points of the bins
uptof = as.integer(partitionf/widthf)+1       #Bin where partition elements fall
uptof[which(uptof==binsf+1)] = binsf          #Correction for last bin
Bmiddlef = getbasismatrix(middleBinsf,Bbasis1,nderiv=0)[,4:53]

## Bins and midpoints
bins = 300                                #Total number of bins
partition = seq(xl,xr,length=bins+1)      #Partition (contains bins+1 elements)
width = partition[2]-partition[1]         #Width of the bin
middleBins = partition[1:bins]+(width/2)  #Middle points of the bins
Bmiddle = getbasismatrix(middleBins,Bbasis1,nderiv=0)[,4:53]
upto = as.integer(tobs/width)+1           #Indicates in which bin t_i falls
upto[which(upto==bins+1)] = bins          #Correction for t.max (bin=301 to 300)

## Penalty matrix and prior precision for beta and gamma coefficients
r = 3                       #Penalty order
D = diag(K)
for (k in 1:r){D = diff(D)} #dim(D) is c(K-r,K)
P = t(D)%*%D                #dim(P) is c(K,K)
P = P+diag(1e-06,K)         #Add a perturbation to make P full rank
zeta = 1e-05                #Precision associated to regression coefficients
a.delta = b.delta = 1e-04   #Parameters for the Gamma hyperprior
nu = 3                      #Parameter for the lambda prior


#### Hessian and gradient
#########################
## Matrix that will be used to compute Block11 of the Hessian
matones = matrix(0,nrow=n,ncol=bins)
for(i in 1:n){matones[i,] = c(rep(1,upto[i]),rep(0,(bins-(upto[i]))))}

## Matrix that will be used to compute Block 22 of the Hessian
matxx = list(); for(i in 1:n){matxx[[i]] = X[i,]%*%t(X[i,])}

## Matrix that will be used to compute the Block 23 of the Hessian
matxz = list(); for(i in 1:n){matxz[[i]] = X[i,]%*%t(W[i,])}

## Matrix that will be used to compute the Block 33 of the Hessian
matzz = list(); for(i in 1:n){matzz[[i]] = W[i,]%*%t(W[i,])}

## Function to apply cumulative sum by row
rowCumSums = function(x){for(i in seq_len(dim(x)[1])){x[i,] = cumsum(x[i,]) };x}

## Computation of the gradient
gradient = function(field){
    theta = field[1:K]; beta = field[(K+1):(K+nbeta)];gamma=field[(K+nbeta+1):H]
    ##
    X.beta = X%*%beta
    W.gamma = W%*%gamma
    ##
    h0.sj = exp(Bmiddle%*%theta)   #Baseline hazard with B-splines at sj
    h0.width = h0.sj*width
    h0.cumsum = cumsum(h0.width)   #Cumulative baseline hazard
    Omega.0 = h0.cumsum[upto]
    prodmat = t(Bmiddle)*matrix(rep(h0.width,K),nrow=K,byrow=TRUE)
    prodmatcum = rowCumSums(prodmat)
    Omega.0k = prodmatcum[,upto]
    ##
    ## Derivative w.r.t the spline coefficients
    term1.vec = as.numeric(t(delta)%*%Bobs)
    term1.mat = t(as.numeric(t(t(delta)*t(exp(W.gamma))))*t(Omega.0k))
    term1 = term1.vec-rowSums(term1.mat)
    ##
    term2 = t(as.numeric(t(exp(X%*%beta+W.gamma)))*
             as.numeric(t(exp(-Omega.0)^(exp(W.gamma))))*t(Omega.0k))
    term2 = rowSums(term2)
    grad.theta = term1-term2
    ##
    ## Derivative w.r.t the beta coefficients
    term1.beta = as.numeric(rowSums(t(X*delta)))
    term2.beta = rowSums(t(X*matrix(rep(exp(X.beta),nbeta),nrow=n,byrow=FALSE)*
                          as.numeric((1-exp(-Omega.0)^(exp(W.gamma))))))
    grad.beta = as.numeric(term1.beta-term2.beta)
    ##
    ## Derivative w.r.t the gamma coefficients
    ci1 = as.numeric((delta-(delta*exp(W.gamma)*Omega.0)-
                  exp(X.beta+W.gamma)*((exp(-Omega.0))^(exp(W.gamma)))*Omega.0))
    ##
    grad.gamma = as.numeric(unname(ci1%*%W))
    ##
    return(c(grad.theta,grad.beta,grad.gamma))
}

## Computation of the Hessian
Hessian = function(field){
    theta = field[1:K]; beta = field[(K+1):(K+nbeta)];gamma=field[(K+nbeta+1):H]
    ## Terms that will be used in the Hessian
    expWgamma = as.numeric(exp(W%*%gamma))
    expXbeta = as.numeric(exp(X%*%beta))
    expXW = as.numeric(exp(X%*%beta+W%*%gamma))
    expX2W = as.numeric(exp(X%*%beta+2*(W%*%gamma)))
    ##
    h0.sj = as.numeric(exp(Bmiddle%*%theta)) #Baseline hazard
    h0.width = h0.sj*width
    h0.cumsum = cumsum(h0.width)             #Cumulative baseline hazard
    Omega.0 = h0.cumsum[upto]
    prodmat = t(Bmiddle)*matrix(rep((h0.sj*width),K),nrow=K,byrow=TRUE)
    prodmatcum = rowCumSums(prodmat)
    Omega.0k = prodmatcum[,upto]
    expexpW = (exp(-Omega.0))^(expWgamma)
    ################## Block 11
    ci1 = (-1)*(delta*expWgamma+expXW*expexpW)
    ci2 = expX2W*expexpW
    term1 = t(Bmiddle)%*%diag((colSums(matones*ci1)*(h0.sj*width)))%*%Bmiddle
    term2 = Omega.0k%*%(t(Omega.0k)*ci2)
    Block11 = term1+term2
    ################ Block 12 (and Block 21)
    ci1 = (-1)*(expXW*expexpW)
    Block12 = unname(t(ci1*t(Omega.0k))%*%X)
    Block21 = t(Block12)
    ############## Block 13 (and Block 31)
    ci1 = (-1)*(delta*expWgamma+expXW*expexpW-expX2W*expexpW*Omega.0)
    Block13 = unname(t(ci1*t(Omega.0k))%*%W)
    Block31 = t(Block13)
    ############## Block 22
    ci1 = (-1)*(expXbeta*(1-expexpW))
    product = mapply("*",matxx,ci1,SIMPLIFY=FALSE)
    Block22 = unname(Reduce("+",product))
    ############# Block 23 (and Block 32)
    ci1 = (-1)*(expXW*expexpW*Omega.0)
    product2 = mapply("*",matxz,ci1,SIMPLIFY=FALSE)
    Block23 = unname(Reduce("+",product2))
    Block32 = t(Block23)
    ########### Block 33
    ci1 = (-1)*(delta*expWgamma*Omega.0+Omega.0*expXW*expexpW-
               Omega.0*Omega.0*expX2W*expexpW)
    product3 = mapply("*",matzz,ci1,SIMPLIFY=FALSE)
    Block33 = unname(Reduce("+",product3))
    ## Inserting the blocks in the Hessian matrix
    Hessianval = matrix(0,nrow=H,ncol=H)
    ##
    Hessianval[1:K,1:K] = Block11
    Hessianval[1:K,((K+1):(K+nbeta))] = Block12
    Hessianval[1:K,((K+nbeta+1):H)] = Block13
    ##
    Hessianval[((K+1):(K+nbeta)),1:K] = Block21
    Hessianval[((K+1):(K+nbeta)),((K+1):(K+nbeta))] = Block22
    Hessianval[((K+1):(K+nbeta)),((K+nbeta+1):(H))] = Block23
    ##
    Hessianval[((K+nbeta+1):H),1:K] = Block31
    Hessianval[((K+nbeta+1):H),((K+1):(K+nbeta))] = Block32
    Hessianval[((K+nbeta+1):H),((K+1+nbeta):H)] = Block33
    ##
    eigenHess = eigen(Hessianval)$values
    ##
    ## Perturbation following Goldfeld et al.(1966)
    if(sum(eigenHess>=0)>=1){
        eigenmax = max(eigenHess)
        Hessianval = Hessianval-((eigenmax+0.01)*diag(1,ncol=H,nrow=H))}
    return(Hessianval)
}

## Log-likelihood function
##########################
loglik = function(field){
    theta = field[1:K]; beta = field[(K+1):(K+nbeta)];gamma=field[(K+nbeta+1):H]
    ##
    W.gamma = W%*%gamma
    X.beta = X%*%beta
    ##
    h0.sj = as.numeric(exp(Bmiddle%*%theta)) #Baseline hazard
    h0.cumsum = cumsum(h0.sj*width)          #Cumulative baseline hazard
    Omega.0 = h0.cumsum[upto]                #Omega_0
    ##
    term1 = as.numeric(delta*(X.beta+W.gamma+Bobs%*%theta-exp(W.gamma)*Omega.0))
    term2 = as.numeric(exp(X.beta)*(1-(exp(-Omega.0)^(exp(W.gamma)))))
    loglikval = sum(term1-term2)
    ##
    return(loglikval)
}

#### Computation of the Posterior and log-posterior functions for lambda
########################################################################

## Prior precision matrix on the latent field (as a function of lambda)
Q = function(lambda){
    Qmat = matrix(0,nrow=H,ncol=H)
    Qmat[1:K,1:K] = lambda*P
    Qmat[(K+1):H,(K+1):H] = zeta*diag(1,(nbeta+ngamma))
    return(Qmat)
}

## Posterior mode of conditional distribution of latent field (Newton-Raphson)
xstar = function(lambda){
    ## Initial conditions
    xstar.init = rep(0,H); tol = 1e-3; dist = 10; counter = 1
    Q.value = Q(lambda)
    ##
    while(dist>tol){
        H.value = Hessian(xstar.init)
        A = Q.value-H.value
        xstar.new = as.numeric(solve(A)%*%(gradient(xstar.init)-
                    H.value%*%xstar.init))
        d = xstar.new-xstar.init
        dist = as.numeric((t(d)%*%A%*%d))
        xstar.init = xstar.new
    }
    return(xstar.new)
}

## Posterior and log-posterior functions for lambda
##
## This p.lamb function returns the log of the lambda posterior density
p.lamb = function(lambda){
  Q.value = Q(lambda)
  eigen.Q = c(svd(P)$d*lambda,rep(zeta,(nbeta+ngamma)))
  ##
  ## Computation of posterior mode without constraint
  xstar.mode = xstar(lambda)
  cov.star = solve(Q.value-Hessian(xstar.mode))
  ## Computation of x.c constrained mode and covariance
  cov.tilde11 = cov.star[K,K]
  cov.tilde21 = cov.tilde12 = cov.star[K,][-K]
  cov.tilde22 = cov.star[-K,-K]
  ##
  cov.xstar.c = cov.tilde22-(1/cov.tilde11)*outer(cov.tilde21,cov.tilde12)
  eigen.cov.xstar.c = svd(cov.xstar.c)$d
  xstar.c.mode = xstar.mode[-K]+cov.tilde21*(1/cov.tilde11)*
                (constr-xstar.mode[K])
  ##
  ## Computation of x.cc
  xstar.cc.mode = c(xstar.c.mode[1:(K-1)],constr,xstar.c.mode[K:(H-1)])
  ##
  outpt = as.numeric(loglik(xstar.cc.mode)-
                    (0.5*t(xstar.cc.mode)%*%Q.value%*%xstar.cc.mode)+
                    0.5*sum(log(eigen.Q))+
                    0.5*sum(log(eigen.cov.xstar.c))+
                    (0.5*nu-1)*log(lambda)-
                    ((0.5*nu+a.delta)*log(b.delta+0.5*nu*lambda)))
  return(outpt)
}


## Plot of posterior of lambda and mass explored
dom.lamb = seq(1,25000,length=300)
pimage.full = lapply(dom.lamb,p.lamb)
pimage.full = exp(unlist(pimage.full)) #Since p.lamb=log(p(lambda|D))
pimage.full = pimage.full/sum(pimage.full) #Normalizing
lamb.ub.idx = tail(which(cumsum(pimage.full)<=0.965),1)
mass.covered = sum(pimage.full[1:lamb.ub.idx])

## Lambda grid
pen.gridlength = 20      #Number of grid points for lambda
pen.grid = seq(1,dom.lamb[lamb.ub.idx],length=pen.gridlength)

## Plot of the lambda posterior
plot(dom.lamb,pimage.full,type="l",col="green",lwd=2,
     ylab="p(lambda|D)",xlab="lambda")
polygon(c(dom.lamb[1:lamb.ub.idx],rev(dom.lamb[1:lamb.ub.idx])),
        c(rep(0,lamb.ub.idx),rev(pimage.full[1:lamb.ub.idx])),col="skyblue")
lines(dom.lamb,pimage.full,type="l",col="green",lwd=2)
rug(pen.grid,col="red",lwd=2)
legend("topright",
c(paste("Mass covered by grid approx:",round(mass.covered*100,2),"%"),
  "Grid","Posterior"),
lty=c(1,1),col=c("gray","red","green"))

pen.grid = seq(1,dom.lamb[lamb.ub.idx],length=pen.gridlength)


## Grid for delta (5 delta gridpoints for each lambda)
m.delta = 5
delta.mat = matrix(0,ncol=m.delta,nrow=length(pen.grid))
delta.shape = (nu/2)+a.delta

lambda.grid.large = rep(pen.grid,c(rep(m.delta,pen.gridlength)))

for(j in 1:length(pen.grid)){

delta.mat[j,] = seq(qgamma(0.025,shape=delta.shape,
              rate=(0.5*nu*pen.grid[j]+b.delta)),
              qgamma(0.975,shape=delta.shape,
              rate=(0.5*nu*pen.grid[j]+b.delta)),length=m.delta)}

## The function below returns the log posterior for lambda, the
## constrained mode at lambda and the diagonal elements of the covariance matrix
p.lamb.mode = function(lambda){
    Q.value = Q(lambda)
    eigen.Q = c(svd(P)$d*lambda,rep(zeta,(nbeta+ngamma)))
    ##
    ## Computation of posterior mode without constraint
    xstar.mode = xstar(lambda)
    cov.star = solve(Q.value-Hessian(xstar.mode))
    ##
    ## Computation of x.c constrained mode and covariance
    cov.tilde11 = cov.star[K,K]
    cov.tilde21 = cov.tilde12 = cov.star[K,][-K]
    cov.tilde22 = cov.star[-K,-K]
    ##
    cov.xstar.c = cov.tilde22-(1/cov.tilde11)*outer(cov.tilde21,cov.tilde12)
    eigen.cov.xstar.c = svd(cov.xstar.c)$d
    xstar.c.mode = xstar.mode[-K]+cov.tilde21*(1/cov.tilde11)*
                  (constr-xstar.mode[K])
    ##
    ## Computation of x.cc
    xstar.cc.mode = c(xstar.c.mode[1:(K-1)],constr,xstar.c.mode[K:(H-1)])
    var.cc = c(diag(cov.xstar.c)[1:(K-1)],0,diag(cov.xstar.c)[K:(H-1)])
    ##
    outpt = as.numeric(loglik(xstar.cc.mode)-
                      (0.5*t(xstar.cc.mode)%*%Q.value%*%xstar.cc.mode)+
                      0.5*sum(log(eigen.Q))+
                      0.5*sum(log(eigen.cov.xstar.c))+
                      (0.5*nu-1)*log(lambda)-
                      ((0.5*nu+a.delta)*log(b.delta+0.5*nu*lambda)))
    return(list(outpt,xstar.cc.mode,var.cc))
}

p.lamb.mode.result = lapply(pen.grid,p.lamb.mode)
pen.image = sapply(p.lamb.mode.result,"[[",1)
pen.image = rep(pen.image,c(rep(m.delta,pen.gridlength)))

delta.grid = as.vector(t(delta.mat))
delta.image = c()

for(i in 1:length(delta.grid)){
    delta.image[i] = dgamma(x=delta.grid[i],shape=delta.shape,
                           rate=(0.5*nu*lambda.grid.large[i]+b.delta))}

Delta.m = rep(apply(delta.mat,1,diff)[1,1:pen.gridlength],each=5)*
    diff(pen.grid)[1]
normal.weights = (pen.image*delta.image*Delta.m)/
    sum(pen.image*delta.image*Delta.m)
grid.number = length(normal.weights)

## Mixture components
mode.mat = t(sapply(p.lamb.mode.result,"[[",2))
var.mat = t(sapply(p.lamb.mode.result,"[[",3))
mode.mat = mode.mat[rep(seq(1,pen.gridlength,by=1),each=m.delta),]
var.mat = var.mat[rep(seq(1,pen.gridlength,by=1),each=m.delta),]

## Mixture mean and variance
mix.mean = colSums(normal.weights*mode.mat)
mix.var = colSums(normal.weights*
                 (((mode.mat-matrix(rep(mix.mean,each=grid.number),ncol=H))^2)+
                    var.mat))

## Posterior distribution of latent field elements
##################################################
post.latent = function(h,domain){
    norm.image = colSums(normal.weights*sapply(domain,dnorm,mean=mode.mat[,h],
                                              sd=sqrt(var.mat[,h])))
    return(norm.image)}

###### 90% and 95% Quantile based intervals (QBI) #####
QBI90 = QBI95 = matrix(0,ncol=5,nrow=2)

## Computation of credible intervals for reg. coefficients
for(k in (K+1):H){
    field.index = k
    dom.lb = mix.mean[field.index]-4*sqrt(mix.var[field.index])
    dom.ub = mix.mean[field.index]+4*sqrt(mix.var[field.index])
    domain = seq(dom.lb,dom.ub,by=0.001)
    ##
    cum.mass = cumsum(post.latent(field.index,domain)*0.001)
    ##
    lb90 = domain[tail(which(cum.mass<=0.05),n=1)]
    ub90 = domain[tail(which(cum.mass<=0.95),n=1)]
    lb95 = domain[tail(which(cum.mass<=0.025),n=1)]
    ub95 = domain[tail(which(cum.mass<=0.975),n=1)]
    ##
    QBI90[,(field.index-K)] = c(lb90,ub90)
    QBI95[,(field.index-K)] = c(lb95,ub95)
}

## Summarizing the estimation results
summary.table = matrix(0,nrow=6,ncol=8)
summary.table[1,1] = ""
summary.table[1,2] = "Parameters"
summary.table[1,3] = "Estimation"
summary.table[1,4] = "QBI 90%"
summary.table[1,5] = "QBI 90%"
summary.table[1,6] = "QBI 95%"
summary.table[1,7] = "QBI 95%"
summary.table[1,8] = "stdev"

#Writing the results in summary table
for(r in 2:6){
    summary.table[r,1] = ""
    if(r==2){summary.table[r,2] = "Intercept"}
    if(r==3){summary.table[r,2] = "Thickness"}
    if(r==4){summary.table[r,2] = "Ulcer"}
    if(r==5){summary.table[r,2] = "Thickness"}
    if(r==6){summary.table[r,2] = "Ulcer"}
    summary.table[r,3] = round(mix.mean[((K-1)+r)],3)
    summary.table[r,4] = round(QBI90[1,(r-1)],3)
    summary.table[r,5] = round(QBI90[2,(r-1)],3)
    summary.table[r,6] = round(QBI95[1,(r-1)],3)
    summary.table[r,7] = round(QBI95[2,(r-1)],3)
    summary.table[r,8] = round(sqrt(mix.var[((K-1)+r)]),3)
}
summary.table = as.data.frame(summary.table[,2:8])
summary.table

toc = proc.time()-tic
toc


