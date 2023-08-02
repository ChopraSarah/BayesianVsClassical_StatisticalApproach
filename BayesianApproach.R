library("sp")

#Plot original data
set.seed(10)
n <- 7
x_orig <- seq(-3, 3, length = n) + rnorm(n, 0, 0.1)
x_orig
## [1] -2.99812538 -2.01842525 -1.13713305 -0.05991677 1.02945451 2.03897943
## [7] 2.87919238
y <- 70 + x_orig + x_orig^3 + rnorm(n, 0, 15)

#Creating the model
n<-length(y)
X <- cbind(rep(1, n), x_orig, x_orig^3)
colnames(X) <- NULL
p<-dim(X)[2]
fit.ls<-lm(y~-1+ X)
beta.0<-fit.ls$coef
nu.0<-1 ; sigma2.0<-sum(fit.ls$res^2)/(n-p)
Sigma.0<- solve(t(X)%*%X)*sigma2.0*n
S<-5000

rmvnorm<-function(n,mu,Sigma)
  
{ # samples from the multivariate normal distribution
  E<-matrix(rnorm(n*length(mu)),n,length(mu))
  t( t(E%*%chol(Sigma)) +c(mu))
}

## some convenient quantites
n<-length(y)
p<-length(beta.0)
iSigma.0<-solve(Sigma.0)
XtX<-t(X)%*%X

## store mcmc samples in these objects
beta.post<-matrix(nrow=S,ncol=p)
sigma2.post<-rep(NA,S)
Sigma.post <- matrix(nrow = S, ncol = p*p)

## starting value
set.seed(1)
sigma2<- var( residuals(lm(y~0+X)) )

for( scan in 1:S) {
  
  #update beta
  V.beta<- solve( iSigma.0 + XtX/sigma2 )
  E.beta<- V.beta%*%( iSigma.0%*%beta.0 + t(X)%*%y/sigma2 )
  beta<-t(rmvnorm(1, E.beta,V.beta) )
  
  #update sigma2
  
  nu.n<- nu.0+n
  ss.n<-nu.0*sigma2.0 + sum( (y-X%*%beta)^2 )
  sigma2<-1/rgamma(1,nu.n/2, ss.n/2)
  
  #save results of this scan
  beta.post[scan,]<-beta
  sigma2.post[scan]<-sigma2
  Sigma.post[scan,] <- c(V.beta)
}

round( apply(beta.post,2,mean), 3)
## [1] 70.796 8.021 0.569
BP <- apply(beta.post, 2, mean)
BP
## [1] 70.7963108 8.0207964 0.5691882
s2P <- mean(sigma2.post)
s2P
## [1] 254.8178
SigmaP <- matrix(apply(Sigma.post, 2, mean), ncol = p, byrow = F)
SigmaP
## [,1] [,2] [,3]
## [1,] 28.3681117 -1.121874 0.2084256
## [2,] -1.1218736 52.612059 -6.8111495
## [3,] 0.2084256 -6.811150 1.0196722
plot(x_orig, y, xlab = "", ylab = "", xlim = c(-10, 10), ylim = c(0, 180), cex = 1, pch = 20, col =
       "black")

xnew <- seq(-10, 10, by = 1)
Xnew <- cbind(rep(1, length(xnew)), xnew, xnew^3)

colnames(Xnew) <- NULL
predP <- Xnew %*% BP
s2P <- sqrt(s2P * diag(Xnew %*% SigmaP %*% t(Xnew)))
ll <- predP - c(qnorm(0.975) * sqrt(s2P))
ul <- predP + c(qnorm(0.975) * sqrt(s2P))
rect(par("usr")[1], par("usr")[3],
     par("usr")[2], par("usr")[4], par("usr")[5],
     col = "green") # Color

polygon(c(xnew, rev(xnew)), c(ll, rev(ul)),
        col = "grey")
lines(xnew, predP, col = "white", lwd =1,lty = "dashed")
lines(xnew, ll, col = "black", lwd = 2, lty = 2)
lines(xnew, ul, col = "black", lwd = 2, lty = 2)
points(x_orig, y, col = "blue",pch = 20)

#Some more random cars
points(6.2, 165, col = "red",pch = 20)
mtext("Road Racing Bayesian", side = 1, line = -24, outer = TRUE,font=2)

set.seed(4)
n <- 7

x_red <- seq(-7, 1, length = n) + rnorm(n, 0, 0.1)
x_red
## [1] -6.9783245 -5.7209159 -4.2442189 -2.9404019 -1.5031049 -0.2644058 0.8718753
y_red <- 10 + x_orig + x_orig^3 + rnorm(n, 0, 5)

points(x_red, y_red, col = "red",pch = 20)
in_or_out = point.in.polygon(x_red,y_red,c(xnew, rev(xnew)),c(ll, rev(ul)))
total_vals = mean(in_or_out)*100
paste(round(total_vals,2) , "% red cars are inside the CI")
## [1] "42.86 % red cars are inside the CI"
#"42.86 % red cars are inside the CI"