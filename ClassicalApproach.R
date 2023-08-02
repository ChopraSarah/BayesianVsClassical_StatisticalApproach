library("sp")

#Plot original data

set.seed(10)
n <- 7
x_orig <- seq(-3, 3, length = n) + rnorm(n, 0, 0.1)
x_orig
## [1] -2.99812538 -2.01842525 -1.13713305 -0.05991677 1.02945451 2.03897943
## [7] 2.87919238
y <- 70 + x_orig + x_orig ^ 3 + rnorm(n, 0, 15)

#Creating the model

X <- cbind(rep(1, n), x_orig, x_orig ^ 3)
colnames(X) <- NULL
X
## [,1] [,2] [,3]
## [1,] 1 -2.99812538 -2.694942e+01

## [2,] 1 -2.01842525 -8.223146e+00
## [3,] 1 -1.13713305 -1.470394e+00
## [4,] 1 -0.05991677 -2.151024e-04
## [5,] 1 1.02945451 1.090992e+00
## [6,] 1 2.03897943 8.476929e+00
## [7,] 1 2.87919238 2.386778e+01
p <- dim(X)[[2]]
p
## [1] 3
B <- solve(t(X) %*% X) %*% t(X) %*% y
B
## [,1]
## [1,] 70.8404944
## [2,] 8.0921003
## [3,] 0.5628859
s2 <- t(y - X %*% B) %*% (y - X %*% B) / (n - p)
s2
## [,1]
## [1,] 186.8761
plot(
  x_orig,
  y,
  xlab = "",
  ylab = "",
  xlim = c(-10, 10),
  ylim = c(0, 180),
  cex = 1,
  pch = 20,
  
  col = "black"
)

xnew <- seq(-10, 10, by = 1)
Xnew <- cbind(rep(1, length(xnew)), xnew, xnew ^ 3)
colnames(Xnew) <- NULL
pred <- Xnew %*% B
ll <- pred - c(qnorm(0.975) * sqrt(s2))
ul <- pred + c(qnorm(0.975) * sqrt(s2))
rect(par("usr")[1],
     par("usr")[3],
     par("usr")[2],
     par("usr")[4],
     par("usr")[5],
     col = "green",
     par(new = TRUE)) # Color
polygon(c(xnew, rev(xnew)), c(ll, rev(ul)),
        col = "grey")
lines(xnew,pred,col = "white",lwd = 1,lty = "dashed")
lines(xnew,ll,col = "black",lwd = 2,lty = 2)
lines(xnew,ul,col = "black",lwd = 2,lty = 2)
points(x_orig, y, col = "blue", pch = 20)

points(6.2, 165, col = "red", pch = 20)
mtext(
  "Road Racing MLE",
  side = 1,
  
  line = -24,
  outer = TRUE,
  font = 2
)
polygon(c(xnew, rev(xnew)), c(ll, rev(ul)),
        col = "grey")

lines(xnew,pred,col = "white",lwd = 1,lty = "dashed")
lines(xnew,ll,col = "black",lwd = 2,lty = 2)
lines(xnew,ul,col = "black",lwd = 2,lty = 2)
points(x_orig, y, col = "blue", pch = 20)
points(6.2, 165, col = "red", pch = 20)
mtext("Road Racing MLE",side = 1, line = -24,outer = TRUE, font = 2)
#Some extra red cars
set.seed(4)
n <- 7
x_red <- seq(-7, 1, length = n) + rnorm(n, 0, 0.1)
x_red
## [1] -6.9783245 -5.7209159 -4.2442189 -2.9404019 -1.5031049 -0.2644058 0.8718753
y_red <- 10 + x_orig + x_orig ^ 3 + rnorm(n, 0, 5)

points(x_red, y_red, col = "red", pch = 20)
in_or_out = point.in.polygon(x_red, y_red, c(xnew, rev(xnew)), c(ll, rev(ul)))
total_vals = mean(in_or_out) * 100
paste(round(total_vals, 2) , "% red cars are inside the CI")
## [1] "28.57 % red cars are inside the CI"

#"28.57 % red cars are inside the CI"