## Not run: 
library(MBA)

##Some useful functions
rmvn <- function(n, mu=0, V = matrix(1)){
  p <- length(mu)
  if(any(is.na(match(dim(V),p)))){stop("Dimension problem!")}
  D <- chol(V)
  t(matrix(rnorm(n*p), ncol=p)%*%D + rep(mu,rep(n,p)))
}

set.seed(1)

##Generate some data
n <- 25 ##number of locations
q <- 2 ##number of outcomes at each location
nltr <- q*(q+1)/2 ##number of triangular elements in the cross-covariance matrix

coords <- cbind(runif(n,0,1), runif(n,0,1))

##Parameters for the bivariate spatial random effects
theta <- rep(3/0.5,q)

A <- matrix(0,q,q)
A[lower.tri(A,TRUE)] <- c(1,-1,0.25)
K <- A%*%t(A)

Psi <- diag(0,q)

C <- mkSpCov(coords, K, Psi, theta, cov.model="exponential")


w.1 <- w[seq(1,length(w),q)]
w.2 <- w[seq(2,length(w),q)]

##Covariate portion of the mean
x.1 <- cbind(1, rnorm(n))
x.2 <- cbind(1, rnorm(n))
x <- mkMvX(list(x.1, x.2))

B.1 <- c(1,-1)
B.2 <- c(-1,1)
B <- c(B.1, B.2)

weight <- 10 ##i.e., trials 
p <- 1/(1+exp(-(x%*%B+w)))
y <- rbinom(n*q, size=rep(weight,n*q), prob=p)

y.1 <- y[seq(1,length(y),q)]
y.2 <- y[seq(2,length(y),q)]

##Call spMvLM
fit <- glm((y/weight)~x-1, weights=rep(weight, n*q), family="binomial")
beta.starting <- coefficients(fit)
beta.tuning <- t(chol(vcov(fit)))

A.starting <- diag(1,q)[lower.tri(diag(1,q), TRUE)]

n.batch <- 100
batch.length <- 50
n.samples <- n.batch*batch.length

starting <- list("beta"=beta.starting, "phi"=rep(3/0.5,q), "A"=A.starting, "w"=0)
tuning <- list("beta"=beta.tuning, "phi"=rep(1,q), "A"=rep(0.1,length(A.starting)),
               "w"=0.5)
priors <- list("beta.Flat", "phi.Unif"=list(rep(3/0.75,q), rep(3/0.25,q)),
               "K.IW"=list(q+1, diag(0.1,q)))

m.1 <- spMvGLM(list(y.1~x.1-1, y.2~x.2-1),
               coords=coords, weights=matrix(weight,n,q),
               starting=starting, tuning=tuning, priors=priors,
               amcmc=list("n.batch"=n.batch,"batch.length"=batch.length,"accept.rate"=0.43),
               cov.model="exponential", n.report=25)

burn.in <- 0.75*n.samples
sub.samps <- burn.in:n.samples

print(summary(window(m.1$p.beta.theta.samples, start=burn.in))$quantiles[,c(3,1,5)])

beta.hat <- t(m.1$p.beta.theta.samples[sub.samps,1:length(B)])
w.hat <- m.1$p.w.samples[,sub.samps]

p.hat <- 1/(1+exp(-(x%*%beta.hat+w.hat)))

y.hat <- apply(p.hat, 2, function(x){rbinom(n*q, size=rep(weight, n*q), prob=p)})

y.hat.mu <- apply(y.hat, 1, mean)

##Unstack to get each response variable fitted values
y.hat.mu.1 <- y.hat.mu[seq(1,length(y.hat.mu),q)]
y.hat.mu.2 <- y.hat.mu[seq(2,length(y.hat.mu),q)]

##Take a look
par(mfrow=c(2,2))
surf <- mba.surf(cbind(coords,y.1),no.X=100, no.Y=100, extend=TRUE)$xyz.est
image(surf, main="Observed y.1 positive trials")
contour(surf, add=TRUE)
points(coords)
zlim <- range(surf[["z"]], na.rm=TRUE)

surf <- mba.surf(cbind(coords,y.hat.mu.1),no.X=100, no.Y=100, extend=TRUE)$xyz.est
image(surf, zlim=zlim, main="Fitted y.1 positive trials")
contour(surf, add=TRUE)
points(coords)

surf <- mba.surf(cbind(coords,y.2),no.X=100, no.Y=100, extend=TRUE)$xyz.est
image(surf, main="Observed y.2 positive trials")
contour(surf, add=TRUE)
points(coords)
zlim <- range(surf[["z"]], na.rm=TRUE)

surf <- mba.surf(cbind(coords,y.hat.mu.2),no.X=100, no.Y=100, extend=TRUE)$xyz.est
image(surf, zlim=zlim, main="Fitted y.2 positive trials")
contour(surf, add=TRUE)
points(coords)

## End(Not run)