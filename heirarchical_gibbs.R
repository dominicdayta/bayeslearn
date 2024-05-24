library(extraDistr)

blood <- data.frame(
  diet_type = c(
    rep("A", times = 4),
    rep("B", times = 6),
    rep("C", times = 6),
    rep("D", times = 8)
  ),
  coagulation = c(
    62, 60, 63, 59, 63, 67, 
    71, 64, 65, 66, 68, 66, 
    71, 67, 68, 68, 56, 62, 
    60, 61, 63, 64, 63, 59
  )
)

# For the group wise averages theta[j] we can take random samples from y[i,j]
# Meanwhile for the overall average mu, we can average out the initial theta[j]
# For tau and sigma, these can be the first parameters drawn during Gibbs

groups <- c("A","B","C","D")
J <- length(groups)
n <- nrow(blood)
maxiter <- 1000

# Initializations for theta and mu
theta <- sapply(groups, function(x){
  this.diet <- blood[blood$diet_type == x, ]$coagulation
  ind <- sample(1:length(this.diet), size = 1)
  return(this.diet[ind])
})

mu <- mean(theta)

# Ready object for saving chain
traces <- list(maxiter)

set.seed(923)
for(i in 1:maxiter){
  
  # draw values for taussq
  tausq.hat <- sum((theta - mu)^2)/(J - 1)
  tausq <- rinvchisq(1, nu = J - 1, tau=theta.hat)
  
  # draw values for sigmasq
  sigmasq.hat <- sapply(groups, function(x){
    this.diet <- blood[blood$diet_type == x, ]$coagulation
    sumsq <- sum((this.diet - theta[x])^2)
    return(sumsq)
  })
  
  sigmasq.hat <- sum(sigmasq.hat)/n
  sigmasq <- rinvchisq(1, nu = n, tau=sigmasq.hat)
  
  # draw values for theta
  theta <- sapply(groups, function(x){
    this.diet <- blood[blood$diet_type == x, ]$coagulation
    nj <- length(this.diet)
    ybar <- mean(this.diet)
    
    thetahat.num <- mu/tausq + ybar * nj/sigmasq
    thetahat.den <- 1/tausq + nj/tausq
    thetahat <- thetahat.num / thetahat.den
    
    V = 1/thetahat.den
    theta <- rnorm(1, thetahat, V)
    return(theta)
  })
  
  # draw values for mu
  muhat <- mean(theta)
  mu <- rnorm(1, muhat, tausq/J)
  
  outputs <- c(i, theta, mu, tausq, sigmasq)
  names(outputs) <- c("iter","theta1","theta2","theta3","theta4","mu","tausq","sigmasq")
  traces[[i]] <- outputs
}

traces <- do.call('rbind', traces)

ggplot(data = data.frame(traces)) +
  geom_line(aes(x = iter, y = theta1))
ggplot(data = data.frame(traces)) +
  geom_line(aes(x = iter, y = theta2))
ggplot(data = data.frame(traces)) +
  geom_line(aes(x = iter, y = theta3))
ggplot(data = data.frame(traces)) +
  geom_line(aes(x = iter, y = theta4))
ggplot(data = data.frame(traces)) +
  geom_line(aes(x = iter, y = mu))
ggplot(data = data.frame(traces)) +
  geom_line(aes(x = iter, y = tausq))
ggplot(data = data.frame(traces)) +
  geom_line(aes(x = iter, y = sigmasq))

cbind(
  mean = apply(traces[200:300,-1], 2, mean),
  median = apply(traces[200:300,-1], 2, median),
  sd = apply(traces[200:300,-1], 2, sd),
  ci_lower = apply(traces[200:300,-1], 2, function(x) quantile(x, probs = c(0.025))),
  ci_upper = apply(traces[200:300,-1], 2, function(x) quantile(x, probs = c(0.975)))
)

