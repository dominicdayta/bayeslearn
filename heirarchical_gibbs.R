library(extraDistr)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(rstan)

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
warmup <- 500

# Initializations for theta and mu
theta <- sapply(groups, function(x){
  this.diet <- blood[blood$diet_type == x, ]$coagulation
  ind <- sample(1:length(this.diet), size = 1)
  return(this.diet[ind])
})

mu <- mean(theta)

# Ready object for saving chain
chains <- list(
  list(maxiter),
  list(maxiter),
  list(maxiter)
)

set.seed(923)
for(i in 1:maxiter){
  
  for(j in 1:3){
    # draw values for taussq
    tausq.hat <- sum((theta - mu)^2)/(J - 1)
    tausq <- rinvchisq(1, nu = J - 1, tau=tausq.hat)
    
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
      thetahat.den <- 1/tausq + nj/sigmasq
      thetahat <- thetahat.num / thetahat.den
      
      V = 1/thetahat.den
      theta <- rnorm(1, thetahat, sqrt(V))
      return(theta)
    })
    
    # draw values for mu
    muhat <- mean(theta)
    mu <- rnorm(1, muhat, sqrt(tausq/J))
    
    outputs <- c(i, theta, mu, tausq, sigmasq)
    names(outputs) <- c("iter","theta1","theta2","theta3","theta4","mu","tausq","sigmasq")
    chains[[j]][[i]] <- outputs
      
  }
  
}

# throw away the first 500 iterations as warmup
chains[[1]] <- chains[[1]][(warmup + 1):maxiter]
chains[[2]] <- chains[[2]][(warmup + 1):maxiter]
chains[[3]] <- chains[[3]][(warmup + 1):maxiter]

trace1 <- do.call('rbind', chains[[1]])
trace2 <- do.call('rbind', chains[[2]])
trace3 <- do.call('rbind', chains[[3]])

# 3 chains for theta 
p1 <- ggplot(data = data.frame(trace1)) +
  geom_line(aes(x = iter, y = theta1), color = "midnightblue") +
  geom_line(data = data.frame(trace2), aes(x = iter, y = theta1), color = "lightblue1") +
  geom_line(data = data.frame(trace3), aes(x = iter, y = theta1), color = "lightblue3")

p2 <- ggplot(data = data.frame(trace1)) +
  geom_line(aes(x = iter, y = theta2), color = "midnightblue") +
  geom_line(data = data.frame(trace2), aes(x = iter, y = theta2), color = "lightblue1") +
  geom_line(data = data.frame(trace3), aes(x = iter, y = theta2), color = "lightblue3")

p3 <- ggplot(data = data.frame(trace1)) +
  geom_line(aes(x = iter, y = theta3), color = "midnightblue") +
  geom_line(data = data.frame(trace2), aes(x = iter, y = theta3), color = "lightblue1") +
  geom_line(data = data.frame(trace3), aes(x = iter, y = theta3), color = "lightblue3")

p4 <- ggplot(data = data.frame(trace1)) +
  geom_line(aes(x = iter, y = theta4), color = "midnightblue") +
  geom_line(data = data.frame(trace2), aes(x = iter, y = theta4), color = "lightblue1") +
  geom_line(data = data.frame(trace3), aes(x = iter, y = theta4), color = "lightblue3")

p5 <- ggplot(data = data.frame(trace1)) +
  geom_line(aes(x = iter, y = mu), color = "midnightblue") +
  geom_line(data = data.frame(trace2), aes(x = iter, y = mu), color = "lightblue1") +
  geom_line(data = data.frame(trace3), aes(x = iter, y = mu), color = "lightblue3")

p6 <- ggplot(data = data.frame(trace1)) +
  geom_line(aes(x = iter, y = tausq), color = "midnightblue") +
  geom_line(data = data.frame(trace2), aes(x = iter, y = tausq), color = "lightblue1") +
  geom_line(data = data.frame(trace3), aes(x = iter, y = tausq), color = "lightblue3")

p7 <- ggplot(data = data.frame(trace1)) +
  geom_line(aes(x = iter, y = sigmasq), color = "midnightblue") +
  geom_line(data = data.frame(trace2), aes(x = iter, y = sigmasq), color = "lightblue1") +
  geom_line(data = data.frame(trace3), aes(x = iter, y = sigmasq), color = "lightblue3")

chains <- ggarrange(p1, p2, p3, p4, p5, p6, p7, nrow = 2, ncol = 4, 
                    common.legend = TRUE, legend = "bottom")

# --- Check if chains have mixed
varChains <- numeric(ncol(trace1[,-1]))
Rhat <- numeric(ncol(trace1[,-1]))

names(varChains) <- colnames(trace1[,-1])
names(Rhat) <- colnames(trace1[,-1])

lasthalf <- round(nrow(trace1)/2 + 1):nrow(trace1)
for(k in 1:ncol(trace1[,-1])){
  chains <- cbind(trace1[lasthalf,k+1], trace2[lasthalf,k+1], trace3[lasthalf,k+1])
  n <- nrow(chains)
  B <- n * var(apply(chains, 2, mean))
  
  W <- mean(apply(chains, 2, var))
  
  varChains[k] <- (n - 1)/n * W + 1/n * B
  Rhat[k] <- sqrt(varChains[k]/W)
}

Rhat

all_traces <- rbind(trace1,trace2,trace3)
cbind(
  mean = apply(all_traces[,-1], 2, mean),
  median = apply(all_traces[,-1], 2, median),
  sd = apply(all_traces[,-1], 2, sd),
  ci_lower = apply(all_traces[,-1], 2, function(x) quantile(x, probs = c(0.025))),
  ci_upper = apply(all_traces[,-1], 2, function(x) quantile(x, probs = c(0.975))),
  Rhat = Rhat
)

# --- Using RStan
N <- 24                     # Number of observations
J <- 4                      # Number of groups

group <- c(
  rep(1, times = 4),
  rep(2, times = 6),
  rep(3, times = 6),
  rep(4, times = 8)
)

y <- as.array(blood$coagulation)

data_list <- list(N = N, J = J, group = group, y = y)

fit <- stan(
  model_code = model,                 # Stan model file
  data = data_list,
  iter = 1000,           # Number of iterations
  chains = 3,            # Number of chains
  warmup = 500,         # Number of warmup samples
  thin = 1               # Thinning interval
)

print(fit)