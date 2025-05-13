library(plotly)

dat <- strsplit(c("0/20","0/20","0/20","0/20","0/20","0/20","0/20","0/19","0/19","0/19",
       "0/19","0/18","0/18","0/17","1/20","1/20","1/20","1/20","1/19","1/19",
       "1/18","1/18","2/25","2/24","2/23","2/20","2/20","2/20","2/20","2/20",
       "2/20","1/10","5/49","2/19","5/46","3/27","2/17","7/49","7/47","3/20",
       "3/20","2/13","9/48","10/50","4/20","4/20","4/20","4/20","4/20","4/20",
       "4/20","10/48","4/19","4/19","4/19","5/22","11/46","12/49","5/20","5/20",
       "6/23","5/19","6/22","6/20","6/20","6/20","16/52","15/47","15/46","9/24"), "/")

y <- as.numeric(sapply(dat, function(x) x[1]))
n <- as.numeric(sapply(dat, function(x) x[2]))

post <- function(a, b){
  lik <- gamma(a + b)/(gamma(a) * gamma(b)) * gamma(a + y) * gamma(b + n - y) / gamma(a + b + n)
  #post <- -1/2 * (a + b)^(-5/2) * prod(lik)
  post <- a * b * (a + b)^(-5/2) * prod(lik)
  return(post)
}

alphas <- seq(0.5, 4.5, by = 0.01)
betas <- seq(5, 25, by = 0.1)

probs <- matrix(0, nrow = length(alphas), ncol = length(betas))
for(i in 1:length(alphas)){
  for(j in 2:length(betas)){
    probs[i,j] <- post(alphas[i], betas[j])    
  }
}

totalprob <- sum(probs)
for(i in 1:length(alphas)){
  for(j in 2:length(betas)){
    probs[i,j] <- probs[i,j] / totalprob
  }
}

plot_ly(x = betas, y = alphas, z = probs, type = "contour")

# Posterior mean
sum(sapply(1:length(betas), function(j){
  sum(probs[,j] * alphas)
}))

sum(sapply(1:length(alphas), function(i){
  sum(probs[i,] * betas)
}))
