c_index <- function(new_model, beta) {
  s0 <- 0
  s1 <- 0
  y0 <- (new_model$X) %*% beta
  y1 <- new_model$Y
  n1 <- nrow(y1) - 1
  
  for (i in 1:n1) {
    if (y1[i, 1] < 0) {
      if (y0[i, 1] < 0) {
        s0 <- s0 + 1
      }
      s1 <- s1 + 1
    }
  }
  return(s0/s1)
}


cross_validation <- function(Y=NULL, status=NULL, X=NULL, seed = 1, n = 200, p = 200, maxit = 50, 
                             lambda = seq(0.01,0.21, 0.02)) {
  if (is.null(Y)) {
    old_model <- gen_AFT(n, p, seed)
  } else {
    n <- nrow(X)
    p <- ncol(X)
    old_model <- list(Y, status, X, n, p)
    names(old_model) <- c('Y', 'status', "X", "n", "p")
    maxit = p - 1
  }
  
  LG_Lasso_min <- 1e6
  lambda_min_LASSO <- 2
  LG_min <- 1e6
  lambda_min <- 2
  k_min <- maxit
  
  new_model <- new_sample(old_model)
  fold <- floor(n/10)
  for (j in 1:length(lambda)) {
    LG_Lasso <- 0
    LG <- rep(0, maxit)
    for (i in 1:10) {
      if (i < 10) {
        out <- c(-((fold*i-fold + 1):(fold*i)))
      } else {
        out <- c(-((fold*i-fold + 1):n))
      }
      old_model1 <- list(matrix(old_model$Y[out,]), matrix(old_model$status[out,]), old_model$X[out,], n-length(out), p)
      names(old_model1) <- c('Y', 'status', "X", "n", "p")
      new_model1 <- new_sample(old_model1)
      temp <- l1.reg(t(new_model1$X), new_model1$Y, lambda = lambda[j]*new_model1$n^2)
      temp1 <- matrix(temp$estimate)
      rank <- Klasso(temp1)
      for (t in 1:maxit) {
        temp2 <- top_k(temp1, rank, t)
        LG[t] <- AFT_Likelihood(new_model, temp2) - AFT_Likelihood(new_model1, temp2) + LG[t]
      }
      LG_Lasso <- AFT_Likelihood(new_model, temp1) - AFT_Likelihood(new_model1, temp1) + LG_Lasso
      cat("-------",j,"-------",i,'\n')
    }
    
    for (t in 1:maxit) {
      if (LG[t] < LG_min) {
        lambda_min <- lambda[j]
        k_min <- t
        LG_min = LG[t]
      }
    }
    
    if (LG_Lasso < LG_Lasso_min) {
      lambda_min_LASSO <- lambda[j]
      LG_Lasso_min = LG_Lasso
    }
  }
  
  beta_LASSO <- l1.reg(t(new_model$X), new_model$Y, lambda = lambda_min_LASSO*new_model$n^2)
  beta_min <- l1.reg(t(new_model$X), new_model$Y, lambda = lambda_min*new_model$n^2)
  rank <- Klasso(matrix(beta_min$estimate))
  beta_min <- top_k(matrix(beta_min$estimate), rank, k_min)
  nz <- nzero(matrix(beta_LASSO$estimate))
  
  out_LASSO <- matrix(c(LG_Lasso_min, lambda_min_LASSO, nz, beta_LASSO$estimate), nrow = 1)
  out_APML0 <- matrix(c(LG_min, lambda_min, k_min, beta_min), nrow = 1)
  out <- list(out_LASSO, out_APML0)
  names(out) <- c("LASSO", "APML0")
  return(out)
}


AFT_Likelihood <- function(new_model, beta) {
  n1 <- new_model$n1
  y0 <-  new_model$Y - new_model$X %*% beta
  s <- 0
  for (i in 1:n1) {
    if (y0[i, 1] < 0) {
      s <- s - y0[i, 1]
    }
  }
  return(s/new_model$n)
}

nzero <- function(beta) {
  n <- length(beta)
  k <- 0
  for (i in 1:n) {
    if (abs(beta[i, 1]) > 1e-5) {
      k <- k + 1
    }
  }
  return(k)
} 

Klasso <- function(beta) {
  n <- length(beta)
  p <- nzero(beta)
  beta <- abs(beta)
  c_beta <- matrix(0, nrow = n, ncol = 1)
  for (i in 1:p) {
    k <- which.max(beta)
    c_beta[k, 1] <- p + 1 -i
    beta[k, 1] <- 0
  }
  c_beta <- p + 1 - c_beta
  return(c_beta)
}

top_k <- function(beta, c_beta, k) {
  n <- length(beta)
  n_beta <- beta 
  for (i in 1:n) {
    if (c_beta[i, 1] <= k) {
      n_beta[i, 1] <- beta[i, 1]
    } else {
      n_beta[i, 1] <- 0
    }
  }
  return(n_beta)
}



##Simulation Data (AFT) Initialized
p <- 200
n <- 200 ##n>15

gen_AFT <- function(n, p, seed) {
  set.seed(seed)
  beta <- (-1)^(seq(1,15,1))*2*exp(seq(0,-14,-1)/15)
  beta <- c(beta, rep(0, p-15))
  beta <- matrix(beta)
  ##K <- diag(0.5,5,5) + matrix(0.5,5,5)
  K <- diag(5)
  K1 <- chol(K)
  
  X <- matrix(rnorm(n*p,0,1), nrow = n, ncol = p)
  X <- cbind(X[1:n, 1:5]%*%K1, X[1:n, 6:10]%*%K1, X[1:n, 11:15]%*%K1, X[1:n, 16:p])
  e <- matrix(rnorm(n,0,1), nrow = n, ncol = 1)
  T <- X %*% beta + e
  a <-  quantile(T, 0.55)
  b <-  quantile(T, 0.65)
  C <- matrix(runif(n, a, b))
  
  theta <- matrix(0, nrow = n, ncol = 1)
  for (i in 1:n) {
    if (T[i, 1] < C[i, 1]) {
      theta[i, 1] <- 1
    }
    T[i, 1] <- min(T[i, 1], C[i, 1])
  }
  
  
  AFT_data <- list(T, theta, X, n, p)
  names(AFT_data) <- c('Y', 'status', "X", "n", "p")
  return(AFT_data)
}

##Simulation Data (AFT) Reformed (new_sample)
new_sample <- function(old_sample) {
  T <- old_sample$Y
  theta <- old_sample$status
  X <- old_sample$X
  n <- old_sample$n
  p <- old_sample$p
  
  n1 <- sum(theta) * (n-1)
  XX <- matrix(0, nrow = (n1 + 1), ncol = p)
  YY <- matrix(0, nrow = (n1 + 1), ncol = 1)
  k <- 1
  for (i in 1:n) {
    for (j in 1:n) {
      if (i != j) {
        if (theta[i] == 1) {
          for (z in 1:p) {
            XX[k, z] <- X[i, z] - X[j, z]
            XX[n1+1, z] <- XX[n1+1, z] - XX[k, z]
          }
          YY[k, 1] <- T[i, 1] - T[j, 1]
          k <- k + 1
        }
      }
    }
  }
  
  YY[n1+1, 1] <- 1e8
  new_data <- list(YY, XX, n, n1, p)
  names(new_data) <- c("Y", "X", "n", "n1", "p")
  return(new_data)
}

##CDLasso Test
library("CDLasso")
test <- gen_AFT(200, 1000, 1)
test1 <- new_sample(test)
write.csv(test$Y,file = "y_input_2.csv")
write.csv(test$status,file = "status_input_2.csv")
write.csv(test$X,file = "x_input_2.csv")

getwd()
temp00 <- l1.reg(t(test1$X), test1$Y, lambda = 0.1*200^2)
temp00$nonzeros
temp01 <- matrix(temp00$estimate)
##rqPen Test
library("rqPen")
temp10 <- matrix(LASSO.fit(YY, XX, tau = 0.5, lambda = 0.1, intercept = FALSE, coef.cutoff = 1e-8))

##Simulation Data (AFT)
n <- 200
p <- 300
maxit <- 50
beta_LASSO <- beta
lambda <- seq(0.01,0.21, 0.02)
LASSO_final <- matrix(NA, nrow = 0, ncol = p+3)
APML0_final <- matrix(NA, nrow = 0, ncol = p+3)
for (z in 1:1) {
  two_cv <- cross_validation(seed = z, n = n, p = p, maxit = maxit)
  test <- two_cv$LASSO
  LASSO_final <- rbind(LASSO_final, two_cv$LASSO)
  APML0_final <- rbind(APML0_final, two_cv$APML0)
  cat("-----------------------------",z)
}

LASSO_final <- out_LASSO
APML0_final <- out_APML0
out_LASSO$LG
out_LASSO$Non_zero
out_APML0$Non_zero
test <- out_LASSO$beta
test2 <- out_APML0$beta

##Real Data Test
setwd("C:/Users/micha/Documents/AFT Model with L0 Regularization/AFT_APML0/AFT_APML0/project")
y_input <- read.csv("y_input_1.csv", header = FALSE)
x_input <- as.matrix(read.csv("x_input_1.csv", header = FALSE))
status_input <- read.csv("status_input_1.csv", header = FALSE)

n <- nrow(y_input)
p <- ncol(x_input)
key <- seq(1, n, 1)
n0 <- floor(n/2)
index_LASSO <- c()
index_APML0 <- c()
LASSO_final1 <- matrix(NA, nrow = 0, ncol = p+3)
APML0_final1 <- matrix(NA, nrow = 0, ncol = p+3)
for (i in 1:10) {
  key_ <- sample(key)
  key0 <- key_[1:n0]
  key1 <- key_[(n0+1):n]
  y0 <- matrix(y_input[key0, 1])
  y1 <- matrix(y_input[key1, 1])
  status0 <- matrix(status_input[key0, 1])
  status1 <- matrix(status_input[key1, 1])
  x0 <- x_input[key0, ]
  x1 <- x_input[key1, ]
  two_cv <- cross_validation(Y = y1, status = status1, X = x1)
  model0 <- list(y0, status0, x0, n0, p)
  names(model0) <- c('Y', 'status', "X", "n", "p")
  model00 <- new_sample(model0)
  index_LASSO <- c(index_LASSO, c_index(model0, matrix(two_cv$LASSO[4:(p+3)])))
  index_APML0 <- c(index_APML0, c_index(model0, matrix(two_cv$APML0[4:(p+3)])))
  LASSO_final1 <- rbind(LASSO_final1, two_cv$LASSO)
  APML0_final1 <- rbind(APML0_final1, two_cv$APML0)
}

old_sample <- model0
new_model$n1
beta <- matrix(two_cv$LASSO[4:(p+3)])
nrow(new_model$X)
new_model$n1