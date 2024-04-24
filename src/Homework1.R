rm(list = ls())
gc()
set.seed(0)
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  plm,
  ggplot2,
  tidyverse,
  fixest,
  tidymodels,
  modelsummary
)

df <- read_csv("./input/MROZ_mini.csv")

##### Question 1 #####

const <- rep(1, nrow(df))
Mat_X <- as.matrix(cbind(const, df[, 1]))
Mat_y <- as.matrix(df[, 3])

# Computing beta_hat
numerical_beta <- solve(t(Mat_X) %*% Mat_X) %*% (t(Mat_X) %*% Mat_y) 

### Q1 Answer ###
print(paste("beta_0: ", numerical_beta[1]))
print(paste("beta_1: ", numerical_beta[2]))


##### Question 2 #####

# Definition of the objective function
compute_ols <- function(theta, df) {
  beta_0 <- theta[1]
  beta_1 <- theta[2]
  J <- 0.0
  for (i in 1:nrow(df)) {
    add_J <- (df$lwage[[i]] - beta_0 - (beta_1*df$educ[[i]]))** 2
    J <- J + add_J
  }
  return(J)
} 

# Set initial value of theta
initial_theta <- c(0, 0)

# Minimizing the objective function by using optim function
result <- optim(par = initial_theta, fn = compute_ols, df = df, method = "BFGS")
result

### Q2 Answer ###
print(paste("Numerical beta_0: ", numerical_beta[1]))
print(paste("Analytical beta_0: ", result$par[1]))
print(paste("Numerical beta_1: ", numerical_beta[2]))
print(paste("Analytical beta_1: ", result$par[2]))


##### Question 3 #####

# compute Asymptotic SE of OLS estimator

# Hayashi p.123 calculating sample mean of S
compute_S_hat <- function(theta, df) {
  S_hat <- matrix(0, ncol = 2, nrow = 2)
  for (i in 1:nrow(df)) {
    x_i_mat <- Mat_X[i, ]
    epsilon_hat <- as.numeric(df$lwage[[i]] - t(x_i_mat) %*% theta)
    add_S_hat <- (epsilon_hat ^ 2) * (x_i_mat %*% t(x_i_mat))
    S_hat <- S_hat + add_S_hat
  }
  S_hat <- (1/nrow(df)) * S_hat
  return(S_hat)
}

S_xx <- (1/nrow(df)) * (t(Mat_X) %*% Mat_X)
S_hat <- compute_S_hat(numerical_beta, df)

# Computing the asymptotic variance estimator
Avar_est <- solve(S_xx) %*% S_hat %*% solve(S_xx)

# Computing the asymptotic SE for beta_0 and beta_1
Asy_std_beta_0 <- sqrt ((1/nrow(df)) * Avar_est[1])
Asy_std_beta_1 <- sqrt ((1/nrow(df)) * Avar_est[4])

### Q3 Answer ###
print(paste("Asymptotic standard error beta 0: ", Asy_std_beta_0))
print(paste("Asymptotic standard error beta 1: ", Asy_std_beta_1))

#For valitation
model <- feols(lwage ~  1   + educ, 
               df, vcov="White"
)
etable(model)

##### Question 4 #####
# See the main body

##### Question 5 #####
# See the main body

##### Question 6 #####

# Define Z as IV
Mat_Z <- as.matrix(cbind(const, df[, 2])) # const + futheduc
# get IV estimator
P_Z = Mat_Z %*% solve(t(Mat_Z) %*% Mat_Z) %*% t(Mat_Z)
numerical_beta_IV <- solve(t(Mat_X) %*% P_Z %*% Mat_X) %*% (t(Mat_X) %*% P_Z %*% Mat_y)
numerical_beta_IV

# Compute asymptotic SE of IV estimator based on Hansen p. 354

# Compute epsilon_hat
compute_epsilon_hat <- function(theta, df) {
  epsilon_hat <- 0
  for (i in 1:nrow(df)) {
    x_i_mat <- Mat_X[i, ]
    z_i_mat <- Mat_Z[i, ]
    add_epsilon_hat <- as.numeric(df$lwage[[i]] - t(x_i_mat) %*% theta)
    add_epsilon_hat <- (add_epsilon_hat) ^ 2
    epsilon_hat <- epsilon_hat + add_epsilon_hat
  }
  epsilon_hat <- (1/nrow(df)) * epsilon_hat
  return(epsilon_hat)
}

Q_xz <- (1/nrow(df)) * (t(Mat_X) %*% Mat_Z)
Q_zx <- (1/nrow(df)) * (t(Mat_Z) %*% Mat_X)
Q_zz <- (1/nrow(df)) * (t(Mat_Z) %*% Mat_Z)
epsilon_hat <- compute_epsilon_hat(numerical_beta_IV, df)

# Computing the asymptotic variance estimator
edge_comp <- solve(Q_xz %*% solve(Q_zz) %*% Q_zx)
Avar_est_IV <- edge_comp * epsilon_hat

# Computing the asymptotic SE for beta_0 and beta_1
Asy_std_beta_IV_0 <- sqrt ((1/nrow(df)) * Avar_est_IV[1])
Asy_std_beta_IV_1 <- sqrt ((1/nrow(df)) * Avar_est_IV[4])

### Q6 Answer ###

# print beta IV
print(paste("Numerical beta_IV_0: ", numerical_beta_IV[1]))
print(paste("Numerical beta_IV_1: ", numerical_beta_IV[2]))

# Print Asy SE for beta
print(paste("Asymptotic standard error beta IV 0: ", Asy_std_beta_IV_0))
print(paste("Asymptotic standard error beta IV 1: ", Asy_std_beta_IV_1))

#For valitation
model <- feols(lwage ~  1  | educ ~ fatheduc, 
               df, vcov="White"
)
etable(model)


kntk_df <- read_csv("./input/data_KinokoTakenoko.csv")
head(kntk_df, 30)

# Set price
kntk_df <- kntk_df %>%
    mutate(p_kino = if_else(occasion == "X1", 200, 0),
           p_kino = if_else(occasion == "X2", 170, p_kino),
           p_kino = if_else(occasion == "X3", 240, p_kino),
           p_kino = if_else(occasion == "X4", 200, p_kino),
           p_kino = if_else(occasion == "X5", 200, p_kino),
           p_take = if_else(occasion == "X1", 200, 0),
           p_take = if_else(occasion == "X2", 200, p_take),
           p_take = if_else(occasion == "X3", 200, p_take),
           p_take = if_else(occasion == "X4", 250, p_take),
           p_take = if_else(occasion == "X5", 180, p_take)
           )

##### Question 2-4 #####
compute_loglikelihood <- function(theta, df) {
  alpha_kino <- theta[1]
  alpha_take <- theta[2]
  beta <- theta[3]
  
  # Calculate probability for i, j, k
  df <- df %>%
    mutate(denominator = 1 + exp(alpha_kino - beta*p_kino) + exp(alpha_take - beta*p_take),
           prob = if_else(choice == 1, exp(alpha_kino - beta*p_kino) / denominator, 0),
           prob = if_else(choice == 2, exp(alpha_take - beta*p_take) / denominator, prob),
           prob = if_else(choice == 0, 1 / denominator, prob))

  N = max(df$id)
  lf <- 0.0
  for (i in 1:N) {
    for (k in 1:5) {
      # get P_{ijk}
      prob <- df %>%
        filter(id == i & occasion == paste0("X", k)) %>%
        pull(prob)
      
      add_lf <- log(prob)
      
      lf <- lf + add_lf
    }
  }
  return(lf)
} 

#For valitation
reg_df <- kntk_df %>% 
  mutate(y = if_else(choice == 1, 1, 0), 
         x = if_else(choice == 1, p_kino, 0),
         x = if_else(choice == 2, p_take, x))

# OLS To get the initial value
model <- feols(y ~  1   + x, 
               reg_df, vcov="White"
)
etable(model)

# Maximizing the Log-Likelihood by using optim function
initial_theta <- c(model$coefficients[1], model$coefficients[1], model$coefficients[2])
start.time <- Sys.time()
MLE_res <- optim(par = initial_theta, fn = compute_loglikelihood, df = kntk_df, method='BFGS', control = list(fnscale = -1))
end.time <- Sys.time()
end.time - start.time

#Answer Q2-4
MLE_res
print(paste("alpha kinoko hat: ", MLE_res$par[1]))
print(paste("alpha takenoko hat: ", MLE_res$par[2]))
print(paste("beta: ", MLE_res$par[3]))

# Robustness check about initial value

#for (i in seq(from = 0, to = 1, by = 0.5)) {
  #initial_theta <- c(i, i, i)
  #print(paste("Initial theta: ", c(i, i, i)))
  #start.time <- Sys.time()
  #MLE_res <- optim(par = initial_theta, fn = compute_loglikelihood, df = kntk_df, method='BFGS', control = list(fnscale = -1))
  #end.time <- Sys.time()
  #end.time - start.time
  #MLE_res
#}

##### Question 3-2 #####
mu <- 2
sigma <- 2
N <- 1000 # num of draw

nu_vec <- rnorm(N, mu, sigma)

E_X2_hat = 0
for (nu_r in nu_vec) {
  E_X2_hat = E_X2_hat + (mu + sigma * nu_r)^2
}
