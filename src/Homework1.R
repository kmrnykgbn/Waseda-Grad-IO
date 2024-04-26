rm(list = ls())
gc()
set.seed(0)
options(digits = 6)
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  plm,
  ggplot2,
  tidyverse,
  fixest,
  knitr,
  kableExtra,
  tidymodels,
  modelsummary,
  ggplot2
)

df <- read_csv("./input/MROZ_mini.csv")

df <- df %>% mutate(kishida = 300)
#For valitation
model <- feols(lwage ~  1  | educ ~ fatheduc, 
               df, vcov="White"
)
etable(model)

##### Question 1 #####

const <- rep(1, nrow(df))
Mat_X <- as.matrix(cbind(const, df[, 1]))
Mat_y <- as.matrix(df[, 3])

# Computing beta_hat
analytical_beta <- solve(t(Mat_X) %*% Mat_X) %*% (t(Mat_X) %*% Mat_y) 

### Q1 Answer ###
print(paste("beta_0: ", analytical_beta[1]))
print(paste("beta_1: ", analytical_beta[2]))

# Create result table
coefs_Q1_1 <- data.frame(
  variable = c("Constant",  "educ", "Num.Obs."),
  OLS_model = c(format(analytical_beta[1], digits = 6),
               format(analytical_beta[2], digits = 6),
               format(nrow(df), digits = 1))
)

# result table
coefs_table_Q1_1  <- kable(coefs_Q1_1, 
                           caption = "Analytical OLS estimation Result",
                           align = c("l", "c"))  %>%
  kable_styling(full_width = FALSE) 
coefs_table_Q1_1 

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
print(paste("Analytical beta_0: ", analytical_beta[1]))
print(paste("Numerical beta_0: ", result$par[1]))
print(paste("Analytical beta_1: ", analytical_beta[2]))
print(paste("Numerical beta_1: ", result$par[2]))

# Create result table
coefs_Q1_2 <- data.frame(
  variable = c("Constant",  "educ", "Num.Obs."),
  Analycal_OLS = c(format(analytical_beta[1], digits = 6),
                    format(analytical_beta[2], digits = 6),
                    format(nrow(df), digits = 1)),
  Numerical_OLS = c(format(result$par[1], digits = 6),
                format(result$par[2], digits = 6),
                format(nrow(df), digits = 1))
)

# result table
coefs_table_Q1_2  <- kable(coefs_Q1_2, 
                           caption = "OLS estimation Result",
                           align = c("l", "c", "c"))  %>%
  kable_styling(full_width = FALSE) 
coefs_table_Q1_2 

##### Question 3 #####

# compute Asymptotic SE of OLS estimator

# calculating sample mean of S
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
S_hat <- compute_S_hat(analytical_beta, df)

# Computing the asymptotic variance estimator
Avar_est <- solve(S_xx) %*% S_hat %*% solve(S_xx)

# Computing the asymptotic SE for beta_0 and beta_1
Asy_std_beta_0 <- sqrt ((1/nrow(df)) * Avar_est[1])
Asy_std_beta_1 <- sqrt ((1/nrow(df)) * Avar_est[4])

### Q3 Answer ###
print(paste("Asymptotic standard error beta 0: ", Asy_std_beta_0))
print(paste("Asymptotic standard error beta 1: ", Asy_std_beta_1))

coefs_Q1_3 <- data.frame(
  variable = c("Constant", "",  "educ", "", "Num.Obs."),
  IV_model = c(analytical_beta[1],
               paste("(", Asy_std_beta_0, ")"),  
               analytical_beta[2], 
               paste("(", Asy_std_beta_1, ")"),
               nrow(df))
)

# result table
coefs_table_Q1_3  <- kable(coefs_Q1_3,
                           format="latex",
                           caption = "Analytical OLS estimation Result",
                           align = c("l", "c"))  %>%
  kable_styling(full_width = FALSE) 
coefs_table_Q1_3 

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
coefs_Q1_6 <- data.frame(
  variable = c("Constant", "",  "educ", "", "Num.Obs."),
  IV_model = c(numerical_beta_IV[1],
            paste("(", Asy_std_beta_IV_0, ")"),  
            numerical_beta_IV[2], 
            paste("(", Asy_std_beta_IV_1, ")"),
            nrow(df))
)

# result table
coefs_table_Q1_6  <- kable(coefs_Q1_6, 
      caption = "IV estimation Result",
      align = c("l", "c"))  %>%
  kable_styling(full_width = FALSE) 

coefs_table_Q1_6 

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
      
      # Sum of d_{ijk} P_{ijk} (Note that we skip calculation when d_{ijk} = 0)
      add_lf <- log(prob)
      lf <- lf + add_lf
    }
  }
  return(lf)
} 

# OLS estimation to get the initial value
reg_df <- kntk_df %>% 
  mutate(y = if_else(choice == 1, 1, 0), 
         x = if_else(choice == 1, p_kino, 0),
         x = if_else(choice == 2, p_take, x))
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

# Create result table
coefs_Q2_4 <- data.frame(
  variable = c("alpha kinoko",  "alpha takenoko", "beta", "Num.Obs."),
  estimates = c(format(MLE_res$par[1], digits = 6),
                   format(MLE_res$par[2], digits = 6),
                   format(MLE_res$par[3], digits = 6),
                   format(nrow(kntk_df), digits = 1))
)

# result table
coefs_table_Q2_4  <- kable(coefs_Q2_4, 
                           caption = "Estimation Result",
                           align = c("l", "c"))  %>%
  kable_styling(full_width = FALSE) 
coefs_table_Q2_4

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
set.seed(0)
mu <- 2
sigma <- sqrt(2)
R <- 1000

# Define Monte Carlo simulation function
compute_E_X2_hat <- function(R) {
  nu_vec <- rnorm(R, 0, 1) # fixed by seed

  E_X2_hat = 0
  for (nu_r in nu_vec) {
    E_X2_hat = E_X2_hat + (mu + sigma * nu_r)^2
  }
  
  E_X2_hat <- (1/ R) *  E_X2_hat
  return(E_X2_hat)
} 

##### Answer Q3-2 #####
res_E_X2_hat <- compute_E_X2_hat(R)
print(paste("estimated E[X^2]: ", res_E_X2_hat))

# Create result table
coefs_Q3_2 <- data.frame(
  variable = c("E[X^2]",  "R"),
  estimates = c(format(res_E_X2_hat, digits = 6),
                format(R, digits = 1))
)

coefs_table_Q3_2  <- kable(coefs_Q3_2, 
                           caption = "Estimation Result",
                           align = c("l", "c"))  %>%
  kable_styling(full_width = FALSE) 
coefs_table_Q3_2

# check dependency on the number of R
res_df = data.frame(matrix(ncol = 2, nrow = 0))
colnames(res_df) <- c("R", "estimated_E_X2")

i = 0
for (r in seq(from = 100, to = 50000, by = 100)) {
  i = i + 1
  
  res_E_X2_hat_r <- compute_E_X2_hat(r)
  
  # Store result
  res_df[i, 1] <- r
  res_df[i, 2] <- res_E_X2_hat_r
}



# Plot R and estimates
ggplot(res_df, aes(x = R, y = estimated_E_X2)) +
  geom_point() +
  labs(x = "R", y = "estimated_E_X2") +
  theme_minimal()