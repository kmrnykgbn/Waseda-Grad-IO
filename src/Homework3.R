rm(list = ls())
gc()
set.seed(1)
options(digits=5) 
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  plm,
  ggplot2,
  tidyverse,
  fixest,
  knitr,
  kableExtra,
  stargazer,
  tidymodels,
  modelsummary,
  ggplot2,
  skimr,
  gmm
)

### Data
orig_df <- read_csv("./input/GMdata.csv")
summary(orig_df)

## Question 1-a

### Create unbalanced and balanced data

# unbalanced panel
unbalanced_df <- orig_df |>
  mutate(d357 = as.factor(if_else(sic3 ==357, 1, 0)),
         yr = as.factor(yr)) |>
  as.data.frame()

# balanced panel
balanced_df <-orig_df |>
  mutate(d357 = as.factor(if_else(sic3 ==357, 1, 0)),
         yr = as.factor(yr)) |>
  group_by(index) |>
  mutate(n = n()) |>
  filter(n == 4) |>
  ungroup() |>
  select(-c(n)) |>
  as.data.frame()

### Summary statistics (unbalanced)
# TODO concat two 

# summary of unbalanced panel
unbalanced_df %>%
  select(ldsal, lemp, ldnpt, ldrst, ldrnd, ldinv) %>%
  stargazer::stargazer(., type = "text")


### Summary statistics (balanced)

# summary of balanced panel
balanced_df  %>%
  select(ldsal, lemp, ldnpt, ldrst, ldrnd, ldinv)  %>%
  stargazer::stargazer(., type = "text")


## Question 1-b and Question 1-c

# TODO coliniality of sic 357 dummy

# eq(0.1) in the balanced panel
bal_01 <- feols(ldsal ~  lemp + ldnpt + ldrst + d357 | yr, 
                balanced_df, vcov="hetero"
)

# eq(0.2) in the balanced panel
bal_02 <- feols(ldsal ~  lemp + ldnpt + ldrst + d357 | yr + index, 
                balanced_df, vcov="hetero"
)

# eq(0.1) in the unbalanced panel
unbal_01 <- feols(ldsal ~  lemp + ldnpt + ldrst + d357 | yr, 
                  unbalanced_df, vcov="hetero"
)

# eq(0.2) in the unbalanced panel
unbal_02 <- feols(ldsal ~  lemp + ldnpt + ldrst + d357 | yr + index, 
                  unbalanced_df, vcov="hetero"
)
# Estimation result
etable(bal_01, bal_02, unbal_01, unbal_02, 
       signif.code=c("***"=0.01,"**"=0.05,"*"=0.10), 
       #style.tex = style.tex("aer"),  tex = TRUE,
       digits=3, 
       digits.stats=3)

## Question 2-a

# First stage
model <- feols(ldsal ~  lemp + yr + yr:d357 + poly(ldnpt, ldrst, ldinv, degree=2), 
               unbalanced_df, vcov="hetero"
)
# Estimation result
etable(model, 
       signif.code=c("***"=0.01,"**"=0.05,"*"=0.10), 
       #style.tex = style.tex("aer"),  tex = TRUE,
       digits=3, 
       digits.stats=3)

## Question 2-b

### First stage
# estimated coef
beta_1 = coef(model)['lemp']
coef_yr78 = coef(model)['yr78']
coef_yr83 = coef(model)['yr83']
coef_yr88 = coef(model)['yr88']
coef_yr73_d357 = coef(model)['yr73:d3571']
coef_yr78_d357 = coef(model)['yr78:d3571']
coef_yr83_d357 = coef(model)['yr83:d3571']
coef_yr88_d357 = coef(model)['yr88:d3571']

# predict phi_it
unbalanced_df　<- unbalanced_df |>
  mutate(d_yr73 = if_else(yr ==73, 1, 0),
         d_yr78 = if_else(yr ==78, 1, 0),
         d_yr83 = if_else(yr ==83, 1, 0),
         d_yr88 = if_else(yr ==88, 1, 0),
         d357 = as.numeric(d357),
  )　|>
  mutate(phi_hat_it = fitted(model)
         - (beta_1 * lemp)
         - (coef_yr78 * d_yr78)
         - (coef_yr83 * d_yr83)
         - (coef_yr88 * d_yr88)
         - (coef_yr73_d357 * d_yr78 * d357)
         - (coef_yr78_d357 * d_yr83 * d357)
         - (coef_yr83_d357 * d_yr88 * d357)
         - (coef_yr88_d357 * d_yr88 * d357)
  ) |>
  group_by(index) |>
  mutate(phi_hat_it_m1 = lag(phi_hat_it, 1),
         ldnpt_m1 = lag(ldnpt, 1),
         ldrst_m1 = lag(ldrst, 1)) |>
  ungroup() |>
  drop_na(phi_hat_it_m1, ldnpt_m1, ldrst_m1)

### Second stage
func_2nd_stage_GMM <- function(theta, df, W) {
  # coefficient
  beta_0 <- theta[1]
  beta_2 <- theta[2]
  beta_3 <- theta[3]
  beta_h <- theta[4]
  beta_h2 <- theta[5]
  
  df <- df |>
    mutate(
      h_it_m1 = beta_h * (phi_hat_it_m1 - beta_0 - beta_2 * ldnpt_m1 
                          - beta_3 * ldrst_m1),
      xi_it = phi_hat_it - beta_0 - beta_2 * ldnpt 
      - beta_3 * ldrst 
      - beta_h * h_it_m1 - beta_h2 * (h_it_m1)^2
    )
  moment_df <-
    cbind(
      df$xi_it * df$ldnpt,
      df$xi_it * df$ldrst
    )
  moment <- apply(moment_df, 2, mean)
  obj_func <- t(moment) %*% W %*% moment
  
  return(obj_func)
}

res_GMM <-
  optim(
    par = c(0.5, 0.5, 0.5, 0.5, 0.5),
    fn = func_2nd_stage_GMM,
    method = "L-BFGS-B",
    df = unbalanced_df,
    W = diag(2)
  )