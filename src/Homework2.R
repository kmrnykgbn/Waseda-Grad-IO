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
  tidymodels,
  modelsummary,
  ggplot2,
  skimr
)


orig_df <- read_csv("./input/WSDR.csv")
head(orig_df, 30)

### Q1 Answer ###
stats <- orig_df %>%
  dplyr::select(move,price, profit, custcoun) %>%
  skimr::skim(.) %>%
  skimr::yank("numeric") %>%
  dplyr::select(skim_variable, mean, sd, p0, p100) %>%
  dplyr::mutate_at(vars(mean, sd, p0, p100), ~round(., 3)) %>%
  kable(format = "latex")
print(stats)

### Q2 Answer ###

# to factor
df <- orig_df |>
  mutate(upc = factor(upc)) |>
  arrange(week, store, upc)

# compute market share
df <- df |>
  dplyr::group_by(store, week) |>
  mutate(M_t = sum(custcoun),
         tot_quant_t = sum(move),
         s_jt = move / M_t,
         s_0t = (M_t- tot_quant_t) / M_t,
         logit_share = log(s_jt/s_0t)) |>
  ungroup() |>
  select(-c(M_t, tot_quant_t))

# compute iV
df <- df |>
  mutate(whole_p_jt = price * (1 - profit*(0.01)))

# OLS estimation in Berry's logit
model1_OLS <- feols(logit_share ~  price + i(upc), 
                    df, vcov="hetero"
)

# IV estimation in Berry's logit
model1_IV <- feols(logit_share ~  i(upc) | price ~ whole_p_jt, 
                   df, vcov="hetero"
)

# First stage
etable(model1_IV, stage = 1, fitstat=~ . + ivfall + ivwaldall.p,
       signif.code=c("***"=0.01,"**"=0.05,"*"=0.10), 
       style.tex = style.tex("aer"),  tex = TRUE,
       digits=3, digits.stats=3)

# Estimation result
etable(model1_OLS, model1_IV, stage = 2, fitstat=~ . + ivfall + ivwaldall.p,
       signif.code=c("***"=0.01,"**"=0.05,"*"=0.10), 
       style.tex = style.tex("aer"),  tex = TRUE,
       digits=3, digits.stats=3)

### Q4 Answer ###

# compute own elasticity for each market
own_elas_jt <- df |>
  select(week, store, upc, descrip, price, move, s_jt, Brand) |>
  mutate(own_elas = model1_IV$coefficients["fit_price"] * price * (1 - s_jt)) 

# take median of own elasticity by upc
own_elas_df <- own_elas_jt |>
  select(upc, descrip, own_elas) |>
  group_by(upc, descrip) |>
  summarize(own_elas = median(own_elas))

# compute cross elasticity for each market
cross_elas_jt <- df |>
  select(week, store, upc, descrip, price, move, s_jt, Brand) |>
  mutate(cross_elas = (-1) * model1_IV$coefficients["fit_price"] * price * s_jt)

# take median of cross elasticity by upc
cross_elas_df <- cross_elas_jt |>
  select(upc, descrip, cross_elas) |>
  group_by(upc, descrip) |>
  summarize(cross_elas = median(cross_elas))

# create elasticity matrix
J <- nrow(cross_elas_df)
elas_mat_med <- matrix( rep(cross_elas_df$cross_elas, J), nrow = J, ncol = J)
diag(elas_mat_med) <- own_elas_df$own_elas
colnames(elas_mat_med) <- rownames(elas_mat_med) <- as.character(cross_elas_df$descrip)
elas_mat_med

### Q5 Answer ###

# compute markup and marginal cost
compute_markup <- function(own_elas_jt, closs_elas_jt,  owner_mat, is_Multi=FALSE) {
  
  unique_week <- unique(own_elas_jt$week)
  unique_store <- unique(own_elas_jt$store)
  
  market_comb <- expand.grid(week = unique_week, store = unique_store)
  
  # store list of markup and mc
  markup_mc_list <- vector("list", nrow(market_comb))
  
  # calc markup and mc for each market
  for (i in 1:nrow(market_comb)) {
    
    week_i <- market_comb$week[i]
    store_i <- market_comb$store[i]
    
    # subset of a particular market
    sub_own_elas <- own_elas_jt |> filter(week == week_i & store == store_i)
    sub_cross_elas <- cross_elas_jt |> filter(week == week_i & store == store_i)
    
    # number of product
    J <- nrow(sub_cross_elas)
    
    # if multi-product bertrand nash
    if (is_Multi) {
      for (j in 1:J) {
        for (k in 1:J) {
          # if j's brand is equal to k's brand, take 1 
          if (sub_cross_elas$Brand[j] == sub_cross_elas$Brand[k]){
            owner_mat[j,k] = 1
          }
        }
      }
    }
    
    # calculate elasticity matrix
    elas_mat <- matrix(rep(sub_cross_elas$cross_elas, J), nrow = J, ncol = J)
    diag(elas_mat) <- sub_own_elas$own_elas
    
    # calculate matrix S(p)
    price_mat <- matrix(rep((1/sub_cross_elas$price), J), nrow = J, ncol = J)
    share_mat <- t(matrix(rep(sub_cross_elas$s_jt, J), nrow = J, ncol = J))
    S_p_mat <- (-1) * elas_mat * price_mat * share_mat
    
    # calculate markup for a paticular market
    D_p_mat <- owner_mat * S_p_mat
    markup <- solve(D_p_mat) %*% sub_cross_elas$s_jt
    mc <- sub_cross_elas$price - markup
    
    return_df <- data.frame(week = market_comb$week[i],
                            store = market_comb$store[i],
                            Brand = sub_cross_elas$Brand, 
                            upc = sub_cross_elas$upc, 
                            descrip = sub_cross_elas$descrip, 
                            price = sub_cross_elas$price, 
                            move = sub_cross_elas$move, 
                            markup = markup, 
                            mc = mc) %>%
      mutate(markup_dev_p = markup/price)
    
    markup_mc_list[[i]] <- return_df
  }
  
  result_df <- do.call(rbind, markup_mc_list)
  return(result_df)
} 

# compute single product Nash equilibrium for each market
J <- nrow(cross_elas_df)
Omega_sin <- diag(J)
markup_sin <- compute_markup(own_elas_jt, cross_elas_jt, Omega_sin)
markup_sin_med <- markup_sin |>
  select(Brand, descrip, markup, markup_dev_p, mc) |>
  group_by(Brand, descrip) |>
  summarize(markup = median(markup), markup_dev_p = median(markup_dev_p), mc = median(mc))

# compute multi product Nash equilibrium for each market
Omega_multi <- matrix(0, J, J)
markup_multi <- compute_markup(own_elas_jt, cross_elas_jt, Omega_multi, is_Multi=TRUE)
markup_multi_med <- markup_multi |>
  select(Brand, descrip, markup, markup_dev_p, mc) |>
  group_by(Brand, descrip) |>
  summarize(markup = median(markup), markup_dev_p = median(markup_dev_p), mc = median(mc))

# joint pricing of all brands
J <- nrow(cross_elas_df)
Omega_joint <- matrix(1, nrow = J, ncol = J)
markup_joint <- compute_markup(own_elas_jt, cross_elas_jt, Omega_joint)
markup_joint_med <- markup_joint |>
  select(Brand, descrip, markup, markup_dev_p, mc) |>
  group_by(Brand, descrip) |>
  summarize(markup = median(markup), markup_dev_p = median(markup_dev_p), mc = median(mc))

res_markup <- left_join(markup_sin_med |> select(Brand, descrip, markup, markup_dev_p),
                        markup_multi_med |> select(Brand, descrip, markup, markup_dev_p),
                        by = c("Brand", "descrip")) |>
  left_join(markup_joint_med |> select(Brand, descrip, markup, markup_dev_p),
            by = c("Brand", "descrip"), suffix = c(".multi", ".joint"))

res_markup <- res_markup %>%
  knitr::kable(format = "latex", booktabs = TRUE, caption = "Estimated markup") %>%
  kableExtra::kable_styling(latex_options = "hold_position")
cat(res_markup)

res_mc <- left_join(markup_sin_med |> select(Brand, descrip, mc),
                    markup_multi_med |> select(Brand, descrip, mc),
                    by = c("Brand", "descrip")) |>
  left_join(markup_joint_med |> select(Brand, descrip, mc),
            by = c("Brand", "descrip"), suffix = c(".multi", ".joint"))

res_mc <- res_mc %>%
  knitr::kable(format = "latex", booktabs = TRUE, caption = "Estimated markup") %>%
  kableExtra::kable_styling(latex_options = "hold_position")
cat(res_mc)

compute_price_eq = function(df, init_p, owner_mat, alpha_hat, is_Multi =TRUE){
  
  # set the threshold
  lambda <- 1e-6
  df$price <- init_p
  distance <- 10000
  while (distance > lambda) {
    
    # compute own elasticity for each market
    own_elas_jt <- df |>
      select(week, store, upc, descrip, price, move, s_jt, Brand, mc) |>
      mutate(own_elas = alpha_hat * price * (1 - s_jt)) 
    
    # compute cross elasticity for each market
    cross_elas_jt <- df |>
      select(week, store, upc, descrip, price, move, s_jt, Brand, mc) |>
      mutate(cross_elas = (-1) * alpha_hat * price * s_jt)
    
    unique_week <- unique(cross_elas_jt$week)
    unique_store <- unique(cross_elas_jt$store)
    
    market_comb <- expand.grid(week = unique_week, store = unique_store)
    
    # store list of p
    p_k_list <- vector("list", nrow(market_comb))
    
    # calc eta_(p) for each market
    for (i in 1:nrow(market_comb)) {
      
      week_i <- market_comb$week[i]
      store_i <- market_comb$store[i]
      
      # subset of a particular market
      sub_own_elas <- own_elas_jt |> filter(week == week_i & store == store_i)
      sub_cross_elas <- cross_elas_jt |> filter(week == week_i & store == store_i)
      
      # number of product
      J <- nrow(sub_cross_elas)
      
      # if multi-product bertrand nash
      if (is_Multi) {
        for (j in 1:J) {
          for (k in 1:J) {
            # if j's brand is equal to k's brand, take 1 
            if (sub_cross_elas$Brand[j] == sub_cross_elas$Brand[k]){
              owner_mat[j,k] = 1
            }
          }
        }
      }
      
      # calculate elasticity matrix
      elas_mat <- matrix(rep(sub_cross_elas$cross_elas, J), nrow = J, ncol = J)
      diag(elas_mat) <- sub_own_elas$own_elas
      
      # calculate matrix S(p)
      price_mat <- matrix(rep((1/sub_cross_elas$price), J), nrow = J, ncol = J)
      share_mat <- t(matrix(rep(sub_cross_elas$s_jt, J), nrow = J, ncol = J))
      S_p_mat <- (-1) * elas_mat * price_mat * share_mat
      
      # calculate p_k for a paticular market
      D_p_mat <- owner_mat * S_p_mat
      kishida <- sub_cross_elas$mc
      p_k_list[[i]] <- sub_cross_elas$mc + solve(D_p_mat) %*% sub_cross_elas$s_jt
    }
    p_k <- do.call(rbind, p_k_list)
    distance <- max(abs(p_k - df$price))
    df$price <- p_k
  }
  return(p_k)
}

# compute price equilibrium
df <- left_join(df, markup_multi |> select(store, week, upc, mc),
                by = c("store", "week", "upc"))
init_p <- rep(0.05, nrow(df))
compute_price_eq(df, init_p, Omega_multi, model1_IV$coefficients["fit_price"])