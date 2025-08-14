pacman::p_load(dplyr, tibble)
set.seed(8)
df <- tibble(
  # Create variables
  y_pre = rnorm(660, 10, 3),
  x = rnorm(660, 40, 15) + 4 * y_pre,
  treat = c(rep(0, 336), rep(1, 324)),
  treat_name = ifelse(treat == 1, "Treatment", "Control"),
  y_post = y_pre + rnorm(660, .5, .1) + treat * 3,
  # Create missing data indicators under diff mechanisms
  miss_mcar = sample(c(0,1), 660, replace = T, prob = c(.75, .25)),
  miss_mar.x = ifelse(percent_rank(x) <= .25, 1, 0),
  miss_mar.x_rand = ifelse(percent_rank(x + rnorm(660, sd = 5)) <= .25, 1, 0),
  miss_mar.x_trt = ifelse(percent_rank(x + 20 * treat) <= .25, 1, 0),
  miss_mar.y_pre = ifelse(percent_rank(y_pre) >= .75, 1, 0),
  miss_mar.y_pre_trt = ifelse(percent_rank(y_pre + 3 * treat) >= .75, 1, 0),
  miss_mnar = ifelse(percent_rank(y_post) <= .25, 1, 0),
  miss_mnar_rand = ifelse(percent_rank(y_post + rnorm(660, sd = 1)) <= .25, 1, 0),
  # Create "observed" vars for plotting, where vals for missing var = 0
  y_plt_mcar = ifelse(miss_mcar == 1, 0, y_post),
  x_plt_mcar = ifelse(miss_mcar == 1, 0, x),
  y_plt_mar.x = ifelse(miss_mar.x == 1, 0, y_post),
  y_plt_mar.x_rand = ifelse(miss_mar.x_rand == 1, 0, y_post),
  x_plt_mar.x = ifelse(miss_mar.x == 1, 0, x),
  y_plt_mar.x_trt = ifelse(miss_mar.x_trt == 1, 0, y_post),
  x_plt_mar.x_trt = ifelse(miss_mar.x_trt == 1, 0, x),
  y_plt_mar.y_pre = ifelse(miss_mar.y_pre == 1, 0, y_post),
  x_plt_mar.y_pre = ifelse(miss_mar.y_pre == 1, 10, x),
  y_plt_mnar = ifelse(miss_mnar == 1, 0, y_post),
  x_plt_mnar = ifelse(miss_mnar == 1, 0, x),
  # Create imputed vars under diff strategies (MAR only)
  y_imp_mean.x = ifelse(
    miss_mar.x == 1, mean(y_post[miss_mar.x == 0]), y_post
  ),
  x_imp_mean.x = ifelse(miss_mar.x == 1, mean(x), x),
  y_imp_mean.x_rand = ifelse(
    miss_mar.x_rand == 1, mean(y_post[miss_mar.x_rand == 0]), y_post
  ),
  x_imp_mean.y_pre = ifelse(
    miss_mar.y_pre == 1, mean(x[miss_mar.y_pre == 0]), x
  ),
  y_imp_mean.x_trt = ifelse(miss_mar.x_trt == 1, mean(y_post), y_post),
  x_imp_mean.x_trt = ifelse(miss_mar.x_trt == 1, mean(x), x),
  x_imp_mean.y_pre_trt = ifelse(
    miss_mar.y_pre_trt == 1, mean(x[miss_mar.y_pre_trt == 0]), x
  )
) |>
  mutate(  # within treat group ranks of variable values
    y_pctile_trt_grp = percent_rank(y_post),
    x_pctile_trt_grp = percent_rank(x),
    y_mean_trt_grp = mean(y_post),
    x_mean_trt_grp = mean(x),
    .by = treat
  ) |>
  mutate(
    y_miss_mnar_treat = case_when(  # for initial extreme example
      treat == 0 & y_pctile_trt_grp <= .34 ~ 1,
      treat == 1 & y_pctile_trt_grp <= .08 ~ 1,
      .default = 0
    ),
    y_imp_mean_trt_grp.x = ifelse(miss_mar.x == 1, y_mean_trt_grp, y_post),
    x_imp_mean_trt_grp.x = ifelse(miss_mar.x == 1, x_mean_trt_grp, x),
    y_imp_mean_trt_grp.x_trt = ifelse(miss_mar.x_trt == 1, y_mean_trt_grp, y_post),
    x_imp_mean_trt_grp.x_trt = ifelse(miss_mar.x_trt == 1, x_mean_trt_grp, x),
  )

preds_mod_x.x <- predict(lm(y_post ~ x, filter(df, miss_mar.x == 0)), df)
preds_mod_trt_x.x <- predict(lm(y_post ~ treat + x, filter(df, miss_mar.x == 0)), df)
preds_mod_sep_regs.x <- c(
  predict(
    lm(y_post ~ x, filter(df, miss_mar.x == 0 & treat == 0)), df[df$treat == 0, ]
  ),
  predict(
    lm(y_post ~ x, filter(df, miss_mar.x == 0 & treat == 1)), df[df$treat == 1, ]
  )
)

preds_mod_x.x_trt <- predict(lm(y_post ~ x, filter(df, miss_mar.x_trt == 0)), df)
preds_mod_trt_x.x_trt <- predict(lm(y_post ~ treat + x, filter(df, miss_mar.x_trt == 0)), df)
preds_mod_sep_regs.x_trt <- c(
  predict(
    lm(y_post ~ x, filter(df, miss_mar.x_trt == 0 & treat == 0)), df[df$treat == 0, ]
  ),
  predict(
    lm(y_post ~ x, filter(df, miss_mar.x_trt == 0 & treat == 1)), df[df$treat == 1, ]
  )
)

df <- df |> 
  mutate(
    y_imp_reg_x.x = ifelse(miss_mar.x == 1, preds_mod_x.x, y_post),
    y_imp_reg_trt_x.x = ifelse(miss_mar.x == 1, preds_mod_trt_x.x, y_post),
    y_imp_reg_sep_regs.x = ifelse(miss_mar.x == 1, preds_mod_sep_regs.x, y_post),
    y_imp_reg_x.x_trt = ifelse(miss_mar.x_trt == 1, preds_mod_x.x_trt, y_post),
    y_imp_reg_trt_x.x_trt = ifelse(miss_mar.x_trt == 1, preds_mod_trt_x.x_trt, y_post),
    y_imp_reg_sep_regs.x_trt = ifelse(miss_mar.x_trt == 1, preds_mod_sep_regs.x_trt, y_post)
  )

# ==== Generate params for adding correct reg lines to plots ===================

# summary(lm(y_post ~ x, df))  # Complete: 0.104463, 0.005375
# summary(lm(y_post ~ treat + x, df))  # RCT Complete
# # (Intercept) 2.225595   0.383151
# #   treat       2.910529   0.181594
# #   x           0.101732   0.004564

## --- Listwise deletion -------------------------------------------------------

# summary(lm(y_post ~ x, filter(df, miss_mar.x_rand == 0)))  # Y: 0.1019, 0.0084
# summary(lm(  # X: 0.0708, 0.0059
#   y_post ~ x, filter(df, miss_mar.y_pre == 0)
# ))
# 
# # RCT
# summary(lm(y_post ~ treat + x, filter(df, miss_mar.x_trt == 0)))  # miss Y
# # (Intercept) 2.096107   0.645322 
# #   treat       3.074620   0.216366 
# #   x           0.102137   0.006821
# summary(lm(  # missing X
#   y_post ~ treat + x, filter(df, miss_mar.y_pre_trt == 0)
# ))
# (Intercept) 4.152437   0.395999
#   treat       1.878091   0.184983
#   x           0.072483   0.004937
# summary(lm(  # missing X & Y, X|Y
#   y_post ~ treat + x, filter(df, miss_mar.x_trt == 0 & miss_mar.y_pre_trt == 0)
# ))
# # (Intercept) 6.007242   0.744251   
# #   treat       1.681257   0.219172   
# #   x           0.052249   0.008169 
# summary(lm(  # missing X & Y, X ind Y
#   y_post ~ treat + x, filter(df, miss_mar.x_trt == 0 & miss_mar.y_pre == 0)
# ))
# # (Intercept) 4.955740   0.677780
# #   treat       2.992607   0.202886
# #   x           0.055940   0.007489

## --- Mean imputation ---------------------------------------------------------

# summary(lm(y_imp_mean.x_rand ~ x, df))  # miss Y: 0.043, 0.005
# summary(lm(y_post ~ x_imp_mean.y_pre, df))  # miss X: 0.071, 0.008
# 
# # RCT
# summary(lm(  # X|Y
#   y_post ~ treat + x_imp_mean.y_pre_trt + miss_mar.y_pre_trt, 
#   filter(df, miss_mar.x_trt == 0)
# ))
# # (Intercept)        6.416058   0.699378
# #   treat              1.464675   0.182905
# #   x_imp_mean         0.048633   0.007802
# 
# summary(lm(  # X ind Y
#   y_post ~ treat + x_imp_mean.y_pre + miss_mar.y_pre, filter(df, miss_mar.x_trt == 0)
# ))
# # (Intercept)      5.449731   0.653460
# #   treat            3.060717   0.168263
# #   x_imp_mean.y_pre 0.053377   0.007323

## --- Reg imputation ----------------------------------------------------------

# summary(lm(y_imp_reg_x.x ~ x, df))  # 0.107462, 0.004749
