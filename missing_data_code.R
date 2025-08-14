#==============================================================================#
# Missing Data Handling Code: FIML, Multiple Imputation & MNAR Sensitivity     #
# by: David Loeb                                                               #
#==============================================================================#

# For explanations of this code see "Missing_Data_Code_and_Explanations.html"

# ==== Setup ===================================================================

if (!require("pak")) install.packages("pak")  # installs pak if needed
pak::pkg_install(  # install/update pkgs if need; does nothing if u have latest
  c("dplyr", "tibble", "lavaan", "blimp-stats/rblimp", "broom", "ggplot2"),
  ask = F  # <- change to `ask = T` if u want to be asked to install/update
)
library(dplyr)

# Simulate data
set.seed(6)
df <- tibble::tibble(  # create data
  id = 1:660,
  income = rnorm(n = 660, mean = 80, sd = 15),
  male = rbinom(n = 660, size = 1, prob = .5),
  test_pre = .1 * income + -.5 * male + rnorm(660, mean = 9, sd = 3),
  treat = c(rep(0, 336), rep(1, 324)),
  test_post = test_pre + .02 * income + 3 * treat + rnorm(660, mean = 3, sd = 1),
  miss_covariates = ifelse(  # create missingness
    percent_rank(.01 * income + rnorm(660, sd = 5)) <= .1, 1, 0
  ),
  miss_test_pre = ifelse(
    percent_rank(treat + rnorm(660, sd = 3)) >= .85, 1, 0
  ),
  miss_test_post = ifelse(
    percent_rank(.2 * test_pre + .5 * treat + rnorm(660)) <= .2, 1, 0
  )
) |> 
  mutate(  # inject missingness
    income = ifelse(miss_covariates == 1, NA, income),
    male = ifelse(miss_covariates == 1, NA, male),
    test_pre = ifelse(miss_test_pre == 1, NA, test_pre),
    test_post = ifelse(miss_test_post == 1, NA, test_post)
  )

# ==== FIML ====================================================================

mod_fiml <- lavaan::sem(  # Run model
  model = "test_post ~ treat + test_pre + income + male",
  data = df,
  missing = "fiml.x"
)
lavaan::summary(mod_fiml)  # Results
broom::tidy(mod_fiml)

# ==== Multiple Imputation =====================================================

## --- Run model ---------------------------------------------------------------

mod_imp <- rblimp::rblimp(
  data = as.data.frame(df),
  ordinal = "treat male",  # binary or ordinal variables
  # nominal = ,  # if you have multinomial variables, specify here
  fixed = "treat",  # variables with no missing data
  model = paste0(  
    "focal.model: test_post ~ treat test_pre income male; ",  # regression model
    "predictor.model: test_pre income male ~ treat;"  # pred imputation mods
  ),
  seed = 6,  # random seed
  burn = 5000,  # number of burn-in iterations (for achieving convergence)
  iter = 5000,  # number of modeling iterations (for results)
  nimps = 20,  # number of imputated data sets to save
  chains = 4  # may want to increase if you have more CPU threads on comp
)

## --- Check Bayesian results --------------------------------------------------

rblimp::output(mod_imp)  # or equivalently: mod_imp@output
rblimp::estimates(mod_imp)
data.frame(mod_imp@estimates)
mod_imp@psr

## --- Extract imputed datasets & save -----------------------------------------

# Extract the datasets from the Blimp model object
df_mult_imp <- purrr::map(
  1:20,
  \(x) mod_imp@imputations[[x]] |> mutate(imp_num = x)  
) |>
  purrr::list_rbind()  # bind all individual datasets into one large df

# Save
saveRDS(df_mult_imp, "imputed_data.rds")
# arrow::write_parquet(df_mult_imp, "imputed_data.parquet")
# readr::write_csv(df_mult_imp, "imputed_data.csv")

## --- Analyze multiply imputed data & get final results -----------------------

df_mult_imp <- readRDS("imputed_data.rds")  # load saved imputed data if need
n_imps <- 20
results_mult_imp <- 1:n_imps |> 
  purrr::map(  # Run individual regression models
    \(x) lm(
      test_post ~ treat + test_pre + income + male, 
      data = filter(df_mult_imp, imp_num == x)
    )
  ) |> 
  purrr::map(  # Put model results into dataframes & include model degrees frdm
    \(x) broom::tidy(x) |> mutate(df = df.residual(x))
  ) |>  
  purrr::imap(\(x, y) mutate(x, imp_num = y)) |>
  purrr::list_rbind() |>
  mutate(variance = std.error^2) |> 
  summarise(  # Compute final results for each term
    .by = term,
    est_final = mean(estimate),
    within_imp_var = mean(variance),
    btwn_imp_var = (sum((estimate - mean(estimate))^2)) / (n_imps - 1),
    total_var = within_imp_var + (1 + 1 / n_imps) * btwn_imp_var,
    se_final = sqrt(total_var),
    
    # Variance from missing data
    prop_var_from_missing = (btwn_imp_var + btwn_imp_var / n_imps) / total_var,
    incrs_var_from_missing = (btwn_imp_var + btwn_imp_var / n_imps) 
      / within_imp_var,
    
    # Degrees of freedom & p-val for hypothesis test
    df_old = (n_imps - 1) / prop_var_from_missing,
    df_complete = mean(df),
    df_observed = (df_complete + 1) / (df_complete + 3) * df_complete 
    * (1 - prop_var_from_missing),
    df_final = (df_old * df_observed) / (df_old + df_observed),
    p = (1 - pt(abs(est_final / se_final), df_final)) * 2
  ) |> 
  relocate(se_final, .after = est_final) |> relocate(p, .after = se_final)

# ==== Sensitivity Analysis: MNAR Pattern Mixture Model ========================

## --- Setup -------------------------------------------------------------------

# Create data
set.seed(6)
df <- tibble::tibble(
  id = 1:660,
  income = rnorm(n = 660, mean = 80, sd = 15),
  treat = c(rep(0, 336), rep(1, 324)),
  test_full = .1 * income + rnorm(660, mean = 50, sd = 10),
  miss_test = ifelse(
    percent_rank(.3 * test_full + 2.5 * (1 - treat) + rnorm(660, sd = 4)) >= .7,
    1, 0
  ),
  test_obs = ifelse(miss_test == 1, NA, test_full)
)

# Estimate "main" model and check results
mod_lwd <- lm(test_obs ~ treat + income, data = df)
broom::tidy(mod_lwd)

## --- Sensitivity modeling ----------------------------------------------------

# Compute sample proportions with observed and missing data
pct_ctrl_miss <- round(mean(df$miss_test[df$treat == 0]), 9)
pct_ctrl_obs <- 1 - pct_ctrl_miss
pct_treat_miss <- round(mean(df$miss_test[df$treat == 1]), 9)
pct_treat_obs <- 1 - pct_treat_miss

# Run sensitivity model
mod_sens <- rblimp::rblimp(
  data = as.data.frame(df),
  ordinal = "treat miss_test",
  fixed = "treat income",
  model = paste0(  # Regression model
    "test_obs ~ 1@int_obs miss_test@int_diff ",  # intercept
      "treat@ate_obs treat*miss_test@ate_diff ",  # treatment effect
      "income; ",  # covariate
    "test_obs ~~ test_obs@resid_var; "  # residual variance
  ),
  parameters = paste0(
    # Pre-specify the missing indicator regression coefficients
    "effect_size_diff_ctrl = .3; ",  # Cohen's D effect size diffs
    "effect_size_diff_treat = .3; ",
    "int_diff = effect_size_diff_ctrl * sqrt(resid_var); ",  # fixed reg coefs
    "ate_diff = effect_size_diff_treat * sqrt(resid_var); ",
    
    # Compute the new ATE under MNAR assumption
    "mean_ctrl_miss = int_obs + int_diff; ",  # missing ctrl group mean
    "mean_treat_obs = int_obs + ate_obs; ",  # observed treat group mean
    "mean_treat_miss = mean_ctrl_miss + ate_obs + ate_diff; ",  # miss treat mean

    "mean_ctrl = int_obs *", pct_ctrl_obs,  # control group
      "+ mean_ctrl_miss *", pct_ctrl_miss, "; ",
    "mean_treat = mean_treat_obs *", pct_treat_obs,  # treatment group
      "+ mean_treat_miss *", pct_treat_miss, "; ",

    "ate_mnar = mean_treat - mean_ctrl; "  # final sensitivity ATE estimate
  ),
  seed = 6,
  burn = 5000,
  iter = 5000
)

mod_sens@estimates  # check results

## --- Testing sensitivity across a range of differences -----------------------

# Set up the modeling function
test_mnar_sensitivity <- function(d_ctrl, d_treat) {
  mod_sens <- rblimp::rblimp(
    data = as.data.frame(df),
    ordinal = "treat miss_test",
    fixed = "treat income",
    model = paste0(
      "test_obs ~ 1@int_obs miss_test@int_diff ",
      "treat@ate_obs treat*miss_test@ate_diff ",
      "income; ", 
      "test_obs ~~ test_obs@resid_var; "
    ),
    parameters = paste0(
      "effect_size_diff_ctrl = ", d_ctrl, "; ",
      "effect_size_diff_treat = ", d_treat, "; ",
      "int_diff = effect_size_diff_ctrl * sqrt(resid_var); ",
      "ate_diff = effect_size_diff_treat * sqrt(resid_var); ",
      
      "mean_ctrl_miss = int_obs + int_diff; ",
      "mean_treat_obs = int_obs + ate_obs; ",
      "mean_treat_miss = mean_ctrl_miss + ate_obs + ate_diff; ",

      "mean_ctrl = int_obs *", pct_ctrl_obs,
      "+ mean_ctrl_miss *", pct_ctrl_miss, "; ",
      "mean_treat = mean_treat_obs *", pct_treat_obs,
      "+ mean_treat_miss *", pct_treat_miss, "; ",

      "ate_mnar = mean_treat - mean_ctrl; "
    ),
    seed = 6,
    burn = 5000,
    iter = 5000
  )
  data.frame(mod_sens@estimates)["Parameter: ate_mnar", ]
}

# Create the effect size df to use as input to the function
eff_sizes <- tidyr::expand_grid(
  d_ctrl = seq(-1, 1, by = .25),
  d_treat = seq(-.5, .5, by = .25)
)

# Run function across each row of the df & combine results into one df
results_sens <- purrr::pmap(eff_sizes, ~ test_mnar_sensitivity(..1, ..2)) |> 
  purrr::list_rbind() |> 
  bind_cols(eff_sizes)

# Add a column specifying the treat-control effect size relation
results_sens <- results_sens |> 
  mutate(
    d_ctrl_plot = case_when(
      d_treat == -.5 ~ d_ctrl - .05,
      d_treat == -.25 ~ d_ctrl - .025,
      d_treat == 0 ~ d_ctrl,
      d_treat == .25 ~ d_ctrl + .025,
      d_treat == .5 ~ d_ctrl + .05
    ),
    d_treat = factor(d_treat, levels = c(.5, .25, 0, -.25, -.5))
  )

# Plot results
library(ggplot2)
ggplot(results_sens, aes(d_ctrl_plot, Estimate, color = d_treat)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = `X2.5.`, ymax = `X97.5.`), width = .05) +
  guides(
    color = guide_legend(
      "ATE Effect Size Diff\nfor People With\nMissing vs Observed\nData"
    )
  ) +
  theme_minimal() +
  theme(
    text = element_text(family = "Roboto", size = 12),
    panel.grid = element_blank(),
    plot.title = element_text(hjust = .5),
    plot.subtitle = element_text(hjust = .5),
    plot.caption = element_text(hjust = 0, size = 8)
  ) +
  geom_hline(yintercept = 0, lty = "dashed") +
  xlab(paste0(
    "Outcome Effect Size Diff for Control Group Members",
    "\nWith Missing vs Observed Data"
  )) +
  ylab("ATE") + 
  labs(
    title = "ATE Estimates Under MNAR Assumptions",
    subtitle = paste0(
      "Results are sensitive to people with missing data scoring",
      "\n.25 SDs higher than those with observed data on average,",
      "\nor having treatment effects .25 SDs smaller on average."
    ),
    caption = paste0(
      "Note: jitter was added to x-axis to ease visualization; all points ",
      "grouped together on\nx-axis fall on the same x-value, i.e. the numbers ",
      "between -1 and 1 in increments of 0.25."
    )
  ) +
  annotate(
    "segment", x = -.04, y = 5, xend = -.01, yend = 3.6751 + .15,
    arrow = arrow(type = "closed", length = unit(.02, "npc"))
  ) +
  annotate(
    "text", x = .2, y = 5.4, label = "Primary model result",
    size = 4, family = "Roboto"
  )
