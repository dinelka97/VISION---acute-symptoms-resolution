## Tasks in progress:

  ## -- (1) Make sure the input dataset is correct (think it's better to directly use the main complete dataset)
  


## -- need to make sure that the RMST and the median values are directly from all the observed observations.


## -- this script uses bootstrap methods to obtain point estimates and inference for RMST and RMST differences
## -- focusses only on symptom resolution here

rm(list = ls())

library(doParallel)
library(foreach)
library(doRNG)
library(data.table)
library(nnet)
library(magrittr)
library(tidyverse)
library(flextable)
library(officer)
library(furrr)
library(future)
library(rlang)
library(progressr)
library(tictoc)
library(boot)
library(rms)
library(survival)
library(flextable)
library(officer)

setwd('G:/projects/cfar/COVID19/VISION/manuscripts/acute_symptoms/')

# import data -------------------------------------------------------------
## -- load already prepped dataset from 'km_rmstUV/R'

df_list_ipcw <- readRDS("data/km/v3/dfs_km.RDS")

  ## -- use only resolution for now (use recovery later if required)
data_ipcw_long <- 
  df_list_ipcw$ipcw_resolution %>%
  select(-c(ipcw_d, prop_Ct, ipcw))

stratum_labels = list(
  sex_birth = c("Female", "Male"),
  age_cat_2_derived = c("18-49", "50-64", "65+"),
  age_cat = c("18-34 yrs", "35-49", "50-64", "65+"),
  pi_cat = c("no prior infections", "1+ prior infections"),
  pi_cat_v2 = c("0", "1", "2+"),
  pi_rec = c("no prior infections", "≤ 6 months", "> 6 months"),
  urban_rural_derived = c("urban", "suburban", "rural"),
  bipoc_race_derived = c("White/Caucasian", "BIPOC"),
  race_new_derived = c("Black/AA", "White/Caucasian", "other"),
  bmi_cat_derived = c("Underweight/Normal", "Overweight", "Obese"),
  vax_rec_derived = c("Unvaccinated", "≤ 6 months", "> 6 months"),
  n_vax_derived = c("0", "1-3", "4-5", "6+"),
  comorb_count_cat = c("None recorded", "1", "2+"),
  comorb_count_cat2 = c("heart, HBP, lung, or diabetes", "other"),
  ethn_hispanic_latinx_derived = c("Not Hispanic/Latinx", "Hispanic/Latinx"),
  smoking_cat_derived_derived = c("Never", "Former", "Current"),
  hybrid_imm = c("Unvaccinated, 0", "Unvaccinated, 1+", "Vaccinated ≤ 6 mo., 0", "Vaccinated ≤ 6 mo., 1+",
                 "Vaccinated > 6 mo., 0", "Vaccinated > 6 mo., 1+"),
  NULL = c("Overall")
)

tau <- 14

bootstrap_rmt_ipcw_overall <- function(n_boot, tau, n_cores = parallel::detectCores() - 1) {
  
  cl = makeCluster(n_cores)
  registerDoParallel(cl)
  registerDoRNG(seed = 123)
  
    ## data.table format required
  dt = 
    data_ipcw_long %>% 
      as.data.table()
  
  # covariates_ipcw <- c("age_derived", "sex_birth", "pi_cat", "urban_rural_derived",
  #                      "bipoc_race_derived", "bmi_cat_derived", "vax_rec_derived",
  #                      "comorb_count_cat", "ethn_hispanic_latinx_derived")
  # 
  # ipcw_form = as.formula(paste("C_t ~", paste(c("rcs(age_derived, c(35, 50, 65))", 
  #                                               covariates_ipcw[-1], "as.factor(tte_start)"), collapse = " + ")))
  
  # Same formula as return to school/work
  ipcw_form = as.formula("C_t ~ rcs(age_derived, c(35, 50, 65)) +
                  as.factor(sex_birth) + as.factor(bipoc_race_derived) +
                  as.factor(ethn_hispanic_latinx_derived) + as.factor(vax_rec_derived) +
                  as.factor(pi_cat) + as.factor(comorb_count_cat) + as.factor(tte_start)")

  rmst_list = foreach(i = 1:n_boot, .combine = 'rbind',
                      .packages = c("data.table", "dplyr", "survival", "splines", "rms", "magrittr")) %dopar% {
  # rmst_list <- list()
  # for(i in seq(n_boot)){
    #resample IDs and keep all rows for each ID
    boot_ids <- sample(unique(dt$id), replace = TRUE)
    boot_index <- data.table(id = boot_ids, replicate = seq_along(boot_ids))
    setkey(dt, id)
    boot_dat <- dt[boot_index, allow.cartesian = TRUE]
    
    # boot_dat %<>%
    #   filter(!if_any(all_of(covariates_ipcw), is.na)) ## -- remove any records that are missing BL covariates
    
    #IPCW
    
    fit_Cd = glm(ipcw_form, data = boot_dat, family="binomial"(link="logit"))
    
    ipcw_denom = predict(fit_Cd, data = boot_dat, type="response")
    
    fit_Cn = boot_dat %>%
      group_by(tte_start) %>% ## CHANGED FROM DINELKA'S TTE_FINAL
      summarise(prop_Ct = mean(C_t, na.rm=T)) %>%
      ungroup()
    
    boot_dat %<>%
      mutate(ipcw_d = ipcw_denom) %>%
      left_join(fit_Cn) %>%
      group_by(replicate) %>%
      mutate(ipcw = cumprod(prop_Ct/ipcw_d),
             ipcw = lag(ipcw, default=1)) %>%
      ungroup()
    
    #KM fit
    
    formula_str = "Surv(time = tte_start, time2 = time, event = event_final) ~ 1" ## -- using only right censored

    fit = survfit(as.formula(formula_str), 
                  data = boot_dat, weights = ipcw, robust = T, id = replicate)
    
    idx = which(fit$time <= tau)
    time_points = c(0, fit$time[idx], tau)
    surv_prob = c(1, fit$surv[idx], tail(fit$surv, 1))
    
    rmst = sum((head(surv_prob, -1) + tail(surv_prob, -1)) / 2 * diff(time_points))
    
    median_tte = unname(summary(fit)$table["median"])
    
    data.frame(replicate = i, rmst = rmst, median_tte = median_tte)
    # rmst_list <- list(rmst_list, data.frame(replicate = i, rmst = rmst))
  }
  
  stopCluster(cl)
  return(rmst_list)
  
}

bootstrap_rmt_ipcw_parallel <- function(group_var, n_boot, tau, n_cores = parallel::detectCores() - 1) {
  
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  registerDoRNG(seed = 123)
  
  dt = 
    data_ipcw_long %>% 
    as.data.table()
  
  # covariates_ipcw <- c("age_derived", "sex_birth", "pi_cat", "urban_rural_derived",
  #                      "bipoc_race_derived", "bmi_cat_derived", "vax_rec_derived", "hybrid_imm",
  #                      "comorb_count_cat", "ethn_hispanic_latinx_derived")
  # 
  # ipcw_form = as.formula(paste("C_t ~", paste(c("rcs(age_derived, c(35, 50, 65))", 
  #                                               covariates_ipcw[-1], "as.factor(tte_start)"), collapse = " + ")))
  # Same formula as return to school/work
  ipcw_form = as.formula("C_t ~ rcs(age_derived, c(35, 50, 65)) +
                  as.factor(sex_birth) + as.factor(bipoc_race_derived) +
                  as.factor(ethn_hispanic_latinx_derived) + as.factor(vax_rec_derived) +
                  as.factor(pi_cat) + as.factor(comorb_count_cat) + as.factor(tte_start)")

  
  #Parallel bootstrap loop
  rmst_list <- foreach(i = 1:n_boot, .combine = rbind, 
                       .packages = c("data.table", "dplyr", "survival", "splines", "rms", "magrittr")) %dopar% {
    boot_ids = sample(unique(dt$id), replace = TRUE)
    boot_index = data.table(id = boot_ids, replicate = seq_along(boot_ids))
    
    setkey(dt, id)
    boot_dat = dt[boot_index, allow.cartesian = TRUE]
    
    # boot_dat %<>%
    #   filter(!if_any(all_of(covariates_ipcw), is.na)) ## -- remove any records that are missing BL covariates
    
    #IPCW
    fit_Cd = glm(ipcw_form, data = boot_dat, family="binomial"(link="logit"))
    
    ipcw_denom = predict(fit_Cd, data=boot_dat, type="response")
    
    
    fit_Cn = boot_dat %>%
      group_by(tte_start) %>% # CHANGED FROM DINELKA'S TTE_RC
      summarise(prop_Ct = mean(C_t, na.rm=T), .groups = "drop") %>%
      ungroup()
    
    boot_dat %<>%
      mutate(ipcw_d = ipcw_denom) %>%
      left_join(fit_Cn) %>%
      group_by(replicate) %>% # CHANGED FROM ID
      mutate(ipcw = cumprod(prop_Ct/ipcw_d),
             ipcw = lag(ipcw, default=1)) %>%
      ungroup()
    
    #KM fit
    formula_str = paste("Surv(time = tte_start, time2 = time, event = event_final) ~", group_var) ## -- using only right censored
    
    fit = survfit(as.formula(formula_str), data = boot_dat, weights = ipcw, robust = T, id = replicate) ## CHECK IF THIS SHOULD BE ID OR REPLICATE

    strata_names = names(fit$strata)
    
    rmst_results <- lapply(seq_along(strata_names), function(j) {
      stratum = strata_names[j]
      start_idx = if (j == 1) 1 else sum(fit$strata[1:(j - 1)]) + 1
      end_idx = sum(fit$strata[1:j])
      
      time_points = fit$time[start_idx:end_idx]
      surv_prob = fit$surv[start_idx:end_idx]
      idx <- which(time_points <= tau)
      time_points = c(0, time_points[idx], tau)
      surv_prob = c(1, surv_prob[idx], tail(surv_prob, 1))
      
      rmst = sum((head(surv_prob, -1) + tail(surv_prob, -1)) / 2 * diff(time_points))
      median_tte = summary(fit)$table[stratum, "median"]
      data.frame(replicate=i, stratum = stratum, rmst = rmst, media_tte = median)
    })

    do.call(rbind, rmst_results)
  }
  
  stopCluster(cl)
  return(rmst_list)
}


bootstrap_rmt_iptw_ipcw_parallel <- function(n_boot, tau, n_cores = parallel::detectCores() - 1) {
  
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  registerDoRNG(seed = 123)
  
  dt = 
    data_ipcw_long %>% 
    as.data.table()
  
  # covariates_ipcw <- c("age_derived", "sex_birth", "pi_cat", "urban_rural_derived",
  #                      "bipoc_race_derived", "bmi_cat_derived", "vax_rec_derived", "hybrid_imm",
  #                      "comorb_count_cat", "ethn_hispanic_latinx_derived")
  
  # ipcw_form = as.formula(paste("C_t ~", paste(c("rcs(age_derived, c(35, 50, 65))", 
                                                # covariates_ipcw[-1], "as.factor(tte_start)"), collapse = " + ")))
  
  # iptw_form = as.formula(paste("hybrid_imm ~", paste(c("rcs(age_derived, c(35, 50, 65))", 
                                                       # covariates_ipcw[!covariates_ipcw %in% c('age_derived', 'hybrid_imm', 'pi_cat', 'vax_rec_derived')]), collapse = " + ")))
  
  # Same formulas as return to school/work
  ipcw_form = as.formula("C_t ~ rcs(age_derived, c(35, 50, 65)) +
                  as.factor(sex_birth) + as.factor(bipoc_race_derived) +
                  as.factor(ethn_hispanic_latinx_derived) + as.factor(vax_rec_derived) +
                  as.factor(pi_cat) + as.factor(comorb_count_cat) + as.factor(tte_start)")
  
  iptw_form = as.formula("hybrid_imm ~ rcs(age_derived, c(35, 50, 65)) + as.factor(sex_birth) + as.factor(comorb_count_cat)")
  
  
  rmst_list <- foreach(i = 1:n_boot, .combine = rbind, 
                       .packages = c("data.table", "dplyr", "survival", "splines", "nnet", "rms", "magrittr")) %dopar% {
     #resample IDs and keep all rows for each ID
     boot_ids <- sample(unique(dt$id), replace = TRUE)
     boot_index <- data.table(id = boot_ids, replicate = seq_along(boot_ids))
     setkey(dt, id)
     setkey(boot_index, id)
     boot_dat <- dt[boot_index, allow.cartesian = TRUE]
     
     # boot_dat %<>%
     #   filter(!if_any(all_of(covariates_ipcw), is.na)) ## -- remove any records that are missing BL covariates
     
     #IPTW
     boot_dat_iptw <- boot_dat %>%
       mutate(hybrid_imm = interaction(vax_rec_derived, pi_cat)) %>%
       group_by(id, replicate) %>%
       slice(1) %>%  # Just 1 row per participant, just using the baseline covariates here
       ungroup()
     
     fit_iptw <- multinom(iptw_form,
                          data = boot_dat_iptw)
     
     iptw_pred <- predict(fit_iptw, data=boot_dat_iptw, type = "probs")
     
     iptw_d <- iptw_pred[cbind(1:nrow(boot_dat_iptw),
                               match(boot_dat_iptw$hybrid_imm, levels(boot_dat_iptw$hybrid_imm)))]
     
     iptw_n <- as.numeric(prop.table(table(boot_dat_iptw$hybrid_imm))[
       match(boot_dat_iptw$hybrid_imm, levels(boot_dat_iptw$hybrid_imm))])
     
     boot_dat_iptw$iptw <- iptw_n / iptw_d
     
     #IPCW
     fit_Cd = glm(ipcw_form, data = boot_dat, family="binomial"(link="logit"))
     
     ipcw_denom = predict(fit_Cd, data=boot_dat, type="response")
     
     fit_Cn = boot_dat %>%
       group_by(tte_start) %>% 
       summarise(prop_Ct = mean(C_t, na.rm=T), .groups = "drop") %>%
       ungroup()
     
     boot_dat %<>%
       mutate(ipcw_d = ipcw_denom) %>%
       left_join(fit_Cn) %>%
       group_by(replicate) %>%
       mutate(ipcw = cumprod(prop_Ct/ipcw_d),
              ipcw = lag(ipcw, default=1)) %>%
       ungroup()
     
     boot_dat <- boot_dat %>%
       left_join(boot_dat_iptw %>%
                   dplyr::select(id, replicate, hybrid_imm, iptw), by=c("id", "replicate")) %>%
       mutate(iptw_ipcw = iptw * ipcw)
     
     #KM fit
     formula_str = paste("Surv(time = tte_start, time2 = time, event = event_final) ~ vax_rec_derived + pi_cat")
     
     fit = survfit(as.formula(formula_str), data = boot_dat, weights = iptw_ipcw, robust = T, id = replicate)
     
     strata_names <- names(fit$strata)
     rmst_results <- lapply(seq_along(strata_names), function(j) {
       stratum <- strata_names[j]
       start_idx <- if (j == 1) 1 else sum(fit$strata[1:(j - 1)]) + 1
       end_idx <- sum(fit$strata[1:j])
       
       time_points <- fit$time[start_idx:end_idx]
       surv_prob <- fit$surv[start_idx:end_idx]
       idx <- which(time_points <= tau)
       time_points <- c(0, time_points[idx], tau)
       surv_prob <- c(1, surv_prob[idx], tail(surv_prob, 1))
       
       rmst <- sum((head(surv_prob, -1) + tail(surv_prob, -1)) / 2 * diff(time_points))
       median_tte <- summary(fit)$table[stratum, "median"]
       
       data.frame(replicate = i, stratum = stratum, rmst = rmst)
     })
     
     do.call(rbind, rmst_results)
  }
  
  stopCluster(cl)
  return(rmst_list)
}

bootstrap_summary_results <- function(rmst_list) {
  has_stratum <- "stratum" %in% colnames(rmst_list)
  
  if (has_stratum) {
    #point estimates and CI
    summary_df <- rmst_list %>%
      group_by(stratum) %>%
      summarise(
        rmst_mean = mean(rmst),
        lower95 = quantile(rmst, 0.025),
        upper95 = quantile(rmst, 0.975),
        .groups = "drop"
      ) %>%
      mutate(across(where(is.numeric), ~ formatC(.x, format = "f", digits = 5)))
    
    #difference estimates and CI
    rmst_wide <- rmst_list %>%
      pivot_wider(names_from = stratum, values_from = rmst, id_cols = replicate)
    
    group_names <- setdiff(names(rmst_wide), "replicate")
    ref_name <- group_names[1]
    
    diff_df <- lapply(group_names, function(g) {
      if (g == ref_name) {
        return(data.frame(stratum = g, diff_mean = "Ref.", diff_lower = "Ref.", diff_upper = "Ref."))
      }
      diff_vec <- rmst_wide[[g]] - rmst_wide[[ref_name]]
      data.frame(
        stratum = g,
        diff_mean = formatC(mean(diff_vec), format = "f", digits = 5),
        diff_lower = formatC(quantile(diff_vec, 0.025), format = "f", digits = 5),
        diff_upper = formatC(quantile(diff_vec, 0.975), format = "f", digits = 5)
      )
    }) %>% bind_rows()
    
    out <- summary_df %>%
      left_join(diff_df, by = "stratum")
    
    #check if ref is in first row -- if not, flip the order
    if (out$diff_mean[1] != "Ref.") {
      out <- out %>%
        arrange(desc(stratum))  # flip order if Ref. is not first
    }
    
  } else {
    #no stratum – overall only
    out <- rmst_list %>%
      summarise(
        rmst_mean = mean(rmst),
        lower95 = quantile(rmst, 0.025),
        upper95 = quantile(rmst, 0.975)
      ) %>%
      mutate(across(everything(), ~ formatC(.x, format = "f", digits = 5)))
  }
  
  return(out)
}


### -- bootstrap results

plan(sequential)

# Set bootstrap settings
boot_reps <- 1e3
day <- 14
output_dir <- "results/rmst/boot_int_seq"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Set up the jobs
bootstrap_jobs <- list(
  overall   = ~bootstrap_summary_results(bootstrap_rmt_ipcw_overall(boot_reps, day)),
  age       = ~bootstrap_summary_results(bootstrap_rmt_ipcw_parallel("age_cat", boot_reps, day)),
  sex       = ~bootstrap_summary_results(bootstrap_rmt_ipcw_parallel("sex_birth", boot_reps, day)),
  race      = ~bootstrap_summary_results(bootstrap_rmt_ipcw_parallel("bipoc_race_derived", boot_reps, day)),
  ethn      = ~bootstrap_summary_results(bootstrap_rmt_ipcw_parallel("ethn_hispanic_latinx_derived", boot_reps, day)),
  bmi       = ~bootstrap_summary_results(bootstrap_rmt_ipcw_parallel("bmi_cat_derived", boot_reps, day)),
  comorb    = ~bootstrap_summary_results(bootstrap_rmt_ipcw_parallel("comorb_count_cat", boot_reps, day)),
  vax       = ~bootstrap_summary_results(bootstrap_rmt_ipcw_parallel("vax_rec_derived", boot_reps, day)),
  prior     = ~bootstrap_summary_results(bootstrap_rmt_ipcw_parallel("pi_cat", boot_reps, day)),
  urban     = ~bootstrap_summary_results(bootstrap_rmt_ipcw_parallel("urban_rural_derived", boot_reps, day)),
  hybrid    = ~bootstrap_summary_results(bootstrap_rmt_iptw_ipcw_parallel(boot_reps, day))
)

# Run and save results sequentially
start <- Sys.time()
for (covariate in names(bootstrap_jobs)) {
  cat("\nRunning:", covariate, "...\n")
  
  tic(paste0("Bootstrap: ", covariate))
  job_fun <- rlang::as_function(bootstrap_jobs[[covariate]])
  result <- job_fun()
  toc()
  
  saveRDS(result, file = file.path(output_dir, paste0(covariate, ".rds")))
  
  # Log completion
  cat(paste(Sys.time(), "-", covariate, "done\n"))
}
end <- Sys.time()
print(paste("Total time: ", end - start))

# ## -- save data frame for Melissa to create forest plot
# 

# List all .RDS files in the folder
rds_files <- list.files(output_dir, pattern = "\\.rds$", full.names = TRUE)

# Read and combine all data frames
combined_df <- rds_files %>%
  map(readRDS) %>%
  bind_rows() %>%
  mutate(
    stratum = ifelse(is.na(stratum), "Overall", stratum)
  )

# Save the combined data frame
saveRDS(combined_df, file = file.path(output_dir, "boot_resolution.rds"))




