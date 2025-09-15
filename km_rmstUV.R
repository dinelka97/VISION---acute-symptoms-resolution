### -- DOCUMENTATION
##### -- author: Dinelka Nanayakkara
##### -- description: Produces KM plots and RMST results for each of the recovery & resolution endpoints
  ### -- also, computes the RMST and converts in table format

  ### -- there is a lot of redundant info in this script that are not needed because we are doing plotting elsewhere


rm(list=ls())

# libraries ---------------------------------------------------------------

library(tidyverse)
library(magrittr)
library(ggplot2)
library(glue)
library(gridExtra)
library(survminer)
library(survival)
library(patchwork)
library(gt)
library(simtrial)
library(survRM2)
library(flextable)
library(officer)
library(haven)
library(cowplot)
library(rms)
library(nnet)

f_remove <- c("dfs_km.RDS", "ipcw_dfs.RDS") ## -- remove files before generating new

# build full paths
rem_path <- file.path(glue("{getwd()}/data/km/v3/"), f_remove)

file.remove(rem_path)

# import data -------------------------------------------------------------

  ## -- tte event data
rds_files <- rds_files <- list.files("data/km/v3", pattern = "^df.*\\.RDS$", full.names = TRUE) ## -- only the df_all_...
df_list <- lapply(rds_files, readRDS)
names(df_list) <- tools::file_path_sans_ext(basename(rds_files))

  ## folder paths
fig_path <- glue("figures/km/v3/")
fig_path_ipcw <- glue("figures/km_ipcw/")

  ## -- baseline data to incorporate covariates

base_dat <- readRDS("data/bl_derived.rds")

## -- SUBSET BASE_DAT TO ONLY RECENT INFECTION COHORT

base_dat %<>%
  filter(inc_answers == 1)

## -- artificial administrative censoring cut point (this was finalized in June 2025 by stats team)
admin_cens <- 18

# km_function -------------------------------------------------------------
  ## -- this function can be used for both crude survival, and with covariates
  ## -- this function provides the KM fit

legend_labels = list(
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
  smk_status = c("Current", "Former", "Never"),
  hybrid_imm = c("Unvaccinated, 0", "Unvaccinated, 1+", "Vaccinated ≤ 6 mo., 0", "Vaccinated ≤ 6 mo., 1+",
                 "Vaccinated > 6 mo., 0", "Vaccinated > 6 mo., 1+"),
  NULL = c("All")
)

### -- INPUTS

## -- RMST Cutoff
cutoff <- 14 ## -- cutoff to compute the RMST

# plot_km <- function(df, km_fit, endpoint, symptoms, covariate=NULL) {
#   
#   legend_title = case_when(
#     is.null(covariate) ~ "",
#     covariate == "sex_birth" ~ "Sex at birth",
#     covariate == "age_cat_2_derived" ~ "Age (years)",
#     covariate == "age_cat" ~ "Age (years)",
#     covariate == "pi_cat" ~ "Number of prior Covid-19 infections",
#     covariate == "pi_cat_v2" ~ "Number of prior Covid-19 infections",
#     covariate == "pi_rec" ~ "Prior infection recency",
#     covariate == "urban_rural_derived" ~ "Urbanicity",
#     covariate == "bipoc_race_derived" ~ "BIPOC status",
#     covariate == "race_new_derived" ~ "Race",
#     covariate == "bmi_cat_derived" ~ "Body Mass Index (BMI)",
#     covariate == "vax_rec_derived" ~ "Vaccination recency",
#     covariate == "n_vax_derived" ~ "Vaccine count (doses)",
#     covariate == "comorb_count_cat" ~ "Comorbidity count",
#     covariate == "comorb_count_cat2" ~ "Comorbidity category",
#     covariate == "ethn_hispanic_latinx_derived" ~ "Ethnicity",
#     covariate == "hybrid_imm" ~ "Vaccination recency, prior infection count"
#   )
#   
#   legend_labels_to_use <- if (is.null(covariate)) {
#     legend_labels[["NULL"]]  # "All" for crude survival without covariates
#   } else {
#     #legend_labels[[covariate]][levels(df[[covariate]])]  # Specific labels based on covariate
#     levels(df[[covariate]])
#   }
#   
#   p = ggsurvplot(
#         km_fit,
#         data = df,
#         fun = "event",
#         title = glue("Acute symptom {endpoint} (using {symptoms} symptoms)"),
#         xlab = "Days since symptom onset",
#         ylab = glue("Probability of acute symptom {endpoint}"),
#         size = 1,                
#         conf.int = TRUE, 
#         conf.int.alpha = 0.35, 
#         risk.table = TRUE,       
#         #risk.table.col = NA,
#         risk.table.height = 0.25,
#         break.time.by = 2,
#         #xlim=c(0, max(df$tte_final) + 1),
#         #ggtheme = theme_bw(),
#         ggtheme = theme_classic(),
#         tables.theme = theme_cleantable(), ## remove if needed
#         palette = c("#F8766D", "#00BFC4", "#C77CFF", "#7CAE00", "#00B0F6", "#E38900"), ## -- to stick to a standard color palette
#         legend.title = legend_title,
#         legend.labs = legend_labels_to_use,
#         risk.table.y.text = FALSE,
#         #risk.table.breaks = seq(0, max(df$tte_final) + 1, by = 2),
#         #risk.table.x.text = TRUE,
#         risk.table.fontsize = 2.25
#   )
#   
#   #t = max(df$tte)
#   #survs = (1-summary(km_fit, times = t)$surv)
#   #labels = paste(round(survs*100), '%', sep='')
#   
#   ## -- get number of categories
#   n_labels <- length(legend_labels_to_use)
#   
#   #legend_rows <- if (n_labels > 3) 2 else 1
#   legend_rows <- 1
# 
#   p$plot = p$plot + 
#     coord_cartesian(ylim = c(0, 1), 
#                     #xlim = c(0, t + 1),
#                     xlim = c(0, 18),
#                     clip = 'on') +
#     theme(#plot.title = element_text(hjust = 0.5, size = 14), 
#           axis.title.x = element_blank(),
#           axis.text.x = element_text(size = 9),
#           axis.title.y = element_text(size = 9),
#           plot.title = element_blank(),
#           legend.position = "top",
#           legend.title = element_text(size = 7)) +
#     guides(color = guide_legend(nrow = legend_rows))
#           
#     #               +
#     # scale_y_continuous(name = NULL, 
#     #                    sec.axis = sec_axis(~., name = NULL, breaks = survs, 
#     #                                        labels = labels)) +
#     # theme(axis.text.y.right = element_text(size = 7.5))
#   
#   
#   p$table = p$table +
#     coord_cartesian(#xlim = c(0, t + 1),
#       xlim = c(0, 18),
#       clip = 'on') +
#     theme(axis.title.y = element_blank(),
#           #axis.text.x = element_text(size = 10),
#           axis.text.x = element_blank(),
#           axis.title.x = element_text(size = 12))
#   
#   
#   return(p)
# }

## ---- NOT USING THIS FUNCTION. USING 'SIMTRIAL' PACKAGE INSTEAD.

## -- FUNCTION TO PRODUCE RMST OUTPUT 

# rmt_init <- function(km_fit, tau){
#   
#   rmt_results = list()
#   
#   if(!is.null(km_fit$strata)){
#     for (i in seq_along(km_fit$strata)) {
#       stratum = names(km_fit$strata)[i]
#       start_idx = if (i == 1) {1} else {sum(km_fit$strata[1:(i - 1)]) + 1}
#       end_idx = sum(km_fit$strata[1:i])
#       
#       end_time = min(km_fit$time[end_idx], cutoff)
#       time_index = which(km_fit$time == end_time)[i] ## -- index from km_fit
#       
#       time_points = km_fit$time[start_idx:time_index]
#       
#       surv_prob = km_fit$surv[start_idx:time_index]
#       time_diff = c(diff(time_points), 0)
#       rmt = round(sum(surv_prob * time_diff), 3) #+ time_points[1]
#       rmt_results[[stratum]] = rmt
#     }
#     for (stratum in names(rmt_results)) {
#       print(glue("Stratum: {sub('.*=', '', stratum)} - rmt up to time {cutoff} = {rmt_results[[stratum]]}, \n"))
#     }
#   }else{
#     
#     end_time = min(max(km_fit$time), cutoff)
#     time_index = which(km_fit$time == end_time) ## -- index from km_fit
#     
#     time_points = km_fit$time[1:time_index]
#     
#     surv_prob = km_fit$surv[1:time_index]
#     time_diff = c(diff(time_points), 0)
#     rmt = round(sum(surv_prob * time_diff), 3) #+ time_points[1]
#     rmt_results[["overall"]] = rmt
#     
#     print(glue("Overall RMT up to time {cutoff} = {rmt_results[['overall']]}, \n"))
#     
#   }
#   return(rmt_results)
# }


### -----

# km_crude_list <- vector("list", length = length(df_list))
# names(km_crude_list) <- substr(names(df_list), 4, nchar(names(df_list)))

covariates <- c("age_cat", "sex_birth", "pi_rec", "pi_cat", "pi_cat_v2", "urban_rural_derived",
                "bipoc_race_derived", "bmi_cat_derived", "vax_rec_derived", "comorb_count_cat",
                "n_vax_derived", "ethn_hispanic_latinx_derived", "hybrid_imm")

covariates_list <- list("age_cat" = "Age (years)", "age_cat_2_derived" = "Age (years)", "age_cat_1_derived" = "Age (years)", 
                        "sex_birth" = "Sex at Birth", "pi_rec" = "Prior Infection Recency", "pi_cat" = "Prior Infections", 
                        "pi_cat_v2" = "Prior Infections", "urban_rural_derived" = "Urbanicity", "bipoc_race_derived" = "Race", 
                        "bmi_cat_derived" = "BMI", "vax_rec_derived" = "Vaccination recency", 
                        "comorb_count_cat" = "Comorbidity Count", "n_vax_derived" = "Vaccination count",
                        "ethn_hispanic_latinx_derived" = "Ethnicity", "hybrid_imm" = "Hybrid Immunity")

endpoints <- c("recovery", "resolution")
symps <- c("all")

rms_comp_results <- vector(mode = "list", length = length(covariates) * length(endpoints)) ## -- comprehensive RMST results
names(rms_comp_results) <- paste(rep(endpoints, each = length(covariates)), covariates, sep = "_")

rms_ipcw_comp_results <- vector(mode = "list", length = length(covariates) * length(endpoints)) ## -- comprehensive RMST results
names(rms_ipcw_comp_results) <- paste(rep(endpoints, each = length(covariates)), covariates, sep = "_")


# Stabilized IPCW ---------------------------------------------------------

  ## -- create a df for IPCW (for the two endpoints separately)

# colnames(df_list[["df_all_symps_resolution"]])

dfs_km <- list()
ipcw_summary <- vector(mode = "list", length = 2); names(ipcw_summary) <- names(df_list)

covariates_ipcw <- c("age_derived", "sex_birth", "pi_cat", "urban_rural_derived",
                "bipoc_race_derived", "bmi_cat_derived", "vax_rec_derived",
                "comorb_count_cat", "ethn_hispanic_latinx_derived")

for (i in names(df_list)){
  
  ### -- get if resolution or recovery
  
  ep = sub(".*_", "", i)
  
  
  ### -- first, for each ID censored for those who only dropped out.
  
  data_ipcw = df_list[[i]]
  
  data_ipcw %<>%
    mutate(ltfu = ifelse(event_rc == 0, 1, 0))
    #mutate(ltfu = ifelse(event_rc == 0 & tte_rc < min(days_since_symptom_onset_v2 + 14 - 1, admin_cens), 1, 0)) - done already
  
  ### -- CODE LONG DATASET FIRST!!!
  
  time_pts <- sort(c(0, unique(data_ipcw$tte_rc)))
  
  data_ipcw_long <- 
    data_ipcw %>%
    mutate(tte_rc = ifelse(tte_rc == 0, 0.0001, tte_rc)) %>%       ## -- if tte_rc == 0, change to 0.0001
    rowwise() %>%
    mutate(
      time = list(time_pts[time_pts > min(time_pts) & time_pts <= tte_rc])  # exclude first (usually 0), keep <= tte_rc
    )  %>% ## -- time should be based on the mid-points too, cause time-varying
    unnest(time) %>%
    group_by(id) %>%
    arrange(id, time) %>%
    mutate(
      event_final = if_else(time == max(time), event_rc, 0),
      tte_start = lag(time),
      tte_start = if_else(is.na(tte_start), 0, tte_start)
    ) %>%
    ungroup() %>%
    mutate(
      tte_final = ifelse(tte_rc == 0, 0.0001, tte_rc), ## -- actually not required but it's okay
      ltfu_t = ifelse(ltfu == 1 & (time == tte_final | tte_final == 0.0001), 1, 0),
      C_t = 1 - ltfu_t
    )

  ipcw_form <- as.formula(paste("C_t ~", paste(c("rcs(age_derived, c(35, 50, 65))", 
                                               covariates_ipcw[-1], "as.factor(time)"), collapse = " + ")))
  
  df_fit_ipcw <- 
    data_ipcw_long %>% 
      left_join(base_dat, by = "id") %>%
      filter(!if_any(all_of(covariates_ipcw), is.na))
  
  fit_Cd <- glm(ipcw_form, data = df_fit_ipcw, family="binomial"(link="logit"))
  
  ipcw_denom <- predict(fit_Cd, data=df_fit_ipcw, type="response")
  
  fit_Cn <- df_fit_ipcw %>%
    group_by(tte_start) %>% ## - change to tte_start
    summarise(prop_Ct = mean(C_t, na.rm=T)) %>%
    ungroup()
  
  df_fit_ipcw %<>%
    mutate(ipcw_d = ipcw_denom) %>%
    left_join(fit_Cn) %>%
    group_by(id) %>%
    mutate(ipcw = cumprod(prop_Ct/ipcw_d),
           ipcw = lag(ipcw, default=1)) %>%
    ungroup()
  
  
  
  ipcw_summary[[i]] = df_fit_ipcw$ipcw
  
  #summary(df_fit_ipcw$ipcw)
  
  dfs_km[[glue("ipcw_{ep}")]] = df_fit_ipcw
  
}

print("IPCW dfs saved")

## -- histogram of IPCW weights

# ipcw_hist <- vector(mode = "list", length = 2); names(ipcw_hist) <- names(dfs_km)
# 
# for(i in names(ipcw_summary)){
# 
#   ipcw_hist[[i]] <- 
#     ggplot(data = data.frame(value = ipcw_summary[[i]]), aes(x = value)) +
#     geom_histogram(color = "white", fill = "steelblue", bins = 50) +  # adjust bins as needed
#     scale_x_continuous(breaks = c(round(min(ipcw_summary[[i]]), 1), 1, 2, 3)) +
#     labs(title = glue("IPW weight distribution: {sub('.*_', '', i)}"),
#          x = "Weight",
#          y = "Frequency") +
#     theme_minimal(base_size = 14) 
#     #theme(plot.title = element_text(hjust = 0.5))
#   
#   ggsave(glue("figures/km_ipcw/hist_ipcw_{sub('.*_', '', i)}_v2.png"), plot = ipcw_hist[[i]],
#          width = 10, height = 6, dpi = 350)
#          
# }
# 
# ## -- create a plot list
# plot_list <- vector(mode = "list", 2)
# names(plot_list) <- endpoints
# 
# for (i in seq_along(plot_list)) {
#   plot_list[[i]] = vector(mode = "list", length = length(covariates) + 1)
#   names(plot_list[[i]]) = c("All", covariates)
# }
# 
# plot_list_ipcw <- vector(mode = "list", 2)
# names(plot_list_ipcw) <- endpoints
# 
# for (i in seq_along(plot_list_ipcw)) {
#   plot_list_ipcw[[i]] = vector(mode = "list", length = length(covariates) + 1)
#   names(plot_list_ipcw[[i]]) = c("All", covariates)
# }
# 
# for(g in covariates){
#   for(s in symps){
#     for(e in endpoints){
#       
#       if (!dir.exists(glue("{fig_path}{g}/"))) {
#         dir.create(glue("{fig_path}{g}/"), recursive = TRUE)
#       }
#       
#       if (!dir.exists(glue("{fig_path_ipcw}{g}/"))) {
#         dir.create(glue("{fig_path_ipcw}{g}/"), recursive = TRUE)
#       }
#       
#       d = glue("df_{s}_symps_{e}")
#       
#         ## -- crude (without accounting for censoring weights)
#       data_km = df_list[[d]] %>% ## -- slice(1) removed because already sliced to have n=1 per ID
#                   left_join(base_dat %>% select(id, g))
#       
#       ## -- this will change depending on if we choose to allow for left-censoring, or if we make the mid-point assumption and have
#       # -- only right censoring
#       #formula_str = paste("Surv(days_since_symptom_onset, tte_final + 0.01, event_final) ~", g)
#       #formula_str = paste("Surv(tte_final + 0.01, event_final) ~", g)
#       #formula_str = paste("Surv(time = left, time2 = right, type = 'interval2') ~", g)
#       formula_str = paste("Surv(time = tte_rc, event = event_rc) ~", g) ## -- using only right censored
#       
#       fit = survfit(as.formula(formula_str), data = data_km)
#       
#       ## -- PLOTTING THE KM-PLOT
#       
#       p = plot_km(data_km, fit, e, s, g) # -- with the risk table
#       
#       combined_plot = (p$plot / p$table) + plot_layout(heights = c(3, 1))
#       combined_plot
#       
#       plot_list[[e]][[g]] = combined_plot
#       
#       #ggsave(glue("{fig_path}{g}/{substr(d, 4, nchar(d))}.png"), plot = combined_plot, width = 10, height = 6, dpi = 350)
#       
#       ## -- INCORPORATING RMST
#       
#       # rms_comp_results[[paste(e, g, sep = "_")]] = rmt_init(fit, cutoff = cutoff) 
#       
#       ## -- RMST using the simtrial package by Merck (seems like)
#       levels = levels(data_km %>% pull(g))
#       ref = levels[1]
#       df_temp = data.frame(level = levels,
#                            rmst = vector(mode = "numeric", length = length(levels)),
#                            rmst_lcl = vector(mode = "numeric", length = length(levels)),
#                            rmst_ucl = vector(mode = "numeric", length = length(levels)),
#                            rmst_diff = vector(mode = "numeric", length = length(levels)),
#                            rmst_diff_lcl = vector(mode = "numeric", length = length(levels)),
#                            rmst_diff_ucl = vector(mode = "numeric", length = length(levels))
#       )
#       
#       for(l in levels[-1]){
#         data = data_km %>% filter(!!sym(g) %in% c(ref, l)) %>% droplevels()
#         out = simtrial:::rmst_two_arm(
#           time_var = data$tte_rc,
#           event_var = data$event_rc,
#           group_var = data %>% pull(g),
#           trunc_time = cutoff,
#           reference = ref
#         )
#         
#         ## -- ref group estimates
#         df_temp[1, "rmst"] = ifelse(which(l == levels) == 2, out$rmst_per_arm[out$rmst_per_arm$group == ref ,"rmst"][1], 
#                                     df_temp[1, "rmst"])
#         df_temp[1, "rmst_lcl"] = ifelse(which(l == levels) == 2, out$rmst_per_arm[out$rmst_per_arm$group == ref,"lcl"][1], 
#                                         df_temp[1, "rmst_lcl"])
#         df_temp[1, "rmst_ucl"] = ifelse(which(l == levels) == 2, out$rmst_per_arm[out$rmst_per_arm$group == ref,"ucl"][1], 
#                                         df_temp[1, "rmst_ucl"])
#         
#         ## -- other group estimates
#         df_temp[which(l == levels), "rmst"] = out$rmst_per_arm[out$rmst_per_arm$group == l,"rmst"]
#         df_temp[which(l == levels), "rmst_lcl"] = out$rmst_per_arm[out$rmst_per_arm$group == l,"lcl"]
#         df_temp[which(l == levels), "rmst_ucl"] = out$rmst_per_arm[out$rmst_per_arm$group == l,"ucl"]
#         
#         df_temp[which(l == levels), "rmst_diff"] = out$rmst_diff[,"rmst_diff"]
#         df_temp[which(l == levels), "rmst_diff_lcl"] = out$rmst_diff[,"lcl"]
#         df_temp[which(l == levels), "rmst_diff_ucl"] = out$rmst_diff[,"ucl"]
#       }
#       
#       rms_comp_results[[paste(e, g, sep = "_")]] = df_temp
#       
#       
#         ## -- ACCOUNTING FOR CENSORING WEIGHTS (IPCW)
#       
#       data_km_ipcw = dfs_km[[d]] %>%
#                     #slice(1) %>%
#                     left_join(base_dat %>% 
#                                 filter(!if_any(all_of(covariates_ipcw), is.na)) %>% 
#                                 select(id, g)) 
#       
#       ## -- this will change depending on if we choose to allow for left-censoring, or if we make the mid-point assumption and have
#       # -- only right censoring
#       #formula_str = paste("Surv(days_since_symptom_onset, tte_final + 0.01, event_final) ~", g)
#       #formula_str = paste("Surv(tte_final + 0.01, event_final) ~", g)
#       #formula_str = paste("Surv(time = left, time2 = right, type = 'interval2') ~", g)
#       formula_str = paste("Surv(time = tte_start, time2 = time, event = event_final) ~", g) ## -- using only right censored
#       
#       fit = survfit(as.formula(formula_str), data = data_km_ipcw, weights = ipcw, robust = T, id = id)
#       
#       ## -- PLOTTING THE KM-PLOT
#       
#       p = plot_km(data_km_ipcw, fit, e, s, g)
#       combined_plot = (p$plot / p$table) + plot_layout(heights = c(3, 1))
#       combined_plot
#       
#       plot_list_ipcw[[e]][[g]] = combined_plot
#       
#       ggsave(glue("{fig_path_ipcw}{g}/{substr(d, 4, nchar(d))}.png"), plot = combined_plot, width = 10, height = 6, dpi = 350)
#       
#       ## -- INCORPORATING RMST
#       
#       # rms_comp_results[[paste(e, g, sep = "_")]] = rmt_init(fit, cutoff = cutoff) 
#       
#       ## -- RMST using the simtrial package by Merck (seems like)
#       # levels = levels(data_km_ipcw %>% pull(g))
#       # ref = levels[1]
#       # df_temp = data.frame(level = levels,
#       #                      rmst = vector(mode = "numeric", length = length(levels)),
#       #                      rmst_lcl = vector(mode = "numeric", length = length(levels)),
#       #                      rmst_ucl = vector(mode = "numeric", length = length(levels)),
#       #                      rmst_diff = vector(mode = "numeric", length = length(levels)),
#       #                      rmst_diff_lcl = vector(mode = "numeric", length = length(levels)),
#       #                      rmst_diff_ucl = vector(mode = "numeric", length = length(levels))
#       # )
#       # 
#       # for(l in levels[-1]){
#       #   data = data_km_ipcw %>% filter(!!sym(g) %in% c(ref, l)) %>% droplevels()
#       #   out = simtrial:::rmst_two_arm(
#       #     time_var = data$tte_rc,
#       #     event_var = data$event_rc,
#       #     group_var = data %>% pull(g),
#       #     trunc_time = cutoff,
#       #     reference = ref
#       #   )
#       #   
#       #   ## -- ref group estimates
#       #   df_temp[1, "rmst"] = ifelse(which(l == levels) == 2, out$rmst_per_arm[out$rmst_per_arm$group == ref ,"rmst"][1], 
#       #                               df_temp[1, "rmst"])
#       #   df_temp[1, "rmst_lcl"] = ifelse(which(l == levels) == 2, out$rmst_per_arm[out$rmst_per_arm$group == ref,"lcl"][1], 
#       #                                   df_temp[1, "rmst_lcl"])
#       #   df_temp[1, "rmst_ucl"] = ifelse(which(l == levels) == 2, out$rmst_per_arm[out$rmst_per_arm$group == ref,"ucl"][1], 
#       #                                   df_temp[1, "rmst_ucl"])
#       #   
#       #   ## -- other group estimates
#       #   df_temp[which(l == levels), "rmst"] = out$rmst_per_arm[out$rmst_per_arm$group == l,"rmst"]
#       #   df_temp[which(l == levels), "rmst_lcl"] = out$rmst_per_arm[out$rmst_per_arm$group == l,"lcl"]
#       #   df_temp[which(l == levels), "rmst_ucl"] = out$rmst_per_arm[out$rmst_per_arm$group == l,"ucl"]
#       #   
#       #   df_temp[which(l == levels), "rmst_diff"] = out$rmst_diff[,"rmst_diff"]
#       #   df_temp[which(l == levels), "rmst_diff_lcl"] = out$rmst_diff[,"lcl"]
#       #   df_temp[which(l == levels), "rmst_diff_ucl"] = out$rmst_diff[,"ucl"]
#       # }
#       # 
#       # rms_ipcw_comp_results[[paste(e, g, sep = "_")]] = df_temp
#     
#     }
#   }
# }
# 
# 
#   ## -- crude survival
# 
# for(s in symps){
#   for(e in endpoints){
#     
#     if (!dir.exists(glue("{fig_path}crude/"))) {
#       dir.create(glue("{fig_path}crude/"), recursive = TRUE)
#     }
#     
#     if (!dir.exists(glue("{fig_path_ipcw}crude/"))) {
#       dir.create(glue("{fig_path_ipcw}crude/"), recursive = TRUE)
#     }
#     
#     d = glue("df_{s}_symps_{e}")
#     data_km = df_list[[d]] %>% group_by(id) %>% slice(1) 
#     
#     #data_km %<>% filter(id %in% data_iptw$id) ## -- comment out later
# 
#     #fit = survfit(as.formula(Surv(days_since_symptom_onset, tte_final + 0.01, event_final) ~ 1), data = data_km)
#     #fit = survfit(as.formula(Surv(time = left, time2 = right, type = "interval2") ~ 1), data = data_km)
#     fit = survfit(as.formula(Surv(time = tte_rc, event = event_rc) ~ 1), data = data_km)
#     
#     
#     p = plot_km(data_km, fit, e, s)
#     combined_plot = (p$plot / p$table) + plot_layout(heights = c(3, 1))
#     
#     plot_list[[e]][["All"]] = combined_plot
#     
#     #ggsave(glue("{fig_path}crude/{substr(d, 4, nchar(d))}.png"), plot = combined_plot, width = 10, height = 6, dpi = 350)
#     
#     ## -- INCORPORATING RMST
#     
#     #rms_comp_results[[paste(e, "overall", sep = "_")]] = rmt_init(fit, cutoff = cutoff) 
#     
#     out = simtrial:::rmst_single_arm(
#       time_var = data_km$tte_rc,
#       event_var = data_km$event_rc,
#       tau = cutoff,
#       group_label = "Crude"
#     )
#     
#     df_temp = data.frame(level = "Overall",
#                          rmst = out$rmst,
#                          rmst_lcl = out$lcl,
#                          rmst_ucl = out$ucl,
#                          rmst_diff = NA, rmst_diff_lcl = NA, rmst_diff_ucl = NA
#                          )
#     
#     rms_comp_results[[paste(e, "overall", sep = "_")]] = df_temp
#     
#     ## -- ACCOUNTING FOR CENSORING WEIGHTS (IPCW)
#     
#     data_km_ipcw =
#       dfs_km[[d]]
#     
#     ## -- this will change depending on if we choose to allow for left-censoring, or if we make the mid-point assumption and have
#     # -- only right censoring
#     #formula_str = paste("Surv(days_since_symptom_onset, tte_final + 0.01, event_final) ~", g)
#     #formula_str = paste("Surv(tte_final + 0.01, event_final) ~", g)
#     #formula_str = paste("Surv(time = left, time2 = right, type = 'interval2') ~", g)
#     
#     formula_str = paste("Surv(time = tte_start, time2 = time, event = event_final) ~ 1") ## -- using only right censored
#     
#     fit = survfit(as.formula(formula_str), data = data_km_ipcw, weights = ipcw, robust = T, id = id)
#     
#     ## -- PLOTTING THE KM-PLOT
#     
#     combined_plot = (p$plot / p$table) + plot_layout(heights = c(3, 1))
#     combined_plot
#     
#     plot_list_ipcw[[e]][["All"]] = combined_plot
# 
#     ggsave(glue("{fig_path_ipcw}crude/{substr(d, 4, nchar(d))}.png"), plot = combined_plot, width = 10, height = 6, dpi = 350)
#     
#     ## -- INCORPORATING RMST
#     
#     #rms_comp_results[[paste(e, "overall", sep = "_")]] = rmt_init(fit, cutoff = cutoff) 
#     
#     out = simtrial:::rmst_single_arm(
#       time_var = data_km_ipcw$tte_rc,
#       event_var = data_km_ipcw$event_rc,
#       tau = cutoff,
#       group_label = "Crude"
#     )
#     
#     df_temp = data.frame(level = "Overall",
#                          rmst = out$rmst,
#                          rmst_lcl = out$lcl,
#                          rmst_ucl = out$ucl,
#                          rmst_diff = NA, rmst_diff_lcl = NA, rmst_diff_ucl = NA
#     )
#     
#     rms_ipcw_comp_results[[paste(e, "overall", sep = "_")]] = df_temp
#   
#   }
# }


# IPCW + IPTW -------------------------------------------------------------

## -- Hybrid Immunity

## IPCW only----
  
# d <- "ipcw_resolution"
# g <- c("vax_rec_derived", "pi_cat")
# 
# data_km_ipcw = 
#   dfs_km[[d]]
  # #slice(1) %>%
  # left_join(base_dat %>% 
  #             filter(!if_any(all_of(covariates_ipcw), is.na)) %>% 
  #             select(id, g)) %>%
    # group_by(id) 
    #mutate(tte_start = 0:(max(tte_rc)-1))

# ipcw_km_hybrid <- survfit(Surv(time = tte_start, time2 = time, event_final) ~ vax_rec_derived + pi_cat, 
#                         data = data_km_ipcw, weights=ipcw, robust=T, id=id)
# 
# # names(ipcw_km_hybrid$strata) <- c("Unvax, no pi", "Unvax, 1+ pi", 
# #                                 "Vax <=6 mo., no pi", "Vax <=6 mo., 1+ pi", 
# #                                 "Vax >6 mo., no pi", "Vax >6 mo., 1+ pi")
# 
# p_ipcw_hybrid <- plot_km(data_km_ipcw, ipcw_km_hybrid, "resolution", "all", covariate = "hybrid_imm")
# 
# p_ipcw_hybrid <- p_ipcw_hybrid[[1]] +
#   labs(x = "Days since symptom onset") +
#   theme(
#     legend.title = element_blank(),
#     plot.title = element_blank(),
#     plot.margin = margin(t = 20, r = 5, b = 5, l = 5),
#     axis.title.x = element_text() ## -- adding this for the x-axis to be displayed. Unsure why it wont ow.
#   )
# 
# p_ipcw_hybrid %>%
#    ggsave(filename = "figures/fig3/resolution_hybrid_imm_ipcw.png", dpi = 300, width = 10, height = 6)


## IPTW+IPCW----

## -- covariates of interest

covariates_iptw <- c("age_derived", "sex_birth", "comorb_count_cat", "hybrid_imm")

## -- prep data

data_iptw = df_list[["df_all_symps_resolution"]]

data_iptw %<>%
  mutate(ltfu = ifelse(event_rc == 0, 1, 0),
         tte_rc = ifelse(tte_rc == 0, 0.0001, tte_rc))  ## -- if tte_rc == 0, change to 0.0001
#mutate(ltfu = ifelse(event_rc == 0 & tte_rc < min(days_since_symptom_onset_v2 + 14 - 1, admin_cens), 1, 0)) - done already

### -- CODE LONG DATASET FIRST!!!

time_pts <- sort(c(0, unique(data_iptw$tte_rc)))

data_iptw_long <- 
  data_iptw %>%
  rowwise() %>%
  mutate(
    time = list(time_pts[time_pts > min(time_pts) & time_pts <= tte_rc])  # exclude first (usually 0), keep <= tte_rc
  )  %>% ## -- time should be based on the mid-points too, cause time-varying
  unnest(time) %>%
  group_by(id) %>%
  arrange(id, time) %>%
  mutate(
    event_final = if_else(time == max(time), event_rc, 0),
    tte_start = lag(time),
    tte_start = if_else(is.na(tte_start), 0, tte_start)
  ) %>%
  ungroup() %>%
  mutate(
    tte_final = ifelse(tte_rc == 0, 0.0001, tte_rc),
    ltfu_t = ifelse(ltfu == 1 & (time == tte_final | tte_final == 0.0001), 1, 0),
    C_t = 1 - ltfu_t
  )

df_fit_iptw <- 
  data_iptw_long %>% 
  left_join(base_dat, by = "id") %>%
  filter(!if_any(all_of(covariates_iptw), is.na)) %>%
  group_by(id) %>%
  slice(1)
  
## Multinomial stabilized IPTW

fit_iptw <- multinom(hybrid_imm ~ rcs(age_derived, c(35, 50, 65)) + as.factor(sex_birth) +
                       as.factor(comorb_count_cat),
                   data = df_fit_iptw) ## -- these are non-time varying BL characteristics

iptw_pred <- predict(fit_iptw, type = "probs")

iptw_d <- iptw_pred[cbind(1:nrow(df_fit_iptw),
                        match(df_fit_iptw$hybrid_imm,
                              levels(df_fit_iptw$hybrid_imm)))]

iptw_n <- as.numeric(prop.table(table(df_fit_iptw$hybrid_imm))
                   [match(df_fit_iptw$hybrid_imm,
                          levels(df_fit_iptw$hybrid_imm))])

iptw_dat <- data.frame(id = df_fit_iptw$id, iptw = iptw_n/iptw_d)

# summary(1/iptw_d)
# summary(iptw_dat$iptw)

dfs_km[["iptw_resolution"]] <-
  dfs_km[["ipcw_resolution"]] %>%
    left_join(iptw_dat) %>%
    mutate(iptw_ipcw = iptw*ipcw)

print("IPTW dfs saved")


# summary(data_km_ipcw$iptw_ipcw)
# 
# iptw_ipcw_km_hybrid <- survfit(Surv(time = tte_start, time2 = time, event_final) ~ 
#                                  vax_rec_derived + pi_cat,
#                              data = data_km_ipcw, weights = iptw_ipcw, robust=T, id=id)
# 
# # names(iptw_ipcw_km_hybrid$strata)<-c("Unvax, no pi", "Unvax, 1+ pi",
# #                                      "Vax <=6 mo., no pi", "Vax <=6 mo., 1+ pi",
# #                                      "Vax >6 mo., no pi", "Vax >6 mo., 1+ pi")
# 
# p_iptw_ipcw_hybrid <- plot_km(data_km_ipcw, iptw_ipcw_km_hybrid, "resolution", "all", covariate = "hybrid_imm")
# 
# plot_wrt = (p_iptw_ipcw_hybrid$plot / p_iptw_ipcw_hybrid$table) + plot_layout(heights = c(3, 1))
# 
# plot_wrt %>%
#   ggsave(filename = "figures/fig3/resolution_hybrid_imm_iptw_ipcw_wrt.png",
#        dpi = 300, width = 10, height = 6)
# 
# p_iptw_ipcw_hybrid <- p_iptw_ipcw_hybrid[[1]] +
#   labs(x = "Days since symptom onset") +
#   theme(
#     legend.title = element_blank(),
#     plot.title = element_blank(),
#     plot.margin = margin(t = 20, r = 5, b = 5, l = 5),
#     axis.title.x = element_text() ## -- adding this for the x-axis to be displayed. Unsure why it wont ow.
#   )
# 
# p_iptw_ipcw_hybrid %>%
#   ggsave(filename = "figures/fig3/resolution_hybrid_imm_iptw_ipcw_wort.png", dpi = 300, width = 10, height = 6)

  ## -- wrt - with risk table, wort - without risk table




##IPTW + IPCW for drop out & treatment----

## -- for now, let's do it only for resolution (but later adapt to recovery as well, if required)

df_do_trt <- df_list[["df_all_symps_resolution"]]  
                ## -- df for dropout and treatment

### Follow-up time----

# check if everyone in the symptom resolution datasets is included in trt_dat data
# setdiff(df_do_trt$id, trt_dat$id)

# symp_res <- function(data.input) { 
#   
#   ## -- this relates to treatment (need trt_ever from before, and create trt_start here)
#   
#   df = 
#     data.input %>%
#     # mutate(
#     #   ind_status = dplyr::recode(school_work_fm_status, "Yes" = "yes", "No" = "no", 
#     #                              .default = "missing") ## -- ?what is this, is this required for me?
#     # ) %>%
#     left_join(trt_dat) %>%
#     mutate(trt_start = ifelse(trt_ever == 1, trt_start, NA))
#   
#   # df = df %>%
#   #   group_by(id) %>%
#   #   arrange(days_since_symptom_onset) %>%
#   #   mutate(
#   #     # Define 2-day consecutive "yes" window
#   #     first_day = first(days_since_symptom_onset),
#   #     first2_yes = (ind_status == "yes") & (lag(ind_status == "yes", default = FALSE)),
#   #     is_first_two_yes = row_number() == 2 & first2_yes & first_day > 1,
#   #     ind_rec2 = case_when(
#   #       is_first_two_yes ~ "interval", 
#   #       first2_yes ~ "yes",
#   #       !first2_yes ~ "no",
#   #       is.na(ind_status) ~ "missing"
#   #     ),
#   #     # censor the day before they missed 3+ days, or the day of their last survey
#   #     ind_cens = slide_lgl(ind_status, ~ all(.x == "missing"), .before = 0, .after = 2, .complete = TRUE),
#   #     t_cens1 = as.numeric(days_since_symptom_onset[which(ind_cens == TRUE)[1]]) - 1,
#   #     t_cens2 = case_when(
#   #       is.na(t_cens1) ~ as.numeric(days_since_symptom_onset)[max(which(ind_status != "missing"), na.rm = TRUE)],
#   #       TRUE ~ NA),
#   #     # admin censoring at day 18
#   #     t_admin = 18,
#   #     # Set event time to the first day in the 2-day window
#   #     event = any(ind_rec2 == "yes"),
#   #     tte_standard = as.numeric(days_since_symptom_onset[which(ind_rec2 == "yes")[1]]) - 1,
#   #     #interval censoring [day 1 since sx onset, study entry]/2
#   #     tte_interval = ifelse(any(ind_rec2 == "interval") & event==TRUE, (1 + first_day) / 2, NA),
#   #     #time to event
#   #     tte = ifelse(!is.na(tte_interval), tte_interval, tte_standard),
#   #     #if the patient took treatment, censor them the day before treatment initiation
#   #     tte_final = min(tte, t_cens1, t_cens2, t_admin, trt_start-1, na.rm = TRUE),
#   #     event_final = case_when(
#   #       event == TRUE & is.na(t_cens1) & tte <= t_admin & (tte<trt_start | is.na(trt_start)) ~ 1,
#   #       event == TRUE & !is.na(t_cens1) & t_cens1 >= tte & tte <= t_admin & (tte<trt_start | is.na(trt_start)) ~ 1,
#   #       TRUE ~ 0
#   #     )
#   #   ) %>%
#   #   ungroup()
#   
#   df_summary = 
#     df 
#     # %>% group_by(id) %>%
#     # slice(1) %>%
#     # dplyr::select(id, event_final, tte_final, t_cens1, t_cens2, tte_interval, tte_standard, tte,
#     #               trt_ever, trt_start, t_admin, first_day)
#   
#   return(df_summary)
# }

df_do_trt %<>%
  filter(tte_rc_trt >= 0) ## -- using tte_rc_trt instead of tte_final

#-- ltfu=censored, t_obs=event time
# school_work_km_wide_trt <- 
#   school_work_time_trt %>%
#     mutate(ltfu = ifelse(event_final == 0 & tte_final < (first_day + 13), 1, 0),
#            t_obs=tte_final)

df_do_trt %<>%
  mutate(ltfu = ifelse(event_rc_trt == 0, 1, 0),
         t_obs = tte_rc_trt)

table(df_do_trt$tte_rc_trt, df_do_trt$ltfu, dnn = c("follow-up time", "ltfu"))

#for event day=0, switch event time by +1 to use the survSplit function

df_do_trt %<>% 
  mutate(tte_final = ifelse(tte_rc_trt == 0, 0.0001, tte_rc_trt)) ## -- think it's better to do this here

df_do_trt_long <- survSplit(Surv(tte_final, event_rc_trt)~., 
                            data = df_do_trt, ## --  
                            cut = c(0, 0.0001, 1, 1.5, 2, 2.5, 3, 3.5, 4:18)) ## -- an easier way to do (keep in mind)

df_do_trt_long %<>%
  group_by(id) %>%
  mutate(n_rows_per_id = n()) %>%
  ungroup() %>%
  mutate(ltfu_t = ifelse(ltfu==1 & (t_obs == tte_final | 
                                      (tte_final == 0.0001 & 
                                         n_rows_per_id==1)), 1, 0), #for ltfu, set ltfu_t=1 on the last follow-up day
         ltfu_t = ifelse(ltfu_t == 1 & 
                           trt_ever == 1, 2, ltfu_t), #keep ltfu_t=1 for drop out, change ltfu_t to 2 when initiating treatment
         C_t = 2-ltfu_t) %>% # C_t=2 if uncensored, C_t=1 if drop out, C_t=0 if censored due to treatment initiation
  dplyr::select(-n_rows_per_id)

###IPCW for treatment----
# base_dat<-base_dat %>%
#   mutate(age_analysis=case_when(age_derived %in% c(18:34) ~ "18-34",
#                                 age_derived %in% c(35:49) ~ "35-49",
#                                 age_derived %in% c(50:64) ~ "50-64",
#                                 age_derived>64 ~ "65+"),
#          prior_inf=case_when(n_prior_infec==0~"0",
#                              n_prior_infec>=1~"1+"),
#          comorb_cat=case_when(rowSums(select(., comorb_var)=="Yes")==0 ~"0 comorbs.",
#                               rowSums(select(., comorb_var)=="Yes")==1 ~ "1 comorbs.",
#                               rowSums(select(., comorb_var)=="Yes")>1 ~ "2+ comorbs."),
#          sex_birth=case_when(sex_birth==1 ~ "Female",
#                              sex_birth==2 ~ "Male"),
#          bipoc_race_derived=case_when(bipoc_race_derived==0 ~ "White/Caucasian",
#                                       bipoc_race_derived==1 ~ "BIPOC"),
#          ethn_hispanic_latinx_derived=case_when(ethn_hispanic_latinx_derived==0 ~ "Not Hispanic/Latinx",
#                                                 ethn_hispanic_latinx_derived==1 ~ "Hispanic/Latinx")) %>%
#   mutate(bipoc_race_derived=factor(bipoc_race_derived, levels=c("White/Caucasian", "BIPOC")),
#          ethn_hispanic_latinx_derived=factor(ethn_hispanic_latinx_derived, levels=c("Not Hispanic/Latinx", "Hispanic/Latinx")))

## -- I THINK MOVE DOWN FROM HERE ON TO OTHER SCRIPT


df_do_trt_long %<>%
  left_join(base_dat)

#change C_t to factor for modeling
df_do_trt_long$C_t <- factor(df_do_trt_long$C_t, levels = c("0", "1", "2")) ## -- change levels to 0,1,2

df_do_trt_long %<>%
  filter(!if_any(all_of(covariates_ipcw), is.na)) ## -- it's the same covariates here as well

fit_Cd_trt <- multinom(C_t ~ rcs(age_derived, c(35, 50, 65))+as.factor(sex_birth)+as.factor(bipoc_race_derived)+
                         as.factor(ethn_hispanic_latinx_derived)+as.factor(vax_rec_derived)+as.factor(pi_cat)+
                         as.factor(comorb_count_cat)+as.factor(tstart), data = df_do_trt_long)

ipcw_pred_trt <- predict(fit_Cd_trt, type = "probs")

#IPCW denominator
# df_do_trt_long %<>%
#   mutate(ipcw_d_trt = ipcw_pred_trt[cbind(1:nrow(df_do_trt_long), 
#                                         match(df_do_trt_long$C_t, levels(df_do_trt_long$C_t)))])

df_do_trt_long %<>%
  mutate(ipcw_d_trt = ipcw_pred_trt[cbind(1:nrow(df_do_trt_long),
                                          match(df_do_trt_long$C_t, levels(df_do_trt_long$C_t)))])


#IPCW numerator
ipcw_marg_probs_trt <- df_do_trt_long %>%
  group_by(tstart, C_t) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(tstart) %>%
  mutate(p=n/sum(n)) %>%
  ungroup()

#IPCW
df_do_trt_long %<>%
  left_join(ipcw_marg_probs_trt, by = c("tstart", "C_t")) %>%
  rename(ipcw_n_trt=p) %>%
  group_by(id) %>%
  mutate(ipcw_trt=cumprod(ipcw_n_trt/ipcw_d_trt),
         ipcw_trt=lag(ipcw_trt, default=1)) %>%
  ungroup()

summary(df_do_trt_long$ipcw_trt)

#truncate extreme weights
low_ipcw <- quantile(df_do_trt_long$ipcw_trt, 0.01)
high_ipcw <- quantile(df_do_trt_long$ipcw_trt, 0.99)

df_do_trt_long %<>%
  mutate(ipcw_trt_trunc = pmin(pmax(ipcw_trt, low_ipcw), high_ipcw)) %>%
  left_join(iptw_dat) %>% ## -- MOVE THIS ALL BELOW HYBRID IMM IPTW + IPCW (w/o trt to get iptw_dat)
  mutate(iptw_ipcw_trt = iptw * ipcw_trt_trunc)

summary(df_do_trt_long$ipcw_trt_trunc)
summary(df_do_trt_long$iptw_ipcw_trt)


  ## -- save df to the list dfs_km

dfs_km[["iptw_trt_resolution"]] = df_do_trt_long

print("IPTW treatment dfs saved")

## -- hybrid plot
# 
# iptw_ipcw_trt_km_hybrid <- survfit(Surv(tstart, tte_final, event_rc_trt) ~ vax_rec_derived + pi_cat, 
#                                  data = df_do_trt_long, weights = iptw_ipcw_trt, robust = T, id=id)
# 
# names(iptw_ipcw_trt_km_hybrid$strata) <- c("Unvax, no pi", "Unvax, 1+ pi", 
#                                          "Vax <=6 mo., no pi", "Vax <=6 mo., 1+ pi", 
#                                          "Vax >6 mo., no pi", "Vax >6 mo., 1+ pi")
# 
# p_iptw_ipcw_trt_hybrid <- plot_km(df_do_trt_long, km_fit = iptw_ipcw_trt_km_hybrid, endpoint = "resolution", 
#                                   symptoms = "all", covariate = "hybrid_imm"
#                                   )
# 
# p_iptw_ipcw_trt_hybrid


## -- SAVE LIST TO LOCAL, to be used in any other scripts (this list will have dfs that include C_t, tte's, and event)

dfs_km %>%
  saveRDS("data/km/v3/dfs_km.RDS") ## -- holds dfs for IPCW (nothing to do with treatment nor IPTW for hybrid imm)

## -- although this list holds for both recovery and resolution the rest (going forward) only include resolution









# Results compilation -----------------------------------------------------

### -- (1) Save all KM plots to km_ipcw_v2 folder
# 
# fig2_all_plots <- plot_list_ipcw[["resolution"]]
# 
# 
# # Create folder if it doesn't exist
# folder_path <- "figures/km_ipcw_v3/"
# 
# if (!dir.exists(folder_path)) {
#   dir.create(folder_path)
# } else {
#   # Delete all files inside it
#   files_to_delete <- list.files(folder_path, full.names = TRUE)
#   file.remove(files_to_delete)
# }
# 
# for (i in seq_along(fig2_all_plots)) {
#   #p = fig2_all_plots[[i]]$plot  # extract the ggplot part
#   ggsave(filename = paste0(folder_path, names(fig2_all_plots)[i], ".png"),
#          plot = fig2_all_plots[[i]],
#          width = 6, height = 4, dpi = 300)
# }
# 
# 
# ### -- (2) KM plots to one combined figure
# 
#   ## -- Figure 2 (b): IPCW KM plots for symptom resolution
# 
# covariates_plot <- c("All", "age_cat", "sex_birth", "bipoc_race_derived", "ethn_hispanic_latinx_derived",
#                      "bmi_cat_derived", "comorb_count_cat", "urban_rural_derived")
# 
# covariates_plot_labels <- c("Overall", "Age(yrs)", "Sex at birth", "Race", "Ethnicity",
#                             "BMI", "Comorbidity count", "Urbanicity")
# 
# plot_list_fig2 <- plot_list_ipcw[["resolution"]][covariates_plot]
# 
# # plot_list_fig2 <- lapply(plot_list_fig2, function(p) {
# # 
# #   # Modify the plot
# #   p = p[[1]] +
# #     labs(x = "Days since symptom onset") +
# #     theme(
# #       legend.title = element_blank(),
# #       plot.title = element_blank(),
# #       plot.margin = margin(t = 20, r = 5, b = 5, l = 5),
# #       axis.title.x = element_text(margin = margin(t = 10))  # Add space above x-axis label
# #     )
# #   
# #   return(p)
# # })
# 
# plot_list_fig2_wrt <- lapply(plot_list_fig2, function(p) {
#   # Extract components
#   surv_plot <- p[[1]] +
#     #labs(x = "Days since symptom onset") +
#     theme(
#       legend.title = element_blank(),
#       plot.title = element_blank(),
#       plot.margin = margin(t = 20, r = 5, b = 5, l = 5),
#       axis.title.x = element_text(margin = margin(t = 10))
#     )
#   
#   risk_table <- p[[2]]
#   
#   # Combine using cowplot with reduced height for table
#   plot_grid(surv_plot, risk_table, ncol = 1, rel_heights = c(0.75, 0.25))  # Adjust rel_heights as needed
# })
# 
# fig2 <- plot_grid(plotlist = plot_list_fig2,
#                            ncol = 4,
#                            labels = glue("({LETTERS[1:length(plot_list_fig2)]}) {covariates_plot_labels}"),
#                            label_size = 12,
#                            label_fontface = "bold",
#                            align = "v")

#fig2

# fig2 %>%
#   ggsave(filename = "figures/fig2/resolution_ipcw_wrt.png", dpi = 300, width = 13, height = 15) ## -- wrt - risk table

# ### -- (3) RMST LIST TO PRESENTABLE TABLE
# 
#   ### -- create function to create table
# 
# create_rmst_table <- function(rms_list, endpoint, ind_ipcw){
#   
#   file_ext = ifelse(ind_ipcw == 1, "ipcw", "crude")
#   
#   if (!dir.exists(glue("results/rmst/html/{file_ext}/"))) {
#     dir.create(glue("results/rmst/html/{file_ext}/"), recursive = TRUE)
#   }
#   
#   if (!dir.exists(glue("results/rmst/word/{file_ext}/"))) {
#     dir.create(glue("results/rmst/word/{file_ext}/"), recursive = TRUE)
#   }
#   
#   out_path_html = glue("results/rmst/html/{file_ext}/")
#   out_path_word = glue("results/rmst/word/{file_ext}/")
# 
#   list_res = rms_list[grepl(paste0("^", endpoint), names(rms_list))]
#   
#   endpoint_label = ifelse(endpoint == "recovery", "recovered/improved", "resolved")
#   
#   char_labels = c(unname(unlist(covariates_list[covariates])), "Overall") ## -- capture characteristic labels from covariates & list
#   
#   # df_rmst <- map_dfr(seq_along(list_res), function(i) {
#   #   item <- list_res[[i]]
#   #   tibble(Row = i, Key = names(item), Value = unlist(item))
#   # })
#   
#   df_rmst = imap_dfr(list_res, ~ mutate(.x, variable = .y))
#   
#   # df_rmst %<>%
#   #   rename(category = Key, RMST = Value) %>%
#   #   group_by(Row) %>%
#   #   mutate(Characteristic = char_labels[Row],
#   #          RMST_diff = case_when(covariate == first(covariate) & Characteristic != "Overall" ~ "Ref",
#   #                                Characteristic != "Overall" ~ sprintf("%.2f", RMST - first(RMST)), 
#   #                                Characteristic == "Overall" ~ "-"),
#   #          RMST = round(RMST, 2)
#   #          ) %>%
#   #   ungroup() %>%
#   #   select(-Row) %>%
#   #   mutate(Stratum = ifelse(Characteristic != "Overall", gsub(".*=", "", covariate), "-")) %>%
#   #   select(-covariate) %>%
#   #   select(Characteristic, Stratum, RMST, RMST_diff)
#   # 
#   
# 
#   df_rmst <-
#     df_rmst %>%
#     relocate(variable) %>%
#     rename(Characteristic = variable)
#   
#   df_rmst$Characteristic <- char_labels[match(df_rmst$Characteristic, unique(df_rmst$Characteristic))]
#   
#   df_rmst %<>%
#     mutate(across(where(is.numeric), ~ round(.x, 2)))
#   
#   
#   ## -- reorder rows so that the Overall row comes to the first
#   df_rmst_ro <- 
#     df_rmst %>%
#     mutate(is_overall = Characteristic == "Overall") %>%
#     arrange(desc(is_overall), Characteristic) %>%
#     select(-is_overall)
#   
#   
#   ## -- est (lcl, ucl) in a nice presentable manner
#   
#   clean_fmt <- function(x) sprintf('%.2f', 0 + round(x, 2))
#   
#   df_rmst_ro %<>%
#     group_by(Characteristic) %>%
#     mutate(
#       rmst_ci = glue("{clean_fmt(rmst)} ({clean_fmt(rmst_lcl)}, {clean_fmt(rmst_ucl)})"),
#       rmst_diff_ci = case_when(
#         Characteristic == "Overall" ~ "-",
#         row_number() == 1 ~ "Ref",
#         TRUE ~ glue("{clean_fmt(rmst_diff)} ({clean_fmt(rmst_diff_lcl)}, {clean_fmt(rmst_diff_ucl)})")
#       )
#     ) %>%
#     ungroup() %>%
#     mutate(level = ifelse(level == "Overall", "-", level)) %>%
#     rename(Stratum = level)
#   
#   ## -- take this to a final dataset
#   
#   ## -- add N to the summary dataframe
#   
#     ## -- create overall df, using all the data points from the primary dataset
#   d = glue("df_all_symps_{endpoint}")
#     
#     ## -- might need to change this based on the IPCW indicator
#   
#   data = 
#     if(ind_ipcw == 0){
#       df_list[[d]] %>% ## -- slice(1) removed because already sliced to have n=1 per ID
#         left_join(base_dat %>% select(id, covariates))
#     }else{
#       df_list_ipcw[[d]] %>%
#         left_join(base_dat %>% 
#                     filter(!if_any(all_of(covariates), is.na)) %>%
#                     select(id, covariates)) %>%
#         group_by(id) %>%
#         slice(1) %>%
#         ungroup()
#     }
#             
#   cat_cols = 
#     data %>% 
#       select(covariates)
#   
#   cat_cols %<>%
#     mutate(across(
#       everything(),
#       ~ as.character(haven::zap_labels(.x))  # remove label + coerce
#     )) %>%
#     select(where(is.character))
#   
#   cat_summary = map_dfr(
#     names(cat_cols),
#     ~ cat_cols %>%
#       count(!!sym(.x), name = "N") %>%
#       rename(Stratum = !!.x) %>%
#       mutate(Variable = .x),
#     .id = NULL
#   ) %>%
#     select(Variable, Stratum, N) %>%
#     filter(!is.na(Stratum))
#   
#   # cat_summary %<>%
#   #   arrange(match(paste(Variable, Stratum), paste(df_rmst_ro$Characteristic, df_rmst_ro$Stratum)))
#   
#   cat_summary %<>%
#     mutate(value = sapply(Variable, function(k) covariates_list[[k]]))
#   
#   df_rmst_ro %<>%
#     left_join(cat_summary %>% select(Variable, Stratum, N, value), by = c("Characteristic" = "value", "Stratum"))
#   
#   df_rmst_ro$N[1] <- nrow(data)
#   
#   df_final =
#     df_rmst_ro %>%
#     select(Characteristic, Stratum, N, rmst_ci, rmst_diff_ci)
#   
#   
#   ## -- save gt table of df_rmst to an html file
#   gt = df_final %>%
#     gt(groupname_col = "Characteristic") %>%
#     tab_options(
#       table.align = "left",
#       column_labels.font.weight = "bold",
#       table.border.top.width = px(1),
#       table.border.bottom.width = px(1),
#       table.border.top.color = "lightgray",
#       table.border.bottom.color = "lightgray"
#     ) %>%
#     cols_align(align = "left", columns = everything()) %>%
#     cols_label(rmst_ci = "RMST (95% CI)", 
#                rmst_diff_ci = "RMST difference (95% CI)") %>%
#     tab_header(title = 
#                  glue("Restricted mean survival time (RMST) and differences for days not {endpoint_label} from COVID-19 during the
#                   first {cutoff} days after symptom onset")
#     )
#     
#     gt %>%
#       gtsave(glue(paste0(out_path_html, "rmst_{endpoint}_{file_ext}_v2.html")))
#     
# 
#   ## -- if user wants it in word (docx) format)
#   
#   df_final %<>%
#     mutate(group_change = lag(Characteristic, default = first(Characteristic)) != Characteristic,
#            row_id = row_number())
#   
#   # Find the row numbers where new groups start (except the first row)
#   group_start_rows = df_final %>%
#     filter(group_change) %>%
#     pull(row_id)
#   
#   df_final %<>%
#     select(-c(group_change, row_id))
#   
#   # Create flextable
#   ft <- 
#     flextable(df_final) %>%
#       set_header_labels(
#         rmst_ci = "RMST (95% CI)",
#         rmst_diff_ci = "RMST difference (95% CI)"
#       ) %>%
#       merge_v(j = "Characteristic") %>%
#       valign(j = "Characteristic", valign = "top") %>%
#       align(align = "center", part = "header") %>%   # Center align just the header
#       align(align = "left", part = "body") %>%       # Left align the body cells
#       fontsize(size = 11, part = "all") %>%
#       font(fontname = "Calibri", part = "all") %>%
#       bold(part = "header") %>%
#       border_remove() %>%
#       hline_top(part = "header", border = fp_border(width = 1.2)) %>%
#       hline(part = "header", border = fp_border(width = 1.2)) %>%
#       hline(i = group_start_rows - 1, border = fp_border(width = 1, color = "black")) %>%
#       autofit() %>%
#       footnote(
#         i = 1, j = "rmst_diff_ci", part = "header",
#         value = as_paragraph("Assuming groups are independent"),
#         ref_symbols = "1"
#       )
#   
#   # Save to Word
#   title_text <- fpar(
#     ftext(glue("Restricted mean survival time (RMST) for days not {endpoint_label} during the first ",
#                "{cutoff} days after symptom onset."), 
#           prop = fp_text(font.family = "Calibri", font.size = 12, bold = TRUE))
#   )
#   
#   read_docx() %>%
#     # body_add_fpar(
#     #   fpar(ftext(
#     #     glue("Restricted mean survival time (RMST) for days not {endpoint_label} during the first ",
#     #          "{cutoff} days after symptom onset."),
#     #     prop = fp_text(font.size = 10)
#     #   ))
#     #   , style = "heading 1") %>%
#     
#     body_add_fpar(
#       title_text,
#       style = "Normal"
#     ) %>%
#     body_add_flextable(ft) %>%
#     print(target = glue(paste0(out_path_word, "rmst_{endpoint}_{file_ext}_v1.docx")))
#   
#   
# }
#   
# ## -- close all word documents
# system2("taskkill", args = c("/IM", "WINWORD.EXE", "/F"))
# 
# create_rmst_table(rms_comp_results, "resolution", ind_ipcw = 0)
# create_rmst_table(rms_comp_results, "recovery", ind_ipcw = 0)
# create_rmst_table(rms_ipcw_comp_results, "recovery", ind_ipcw = 1)
# create_rmst_table(rms_ipcw_comp_results, "resolution", ind_ipcw = 1)
# 
# 
# shell.exec(normalizePath(glue("G:/projects/cfar/COVID19/VISION/manuscripts/acute_symptoms/results/rmst/",
#                               "word/crude/rmst_resolution_crude_v1.docx")))
# 
# shell.exec(normalizePath(glue("G:/projects/cfar/COVID19/VISION/manuscripts/acute_symptoms/results/rmst/",
#                               "word/ipcw/rmst_resolution_ipcw_v1.docx")))
# 
# 












