rm(list=ls())
# libraries ---------------------------------------------------------------

library(magrittr)
library(tidyverse)
library(scales)
library(glue)

# import data -------------------------------------------------------------

# file_date <- "2025-01-17_1723"
file_date <- "2025-09-12_1439"
source(glue("FischerVISIONStudy-AcuteSymptomEndpoint_R_{file_date}.r"))

load("../../manuscripts/baseline/data/baseline_dataset.Rda") #base_dat

# pre-processing ----------------------------------------------------------
  
  ## !! - make adjustment for the new BL date (consent date to calculate days since sx onset)

  ## -- only include verified participants

verified <- readRDS("G:/projects/cfar/COVID19/VISION/data/derv/verified_participants_baseline_Dinelka.RDS")
df_ver <- subset(data, id %in% verified$id)

df_ver %<>% left_join(base_dat %>% dplyr::select(id, days_since_symptom_onset, inc_answers), by = "id")

  ## -- remove 'general_arm_1' & only days 1-14; only days_since_symp <= 7

filt_days_since_so <- 7

## -- artificial administrative censoring cut point (this was finalized in June 2025 by stats team)
admin_cens <- 18

df_ver %<>% 
  filter(redcap_event_name.factor %in% c("Baseline Day 1", paste0("Day ", 1:14)) &
                     redcap_event_name.factor != "General" &
           days_since_symptom_onset <= filt_days_since_so &
           inc_answers == 1) ## inc_answers - ensures only recent infection cohort, excludes rebounders

df_ver %<>%
  mutate(days_since_symptom_onset_v2 = days_since_symptom_onset + 
           as.numeric(gsub("\\D", "", redcap_event_name.factor)) - 1)

  ## -- some symptoms were asked based on an initial question. Need to manually impute.
    ### -- to give more context, some don't get to answer no to certain symptoms because they answered
    ### no to an initial question. We have to manually change these to "No".

vars_to_update <- setdiff(grep("\\.factor$", colnames(df_ver), value = TRUE), 
          c("redcap_event_name.factor", "redcap_repeat_instrument.factor", "sym_yn.factor"))

df_ver %<>%
  mutate(across(vars_to_update, ~ ifelse(sym_yn.factor == "No" & is.na(.), "No", as.character(.)))) %>%
  mutate(across(vars_to_update, as.factor))

    ### --- now for each specific sub-type of symptom

symp_types <- c("flu", "nvd", "heent", "mental", "other")
flu_symps <- c("fever", "sob_rest", "sob_activity", "cough", "chills", "headache", "sore_throat",
               "pain", "weakness", "joint", "fatigue")
nvd_symps <- c("nausea", "vomit", "diarrhea")
heent_symps <- c("runny_nose", "taste_none", "smell_none", "ringing", "vision")
mental_symps <- c("concentration", "memory", "dizzy")
other_symps <- c("tachycardia", "chest_pain", "rash", "tingling")

symps_list <- list(flu_symps, nvd_symps, heent_symps, mental_symps, other_symps)
names(symps_list) <- symp_types

all_symps <- paste0(unname(unlist(symps_list)), ".factor")

for(symp in symp_types){
  df_ver <- df_ver %>%
    mutate(across(paste0(symps_list[[symp]], ".factor"), ~ ifelse(!!sym(paste0(symp, "_acute_yn.factor")) == "No" & is.na(.), "No", as.character(.)))) %>%
    mutate(across(paste0(symps_list[[symp]], ".factor"), as.factor))
}



  ## -- we will skip handling this for now (doubt we have a high percentage of missingness)


## -- expanding the time points, to ensure that all records have times starting from 0 to their possible maximum

df_ver %<>%
  mutate(days_since_symptom_onset_v2 = as.numeric(days_since_symptom_onset_v2))

df_limits <- df_ver %>%
  group_by(id) %>%
  summarise(
    min_time = first(days_since_symptom_onset_v2),
    #min_time = 0,
    max_time = min(first(days_since_symptom_onset_v2) + 14 - 1, admin_cens),
    .groups = "drop"
  )

## -- summary of days since symptom onset

df_limits %>%
  filter(min_time == 0) %>%
  summarise(prop = percent(n()/nrow(df_limits), accuracy=0.01))

    ## -- such a small percentage of individuals join VISION on the exact day of symptom onset

df_ver %>%
  group_by(id) %>%
  slice(1) %>%
  ggplot(aes(x = days_since_symptom_onset_v2)) +
  geom_histogram()

all_times <- df_limits %>%
  rowwise() %>%
  mutate(days_since_symptom_onset_v2 = list(seq(min_time, max_time))) %>%
  #mutate(days_since_symptom_onset_v2 = list(seq(0, max_time))) %>%
  unnest(days_since_symptom_onset_v2)  # Expand time sequence for each ID

# Merge with the original dataframe
df_preproc <- all_times %>%
  left_join(df_ver, by = c("id", "days_since_symptom_onset_v2"))

## -- how many have days_since_sx_onset > 1 and symptoms resolve on either day 1 (study entry) or 2


# missing data ------------------------------------------------------------

df_miss_symps <- 
  df_preproc %>%
  mutate(any_na = rowSums(is.na(dplyr::select(., all_symps))) > 0) %>%
  filter(any_na == TRUE)

nrow(df_miss_symps) / nrow(df_preproc) * 100



# export data -------------------------------------------------------------

saveRDS(df_preproc, "data/df_preproc.RDS")
saveRDS(df_preproc, "data/df_preproc_v2.RDS")



# export only symptom data ------------------------------------------------

df_preproc %>%
  dplyr::select(id, redcap_event_name.factor,
                days_since_symptom_onset_v2, all_symps) %>%
  saveRDS("data/symptoms.RDS")








