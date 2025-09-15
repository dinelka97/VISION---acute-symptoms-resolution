# Description -------------------------------------------------------------
  
  ### -- use this script to pre-process the BL covariates that are to be used in downstream analysis. Need to convert some 
  ### -- to categorical, introduce new variables, etc.
  ### -- also used to create table 1 specific for symptom resolution


# Load libraries ----------------------------------------------------------

library(tidyverse)
library(magrittr)
library(tableone)
library(gtsummary)
library(flextable)
library(officer)
library(gt)

# Load Data ---------------------------------------------------------------

## -- baseline data to incorporate covariates
## -- load follow-up treatment data, as treatment is a time-varying covariate

load("../baseline/data/baseline_dataset.Rda")
#base_dat <- read_csv("../baseline/data/baseline_dataset.csv")
trt_data <- readRDS("data/trt_date.RDS")

## -- load tte event data to create table 1 (need to subset to tte ids only)

## -- tte event data
rds_files <- list.files("data/km/v3", pattern = "\\.RDS$", full.names = TRUE)
df_list <- lapply(rds_files, readRDS)
names(df_list) <- tools::file_path_sans_ext(basename(rds_files))


# Covariates --------------------------------------------------------------

## -- 4 age categories instead of 3
base_dat$age_cat <- cut(base_dat$age_derived,
                        breaks = c(17, 34, 49, 64, Inf),
                        labels = c("18-34", "35-49", "50-64", "65+"),
                        right = TRUE)


## -- add labels for sex at birth
base_dat$sex_birth <- factor(base_dat$sex_birth,
                             levels = c(1, 2),
                             labels = c("Female", "Male"))


## -- prior infection (binary)
base_dat$pi_cat <- factor(base_dat$covid_past_derived, 
                          levels = c(0, 1), 
                          labels = c("no prior infections", "1+ prior infections"))

## -- prior infection recency (three groups, similar to vaccination recency)
base_dat$pi_rec <- factor(base_dat$covid_past5_derived,
                          levels = c("None", "≤ 6 months", "> 6 months"),
                          labels = c("No prior infections", "≤ 6 months", "> 6 months"))


base_dat$bipoc_race_derived <- factor(base_dat$bipoc_race_derived, 
                                      levels = c(0, 1), 
                                      labels = c("White/Caucasian", "BIPOC"))

base_dat$race_new_derived <- factor(base_dat$bipoc_race_derived, 
                                    levels = c(0, 1, 2), 
                                    labels = c("Black/AA", "White/Caucasian", "other"))
base_dat$ethn_hispanic_latinx_derived <-
  factor(base_dat$ethn_hispanic_latinx_derived,
         levels = c(0, 1),
         labels = c("Not Hispanic/Latinx", "Hispanic/Latinx"))

base_dat %<>%
  mutate(n_vax_derived = factor(case_when(
    vax_number == 0 ~ "0",
    vax_number %in% c(1,2,3) ~ "1-3",
    vax_number %in% c(4,5) ~ "4-5",
    as.numeric(vax_number) >= 6 ~ "6+"
  )))

base_dat %<>%
  mutate(across(starts_with("med_history_v2___") & !ends_with(c("77", "77_neg")), ~ 
                  ifelse(. == "Yes", 1, ifelse(. == "No", 0, NA)))) %>%
  mutate(comorb_count = rowSums(select(., starts_with("med_history_v2___") & 
                                          !ends_with(c("77", "77_neg"))), na.rm = TRUE),
         comorb_count_cat = case_when(
           comorb_count == 0 ~ "0",
           comorb_count == 1 ~ "1",
           comorb_count >= 2 ~ "2+"
         ),
         comorb_count_cat = factor(comorb_count_cat,
                                   levels = c("0", "1", "2+"),
                                   labels = c("0", "1", "2+")))


base_dat %<>%
  mutate(
    comorb_count_cat2 = case_when(
      if_any(starts_with("med_history_v2___") & ends_with(paste(2:5)), ~ .x == 1) ~ "heart, HBP, lung, or diabetes",
      TRUE ~ "other"),
    comorb_count_cat2 = factor(comorb_count_cat2)
  )


## -- prior infections to deeper categories (collapse to 0, 1, 2+)


base_dat %<>%
  mutate(pi_cat_v2 = case_when(
    covid_past_derived == 0 ~ "0",
    covid_past1_derived == 1 ~ "1",
    covid_past1_derived >= 2 ~ "2+",
    #covid_past1_derived >= 3 ~ "3+"
  ),
  pi_cat_v2 = factor(pi_cat_v2))


## -- change labels for VAX_REC

base_dat$vax_rec_derived <- factor(base_dat$vax_rec_derived,
                                   #levels = c("Unvaccinated", "<= 6 months", "> 6 months"),
                                   levels = c("<= 6 months", "> 6 months", "Unvaccinated"),
                                   #labels = c("Unvaccinated", "≤ 6 months", "> 6 months"),
                                   labels = c("≤ 6 months", "> 6 months", "Unvaccinated"))

    ## -- change order to "<= 6 months", "> 6 months", "Unvaccinated"


## -- it is not required to create a new variable for smoking because Jess already has done it.
## -- but, I will go ahead and create the levels here.

base_dat$smoking_cat_derived_derived <- factor(base_dat$smoking_cat_derived_derived,
                                               levels = c(0, 1, 2),
                                               labels = c("Never", "Former", "Current"))


## -- create a variable for HYBRID IMMUNITY

base_dat %<>%
  mutate(hybrid_imm = 
           factor(
             interaction(vax_rec_derived, pi_cat, sep = "_"),
             levels = c("Unvaccinated_no prior infections", "Unvaccinated_1+ prior infections", "≤ 6 months_no prior infections",
                        "≤ 6 months_1+ prior infections", "> 6 months_no prior infections", "> 6 months_1+ prior infections"),
             labels = c("Unvaccinated, 0", "Unvaccinated, 1+", "Vaccinated ≤ 6 mo., 0", "Vaccinated ≤ 6 mo., 1+",
                        "Vaccinated > 6 mo., 0", "Vaccinated > 6 mo., 1+"))
  )

## -- SUBSET BASE_DAT TO ONLY RECENT INFECTION COHORT

base_dat %<>%
  filter(inc_answers == 1)


### -- EXTRACT FINAL BL DATASET TO BE USED DOWNSTREAM

base_dat %>% write.csv(file = "data/bl_derived.csv", row.names = FALSE)
base_dat %>% saveRDS(file = "data/bl_derived.rds")


##### --- FOLLOW-UP COVARIATES (TREATMENT)

trt_data %>% head()




# CREATE TABLE 1 ----------------------------------------------------------

covariates_tab1 <- c("age_derived", "age_cat", "sex_birth", "bipoc_race_derived", "ethn_hispanic_latinx_derived",
                     "bmi_cat_derived", "comorb_count_cat", "vax_rec_derived", "pi_cat", "urban_rural_derived",
                     "smoking_cat_derived_derived")

df_tab1 <- 
  base_dat %>%
    filter(id %in% df_list[["df_all_symps_resolution"]]$id) %>%
    select(c(id, all_of(covariates_tab1)))

df_tab1 <- df_tab1 %>%
  mutate(across(everything(), ~ `attr<-`(.x, "label", NULL))) ## -- remove labels

table1 <- 
  df_tab1 %>%
    select(all_of(covariates_tab1)) %>%
    tbl_summary(
      type = c(all_dichotomous() ~ "categorical", all_continuous() ~ "continuous2"),
      statistic = list(
        all_continuous() ~ c("{median} ({p25}, {p75})"),
        all_categorical() ~ "{n} ({p}%)"
      ),
      missing = "ifany",                 
      percent = "cell",
      missing_text = "Missing",
      label = list(sex_birth ~ "Sex at birth",
                   bipoc_race_derived ~ "BIPOC race",
                   ethn_hispanic_latinx_derived ~ "Hipanic/Latinx ethnicity",
                   bmi_cat_derived ~ "BMI",
                   comorb_count_cat ~ "Comorbidity count",
                   vax_rec_derived ~ "Vaccination recency",
                   pi_cat ~ "Prior infections",
                   urban_rural_derived ~ "Resdential urbanicity",
                   smoking_cat_derived_derived ~ "Smoking status")
    ) %>%
    # modify_variable_label(
    #   age_derived = "Median (Q1, Q3)",
    #   age_cat = "Age"
    # ) %>%
    modify_header(all_stat_cols(FALSE) ~ "**{level}**<br>N = {n}") %>%
    add_variable_group_header(
      header = "Age",  # The desired group header
      variables = c(age_derived, age_cat) # The variables to group
    ) %>%
    modify_table_body(
      ~.x %>%
        dplyr::filter(! (variable %in% c("age_derived", "age_cat") & row_type == "label"))
    ) %>%
    bold_labels() %>%
    modify_table_styling(
      rows = label == "Age", # Select the "Age" row based on its label
      columns = label,           # Apply the styling to the "label" column
      text_format = "bold"       # Set the text format to "bold"
    )

table1

tbl_stack(list(table1)) %>% 
  as_gt() %>%
  #tab_header(title = md("## Baseline characteristics of VISION analysis cohort for acute symptom resolution")) %>% 
  #set_table_properties(width = 1, layout = "auto fit") %>% # Set full width and auto fit %>%
  gt::tab_options(table.width = gt::pct(100)) %>%
  #gt::cols_width(label~px(300), everything()~px(120)) %>%
  gtsave("tables/table1_resolution.docx") ## -- need to auto fit to window manually from word




