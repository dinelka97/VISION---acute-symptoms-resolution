### -- this is the upto date version that is to be used.

rm(list=ls())
# libraries ---------------------------------------------------------------

library(tidyverse)
library(magrittr)
library(ggplot2)
library(glue)
library(gridExtra)
library(slider)

# import data -------------------------------------------------------------

df <- readRDS("data/df_preproc_v2.RDS")

  ## -- import treatment data from Ning's code, because 

## -- this is directly sourced in case NZ makes any changes to this script it is automatically reflected here.
source("G:/projects/cfar/COVID19/VISION/code/create_tlf/Analysis/return_to_health/treatment_acute_paper.R") #trt_dat

  ## -- join with the df

df %<>%
  left_join(trt_dat, by = "id") %>%
  mutate(trt_start = ifelse(trt_ever == 1, trt_start, NA))

#df <- df %>% group_by(id) %>% filter(!is.na(days_since_symptom_onset))

symp_types <- c("flu", "nvd", "heent", "mental", "other")
flu_symps <- c("fever", "sob_rest", "sob_activity", "cough", "chills", "headache", "sore_throat",
               "pain", "weakness", "joint", "fatigue")
nvd_symps <- c("nausea", "vomit", "diarrhea")
heent_symps <- c("runny_nose", "taste_none", "smell_none", "ringing", "vision")
mental_symps <- c("concentration", "memory", "dizzy")
other_symps <- c("tachycardia", "chest_pain", "rash", "tingling")

activ2_symps <- c("fever", "cough", "sore_throat", "pain", "fatigue", "headache", "chills",
                  "sob_rest", "sob_activity", "runny_nose", "nausea", "vomit", "diarrhea") ## -- active2 combine sob. VISION combines runny & stuffy nose

symps_list <- list(flu_symps, nvd_symps, heent_symps, mental_symps, other_symps)
names(symps_list) <- symp_types

all_symps <- paste0(unname(unlist(symps_list)), ".factor")
activ2_symps <- paste0(activ2_symps, ".factor")

## -- censoring due to missingness
cens_miss_thres <- 3 # - if 3+ days of missingness, censor at last available diary

# creating endpoint variables ---------------------------------------------

  ## -- I have used the baseline date (arm1_dt) as the starting point of time to recovery
  ## -- a new outcome variable is created for the time to acute symptom recovery (in days)


# create function to obtain the time to event (tte) and the event status (event) ------------------------------------------

  ## -- the function for recovery is a tad bit slow because we are manipulating 26 symptom columns.
  ## -- need to look into it a little more and make it more efficient.
  ## -- takes about 3.5 minutes for it to run.

create_df <- function(df_init, symptoms, endpoint, symptom_ind = "complete"){
  
  #### -- TWO DIFFERENT ALGORITHMS FOR RECOVERY/IMPROVEMENT & RESOLUTION
  
  #### -- TIME THIS, CAUSE IT IS NOT TOO EFFICIENT
  
  start_time = Sys.time()
  
  ##### -- RECOVERY/IMPROVEMENT
  
  if(endpoint == "recovery"){

    ## -- change the levels of the symptom factor variables
    
    severity_levels <- c("No", "Mild", "Moderate", "Severe")
    
    df_init %<>%
      mutate(across(
        .cols = symptoms,  # Replace with your actual column names
        .fns = ~ factor(.x, levels = severity_levels, ordered = TRUE)
      ))
    
    ## -- indicator for those who had no symptoms at a given day
    
    df_temp =
      df_init %>%
        ungroup() %>% ## -- the dataset it comes from has been grouped by ID
        mutate(ind_bl =
                 case_when(#!left_cens ~ "FALSE",
                   apply(df_init[, symptoms], 1, function(x) any(is.na(x))) ~ NA,
                   endpoint %in% c("resolution", "recovery") ~ 
                     as.character(apply(df_init[, symptoms], 1, function(x) all(x == "No")))
                   )
        )
               
    ## -- identifying those IDs who had symptom resolution/improvement by study entry. For these tte is computed separately
    
    bl_ep_ids = 
      df_temp %>%
        group_by(id) %>%
        slice_head(n = 2) %>%
        summarise(keep = all(ind_bl == TRUE), .groups = "drop") %>%
        filter(keep) %>%
        pull(id)
    
    df_sub =
      df_temp %>%
        filter(!id %in% bl_ep_ids)
    
    ## -- for those in this, improvement will be met based if following met for 2+ consecutive days
      ## -- 1. Mild to No, 2. Moderate to Mild/No, 3. Severe to Mild/No, or 4. No to no
    
      ## -- try checking this on one symptom (and then apply to other symptoms)
      ## -- helper function to get improvement for one symptom
    get_improvement <- function(vec, id, timepoint) {
      # Convert to ordered factor (if not already)
      
      lag_vec = ave(vec, id, FUN = function(x) dplyr::lag(x))
      
      case_when(
        is.na(lag_vec) ~ NA_character_,
        lag_vec == "Severe" & vec %in% c("No", "Mild") ~ "Improved",
        (lag_vec == "Severe" & vec %in% c("Severe", "Moderate")) |
          (lag_vec == vec) ~ "Unchanged",
        lag_vec > vec ~ "Improved",
        lag_vec < vec ~ "Worsened"
      )
    }
    
      ## -- helper function to get recovery status for each symptom
    # Recovery status helper
    
    get_recovery <- function(imp_vec, sev_vec, id) {
      ave(seq_along(imp_vec), id, FUN = function(i) {
        x = imp_vec[i]
        sev = sev_vec[i]
        lead1 = dplyr::lead(x, 1)
        
        dplyr::case_when(
          is.na(x) ~ NA_character_,
          x == "Improved" & lead1 %in% c("Improved", "Unchanged") ~ "Recovered",
          sev == "No" & lead1 == "Unchanged" ~ "Recovered",
          TRUE ~ "Not Recovered"
        )
      })
    }
    
      ## -- get improvement indicators for all the symptoms
    df_sub %<>%
      bind_cols(
        map_dfc(symptoms, function(col) {
          imp_col <- paste0("imp_", col)
          tibble(!!imp_col := get_improvement(df_sub[[col]], df_sub$id, df_sub$redcap_event_name.factor))
        })
      )
    
    recovery_cols <- map_dfc(symptoms, function(symptom) {
      imp_col = paste0("imp_", symptom)
      recov_col = paste0("recov_", symptom)
      
      tibble(!!recov_col := get_recovery(
        imp_vec = df_sub[[imp_col]],
        sev_vec = df_sub[[symptom]],
        id = df_sub$id
      ))
    })
    
    
    df_sub <- bind_cols(df_sub, recovery_cols)
    
    ## -- indicator column if all symptoms have recovered
    
    df_sub %<>%
      rowwise() %>%
      mutate(
        all_recovered = {
          vals = c_across(starts_with("recov_"))
          if (all(is.na(vals))) {
            NA
          } else {
            all(vals == "Recovered", na.rm = TRUE)*1
          }
        }
      ) %>%
      ungroup()
    
    
    ## -- use view to do a quick check (comment after done)
    # selected_cols = 
    #   df_sub %>% 
    #     select(matches("(_)?(headache|cough|chills)\\.factor$")) %>%
    #     colnames()
    # 
    # View(
    #   df_sub %>%
    #     select(id, selected_cols, all_recovered)
    # )
    
    ## -- now get the tte, and combine with the other IDs who experienced recovery/resolution at BL.
    
     df_sub %<>%
       group_by(id) %>%
       mutate(
         tte_temp = 
           if (any(all_recovered == 1, na.rm = TRUE)) {
            days_since_symptom_onset_v2[which(all_recovered == 1)[1]]
           } else {
            NA ## -- account for right censoring later
         }
       ) %>%
       ungroup()
     
     # View(
     #   df_sub %>%
     #     select(id, selected_cols, all_recovered, tte)
     # )
     
     ## -- prior to taking this dataset to the next step, let's keep only the all_resolved, and 
      # -- tte columns (this won't be the final tte)
     
     symps = sub("\\.factor$", "", symptoms)
     prefixes <- c("imp", "recov")  
  
     # collapse prefixes and symptoms to regex patterns
     prefix_pattern <- paste(prefixes, collapse = "|")
     symp_pattern <- paste(symps, collapse = "|")
     
     df_sub %<>%
         select(-matches(paste0("^(", prefix_pattern, ")_(", symp_pattern, ")\\.factor$")), -ind_bl)
     
     
     ### -- account for censoring due to missingness
     
     
     df_sub %<>%
       mutate(
         na_ind = apply(df_sub[, symptoms], 1, function(x) any(is.na(x)))
       ) %>%
       group_by(id) %>%
       mutate(
         ind_cens = slide_lgl(na_ind, ~ all(.x == TRUE), .before = 0, .after = cens_miss_thres - 1, .complete = TRUE),
         t_cens1 = as.numeric(days_since_symptom_onset_v2[which(ind_cens == TRUE)[1]]) - 1, ## time to diary right before threshold
         t_cens2 = case_when(is.na(t_cens1) ~ as.numeric(days_since_symptom_onset_v2)[max(which(!na_ind), na.rm = TRUE)],
                             TRUE ~ NA),
         event_temp = any(all_recovered == 1),
         tte = min(tte_temp, t_cens1, t_cens2, na.rm = TRUE),
         event = 
           case_when(event_temp == TRUE & is.na(t_cens1) ~ 1,
                     event_temp == TRUE & !is.na(t_cens1) & t_cens1 >= tte ~ 1, 
                     TRUE ~ 0)
      ) 
     
     ## -- save censoring due to missingness percentage
     
     n_miss_cens = df_sub %>%
       filter(ind_cens) %>%
       distinct(id) %>%
       nrow()
     
     miss_cens_prop = n_miss_cens / (df_sub %>% distinct(id) %>% nrow())
     
     line <- paste0(endpoint, ": ", miss_cens_prop)
     
     write(paste0(line, "\n"), file = "cens_miss.txt", append = TRUE)
     
     df_sub %<>%
        select(-c(all_recovered, tte_temp, event_temp, na_ind, ind_cens, t_cens1, t_cens2))
     
     
     ### -- prepare the tte variable for the set which had recovery at BL
     df_sub %<>%
         bind_rows(
           df_init %>%
             filter(id %in% bl_ep_ids) %>%
             group_by(id) %>%
             mutate(tte = unique(min_time),
                    event = 1) %>%
             ungroup()
         )
     
     ## -- keep only needed variable columns
     
     df_sub %<>%
       select(id, symptoms, days_since_symptom_onset, days_since_symptom_onset_v2, tte, trt_start, trt_ever, event) %>%
       arrange(days_since_symptom_onset, desc(event), tte, id, days_since_symptom_onset_v2)
     
     
     
  } else{ ## -- this is for symptom resolution. Code is different because recovery was done after and thinking was different (DN)
    
    df_sub = 
      df_init %>%
      ungroup() %>%
      mutate(ind_rec = 
               case_when(#!left_cens ~ "FALSE",
                 apply(df_init[, symptoms], 1, function(x) any(is.na(x))) ~ NA,
                 endpoint == "resolution" ~ as.character(apply(df_init[, symptoms], 1, function(x) all(x == "No"))))
      ) %>%
      mutate(
        na_ind = apply(df_init[, symptoms], 1, function(x) any(is.na(x))),
        ind_rec = recode(ind_rec, "TRUE" = "yes", "FALSE" = "no", .missing = "missing"))
    
    df_sub %<>%
      group_by(id) %>%
      mutate(ind_rec2 = lag(ind_rec == "yes", default = FALSE) & ind_rec == "yes",
             ind_rec2 = case_when(ind_rec2 == TRUE ~ "yes",
                                  ind_rec2 == FALSE ~ "no",
                                  is.na(ind_rec2) ~ "missing"),
             ind_cens = slide_lgl(ind_rec, ~ all(.x == "missing"), .before = 0, .after = 3, .complete = TRUE), ## 4 cons. missing diaries
             t_cens1 = as.numeric(days_since_symptom_onset_v2[which(ind_cens == TRUE)[1]]) - 1, ## time to diary right before 4 cons.
             t_cens2 = case_when(is.na(t_cens1) ~ as.numeric(days_since_symptom_onset_v2)[max(which(ind_rec != "missing"), na.rm = TRUE)],
                                 TRUE ~ NA),
             event_temp = any(ind_rec2 == "yes"),
             tte_temp = as.numeric(days_since_symptom_onset_v2[which(ind_rec2 == "yes")[1]-1]), ## -- pick the first day from the two consecutive days of recovery/resolution
             tte = min(tte_temp, t_cens1, t_cens2, na.rm = TRUE),
             event = 
               case_when(event_temp == TRUE & is.na(t_cens1) ~ 1,
                         event_temp == TRUE & !is.na(t_cens1) & t_cens1 >= tte ~ 1, 
                         TRUE ~ 0)) %>%
      ungroup()
    
    ## -- save censoring due to missingness percentage
    
    n_miss_cens = df_sub %>%
      filter(ind_cens) %>%
      distinct(id) %>%
      nrow()
    
    miss_cens_prop = n_miss_cens / (df_sub %>% distinct(id) %>% nrow())
    
    line <- paste0(endpoint, ": ", miss_cens_prop)
    
    write(paste0(line, "\n"), file = "cens_miss.txt", append = TRUE)
    
    
    df_sub %<>%
      select(id, symptoms, days_since_symptom_onset, days_since_symptom_onset_v2, trt_start, trt_ever, tte, event) %>%
      arrange(days_since_symptom_onset, desc(event), tte, id, days_since_symptom_onset_v2)
    
  }
   
  
  ### -- adjustments due to censoring. The final dataset will include the times for both left censoring (left and right times),
    # -- and also if we want to do all right censored methods, a mid-point based on the left and right extremes.
  
  df_sub %<>%
    group_by(id) %>%
    mutate(left = 
             case_when(
               tte == days_since_symptom_onset & event == 1 ~ 1,
               TRUE ~ tte
             ),
           right = 
             case_when(
               tte == days_since_symptom_onset & event == 1 ~ as.numeric(first(days_since_symptom_onset)),
               tte != days_since_symptom_onset & event == 1 ~ tte,
               event == 0 ~ Inf,
             ),
           event_ic = ## -- event_ic: interval censoring event time without considering any dropout 
             case_when(
               tte == days_since_symptom_onset & event == 1 & left != right ~ 3,
               TRUE ~ event
             )
    ) %>%
    ungroup()
  
  ## -- for those with left-censoring, pick mid-point as event date (can use this with only right-censoring then)
   # -- if this approach is considered, then for those left-censored previously will now be considered to have an event
  
  df_sub %<>%
    group_by(id) %>%
    arrange(id, days_since_symptom_onset_v2) %>%
    mutate(tte_rc = ## -- tte_rc: time to event prior to censoring out the treated population 
             case_when(
               tte == days_since_symptom_onset & event == 1 ~ (left + right)/2, 
               TRUE ~ tte
             ),
           event_rc =  ## -- event_rc: right censored event times prior to considering any censoring out the treated population
             case_when(
               tte == days_since_symptom_onset & event == 1 ~ 1,
               TRUE ~ event
             ),
           ## -- need to code up to account for censoring from treatment (event_final is one that accounts for treatment also)
           
           event_rc_trt = ## -- also considering treatment, and subsetting to only the untreated population
             case_when(
               event_rc == 1 & (tte_rc < trt_start | is.na(trt_start)) ~ 1,
               #event == 1 & (tte_rc < trt_start | is.na(trt_start)) ~ 1,
               TRUE ~ 0
             ),
           
           tte_rc_trt = ## -- similar logic to event_rc_trt, look at the minimum between tte_rc and now, trt start date
             min(tte_rc, trt_start - 1, na.rm = TRUE)
  
           ) %>%
    ungroup()
  
    ## -- all seems to be fine. Let's move df_sub to be df_final
    df_final = 
      df_sub %>%
        group_by(id) %>%
        slice(1) %>% 
        ungroup()
        
    
    end_time = Sys.time()
    
    print(glue("Time (mins): {as.numeric(end_time - start_time, units = 'mins')}"))
    
    return(df_final)
    
}
    
      ## -- quick recap (user can use whatever tte and event, depending on the method used)
      ## -- tte: time to event (this includes left-censoring events but does not provide enough info for it)
      ## -- event: 0 - right-censored, 1 - event/recovery/resolution, 3 - interval censored
      ## -- left, right: the left and right ends of the interval. Equal if event observed
      ## -- tte_rc: if we assume not interval censoring, and for those interval censored the tte is the midpoint of the interval
      ## -- event_rc: this is binary. 0 - all right-censored, 1 - event (left-censored from before will join this category)
  
  
    
    ### -- the function 'plot_func' is incomplete.
    
# plot_func <- function(df, plot_range){
#   
#   if((df_final %>% filter(is.na(tte_rc)) %>% nrow()) > 0){
#     stop("Missing (NA) values found in the right-censored tte")
#   }
#   
#   if(
#     df_final %>%
#     group_by(id) %>%
#     summarise(all_same = (n_distinct(tte_rc) == 1) & (n_distinct(event_rc) == 1)) %>%
#     filter(all_same == FALSE) %>%
#     nrow() > 0){
#     stop("Different tte or event values found for a given ID")
#   }
#   
#   df_plot.id <- df_plot %>%
#     distinct(id) %>%
#     mutate(newid = row_number())
#   
#   df_plot %<>% 
#     left_join(df_plot.id, by = "id")
#   
#   # visualizing time to event (recovery & resolution) --------------------------------------------
#   
#   idrange = df_plot.id$id[plot_range]
#   
#   p <- ggplot(df_plot %>% filter(id %in% idrange), aes(x = days_since_symptom_onset_v2, y = factor(newid), fill = factor(ind_rec, levels = c("yes", "no", "missing")))) +
#     geom_tile(color = "white") +
#     scale_fill_manual(values = c("yes"="#00A5AD", "no"="#EF426F", "missing"="grey"), labels = c("yes", "no", "missing")) +
#     labs(x = "Days since symptom onset", y = "ID", fill = glue("Acute symptom {endpoint} indicator")) +
#     geom_point(aes(shape = factor(event_final, levels = c("0", "1")), x = as.numeric(tte_final), y = factor(newid)), color = "black", size = 2) +
#     scale_shape_manual(values = c("1" = 16, "0" = 1), name="Outcome", labels=c("Censor", "Event")) +
#     scale_y_discrete(labels = df_plot %>%
#                        filter(id %in% idrange) %>%
#                        pull(id) %>%
#                        unique()) +
#     guides(
#       fill = guide_legend(override.aes = list(shape = NA)),
#       shape = guide_legend(order = 2)
#     ) +
#     theme_minimal()
#   
#   output = list(plot = p, sub_data = df_plot %>% filter(id %in% idrange),
#                 all_data = df_plot %>% select(-all_of(symptoms)))
#   
#   return(output)
#   
# }

  
  

# export data -------------------------------------------------------------

  ## -- first gather all datasets with tte, day of entry, and event

symps <- list(all_symps, activ2_symps); names(symps) <- c("all_symps", "activ2_symps")
endpoints <- c("recovery", "resolution")
df_out_list <- vector("list", length = length(symps) * length(endpoints))
names(df_out_list) <- apply(expand.grid(names(symps), endpoints), 1, paste, collapse = "_")

  ## -- folder v2 contains with time points starting from 0 (no need to consider an extended KM plot here)

#for(s in names(symps)){
for(s in "all_symps"){
  for(e in endpoints){
    #x = acute_symptom(df, symps[[s]], e, plot_range = 1:250)
    #df_out_list[[glue(s, "_", e)]] = x[["all_data"]]
    
    x = create_df(df, symps[[s]], e)
  
    #saveRDS(x[["all_data"]], glue("data/km/v2/df_{s}_{e}.RDS"))
    saveRDS(x, glue("data/km/v3/df_{s}_{e}.RDS"))
    
    print(glue("Endpoint derivation complete: {s}, {e}"))
  }
}


rm(list=ls())


  
  
  
  
  
  
