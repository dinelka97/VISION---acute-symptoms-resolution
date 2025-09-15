# Code for creating KM plots for symptom resolution
# Just symptom resolution, not recovery
# Using ipcw data instead of crude
  ## -- datasets used in this script are derived in km_rmstUV.R, and previous dependent scripts.
  ## -- only the datasets have been brought to this script to keep this script at its minimum
  ## -- a list generated in km_rmstUV.R holds all required dataframes to run this script

rm(list = ls())

library(ggplot2)
library(tidyverse)
library(RColorBrewer)
library(survminer)
library(survival)
library(patchwork)
library(officer)
library(rms)
library(nnet)
library(png)
library(grid)
library(grImport2)
library(jpeg)
library(cowplot)
library(glue)
library(officer)

setwd("G:/projects/cfar/COVID19/VISION/manuscripts/acute_symptoms")

## folder paths
fig_path_ipcw <- "figures/resolution_final/"

## -- load already prepped dataset from 'km_rmstUV/R'

dfs_km <- readRDS("data/km/v3/dfs_km.RDS") # from km_rmstUV.R

covariates <- c(
  'age_cat',
  'sex_birth',
  'bipoc_race_derived',
  'ethn_hispanic_latinx_derived',
  "bmi_cat_derived",
  "comorb_count_cat",
  'vax_rec_derived',
  'pi_cat',
  "urban_rural_derived",
  "smoking_cat_derived_derived" ## -- using the variable derived by Jess
)
  
title_map <- list(
  "overall" = 'Overall',
  "age_cat" = "Age (years)",
  "sex_birth" = "Sex at Birth", 
  "bipoc_race_derived" = "Race", 
  "ethn_hispanic_latinx_derived" = "Ethnicity", 
  "urban_rural_derived" = "Urbanicity",
  "bmi_cat_derived" = "BMI", 
  "pi_cat" = "Prior Infections",
  "comorb_count_cat" = "Comorbidity",
  "vax_rec_derived" = "Vaccination recency",
  "smoking_cat_derived_derived" = "Smoking status"
)


kmplot <- function(km_fit, km_dat, plot_title, table.font = 4.2) {
  
  if(!is.null(km_fit$strata)){
    km_fit$strata <- setNames(
      km_fit$strata, gsub(".*?=", "", gsub(", *.*?=", ", ", names(km_fit$strata)))
    )
  }
  
  p = ggsurvplot(
    km_fit,
    data = km_dat,
    fun = "event",
    title = plot_title,
    xlab = "Days since symptom onset",
    ylab = "Probability of acute symptom resolution",
    size = 1,                
    conf.int = TRUE, 
    conf.int.alpha = 0.35, 
    risk.table = TRUE,       
    #risk.table.col = NA,
    risk.table.height = 0.25,
    break.time.by = 2,
    #xlim=c(0, max(df$tte_final) + 1),
    #ggtheme = theme_bw(),
    # ggtheme = theme_classic(base_size=9),
    ggtheme = theme_bw(base_size=9),
    tables.theme = theme_cleantable(), ## remove if needed
    # tables.theme = theme(element_text(size = 9, lineheight = 1.8)),
    palette=brewer.pal(
      n = length(names(km_fit$strata)),
      name = "Set2"),
    # legend.labs = legend_labels_to_use,
    risk.table.y.text = FALSE,
    #risk.table.breaks = seq(0, max(df$tte_final) + 1, by = 2),
    #risk.table.x.text = TRUE,
    risk.table.fontsize = table.font, #2.25,
    legend.title=""
  )
  
  ## -- get number of categories
  n_labels <- length(names(km_fit$strata))
  
  legend_rows <- if (n_labels > 4) 2 else 1
  
  p$plot = p$plot + 
    coord_cartesian(ylim = c(0, 1), 
                    #xlim = c(0, t + 1),
                    xlim = c(0, 18),
                    clip = 'on') +
    theme(
      plot.title = element_text(#hjust = 0.5, 
                                size = 14, face = "bold"), 
      axis.title.x = element_blank(),
      axis.text.x = element_text(size = 12),
      axis.title.y = element_text(size = 12),
      axis.text.y = element_text(size = 12),
      legend.position = "top",          # inside plot (x, y in 0–1)
      #legend.justification = c("right", "top"),
      legend.title = element_text(size = 7),
      legend.text  = element_text(size = 11.5)) +
    guides(color = guide_legend(nrow = legend_rows))

  p$table = p$table +
    coord_cartesian(#xlim = c(0, t + 1),
      xlim = c(0, 18),
      clip = 'on') +
    theme(axis.title.y = element_blank(),
          #axis.text.x = element_text(size = 10),
          axis.text.x = element_blank(),
          axis.title.x = element_text(size = 12),
          plot.title = element_text(size = 14, face = "bold"))
  
  combined_plot <- 
    (p$plot / p$table) + 
      plot_layout(heights = c(3, 2)) &
      theme(plot.margin = margin(5, 2, 10, 5))

  return(combined_plot)
}



# IPCW - univariate  ----------------------------------------------

data_km_ipcw <- dfs_km[["ipcw_resolution"]] 

for(cov in c('overall', covariates)){
  if(cov == 'overall'){
    formula_str <- "Surv(time = tte_start, time2 = time, event = event_final) ~ 1"
  } else {
    formula_str <- paste("Surv(time = tte_start, time2 = time, event = event_final) ~", cov)
  }
  
  fit <- survfit(as.formula(formula_str), data = data_km_ipcw, weights = ipcw, robust = T, id = id)
  
  title <- title_map[[cov]]
  p <- kmplot(fit, data_km_ipcw, title)
  
  # Save to jpg
  filename <- sprintf("%s_resolution_ipcw.jpg", cov)
  ggsave(paste0(fig_path_ipcw, filename), width = 5.3, height = 5.3, units = "in", dpi = 300)
  #print(p)  
  #dev.off()
}


## Hybrid immunity: IPCW + IPTW -------------------------------------------------------------

data_iptw <- dfs_km[["iptw_resolution"]]

iptw_ipcw_km_hybrid <- survfit(Surv(time = tte_start, time2 = time, event_final) ~ 
                                 vax_rec_derived + pi_cat,
                               data = data_iptw, weights = iptw_ipcw, robust=T, id=id)

names(iptw_ipcw_km_hybrid$strata) <- c("Unvax, no pi", "Unvax, 1+ pi", "Vax <=6 mo., no pi", "Vax <=6 mo., 1+ pi",
                                       "Vax >6 mo., no pi", "Vax >6 mo., 1+ pi")

iptw_p <- kmplot(km_fit = iptw_ipcw_km_hybrid, km_dat = data_km_ipcw, plot_title = "Hybrid Immunity")


# Save to jpg
iptw_filename <- "hybrid_imm_resolution_ipcw.jpg"
ggsave(paste0(fig_path_ipcw, iptw_filename), width = 5.3, height = 5.3, units = "in", dpi = 600)
# print(iptw_p)  
# dev.off()


### -- paneling generated figures

  ## -- PRIMARY PAPER

#   ## (a) panel including smoking
# 
# covs_main_plot <- c("overall", "age", "sex", "comorb", "vax", "pi", "smoking") ## -- includes just the prefixes
# 
# ind_plots <- list.files(fig_path_ipcw, pattern = "\\.png$", full.names = TRUE)
# ind_plots <- sapply(covs_main_plot, function(prefix) { grep(prefix, ind_plots, value = TRUE)})
# 
# ## -- read them in as grobs
# 
# plots <- lapply(ind_plots, function(f) {
#   img = readPNG(f)
#   rasterGrob(img, interpolate = TRUE)
# })
# 
# plot_check <- plot_grid(plotlist = plots, 
#                         ncol = 2,
#                         labels = paste0("(", LETTERS[1:length(plots)], ")"),
#                         label_size = 6,
#                         #label_x = 0.02,                      # label position
#                         label_y = 0.95
# )
# 
# ggsave(glue(fig_path_ipcw, "/primary/primary_wsmoking.png"), plot = plot_check, height = 6, width = 8, dpi = 900)
# #ggsave(glue(fig_path_ipcw, "/primary/panel_wsmoking.pdf"), plot = plot_check, height = 6, width = 8, dpi = 600)


## (b) panel excluding smoking

covs_main_plot <- c("overall", "age", "sex", "vax", "pi", "hybrid_imm") ## -- includes just the prefixes

ind_plots <- list.files(fig_path_ipcw, pattern = "\\.jpg$", full.names = TRUE)
ind_plots <- sapply(covs_main_plot, function(prefix) { grep(prefix, ind_plots, value = TRUE)})

## -- read them in as grobs

plots <- lapply(ind_plots, function(f) {
  img = readJPEG(f)
  rasterGrob(img, interpolate = TRUE)
})

plot <- plot_grid(plotlist = plots, 
                        ncol = 2,
                        labels = paste0("(", LETTERS[1:length(plots)], ")"),
                        label_size = 7.5,
                        label_x = 0.02,                      # label position
                        label_y = 0.99,
                        align = "hv"
)

ggsave(glue(fig_path_ipcw, "/primary/primary_v1.jpg"), plot = plot, height = 11, width = 8, dpi = 600)
#ggsave(glue(fig_path_ipcw, "/primary/primary_wosmoking.pdf"), plot = plot_check, height = 4, width = 8, dpi = 600)

doc_path <- "figures/resolution_final/primary/resolution_main.docx"

# open doc if exists, ow create new file
if (file.exists(doc_path)) {
  doc = read_docx(doc_path)
} else {
  doc = read_docx()
}

doc <- doc %>%
  body_add_img(src = glue(fig_path_ipcw, "/primary/primary_v1.jpg"),
               width = NA,   # typical width for Word with 1in margins
               height = NA)   # auto scale height

# Save back
print(doc, target = doc_path)

print("Saved paneled figure for main")


  ## SUPPLEMENT

  ## -- (a) with smoking included

covs_supp_plot <- c("bipoc", "ethn", "bmi", "comorb", "smoking", "urban") ## -- includes just the prefixes

ind_plots <- list.files(fig_path_ipcw, pattern = "\\.jpg$", full.names = TRUE)
ind_plots <- sapply(covs_supp_plot, function(prefix) { grep(prefix, ind_plots, value = TRUE)})

## -- read them in as grobs

plots <- lapply(ind_plots, function(f) {
  img = readJPEG(f)
  rasterGrob(img, interpolate = TRUE)
})

plot <- plot_grid(plotlist = plots, 
                        ncol = 2,
                        labels = paste0("(", LETTERS[1:length(plots)], ")"),
                        label_size = 7.5,
                        label_x = 0.02,                      # label position
                        label_y = 0.99,
                        align = "hv"
)

ggsave(glue(fig_path_ipcw, "/supplement/supp_v1.jpg"), plot = plot, height = 11, width = 8, dpi = 600)

print("Saved paneled figure for supplement")


## -- (b) with smoking excluded
# 
# covs_main_plot <- c("") ## -- includes just the prefixes (confirm the order)
# 
# ind_plots <- list.files(fig_path_ipcw, pattern = "\\.png$", full.names = TRUE)
# ind_plots <- sapply(covs_main_plot, function(prefix) { grep(prefix, ind_plots, value = TRUE)})
# 
# ## -- read them in as grobs
# 
# plots <- lapply(ind_plots, function(f) {
#   img = readPNG(f)
#   rasterGrob(img, interpolate = TRUE)
# })
# 
# plot_check <- plot_grid(plotlist = plots, 
#                         ncol = 2,
#                         labels = paste0("(", LETTERS[1:length(plots)], ")"),
#                         label_size = 6,
#                         #label_x = 0.01,                      # label position
#                         label_y = 0.95
# )
# 
# ggsave(glue(fig_path_ipcw, "/supplement/supp_wosmoking.png"), plot = plot_check, height = 4, width = 8, dpi = 600)


## Hybrid immunity: IPCW + IPTW (w/ dropout and treatment censoring) -------------------

data_iptw_trt <- dfs_km[["iptw_trt_resolution"]]

iptw_ipcw_trt_km_hybrid <- survfit(Surv(time = tstart, time2 = tte_final, event = event_rc_trt) ~  ## -- tstart, tte_final from prev der
                                 vax_rec_derived + pi_cat,
                               data = data_iptw_trt, weights = iptw_ipcw_trt, robust=T, id=id)

names(iptw_ipcw_trt_km_hybrid$strata) <- c("Unvax, no pi", "Unvax, 1+ pi", "Vax <=6 mo., no pi", "Vax <=6 mo., 1+ pi",
                                            "Vax >6 mo., no pi", "Vax >6 mo., 1+ pi")

iptw_p <- kmplot(km_fit = iptw_ipcw_trt_km_hybrid, km_dat = data_iptw_trt, plot_title = "Hybrid Immunity")


# Save fig
iptw_trt_filename <- "hybrid_imm_resolution_iptw_trt.jpg"
ggsave(filename = paste0(fig_path_ipcw, iptw_trt_filename), plot = iptw_p, 
       width = 5.3, height = 5.3, units = "in", dpi = 600)
#print(iptw_p)  
#dev.off()

print("Save hybrid immunity plot for untreated scenario")



# Combine to a word document

# png_files <- list.files(fig_path_ipcw, pattern = "\\.png$", full.names = TRUE)
# # Reorder so "overall_resolution_ipcw.png" is first
# target_file <- "overall_resolution_ipcw.png"
# first_file <- png_files[basename(png_files) == target_file]
# other_files <- png_files[basename(png_files) != target_file]
# png_files <- c(first_file, other_files)

# Reorder
# file_prefixes <- c("overall", covariates, "hybrid_imm")
# png_files <- paste0(fig_path_ipcw, file_prefixes, "_resolution_ipcw.png")
# 
# # Create a new Word document
# doc <- read_docx()
# 
# # Loop through PNG files and add each to the document
# for (img in png_files) {
#   doc <- doc %>%
#     body_add_img(src = img, width = 6, height = 4) %>%
#     body_add_par("") %>% # spacer
#     body_add_par("")  # spacer
# }
# 
# # Save the document
# print(doc, target = paste(fig_path_ipcw, "sym_resolution_km_plots.docx"))








