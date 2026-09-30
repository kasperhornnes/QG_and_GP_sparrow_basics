# =============================================================================
#                            Setup
# =============================================================================
# ---- Packages ----
library(INLA)
library(dplyr)
library(data.table)
library(MCMCglmm)
library(HDInterval)
library(ggplot2)
library(tidyr)


#load utility functions and separate r files
source("r/func.R")
source("r/01_prepare_data.R")
source("r/02_make_pcs.R")
source("r/02_make_grm.R")
source("r/03_bpcrr.R")
source("r/03_gblup.R")
source("r/04_evaluation.R")
source("r/05_plotting.R")

#For reproducability
set.seed(2)


# =============================================================================
#                            Analysis settings
# =============================================================================



# ---- Subsetting the phenotype data ----

# response_colname <- "body_mass"
# response <- "mass"

# response_colname <- "thr_tarsus"
# response <- "tarsus"

response_colname <- "thr_wing"
response <- "wing"

####### Helgeland islands:
isls <- c(20, 22, 23, 24, 26, 27, 28, 33, 331, 332, 34, 35, 38)
sys_name <- "helgeland"
####### Southern islands:
# isls <- c(60, 61, 63, 67, 68),
# sys_name <- "southern"
####### Helgeland and southern together:
# isls <- c(c(20, 22, 23, 24, 26, 27, 28, 33, 331, 332, 34, 35, 38),
#           c(60, 61, 63, 67, 68))
# sys_name <- "all"

# Development subset, NULL to have all observations
testing_n <- NULL

#number of pcs in the BPCRR model, and number of pcs available from the PCA
n_pcs_available <- 1000
n_pcs_model <- 800

#prior for the PCs in the BPCRR model, either "inla_default" or "fixed"
pc_prior <- "inla_default"
varA_prior = NULL #(bodymass: 5.2 * 0.3)

#cutoff year for training, and the years to evaluate after the cutoff
cutoff_year <- 2005
target_years <- 2006:2020

# =============================================================================
#                            Cache / results
# =============================================================================

cache_dir <- "results/cache"
model_dir <- "results/models"

dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)

# Set TRUE when you deliberately want to rerun everything
force_refit <- FALSE



testing_tag <- if (is.null(testing_n)) {
  "all"
} else {
  paste0("n", testing_n)
}


# ---- File paths ----

pheno_file <- "data/AdultMorphology_20241117.csv"

orig_geno_files <- paste0(
  "data/combined_200k_70k_sparrow_genotype_data/",
  "combined_200k_70k_helgeland_south_corrected_snpfiltered_2024-02-05",
  c(".map", ".ped", ".fam", ".bim", ".bed")
)


# =============================================================================
#                            Prepare data
# =============================================================================
prepared_file <- file.path(
  cache_dir,
  paste0("prepared_", sys_name, "_", response, "_", testing_tag, ".rds")
)

prepared <- load_or_run(
  file = prepared_file,
  force = force_refit,
  fun = function() {

    prepare_data(
      pheno_file = pheno_file,
      orig_geno_files = orig_geno_files,
      islands = isls,
      sys_name = sys_name,
      froh_file = "data/FROH2.5_helgeland.txt",
      response_colname = response_colname,
      response = response,
      testing = testing_n
    )
  }
)

pheno_data <- prepared$pheno_data

# =============================================================================
#                      Create GRM for GBLUP and PCs for BPCRR
# =============================================================================



pcs_file <- file.path(
  cache_dir,
  paste0("pcs_", sys_name, "_", response, "_", testing_tag, "_nPCs_max=", n_pcs_available, ".rds")
)

pcs <- load_or_run(
  file = pcs_file,
  force = force_refit,
  fun = function() {

    make_bpcrr_pcs(
      prepared = prepared,
      n_pcs = n_pcs_available,
      out_dir = file.path(
        "data",
        "bpcrr_pca",
        paste0(sys_name, "_", response)
      )
    )
  }
)

grm_file <- file.path(
  cache_dir,
  paste0("grm_", sys_name, "_", response, "_", testing_tag, "LDpruned.rds")
)

grm <- load_or_run(
  file = grm_file,
  force = force_refit,
  fun = function() {

    make_analysis_grm(
      prepared = prepared,
      bfile = pcs$pruned_bfile,
    )
  }
)


# =============================================================================
#         Fitting the genomic animal model and the BPCRR, train before 2005
# =============================================================================



base_train_years <- sort(
  unique(
    pheno_data$year_num[
      pheno_data$year_num <= cutoff_year
    ]
  )
)

gblup_file <- file.path(
  model_dir,
  paste0("gblup_", sys_name, "_", response,
         "_", testing_tag, "_cutoff_", cutoff_year
         , "LDpruned.rds")
)

fit_gblup <- load_or_run(
  file = gblup_file,
  force = force_refit,
  fun = function() {

    fit_gblup <- fit_temporal_gblup(
      pheno_data = pheno_data,
      inverse_relatedness_matrix = grm$inv_grm,
      train_years = base_train_years

    )
  }
)




bpcrr_file <- file.path(
  model_dir,
  paste0(
    "bpcrr_", sys_name, "_", response,
    "_", testing_tag, "_cutoff_", cutoff_year,
    "_prior_", pc_prior,
    "_nPCs=", n_pcs_model, ".rds"))

fit_bpcrr <- load_or_run(
  file = bpcrr_file,
  force = force_refit,
  fun = function() {

    fit_bpcrr <- fit_temporal_bpcrr(
      pheno_data = pheno_data,
      pcs = pcs,
      train_years = base_train_years,
      n_pcs = n_pcs_model,
      pc_prior = pc_prior,
      varA_prior = varA_prior
    )
  }
)
















# =============================================================================
#              Evaluate all available years after cutoff
# =============================================================================


eval_gblup <- evaluate_target_years(
  fit = fit_gblup,
  pheno_data = pheno_data,
  target_years = target_years,
  method = "GBLUP"
)

eval_bpcrr <- evaluate_target_years(
  fit = fit_bpcrr,
  pheno_data = pheno_data,
  target_years = target_years,
  method = "BPCRR"
)

evaluation_all_years <- dplyr::bind_rows(
  eval_gblup$summary,
  eval_bpcrr$summary
)


# =============================================================================
#             Calculating relatedness between training and test sets
# =============================================================================

relatedness_by_year <- calculate_temporal_relatedness(
  pheno_data = pheno_data,
  grm = grm$grm,
  cutoff_year = cutoff_year
)

relatedness_by_year


# =============================================================================
#                       Plotting
# =============================================================================

plot_unseen <- plot_unseen_accuracy(
  evaluation_data = evaluation_all_years,
  min_year = min(target_years),
  max_year = max(target_years)
)

plot_unseen

plot_relatedness <- plot_relatedness_over_time(
  relatedness_by_year,
  min_year = min(target_years),
  max_year = max(target_years)
)

plot_relatedness

plot_rel_hist <- plot_relatedness_histogram_compare(
  pheno_data = pheno_data,
  grm = grm$grm,
  train_years = base_train_years,
  compare_years = c(2006, 2016),
  unseen_only = TRUE,
  bins = 200
)

plot_rel_hist




#testing the corralation between running the model with the test phenotypes as NA, versus removing test rows
# =============================================================================
#                   Comparing full fit vs train only
# =============================================================================


#
#
#
# source("r/test_full_vs_train_only.R")
#
# gblup_train <- fit_temporal_gblup(
#   pheno_data = pheno_data,
#   inverse_relatedness_matrix = grm$inv_grm,
#   cutoff_year = 2005,
#   verbose = TRUE
# )
#
# gblup_full <- fit_gblup_full_na(
#   pheno_data = pheno_data,
#   inverse_relatedness_matrix = grm$inv_grm,
#   cutoff_year = 2005,
#   verbose = TRUE
# )
#
#
# bpcrr_train <- fit_temporal_bpcrr(
#   pheno_data = pheno_data,
#   pcs = pcs,
#   cutoff_year = 2005,
#   n_pcs = 100,
#   verbose = TRUE
# )
#
# bpcrr_full <- fit_bpcrr_full_na(
#   pheno_data = pheno_data,
#   pcs = pcs,
#   cutoff_year = 2005,
#   n_pcs = 100,
#   verbose = TRUE
# )
#
# compare_bv_fits(
#   gblup_full,
#   gblup_train
# )
#
# compare_bpcrr_bv_fits(
#   bpcrr_full,
#   bpcrr_train,
#   pcs
# )
#
# compare_target_bv(
#   gblup_full,
#   gblup_train,
#   pheno_data,
#   target_year = 2006
# )
#
# compare_target_bv(
#   bpcrr_full,
#   bpcrr_train,
#   pheno_data,
#   target_year = 2006
# )
#
# gblup_horizon_comparison <- compare_bv_horizons(
#   fit_full = gblup_full,
#   fit_train = gblup_train,
#   pheno_data = pheno_data,
#   cutoff_year = 2005
# )
#
# gblup_horizon_comparison
