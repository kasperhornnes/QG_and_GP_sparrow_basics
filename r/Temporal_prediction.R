# =============================================================================
#                            Setup
# =============================================================================
# ---- Packages ----
library(INLA)
library(dplyr)
library(data.table)
library(MCMCglmm)
library(HDInterval)

#load utility functions and separate r files
source("r/func.R")
source("r/01_prepare_data.R")
source("r/02_make_pcs.R")
source("r/02_make_grm.R")
source("r/03_bpcrr.R")
source("r/03_gblup.R")
source("r/04_evaluation.R")

#For reproducability
set.seed(2)


# ---- Subsetting the phenotype data ----
##### Just body mass measurements
response_colname <- "body_mass"
response <- "mass"
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
prepared <- prepare_data(
  pheno_file = pheno_file,
  orig_geno_files = orig_geno_files,
  islands = isls,
  sys_name = sys_name,
  froh_file = "data/FROH2.5_helgeland.txt",
  response_colname = response_colname,
  response = response,
  testing = testing_n
)

pheno_data <- prepared$pheno_data


# =============================================================================
#                      Create GRM for GBLUP and PCs for BPCRR
# =============================================================================

grm <- make_analysis_grm(
  prepared = prepared
)


pcs <- make_bpcrr_pcs(
  prepared = prepared,
  n_pcs = 1000
)


# =============================================================================
#                 Fitting the genomic animal model and the BPCRR
# =============================================================================


fit_gblup <- fit_temporal_gblup(
  pheno_data = pheno_data,
  inverse_relatedness_matrix = grm$inv_grm,
  cutoff_year = 2005,
)

fit_bpcrr <- fit_temporal_bpcrr(
  pheno_data = pheno_data,
  pcs = pcs,
  cutoff_year = 2005,
  n_pcs = 650
)



# =============================================================================
#                   Prediction and evaluation
# =============================================================================

eval_gblup <- evaluate_horizons(
  fit = fit_gblup,
  pheno_data = pheno_data,
  horizons = c(1, 2, 3, 5, 10),
  method = "GBLUP"
)

eval_bpcrr <- evaluate_horizons(
  fit = fit_bpcrr,
  pheno_data = pheno_data,
  horizons = c(1, 2, 3, 5, 10),
  method = "BPCRR"
)

evaluation_summary <- dplyr::bind_rows(
  eval_gblup$summary,
  eval_bpcrr$summary
)

evaluation_summary





# =============================================================================
#                   Comparing full fit vs train only
# =============================================================================


source("r/test_full_vs_train_only.R")

gblup_train <- fit_temporal_gblup(
  pheno_data = pheno_data,
  inverse_relatedness_matrix = grm$inv_grm,
  cutoff_year = 2005,
  verbose = TRUE
)

gblup_full <- fit_gblup_full_na(
  pheno_data = pheno_data,
  inverse_relatedness_matrix = grm$inv_grm,
  cutoff_year = 2005,
  verbose = TRUE
)


bpcrr_train <- fit_temporal_bpcrr(
  pheno_data = pheno_data,
  pcs = pcs,
  cutoff_year = 2005,
  n_pcs = 100,
  verbose = TRUE
)

bpcrr_full <- fit_bpcrr_full_na(
  pheno_data = pheno_data,
  pcs = pcs,
  cutoff_year = 2005,
  n_pcs = 100,
  verbose = TRUE
)

compare_bv_fits(
  gblup_full,
  gblup_train
)

compare_bpcrr_bv_fits(
  bpcrr_full,
  bpcrr_train,
  pcs
)

compare_target_bv(
  gblup_full,
  gblup_train,
  pheno_data,
  target_year = 2006
)

compare_target_bv(
  bpcrr_full,
  bpcrr_train,
  pheno_data,
  target_year = 2006
)

gblup_horizon_comparison <- compare_bv_horizons(
  fit_full = gblup_full,
  fit_train = gblup_train,
  pheno_data = pheno_data,
  cutoff_year = 2005
)

gblup_horizon_comparison
