# =============================================================================
# Introductory script with basics for house sparrow quant. gen.
# Author: Kenneth Aase
# Expand on or change the code for your own workflow/analysis,
# this is just meant to be a simple and instructive starting point.
# ====================================s=========================================
# =============================================================================
#                            Part 0: setup
# =============================================================================
# ---- Packages necesarry to run this script ----
library(INLA) # Cannot be installed from CRAN, for install instructions for your
# system see: https://www.r-inla.org/download-install
library(dplyr)
library(data.table)
library(MCMCglmm)
library(HDInterval)
# ---- Data paths ----
# File path to morphological phenotype data
pheno_file <- "data/AdultMorphology_20241117.csv"
# Check if file exist
stopifnot(file.exists(pheno_file))

# File paths to genomic data
orig_geno_files <-
  paste0("data/combined_200k_70k_sparrow_genotype_data/",
         "combined_200k_70k_helgeland_south_corrected_snpfiltered_2024-02-05",
         c(".map", ".ped", ".fam", ".bim", ".bed"))
# Check if files exist
for (file in orig_geno_files) {
  stopifnot(file.exists(file))
}
# ---- File path to plink program here ----
# The script relies on external calls to the very useful PLINK 1.9 program for
# genomic data handling. Download the version matching your system from here:
# https://www.cog-genomics.org/plink/1.9/
# And make a path to the program file:
get_plink_path <- function() {
  "PLINK/plink"
}
# Check if program file exists
stopifnot(file.exists(get_plink_path()))

# ----  Load some utility functions ----
# Some details are deliberately hidden away inside these functions, so if you
# want to make changes somewhere, it would be wise to check inside these
source("r/func.R")
set.seed(2) # Set seed for reproducibility
# =============================================================================
#                     Part 1: Quality Control (QC)
# =============================================================================

# Vector of ID codes (ring numbers) for genotyped individuals
genotyped_inds <- get_genotyped_inds(fam_file = orig_geno_files[3], sel = 1)

# Basic quality control of the genomic data
# Remove SNPs and individuals with too much mssing data,
# and SNPs with too little genetic variation (maf)
qc_filters <- list(genorate_ind = 0.05,
                   genorate_snp = 0.1,
                   maf = 0.01)
# See inside the do_qc function for more details.
# For more filtering options, see https://www.cog-genomics.org/plink/1.9/filter
qc_overall <- do_qc(fam_file = orig_geno_files[3],
                    ncores = 8,
                    mem = 8 * 6000,
                    qc_filt = qc_filters,
                    keep_inds = genotyped_inds,
                    sys = "",
                    resp = "overall")
genotyped_inds_qc <- get_genotyped_inds(fam_file = qc_overall[2], sel = 1)

# =============================================================================
#                     Part 2: Genomic animal model
# =============================================================================
# ---- Here is a simple genomic animal model for body mass fitted in INLA ----

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

# ---- Create data frame of chosen phenotype data ----
pheno_data <- pheno_wrangle(
  filepath = pheno_file,
  genotyped_inds = genotyped_inds_qc,
  islands = isls,
  y_col_name = response_colname,
  # Use this to include only e.g. 250 observations, for fast testing of code
  # (can cause numerical problems if set too low):
  testing = 250
  # Use this to include all observations:
  # testing = NULL
)
# Inspect:
head(pheno_data)

# ---- Creating the Genomic Relatedness Matrix (GRM) ----
# Make new genomic data files that include only relevant phenotyped individuals
geno_files <- do_qc(fam_file = qc_overall[2],
                    ncores = 8,
                    mem = 8 * 6000,
                    qc_filt = qc_filters,
                    keep_inds = unique(pheno_data$ringnr),
                    sys = sys_name,
                    resp = response)
# Create files for GRM
grm_files <- make_raw_grm(analysis_inds = unique(pheno_data$ringnr),
                          bfile = gsub(".{4}$", "", geno_files[2]),
                          frq_file = geno_files[5],
                          genorate_ind = qc_filters$genorate_ind,
                          genorate_snp = qc_filters$genorate_snp,
                          ncores = 8,
                          mem = 8 * 6000,
                          maf = qc_filters$maf,
                          response = response,
                          geno_set = paste0(sys_name, "_70K"))
# Load data files GRM, and create its inverse
grm_obj <- compute_grm_obj(frq_file = geno_files[5],
                           rel_file = grm_files[2],
                           id_file = grm_files[1],
                           bim_file = grm_files[3],
                           pheno_data = pheno_data,
                           id_col = 1)
# Inspect:
grm_obj$grm[1:6, 1:6]
grm_obj$inv_grm[1:6, 1:6]

# ---- Vector of terms we want to include in our INLA model ----
gam_effect_strs <-
  c("1", # Intercept
    "sex", # Fixed effect of sex
    "month", # Fixed effect of month
    "age", # Fixed effect of age
    # Random effect accounting for hatch year:
    "f(hatch_year, model = \"iid\", hyper = prior$hyperpar_var)",
    # Random effect accounting for island
    "f(island, model = \"iid\", hyper = prior$hyperpar_var)",
    # Random effect accounting for identity
    "f(id2, model = \"iid\", hyper = prior$hyperpar_var)",
    # Random effect accounting for measurement session
    "f(day_session, model = \"iid\", hyper = prior$hyperpar_var)",
    # Random effect for breeding values:
    paste0("f(id1, values = as.numeric(colnames(inverse_relatedness_matrix)),",
           "model = \"generic0\", hyper = prior$hyperpar_var, constr = FALSE,",
           " Cmatrix = inverse_relatedness_matrix)"))
# For models with Gaussian likelihood, INLA automatically includes a residual

# ---- Priors ----
# Specify priors for random effect variances. We use a Penalized Complexity (PC)
# prior.
prior <- (pheno_data %>%
            getElement("y") %>%
            var() %>%
            make_prior(pc_prec_upper_var = . / 3,
                       var_init = . / 5,
                       tau = 0.05))
# By default INLA uses N(0, 10^4) priors for fixed effect coefficients

# ---- Fit the Genomic Animal Model (GAM) ----
gam <- run_gp(pheno_data = pheno_data,
              train_inds = unique(pheno_data$ringnr), # Use all inds to fit GAM
              test_inds = character(0), # No test set in the GAM
              inverse_relatedness_matrix = grm_obj$inv_grm,
              effects_vec = gam_effect_strs,
              prior = prior,
              # Output detailed information from INLA while fitting or not:
              verbose = TRUE,
              # This option drastically increases size of model object, but is
              # necessary to be able to generate joint posterior samples from
              # the latent field:
              control.compute.config = TRUE)

# Posterior statistics for fixed effects:
gam$summary.fixed
# INLA works with precision (1/variance) rather than variance.
# Posterior statistics for precisions:
gam$summary.hyperpar

# To find posterior statistics for variance components we need to do a transform
hyp_vars <- do.call(rbind,
                    lapply(gam$marginals.hyperpar,
                           inla_posterior_variances))
row.names(hyp_vars) <- gsub(pattern = "Precision",
                            replacement = "Variance",
                            x = row.names(hyp_vars))
hyp_vars
# Note: "Variance for the Gaussian observations" refers to the residual term
# =============================================================================
#             Part 3: Genomic prediction with Bayesian G-BLUP
# =============================================================================
# ---- We predict breeding values for genotyped but non-phenotyped inds ----

# Create folds for a cross validation (CV)
cv_test_sets <- make_cv_test_sets(analysis_inds = unique(pheno_data$ringnr),
                                  num_folds = 10)

# New prior (not using information from test set)
prior_cv <- (pheno_data %>%
               dplyr::filter(!ringnr %in% cv_test_sets[[1]]) %>%
               getElement("y") %>%
               var() %>%
               make_prior(pc_prec_upper_var = . / 3,
                          var_init = . / 5,
                          tau = 0.05))

# Run a genomic prediction for one of the folds in the CV. You can make a
# for-loop do the full CV, but it will be very slow if using the full data set.
# I recommend using a computing cluster to do it in parallel.
gpgrm_cv <- run_gp(pheno_data = pheno_data,
                   train_inds = setdiff(unique(pheno_data$ringnr),
                                        cv_test_sets[[1]]),
                   test_inds = cv_test_sets[[1]],
                   inverse_relatedness_matrix = grm_obj$inv_grm,
                   effects_vec = gam_effect_strs,
                   prior = prior_cv,
                   verbose = TRUE,
                   control.compute.config = TRUE)

# ---- # Extracting results from the INLA object ----
n <- length(unique(pheno_data$ringnr))
# Breeding value posterior statistics
bv_stats <- gpgrm_cv$summary.random$id1
# Predicted breeding values (posterior means)
pred_bv <- bv_stats$mean[order(bv_stats$ID)][1:n]
# Posterior standard deviations
pred_bv_sd <- bv_stats$sd[order(bv_stats$ID)][1:n]
# Breeding values repeated to match repeated measurements
pheno_data$pred_bv_rep <- pred_bv[pheno_data$id1]
# Predicted phenotypes (posterior means)
pred_pheno <- gpgrm_cv$summary.fitted.values$mean
pred_pheno_sd <- gpgrm_cv$summary.fitted.values$sd

# =============================================================================
#                     Part 4: Posterior sampling
# =============================================================================
# ---- Generate joint posterior samples from INLA model ----
samp <- INLA::inla.posterior.sample(n = 1e4, # number of samples
                                    result = gam, # INLA model
                                    add.names = FALSE)
# Example:
# ---- Computing VA in a cohort, e.g. sparrows born in 1997 on Hestmannøy ----
pheno_data %>%
  dplyr::filter(hatch_year == 2004 & first_island == 27) %>%
  getElement("id1") %>%
  unique() %>%
  as.numeric() ->
  rows
cohort_va_samps <- INLA::inla.posterior.sample.eval(
  samples = samp,
  rows = rows,
  fun = function(..., rows) {
    id1 %>% # Breeding values
      `[`(rows) %>% # Subset
      var() # Compute VA
  })
# Posterior statistics for cohort VA
extract_stats(cohort_va_samps[1, ])
