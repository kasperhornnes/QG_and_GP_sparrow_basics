# =============================================================================
# 03_gblup.R
#
# Functions for temporal genomic prediction using GBLUP.
#
# This file:
#   - Defines the GBLUP model
#   - Creates temporal train/test splits
#   - Creates training-data priors
#   - Fits one temporal GBLUP scenario
#
# Evaluation is handled elsewhere.
# =============================================================================


make_gblup_effects <- function() {

  c(
    "1",

    # Fixed effects
    "sex",
    "month",
    "age",
    "f_roh",

    # Random environmental effects
    "f(hatch_year, model = \"iid\", hyper = prior$hyperpar_var)",
    "f(island, model = \"iid\", hyper = prior$hyperpar_var)",
    "f(id2, model = \"iid\", hyper = prior$hyperpar_var)",
    "f(day_session, model = \"iid\", hyper = prior$hyperpar_var)",

    # Genomic breeding value
    paste0(
      "f(id1, ",
      "values = as.numeric(colnames(inverse_relatedness_matrix)), ",
      "model = \"generic0\", ",
      "hyper = prior$hyperpar_var, ",
      "constr = FALSE, ",
      "Cmatrix = inverse_relatedness_matrix)"
    )
  )
}

#making temporal trainging and testing splits


make_temporal_prior <- function(
    pheno_data,
    train_rows
) {

  training_variance <- var(
    pheno_data$y[train_rows],
    na.rm = TRUE
  )

  make_prior(
    pc_prec_upper_var = training_variance / 3,
    var_init = training_variance / 5,
    tau = 0.05
  )
}

fit_temporal_gblup <- function(
    pheno_data,
    inverse_relatedness_matrix,
    train_years,
    verbose = TRUE,
    control.compute.config = FALSE
) {

  # ---------------------------------------------------------------------------
  # 1. Define training data
  # ---------------------------------------------------------------------------

  split <- make_training_split(
    pheno_data = pheno_data,
    train_years = train_years
  )


  # ---------------------------------------------------------------------------
  # 2. Prior based only on training phenotypes
  # ---------------------------------------------------------------------------

  prior <- make_temporal_prior(
    pheno_data = pheno_data,
    train_rows = split$train_rows
  )


  # ---------------------------------------------------------------------------
  # 3. Model effects
  # ---------------------------------------------------------------------------

  effects_vec <- make_gblup_effects()


  # ---------------------------------------------------------------------------
  # 4. Fit GBLUP
  # ---------------------------------------------------------------------------

  fit <- run_gp_training(
    pheno_data = pheno_data,
    train_rows = split$train_rows,
    inverse_relatedness_matrix = inverse_relatedness_matrix,
    effects_vec = effects_vec,
    prior = prior,
    verbose = verbose,
    control.compute.config = control.compute.config
  )


  # ---------------------------------------------------------------------------
  # 5. Store information about the experiment
  # ---------------------------------------------------------------------------

  fit$train_rows <- split$train_rows
  fit$train_years <- split$train_years
  fit$training_birds <- split$training_birds

  fit$n_training_obs <- split$n_training_obs
  fit$n_training_birds <- split$n_training_birds

  fit$last_training_year <- max(
    split$train_years
  )

  fit$method <- "GBLUP"


  fit
}
