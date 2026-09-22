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

make_temporal_split <- function(
    pheno_data,
    cutoff_year
) {

  train_rows <- pheno_data$year_num <= cutoff_year

  # Basic checks
  if (sum(train_rows) == 0) {
    stop("No training observations available for cutoff year ", cutoff_year)
  }

  list(
    cutoff_year = cutoff_year,
    train_rows = train_rows
  )
}

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
    cutoff_year,
    verbose = TRUE,
    control.compute.config = FALSE
) {

  # ---------------------------------------------------------------------------
  # 1. Temporal split
  # ---------------------------------------------------------------------------

  split <- make_temporal_split(
    pheno_data = pheno_data,
    cutoff_year = cutoff_year
  )


  # ---------------------------------------------------------------------------
  # 2. Prior based only on training data
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

  fit <- run_gp_cutoff(
    pheno_data = pheno_data,
    cutoff_year = cutoff_year,
    inverse_relatedness_matrix = inverse_relatedness_matrix,
    effects_vec = effects_vec,
    prior = prior,
    verbose = verbose,
    control.compute.config = control.compute.config
  )


  # ---------------------------------------------------------------------------
  # 5. Store information about the experiment
  # ---------------------------------------------------------------------------

  fit$cutoff_year <- cutoff_year

  fit
}

