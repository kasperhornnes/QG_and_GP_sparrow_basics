# =============================================================================
# 03_bpcrr.R
#
# Functions for temporal genomic prediction using BPCRR.
#
# This file:
#   - Defines the BPCRR model
#   - Creates temporal train/test splits
#   - Creates training-data priors
#   - Fits one temporal BPCRR scenario
#
# Evaluation is handled elsewhere.
# =============================================================================


make_bpcrr_effects <- function() {

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

    # Genomic breeding value using PCs
    paste0(
      "f(id1, ",
      "model = \"z\", ",
      "Z = Z_pc, ",
      "hyper = list(",
      "prec = list(",
      "initial = log(1 / u_prior_var), ",
      "fixed = TRUE",
      ")",
      "))"
    )
  )
}


# =============================================================================
# Temporal training split
# =============================================================================

make_temporal_split_bpcrr <- function(
    pheno_data,
    cutoff_year
) {

  train_rows <- pheno_data$year_num <= cutoff_year

  if (sum(train_rows) == 0) {
    stop(
      "No training observations available for cutoff year ",
      cutoff_year
    )
  }

  list(
    cutoff_year = cutoff_year,
    train_rows = train_rows
  )
}


# =============================================================================
# Prior based on training phenotypes
# =============================================================================

make_temporal_prior_bpcrr <- function(
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


# =============================================================================
# Fit BPCRR for one cutoff year
# =============================================================================

fit_temporal_bpcrr <- function(
    pheno_data,
    pcs,
    cutoff_year,
    n_pcs = 1000,
    varA_prior = 5.2 * 0.3,
    verbose = TRUE,
    control.compute.config = FALSE,
    num.threads = 8
) {

  # ---------------------------------------------------------------------------
  # 1. Temporal split
  # ---------------------------------------------------------------------------

  split <- make_temporal_split_bpcrr(
    pheno_data = pheno_data,
    cutoff_year = cutoff_year
  )


  # ---------------------------------------------------------------------------
  # 2. Prior based only on training data
  # ---------------------------------------------------------------------------

  prior <- make_temporal_prior_bpcrr(
    pheno_data = pheno_data,
    train_rows = split$train_rows
  )


  # ---------------------------------------------------------------------------
  # 3. Select PCs
  # ---------------------------------------------------------------------------

  if (n_pcs > ncol(pcs$Z)) {
    stop(
      "Requested ", n_pcs,
      " PCs, but only ",
      ncol(pcs$Z),
      " are available."
    )
  }

  Z_pc <- pcs$Z[
    ,
    seq_len(n_pcs),
    drop = FALSE
  ]


  # ---------------------------------------------------------------------------
  # 4. Fixed ridge prior for PC effects
  # ---------------------------------------------------------------------------

  pc_variances <- apply(
    Z_pc,
    2,
    var
  )

  total_pc_variance <- sum(
    pc_variances
  )

  u_prior_var <-
    varA_prior / total_pc_variance


  # ---------------------------------------------------------------------------
  # 5. Model effects
  # ---------------------------------------------------------------------------

  effects_vec <- make_bpcrr_effects()


  # ---------------------------------------------------------------------------
  # 6. Training phenotype data only
  # ---------------------------------------------------------------------------

  data_model <- pheno_data[
    split$train_rows,
    ,
    drop = FALSE
  ]

  data_model$y_na <- data_model$y


  # ---------------------------------------------------------------------------
  # 7. Fit BPCRR
  # ---------------------------------------------------------------------------

  inla_formula <- stats::reformulate(
    effects_vec,
    response = "y_na"
  )

  model <- INLA::inla(
    inla_formula,
    family = "gaussian",
    data = data_model,
    verbose = verbose,

    control.compute = list(
      config = control.compute.config
    ),

    control.family = list(
      hyper = prior$hyperpar_var
    ),

    num.threads = num.threads
  )


  # ---------------------------------------------------------------------------
  # 8. Store information about the experiment
  # ---------------------------------------------------------------------------

  fit <- list(
    model = model,
    data = data_model,
    cutoff_year = cutoff_year,
    n_pcs = n_pcs,
    varA_prior = varA_prior,
    u_prior_var = u_prior_var,
    Z_pc = Z_pc,
    method = "BPCRR"
  )

  fit
}
