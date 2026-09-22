
#testing "full" model vs train only model

fit_gblup_full_na <- function(
    pheno_data,
    inverse_relatedness_matrix,
    cutoff_year,
    verbose = TRUE
) {

  split <- make_temporal_split(
    pheno_data,
    cutoff_year
  )

  prior <- make_temporal_prior(
    pheno_data,
    split$train_rows
  )

  effects_vec <- make_gblup_effects()

  # Keep ALL rows
  data_model <- pheno_data

  # But hide future phenotypes
  data_model$y_na <- data_model$y
  data_model$y_na[!split$train_rows] <- NA

  formula <- stats::reformulate(
    effects_vec,
    response = "y_na"
  )

  model <- INLA::inla(
    formula,
    family = "gaussian",
    data = data_model,
    verbose = verbose,
    control.compute = list(config = FALSE),
    control.family = list(
      hyper = prior$hyperpar_var
    )
  )

  list(
    model = model,
    data = data_model,
    cutoff_year = cutoff_year
  )
}

fit_bpcrr_full_na <- function(
    pheno_data,
    pcs,
    cutoff_year,
    n_pcs = 100,
    varA_prior = 5.2 * 0.3,
    verbose = TRUE,
    num.threads = 8
) {

  split <- make_temporal_split_bpcrr(
    pheno_data,
    cutoff_year
  )

  prior <- make_temporal_prior_bpcrr(
    pheno_data,
    split$train_rows
  )

  # Same PCs as train-only model
  Z_pc <- pcs$Z[
    ,
    seq_len(n_pcs),
    drop = FALSE
  ]

  pc_variances <- apply(
    Z_pc,
    2,
    var
  )

  u_prior_var <-
    varA_prior / sum(pc_variances)

  effects_vec <- make_bpcrr_effects()

  # Keep all rows
  data_model <- pheno_data

  # Hide future phenotypes
  data_model$y_na <- data_model$y
  data_model$y_na[!split$train_rows] <- NA

  formula <- stats::reformulate(
    effects_vec,
    response = "y_na"
  )

  model <- INLA::inla(
    formula,
    family = "gaussian",
    data = data_model,
    verbose = verbose,
    control.compute = list(config = FALSE),
    control.family = list(
      hyper = prior$hyperpar_var
    ),
    num.threads = num.threads
  )

  list(
    model = model,
    data = data_model,
    cutoff_year = cutoff_year,
    n_pcs = n_pcs
  )
}

compare_bv_fits <- function(
    fit_full,
    fit_train
) {

  bv_full <- fit_full$model$summary.random$id1
  bv_train <- fit_train$model$summary.random$id1

  comparison <- merge(
    bv_full[, c("ID", "mean")],
    bv_train[, c("ID", "mean")],
    by = "ID",
    suffixes = c("_full", "_train")
  )

  diff <- comparison$mean_full -
    comparison$mean_train

  data.frame(
    n = nrow(comparison),

    correlation = cor(
      comparison$mean_full,
      comparison$mean_train
    ),

    mean_abs_difference =
      mean(abs(diff)),

    max_abs_difference =
      max(abs(diff))
  )
}

compare_target_bv <- function(
    fit_full,
    fit_train,
    pheno_data,
    target_year
) {

  bv_full <- fit_full$model$summary.random$id1
  bv_train <- fit_train$model$summary.random$id1

  comparison <- merge(
    bv_full[, c("ID", "mean")],
    bv_train[, c("ID", "mean")],
    by = "ID",
    suffixes = c("_full", "_train")
  )

  target_ids <- unique(
    pheno_data$id1[
      pheno_data$year_num == target_year
    ]
  )

  comparison <- comparison[
    comparison$ID %in% target_ids,
  ]

  diff <- comparison$mean_full -
    comparison$mean_train

  data.frame(
    target_year = target_year,
    n = nrow(comparison),

    correlation = cor(
      comparison$mean_full,
      comparison$mean_train
    ),

    mean_abs_difference =
      mean(abs(diff)),

    max_abs_difference =
      max(abs(diff))
  )
}

compare_bpcrr_bv_fits <- function(
    fit_full,
    fit_train,
    pcs
) {

  n_animals <- nrow(pcs$Z)

  bv_full <- fit_full$model$summary.random$id1[
    seq_len(n_animals),
  ]

  bv_train <- fit_train$model$summary.random$id1[
    seq_len(n_animals),
  ]

  comparison <- data.frame(
    id1 = seq_len(n_animals),
    mean_full = bv_full$mean,
    mean_train = bv_train$mean
  )

  diff <- comparison$mean_full -
    comparison$mean_train

  data.frame(
    n = n_animals,

    correlation = cor(
      comparison$mean_full,
      comparison$mean_train
    ),

    mean_abs_difference =
      mean(abs(diff)),

    max_abs_difference =
      max(abs(diff))
  )
}

compare_bv_horizons <- function(
    fit_full,
    fit_train,
    pheno_data,
    cutoff_year,
    horizons = c(1, 2, 3, 5, 10)
) {

  bv_full <- fit_full$model$summary.random$id1
  bv_train <- fit_train$model$summary.random$id1

  results <- lapply(horizons, function(h) {

    target_year <- cutoff_year + h

    target_ids <- unique(
      pheno_data$id1[
        pheno_data$year_num == target_year
      ]
    )

    # Only IDs available in both fits
    target_ids <- target_ids[
      target_ids <= nrow(bv_full) &
        target_ids <= nrow(bv_train)
    ]

    mean_diff <-
      bv_full$mean[target_ids] -
      bv_train$mean[target_ids]

    sd_diff <-
      bv_full$sd[target_ids] -
      bv_train$sd[target_ids]

    data.frame(
      horizon = h,
      target_year = target_year,
      n = length(target_ids),

      # Breeding-value estimates
      bv_correlation = cor(
        bv_full$mean[target_ids],
        bv_train$mean[target_ids]
      ),

      bv_mean_abs_diff = mean(
        abs(mean_diff)
      ),

      # Breeding-value uncertainty
      sd_correlation = cor(
        bv_full$sd[target_ids],
        bv_train$sd[target_ids]
      ),

      sd_mean_abs_diff = mean(
        abs(sd_diff)
      )
    )
  })

  do.call(
    rbind,
    results
  )
}
