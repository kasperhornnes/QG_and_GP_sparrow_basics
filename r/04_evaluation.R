# =============================================================================
# 04_evaluation.R
#
# Evaluation functions for temporal genomic prediction.
#
# Works for:
#   - GBLUP
#   - BPCRR
#
# Prediction accuracy is evaluated as the correlation between:
#
#   mean observed phenotype per bird in the target year
#
# and
#
#   predicted breeding value
#
# =============================================================================


# =============================================================================
# Extract breeding values
# =============================================================================

extract_breeding_values <- function(
    fit,
    method = NULL
) {

  # Try to identify model automatically
  if (is.null(method)) {

    if (!is.null(fit$method)) {

      method <- toupper(fit$method)

    } else if (!is.null(fit$Z_pc)) {

      method <- "BPCRR"

    } else {

      method <- "GBLUP"
    }
  }


  # ---------------------------------------------------------------------------
  # GBLUP
  # ---------------------------------------------------------------------------

  if (method == "GBLUP") {

    bv <- fit$model$summary.random$id1

    return(
      data.frame(
        id1 = as.numeric(bv$ID),
        bv = bv$mean,
        bv_sd = bv$sd
      )
    )
  }


  # ---------------------------------------------------------------------------
  # BPCRR
  # ---------------------------------------------------------------------------

  if (method == "BPCRR") {

    bv <- fit$model$summary.random$id1

    # Important:
    # model = "z" contains:
    #
    # first N rows  = individual breeding values
    # final K rows  = PC effects

    n_animals <- nrow(
      fit$Z_pc
    )

    bv <- bv[
      seq_len(n_animals),
      ,
      drop = FALSE
    ]

    return(
      data.frame(
        id1 = seq_len(n_animals),
        bv = bv$mean,
        bv_sd = bv$sd
      )
    )
  }


  stop(
    "method must be GBLUP or BPCRR"
  )
}


# =============================================================================
# Evaluate one temporal horizon
# =============================================================================

#help function to calculate correlation while handling missing values and small sample sizes
safe_cor <- function(x, y) {

  keep <- complete.cases(x, y)

  # Need at least 2 observations to calculate a correlation
  if (sum(keep) < 2) {
    return(NA_real_)
  }

  cor(
    x[keep],
    y[keep]
  )
}

evaluate_temporal_model <- function(
    fit,
    pheno_data,
    target_year,
    method = NULL
) {

  # ---------------------------------------------------------------------------
  # 1. Determine model type
  # ---------------------------------------------------------------------------

  if (is.null(method)) {

    if (!is.null(fit$method)) {

      method <- toupper(
        fit$method
      )

    } else if (!is.null(fit$Z_pc)) {

      method <- "BPCRR"

    } else {

      method <- "GBLUP"
    }
  }


  # ---------------------------------------------------------------------------
  # 2. Check temporal ordering
  # ---------------------------------------------------------------------------

  if (is.null(fit$train_years)) {
    stop(
      "Fit does not contain train_years. ",
      "The model may have been fitted with the old code."
    )
  }

  if (any(fit$train_years >= target_year)) {
    stop(
      "For prospective temporal prediction, all training years ",
      "must be earlier than the target year."
    )
  }


  # ---------------------------------------------------------------------------
  # 3. Birds actually used for training
  # ---------------------------------------------------------------------------

  training_birds <- fit$training_birds


  # ---------------------------------------------------------------------------
  # 4. Target-year phenotype
  # ---------------------------------------------------------------------------

  target_data <-
    pheno_data |>
    dplyr::filter(
      year_num == target_year
    ) |>
    dplyr::group_by(
      ringnr,
      id1
    ) |>
    dplyr::summarise(
      observed = mean(
        y,
        na.rm = TRUE
      ),

      n_observations = dplyr::n(),

      .groups = "drop"
    )


  if (nrow(target_data) == 0) {

    stop(
      "No phenotype observations for target year ",
      target_year
    )
  }


  # ---------------------------------------------------------------------------
  # 5. Extract predicted breeding values
  # ---------------------------------------------------------------------------

  bv <- extract_breeding_values(
    fit = fit,
    method = method
  )


  # ---------------------------------------------------------------------------
  # 6. Match phenotype and breeding value
  # ---------------------------------------------------------------------------

  predictions <-
    target_data |>
    dplyr::left_join(
      bv,
      by = "id1"
    ) |>
    dplyr::mutate(

      seen_before =
        ringnr %in% training_birds,

      method = method,

      target_year = target_year,

      last_training_year =
        fit$last_training_year,

      years_since_last_training =
        target_year -
        fit$last_training_year
    )


  if (anyNA(predictions$bv)) {

    warning(
      sum(is.na(predictions$bv)),
      " target birds have no predicted breeding value."
    )
  }


  # ---------------------------------------------------------------------------
  # 7. Prediction accuracy
  # ---------------------------------------------------------------------------

  accuracy_overall <- safe_cor(
    predictions$observed,
    predictions$bv
  )


  accuracy_seen <- safe_cor(
    predictions$observed[
      predictions$seen_before
    ],
    predictions$bv[
      predictions$seen_before
    ]
  )


  accuracy_unseen <- safe_cor(
    predictions$observed[
      !predictions$seen_before
    ],
    predictions$bv[
      !predictions$seen_before
    ]
  )


  # ---------------------------------------------------------------------------
  # 8. Summary
  # ---------------------------------------------------------------------------

  summary <- data.frame(

    method = method,

    target_year = target_year,

    last_training_year =
      fit$last_training_year,

    years_since_last_training =
      target_year -
      fit$last_training_year,

    n_training_obs =
      fit$n_training_obs,

    n_training_birds =
      fit$n_training_birds,

    n_target =
      nrow(predictions),

    n_seen = sum(
      predictions$seen_before
    ),

    n_unseen = sum(
      !predictions$seen_before
    ),

    accuracy_overall =
      accuracy_overall,

    accuracy_seen =
      accuracy_seen,

    accuracy_unseen =
      accuracy_unseen,

    mean_bv_sd = mean(
      predictions$bv_sd,
      na.rm = TRUE
    )
  )


  list(
    summary = summary,
    predictions = predictions
  )
}

# =============================================================================
# Evaluate several horizons
# =============================================================================
evaluate_target_years <- function(
    fit,
    pheno_data,
    target_years,
    method = NULL
) {

  evaluations <- lapply(
    target_years,
    function(year) {

      evaluate_temporal_model(
        fit = fit,
        pheno_data = pheno_data,
        target_year = year,
        method = method
      )
    }
  )


  summary <- dplyr::bind_rows(
    lapply(
      evaluations,
      function(x) x$summary
    )
  )


  predictions <- dplyr::bind_rows(
    lapply(
      evaluations,
      function(x) x$predictions
    )
  )


  list(
    summary = summary,
    predictions = predictions
  )
}



calculate_temporal_relatedness <- function(
    pheno_data,
    grm,
    cutoff_year
) {

  # ---------------------------------------------------------------------------
  # Training individuals
  # ---------------------------------------------------------------------------

  training_ids <- unique(
    pheno_data$id1[
      pheno_data$year_num <= cutoff_year
    ]
  )

  training_birds <- unique(
    pheno_data$ringnr[
      pheno_data$year_num <= cutoff_year
    ]
  )


  # ---------------------------------------------------------------------------
  # Future years
  # ---------------------------------------------------------------------------

  future_years <- sort(
    unique(
      pheno_data$year_num[
        pheno_data$year_num > cutoff_year
      ]
    )
  )


  # ---------------------------------------------------------------------------
  # Calculate train-target relatedness for each year
  # ---------------------------------------------------------------------------

  results <- lapply(
    future_years,
    function(target_year) {

      target_data <- pheno_data |>
        dplyr::filter(
          year_num == target_year,
          !(ringnr %in% training_birds)
        ) |>
        dplyr::distinct(
          ringnr,
          id1
        )

      target_ids <- target_data$id1


      # Need target birds
      if (length(target_ids) == 0) {

        return(
          data.frame(
            target_year = target_year,
            n_target = 0,
            mean_relatedness = NA_real_,
            sd_relatedness = NA_real_,
            var_relatedness = NA_real_,
            precision_relatedness = NA_real_
          )
        )
      }


      # GRM row/column names are id1
      train_idx <- match(
        as.character(training_ids),
        rownames(grm)
      )

      target_idx <- match(
        as.character(target_ids),
        colnames(grm)
      )


      if (anyNA(train_idx)) {
        stop("Some training individuals are missing from GRM.")
      }

      if (anyNA(target_idx)) {
        stop(
          "Some target individuals are missing from GRM in year ",
          target_year
        )
      }


      # Rectangle of relationships:
      #
      # training individuals x target individuals
      relatedness <- grm[
        train_idx,
        target_idx,
        drop = FALSE
      ]

      relatedness_values <- as.vector(
        relatedness
      )


      data.frame(
        target_year = target_year,

        n_target = length(target_ids),

        mean_relatedness = mean(
          relatedness_values,
          na.rm = TRUE
        ),

        sd_relatedness = sd(
          relatedness_values,
          na.rm = TRUE
        ),

        var_relatedness = var(
          relatedness_values,
          na.rm = TRUE
        ),

        precision_relatedness =
          1 / var(
            relatedness_values,
            na.rm = TRUE
          )
      )
    }
  )


  dplyr::bind_rows(results)
}


extract_pairwise_relatedness <- function(
    pheno_data,
    grm,
    train_years,
    target_year,
    unseen_only = TRUE
) {

  # ---------------------------------------------------------------------------
  # Training birds and IDs
  # ---------------------------------------------------------------------------

  split <- make_training_split(
    pheno_data = pheno_data,
    train_years = train_years
  )

  training_ids <- unique(
    pheno_data$id1[
      split$train_rows
    ]
  )

  training_birds <- split$training_birds


  # ---------------------------------------------------------------------------
  # Target birds and IDs
  # ---------------------------------------------------------------------------

  target_data <- pheno_data |>
    dplyr::filter(
      year_num == target_year
    ) |>
    dplyr::distinct(
      ringnr,
      id1
    )

  if (unseen_only) {
    target_data <- target_data |>
      dplyr::filter(
        !(ringnr %in% training_birds)
      )
  }

  if (nrow(target_data) == 0) {
    stop(
      "No target birds found for year ", target_year
    )
  }


  # ---------------------------------------------------------------------------
  # Match to GRM rows/columns
  # ---------------------------------------------------------------------------

  train_idx <- match(
    as.character(training_ids),
    rownames(grm)
  )

  target_idx <- match(
    as.character(target_data$id1),
    colnames(grm)
  )

  if (anyNA(train_idx)) {
    stop("Some training IDs are missing from the GRM.")
  }

  if (anyNA(target_idx)) {
    stop("Some target IDs are missing from the GRM.")
  }


  # ---------------------------------------------------------------------------
  # Extract all pairwise training-target relatedness values
  # ---------------------------------------------------------------------------

  relatedness_matrix <- grm[
    train_idx,
    target_idx,
    drop = FALSE
  ]

  data.frame(
    target_year = target_year,
    relatedness = as.vector(relatedness_matrix)
  )
}
