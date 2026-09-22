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
    horizon,
    method = NULL
) {

  cutoff_year <- fit$cutoff_year

  target_year <-
    cutoff_year + horizon


  # ---------------------------------------------------------------------------
  # 1. Determine model type
  # ---------------------------------------------------------------------------

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
  # 2. Birds used for training
  # ---------------------------------------------------------------------------

  training_birds <- unique(
    pheno_data$ringnr[
      pheno_data$year_num <= cutoff_year
    ]
  )


  # ---------------------------------------------------------------------------
  # 3. Target-year phenotypes
  #
  # Average repeated measurements within bird.
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
  # 4. Extract predicted breeding values
  # ---------------------------------------------------------------------------

  bv <- extract_breeding_values(
    fit = fit,
    method = method
  )


  # ---------------------------------------------------------------------------
  # 5. Match phenotype and BV
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

      cutoff_year = cutoff_year,

      target_year = target_year,

      horizon = horizon
    )


  # Check prediction availability
  if (anyNA(predictions$bv)) {

    warning(
      sum(is.na(predictions$bv)),
      " target birds have no predicted breeding value."
    )
  }


  # ---------------------------------------------------------------------------
  # 6. Prediction accuracy
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
  # 7. Summary
  # ---------------------------------------------------------------------------

  summary <- data.frame(

    method = method,

    cutoff_year = cutoff_year,

    horizon = horizon,

    target_year = target_year,

    n_target = nrow(predictions),

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

evaluate_horizons <- function(
    fit,
    pheno_data,
    horizons = c(1, 2, 3, 5, 10),
    method = NULL
) {

  evaluations <- lapply(
    horizons,
    function(h) {

      evaluate_temporal_model(
        fit = fit,
        pheno_data = pheno_data,
        horizon = h,
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
