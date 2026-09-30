# Plotting functions for temporal genomic prediction
plot_unseen_accuracy <- function(
    evaluation_data,
    min_year = NULL,
    max_year = NULL
) {

  plot_data <- evaluation_data |>
    dplyr::filter(
      !is.na(accuracy_unseen)
    )

  if (!is.null(min_year)) {
    plot_data <- plot_data |>
      dplyr::filter(target_year >= min_year)
  }

  if (!is.null(max_year)) {
    plot_data <- plot_data |>
      dplyr::filter(target_year <= max_year)
  }

  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = target_year,
      y = accuracy_unseen,
      color = method
    )
  ) +
    ggplot2::geom_hline(
      yintercept = 0,
      linetype = "dashed",
      linewidth = 0.5,
      color = "grey60"
    ) +
    ggplot2::geom_line(
      linewidth = 1
    ) +
    ggplot2::geom_point(
      size = 2.8
    ) +
    ggplot2::scale_x_continuous(
      breaks = sort(
        unique(plot_data$target_year)
      )
    ) +
    ggplot2::scale_color_brewer(
      palette = "Dark2"
    ) +
    ggplot2::labs(
      title = "Temporal genomic prediction accuracy",
      subtitle = "Birds not observed in the training period",
      x = "Prediction year",
      y = "Prediction accuracy (Pearson correlation)",
      color = "Method"
    ) +
    ggplot2::theme_minimal(
      base_size = 13
    ) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "top",
      axis.text.x = ggplot2::element_text(
        angle = 45,
        hjust = 1
      ),
      plot.title = ggplot2::element_text(
        face = "bold"
      )
    )
}

plot_relatedness_over_time <- function(
    relatedness_data,
    min_year = NULL,
    max_year = NULL
) {

  # ---------------------------------------------------------------------------
  # Select year range
  # ---------------------------------------------------------------------------

  plot_data <- relatedness_data

  if (!is.null(min_year)) {
    plot_data <- plot_data |>
      dplyr::filter(target_year >= min_year)
  }

  if (!is.null(max_year)) {
    plot_data <- plot_data |>
      dplyr::filter(target_year <= max_year)
  }


  # ---------------------------------------------------------------------------
  # Reshape for plotting
  # ---------------------------------------------------------------------------

  plot_data <- plot_data |>
    dplyr::select(
      target_year,
      mean_relatedness,
      precision_relatedness
    ) |>
    tidyr::pivot_longer(
      cols = c(
        mean_relatedness,
        precision_relatedness
      ),
      names_to = "metric",
      values_to = "value"
    ) |>
    dplyr::mutate(
      metric = dplyr::recode(
        metric,
        mean_relatedness = "Mean relatedness",
        precision_relatedness = "Precision"
      )
    )


  # ---------------------------------------------------------------------------
  # Plot
  # ---------------------------------------------------------------------------

  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = target_year,
      y = value
    )
  ) +
    ggplot2::geom_line(
      linewidth = 1
    ) +
    ggplot2::geom_point(
      size = 2.5
    ) +
    ggplot2::facet_wrap(
      ~ metric,
      scales = "free_y",
      ncol = 1
    ) +
    ggplot2::scale_x_continuous(
      breaks = sort(
        unique(plot_data$target_year)
      )
    ) +
    ggplot2::labs(
      title = "Genomic relatedness to the training population over time",
      subtitle = "Unseen birds relative to birds in the training period",
      x = "Target year",
      y = NULL
    ) +
    ggplot2::theme_minimal(
      base_size = 13
    ) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(
        face = "bold"
      ),
      plot.title = ggplot2::element_text(
        face = "bold"
      ),
      axis.text.x = ggplot2::element_text(
        angle = 45,
        hjust = 1
      )
    )
}

plot_relatedness_histogram_compare <- function(
    pheno_data,
    grm,
    train_years,
    compare_years = c(2006, 2016),
    unseen_only = TRUE,
    bins = 50
) {

  plot_data <- dplyr::bind_rows(
    lapply(
      compare_years,
      function(year) {
        extract_pairwise_relatedness(
          pheno_data = pheno_data,
          grm = grm,
          train_years = train_years,
          target_year = year,
          unseen_only = unseen_only
        )
      }
    )
  )

  summary_data <- plot_data |>
    dplyr::group_by(target_year) |>
    dplyr::summarise(
      mean_relatedness = mean(relatedness, na.rm = TRUE),
      sd_relatedness = sd(relatedness, na.rm = TRUE),
      var_relatedness = var(relatedness, na.rm = TRUE),
      precision_relatedness = 1 / var_relatedness,
      .groups = "drop"
    ) |>
    dplyr::mutate(
      label = paste0(
        "Mean = ", round(mean_relatedness, 4), "\n",
        "SD = ", round(sd_relatedness, 4), "\n",
        "Precision = ", round(precision_relatedness, 1)
      )
    )

  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = relatedness,
      y = after_stat(density)
    )
  ) +
    ggplot2::geom_histogram(
      bins = bins,
      fill = "grey75",
      color = "white"
    ) +
    ggplot2::geom_vline(
      data = summary_data,
      ggplot2::aes(
        xintercept = mean_relatedness
      ),
      color = "red",
      linewidth = 0.8,
      inherit.aes = FALSE
    ) +
    ggplot2::geom_text(
      data = summary_data,
      ggplot2::aes(
        x = Inf,
        y = Inf,
        label = label
      ),
      hjust = 1.05,
      vjust = 1.1,
      size = 3.5,
      inherit.aes = FALSE
    ) +
    ggplot2::facet_wrap(
      ~ target_year,
      ncol = 1,
      scales = "fixed"
    ) +
    ggplot2::labs(
      title = "Distribution of training-target genomic relatedness",
      subtitle = "Each panel shows all pairwise GRM values between training birds and unseen target-year birds",
      x = "Pairwise genomic relatedness",
      y = "Density"
    ) +
    ggplot2::theme_minimal(
      base_size = 13
    ) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(face = "bold")
    ) +
    ggplot2::coord_cartesian(
      xlim = c(-0.12, 0.12))
}
