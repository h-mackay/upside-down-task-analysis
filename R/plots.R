# Figures: residual diagnostics; bar figures; rmcorr scatterplots

# Diagnostics ----
# Q-Q + residuals vs fitted for one DV (and phase); returns residual skew
check_residuals <- function(
    model_frame,
    dv,
    within,
    id = "participant_id",
    label = dv,
    phase = NULL,
    out_dir,
    prefix,
    width = 6,
    height = 5,
    dpi = 150
) {

  stopifnot(!anyNA(model_frame[[dv]]))

  # Participant as fixed effect: approximates repeated-measures ANOVA residuals
  formula_rhs <- paste(id, "+", paste(within, collapse = " * "))
  lm_formula <- as.formula(paste(dv, "~", formula_rhs))
  lm_fit <- lm(lm_formula, data = model_frame)

  residual_data <- tibble::tibble(
    fitted = fitted(lm_fit),
    resid = residuals(lm_fit),
    std_resid = rstandard(lm_fit)
  )

  plot_title <- if (is.null(phase)) {
    label
  } else {
    paste0(stringr::str_to_sentence(phase), ": ", label)
  }

  qq_plot <- ggplot2::ggplot(residual_data, ggplot2::aes(sample = std_resid)) +
    ggplot2::stat_qq(alpha = 0.6) +
    ggplot2::stat_qq_line() +
    ggplot2::labs(
      title = plot_title,
      x = "Theoretical quantiles",
      y = "Residuals (standardized)"
    ) +
    ggplot2::theme_bw()

  rf_plot <- ggplot2::ggplot(residual_data, ggplot2::aes(x = fitted, y = resid)) +
    ggplot2::geom_point(alpha = 0.6) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
    ggplot2::geom_smooth(
      method = "loess",
      formula = y ~ x,
      se = FALSE,
      colour = "steelblue"
    ) +
    ggplot2::labs(
      title = plot_title,
      x = "Fitted values",
      y = "Residuals"
    ) +
    ggplot2::theme_bw()

  file_stem <- paste(c("diag", prefix, dv, phase), collapse = "_")
  qq_path <- file.path(out_dir, paste0(file_stem, "_qq.png"))
  rf_path <- file.path(out_dir, paste0(file_stem, "_resid_fitted.png"))
  ggplot2::ggsave(qq_path, qq_plot, width = width, height = height, dpi = dpi)
  ggplot2::ggsave(rf_path, rf_plot, width = width, height = height, dpi = dpi)

  tibble::tibble(
    dv = dv,
    phase = phase,
    n = dplyr::n_distinct(model_frame[[id]]),
    skew_resid = mean(scale(residual_data$resid)^3)  # skewness
  )
}

# Theme ----
# Shared theme for presentation figures (diagnostics use theme_bw)
theme_upsidedown <- function(base_size = 12) {
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold"),
      legend.position = "bottom",
      plot.caption = ggplot2::element_text(size = base_size - 3, colour = "grey30")
    )
}

# Eye-tracking figures ----
# by_aoi = TRUE: bars by AOI, panels by position; FALSE: proportions vs .5
make_bar_figure <- function(
    data,
    dv,
    phase,
    by_aoi = TRUE,
    show_points = FALSE,
    cfg = load_config()
) {
  stopifnot(phase %in% levels(data$phase))
  phase_data <- droplevels(dplyr::filter(data, phase == .env$phase))

  dv_specs <- c(cfg$dvs_aoi, cfg$dvs_prop)
  dv_spec <- purrr::detect(dv_specs, ~ .x$name == dv)
  y_lab <- if (is.null(dv_spec)) dv else dv_spec$label

  conditional <- grepl("tff|fix_dur", dv)  # only defined when aoi was fixated
  caption <- paste0(
    "Error bars = within-subject (Cousineau-Morey) SE.",
    if (conditional) " Conditional means: cells with no fixations excluded." else "",
    if (!by_aoi) " Dashed line = chance (.5)." else ""
  )

  if (by_aoi) {
    bar_data <- cousineau_morey_se(phase_data, dv, change_type, position, aoi)
    fill_var <- "aoi"
    legend_title <- "AOI"
    fill_colours <- c(top = "white", bottom = "grey60")
  } else {
    bar_data <- cousineau_morey_se(phase_data, dv, change_type, position)
    fill_var <- "position"
    legend_title <- "Position"
    fill_colours <- c(upright = "white", inverted = "grey60")
  }

  dodge <- ggplot2::position_dodge(width = 0.75)
  fig <- ggplot2::ggplot(
    bar_data,
    ggplot2::aes(x = change_type, y = mean, fill = .data[[fill_var]])
  ) +
    ggplot2::geom_col(
      position = dodge,
      width = 0.66,
      colour = "black",
      linewidth = 0.4
    )

  if (!by_aoi) {
    fig <- fig + ggplot2::geom_hline(
      yintercept = 0.5,
      linetype = "dashed",
      colour = "grey40",
      linewidth = 0.5
    )
  }

  if (show_points) {
    fig <- fig + ggplot2::geom_point(
      data = phase_data,
      ggplot2::aes(y = .data[[dv]], fill = .data[[fill_var]]),
      position = ggplot2::position_jitterdodge(
        jitter.width = 0.16,
        dodge.width = 0.75,
        seed = cfg$seed  # same jitter every run
      ),
      size = 1.1,
      alpha = 0.3,
      shape = 16,
      colour = "black",
      na.rm = TRUE
    )
  }

  fig <- fig +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = mean - se_within, ymax = mean + se_within),
      position = dodge,
      width = 0.18,
      linewidth = 0.5
    ) +
    ggplot2::scale_fill_manual(
      name = legend_title,
      values = fill_colours,
      labels = stringr::str_to_title
    ) +
    ggplot2::scale_x_discrete(
      labels = stringr::str_to_title
    ) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0, 0.05)),
      limits = if (!by_aoi) c(0, 1) else c(0, NA)
    ) +
    ggplot2::labs(
      title = paste0(stringr::str_to_sentence(phase), ": ", y_lab),
      x = "Change type",
      y = y_lab,
      caption = caption
    ) +
    theme_upsidedown()

  if (by_aoi) {
    fig <- fig + ggplot2::facet_wrap(
      ggplot2::vars(position),
      labeller = ggplot2::as_labeller(stringr::str_to_title)
    )
  }
  fig
}

# Behavioural figures ----
make_behav_figure <- function(
    data,
    dv,
    by_object = FALSE,
    cfg = load_config()
) {
  dv_spec <- purrr::detect(cfg$dvs_behav, ~ .x$name == dv)
  y_lab <- if (is.null(dv_spec)) dv else dv_spec$label

  if (by_object) {
    bar_data <- cousineau_morey_se(data, dv, object, change_type, position)
  } else {
    # Average over objects within each participant
    subject_means <- data |>
      dplyr::group_by(participant_id, change_type, position) |>
      dplyr::summarise(!!dv := mean(.data[[dv]]), .groups = "drop")
    bar_data <- cousineau_morey_se(subject_means, dv, change_type, position)
  }

  dodge <- ggplot2::position_dodge(width = 0.75)
  fig <- ggplot2::ggplot(
    bar_data,
    ggplot2::aes(x = change_type, y = mean, fill = position)
  ) +
    ggplot2::geom_col(
      position = dodge,
      width = 0.66,
      colour = "black",
      linewidth = 0.4
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = mean - se_within, ymax = mean + se_within),
      position = dodge,
      width = 0.18,
      linewidth = 0.5
    ) +
    ggplot2::scale_fill_manual(
      name = "Position",
      values = c(upright = "white", inverted = "grey60"),
      labels = stringr::str_to_title
    ) +
    ggplot2::scale_x_discrete(
      labels = stringr::str_to_title
    ) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0, 0.05)),
      limits = if (dv == "accuracy") c(0, 1) else c(0, NA)
    ) +
    ggplot2::labs(
      x = "Change type",
      y = y_lab,
      caption = "Error bars = within-subject (Cousineau-Morey) SE."
    ) +
    theme_upsidedown()
  if (by_object) {
    fig <- fig + ggplot2::facet_wrap(
      ggplot2::vars(object),
      labeller = ggplot2::as_labeller(stringr::str_to_title)
    )
  }
  fig
}

# rmcorr figures ----
make_rmcorr_figure <- function(
    data,
    predictor,
    outcome = "ies",
    id = "participant_id",
    stats = NULL,
    title = NULL,
    x_lab = predictor,
    y_lab = outcome
) {
  model_frame <- data.frame(
    pid = factor(data[[id]]),
    x = data[[predictor]],
    y = data[[outcome]]
  )
  model_frame <- tidyr::drop_na(model_frame)
  model_frame <- droplevels(model_frame)

  # Same model as rmcorr (error and df used for ribbon)
  lm_fit <- lm(y ~ pid + x, data = model_frame)

  subtitle <- if (!is.null(stats)) {
    sprintf(
      "r(%d) = %s, p %s, 95%% CI [%s, %s], n = %d",
      stats$df,
      fmt_r(stats$r),
      paste0(if (stats$p < .001) "" else "= ", fmt_p(stats$p)),  # fmt_p auto adds "<"
      fmt_r(stats$ci_low),
      fmt_r(stats$ci_high),
      stats$n
    )
  } else NULL

  # Remove each participant's mean; add back grand mean (rmcorr fits)
  model_frame$x_adj <- model_frame$x -
    ave(model_frame$x, model_frame$pid) +
    mean(model_frame$x)
  model_frame$y_adj <- model_frame$y -
    ave(model_frame$y, model_frame$pid) +
    mean(model_frame$y)
  adj_fit <- lm(y_adj ~ x_adj, data = model_frame)

  grid_x <- data.frame(
    x_adj = seq(min(model_frame$x_adj), max(model_frame$x_adj), length.out = 80)
  )
  predicted <- predict(adj_fit, newdata = grid_x, se.fit = TRUE)
  # Rescale adjusted fit's CI to rmcorr model's error and df
  se_scale <- sigma(lm_fit) / sigma(adj_fit)
  t_crit <- qt(0.975, df.residual(lm_fit))
  ribbon <- data.frame(
    x_adj = grid_x$x_adj,
    fit = predicted$fit,
    lo = predicted$fit - t_crit * predicted$se.fit * se_scale,
    hi = predicted$fit + t_crit * predicted$se.fit * se_scale
  )

  fig <- ggplot2::ggplot(model_frame, ggplot2::aes(x = x_adj, y = y_adj)) +
    ggplot2::geom_ribbon(
      data = ribbon,
      ggplot2::aes(x = x_adj, ymin = lo, ymax = hi),
      inherit.aes = FALSE,
      fill = "grey75",
      alpha = 0.5
    ) +
    ggplot2::geom_point(size = 1.5, colour = "black", alpha = 0.55) +
    ggplot2::geom_line(
      data = ribbon,
      ggplot2::aes(x = x_adj, y = fit),
      inherit.aes = FALSE,
      linewidth = 1.1,
      colour = "royalblue"
    ) +
    ggplot2::labs(
      caption = paste(
        "Within-participant adjusted. Line = rmcorr common slope;",
        "ribbon = 95% CI."
      )
    )
  fig + ggplot2::labs(title = title, subtitle = subtitle, x = x_lab, y = y_lab) +
    theme_upsidedown()
}
