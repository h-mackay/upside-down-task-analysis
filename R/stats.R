# ANOVAs + effect sizes; follow-up tests; repeated-measures correlations

# ANOVA ----
# Repeated-measures ANOVA: Type III; Greenhouse-Geisser correction; partial eta squared
run_anova <- function(
    data,
    dv,
    within,
    id = "participant_id",
    phase = NULL,
    type = 3,
    correction = "GG",
    es = "pes",
    ci = 0.90
) {

  if (!is.null(phase)) {
    data <- droplevels(dplyr::filter(data, phase == .env$phase))
  }

  # Participants missing any cell = dropped (recorded in dropped_ids)
  na_by_id <- tapply(data[[dv]], data[[id]], function(x) any(is.na(x)))

  kept_ids <- names(na_by_id)[!na_by_id]
  dropped_ids <- names(na_by_id)[na_by_id]
  model_frame <- droplevels(dplyr::filter(data, .data[[id]] %in% kept_ids))
  n <- length(kept_ids)

  expected_cells <- prod(lengths(purrr::map(within, ~ levels(model_frame[[.x]]))))
  stopifnot(
    all(table(model_frame[[id]]) == expected_cells),
    !anyDuplicated(model_frame[c(id, within)])
  )

  aov_fit <- afex::aov_ez(
    id = id,
    dv = dv,
    data = model_frame,
    within = within,
    type = type,
    anova_table = list(correction = correction, es = es)
  )
  anova_gg <- as.data.frame(aov_fit$anova_table)

  # Partial eta squared + 90% CI (for .05 F test)
  anova_uncorr <- as.data.frame(anova(aov_fit, correction = "none"))
  es_ci <- effectsize::F_to_eta2(
    f = anova_gg[["F"]],
    df = anova_uncorr[["num Df"]],
    df_error = anova_uncorr[["den Df"]],
    ci = ci,
    alternative = "two.sided"  # default = one-sided
  )

  # Reported degrees of freedom = corrected (to match corrected p)
  anova_effects <- tibble::tibble(
    effect = rownames(anova_gg),
    num_df = anova_gg[["num Df"]],
    den_df = anova_gg[["den Df"]],
    f = anova_gg[["F"]],
    p = anova_gg[["Pr(>F)"]],  # corrected
    pes = anova_gg[[es]],  # same with/without correction
    pes_low = es_ci[["CI_low"]],
    pes_high = es_ci[["CI_high"]]
  )

  # Sphericity: Mauchly's test + Greenhouse-Geisser / Huynh-Feldt epsilons
  # Hides "HF eps > 1" warning; only GG correction used
  sph_raw <- suppressWarnings(summary(aov_fit$Anova))
  no_mauchly <- tibble::tibble(
    effect = character(),
    mauchly_w = double(),
    mauchly_p = double(),
    gg_eps = double(),
    hf_eps = double()
  )
  mauchly <- sph_raw$sphericity.tests
  sphericity <-
    if (is.null(mauchly) || nrow(mauchly) == 0L) {  # 2-level factors (no test)
      no_mauchly
    } else {
      epsilons <- sph_raw$pval.adjustments
      epsilons <- epsilons[rownames(mauchly), , drop = FALSE]
      tibble::tibble(
        effect = rownames(mauchly),
        mauchly_w = mauchly[, "Test statistic"],
        mauchly_p = mauchly[, "p-value"],
        gg_eps = epsilons[, "GG eps"],
        hf_eps = epsilons[, "HF eps"]
      )
    }

  list(
    dv = dv,
    phase = phase,
    within = within,
    n = n,
    dropped_ids = dropped_ids,
    model_frame = model_frame,
    aov_fit = aov_fit,
    anova_effects = anova_effects,
    sphericity = sphericity
  )
}

# Follow-up tests ----
# Pairwise tests within each level of `by`; each uses own error term (df = n - 1)
run_simple_effects <- function(
    anova_result,
    compare,
    by,
    adjust = "bonferroni",
    level = 0.95,
    model = "multivariate"
) {

  emm <- emmeans::emmeans(
    anova_result$aov_fit,
    specs = as.formula(paste0("~", compare, "|", by)),
    model = model
  )
  contrast_family <- rbind(pairs(emm))
  s <- summary(contrast_family, adjust = adjust, infer = TRUE, level = level)

  k_expected <- nlevels(anova_result$model_frame[[by]]) *
    choose(nlevels(anova_result$model_frame[[compare]]), 2)
  stopifnot(
    nrow(s) == k_expected,
    all(s$df == anova_result$n - 1)
  )

  tibble::tibble(
    dv = anova_result$dv,
    phase = anova_result$phase,
    by = by,
    by_level = as.character(s[[by]]),
    contrast = as.character(s$contrast),
    estimate = s$estimate,
    se = s$SE,
    df = s$df,
    ci_low = s$lower.CL,
    ci_high = s$upper.CL,
    t = s$t.ratio,
    p_adj = s$p.value,
    k = k_expected,
    adjust = adjust
  )
}

# Each cell mean vs fixed value (e.g. proportions vs .5); df = n - 1
run_cells_vs_null <- function(
    anova_result,
    null = 0.5,
    adjust = "bonferroni",
    level = 0.95,
    model = "multivariate"
) {

  emm_formula <- as.formula(paste0("~", paste(anova_result$within, collapse = " * ")))
  emm <- emmeans::emmeans(anova_result$aov_fit, specs = emm_formula, model = model)
  s <- summary(emm, infer = TRUE, null = null, adjust = adjust, level = level)

  k_expected <- prod(
    purrr::map_int(anova_result$within, ~ nlevels(anova_result$model_frame[[.x]]))
  )
  stopifnot(
    nrow(s) == k_expected,
    all(s$df == anova_result$n - 1)
  )

  cells <- tibble::as_tibble(s[anova_result$within])
  cells <- dplyr::mutate(
    cells,
    dv = anova_result$dv,
    phase = anova_result$phase,
    .before = 1
  )

  dplyr::mutate(
    cells,
    emmean = s$emmean,
    se = s$SE,
    df = s$df,
    ci_low = s$lower.CL,
    ci_high = s$upper.CL,
    null = null,
    t = s$t.ratio,
    p_adj = s$p.value,
    k = k_expected,
    adjust = adjust
  )
}

# Repeated-measures correlation ----
run_rmcorr <- function(
    data,
    outcome,
    predictor,
    id = "participant_id",
    min_obs = 3,
    ci = 0.95
) {

  model_frame <- data.frame(
    pid = factor(data[[id]]),
    y = data[[outcome]],
    x = data[[predictor]]
  )
  model_frame <- tidyr::drop_na(model_frame)

  n_by_id <- table(model_frame$pid)
  # Participants with < min_obs complete pairs = dropped
  kept_ids <- names(n_by_id)[n_by_id >= min_obs]
  dropped_ids <- setdiff(unique(as.character(data[[id]])), kept_ids)
  model_frame <- droplevels(dplyr::filter(model_frame, pid %in% kept_ids))
  stopifnot(length(kept_ids) >= 2)

  rmcorr_fit <- rmcorr::rmcorr(
    participant = pid,
    measure1 = y,
    measure2 = x,
    dataset = model_frame,
    CI.level = ci
  )

  tibble::tibble(
    outcome = outcome,
    predictor = predictor,
    n = length(kept_ids),
    n_obs = nrow(model_frame),
    dropped = paste(dropped_ids, collapse = ";"),
    r = rmcorr_fit$r,
    df = rmcorr_fit$df,
    p = rmcorr_fit$p,
    ci_low = rmcorr_fit$CI[1],
    ci_high = rmcorr_fit$CI[2]
  )
}
