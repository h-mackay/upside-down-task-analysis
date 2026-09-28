# et_data + behav_data -> figures

suppressPackageStartupMessages(library(tidyverse))
source(here::here("R", "utils.R"))
source(here::here("R", "aggregate.R"))
source(here::here("R", "plots.R"))
cfg <- load_config()

et_data <- readRDS(here::here(cfg$paths$et_data))

aoi_obj <- set_factor_levels(et_data$aoi, cfg)
prop_obj <- set_factor_levels(et_data$prop, cfg)
collapsed <- collapse_objects(aoi_obj, prop_obj)

# ET single figures ----
# One per DV x phase
et_dv_specs <- bind_rows(
  tibble(dv = map_chr(cfg$dvs_aoi, "name"), dv_type = "aoi"),
  tibble(dv = map_chr(cfg$dvs_prop, "name"), dv_type = "prop")
)
et_figure_combos <- crossing(et_dv_specs, phase = c("test", "target"))

et_figure_dir <- file.path(here::here(cfg$paths$out_figures), "et")
for (subdir in c("singles", "sheets/bars", "sheets/points")) {
  dir.create(file.path(et_figure_dir, subdir), recursive = TRUE, showWarnings = FALSE)
}

pwalk(et_figure_combos, function(dv, dv_type, phase) {
  fig <- make_bar_figure(
    data = if (dv_type == "aoi") collapsed$aoi else collapsed$prop,
    dv = dv,
    phase = phase,
    by_aoi = dv_type == "aoi",
    cfg = cfg
  )
  fname <- paste0("fig_et_", dv, "_", phase, ".png")
  ggsave(
    file.path(et_figure_dir, "singles", fname),
    fig,
    width = if (dv_type == "aoi") 7.5 else 5.5,
    height = 5,
    dpi = 300
  )
})

# ET sheets ----
# Four panels each; by_aoi = panels split by AOI; with and without participant points
aoi_metrics <- c(
  "n_fix",
  "fix_dur_ms",
  "dwell_ms",
  "n_sacc"
)
top_metrics <- c(
  "tff_ms",
  "prop_top_dwell",
  "prop_top_fix",
  "prop_top_sacc"
)

sheet_specs <- tribble(
  ~sheet, ~dvs, ~by_aoi,
  "aoi_metrics", aoi_metrics, c(TRUE, TRUE, TRUE, TRUE),
  "top_metrics", top_metrics, c(TRUE, FALSE, FALSE, FALSE)
)

et_sheet_combos <- expand_grid(
  sheet_specs,
  phase = c("test", "target"),
  points = c(FALSE, TRUE)
)
pwalk(et_sheet_combos, function(sheet, dvs, by_aoi, phase, points) {
  panels <- map2(dvs, by_aoi, function(dv, panel_by_aoi) {
    make_bar_figure(
      data = if (panel_by_aoi) collapsed$aoi else collapsed$prop,
      dv = dv,
      phase = phase,
      by_aoi = panel_by_aoi,
      show_points = points,
      cfg = cfg
    ) + labs(title = NULL, caption = NULL)
  })
  sheet_fig <- patchwork::wrap_plots(panels, ncol = 2, guides = "collect") +
    patchwork::plot_annotation(
      tag_levels = "A",
      title = paste0(
        str_to_sentence(phase),
        ": ",
        if (sheet == "aoi_metrics") "AOI metrics" else "Top-AOI metrics"
      ),
      caption = paste(
        "Error bars = within-subject (Cousineau-Morey) SE.",
        if (sheet == "aoi_metrics") {
          "Mean fixation duration: cells with no fixations excluded."
        } else {
          "Time to first fixation: cells with no fixations excluded. Dashed line = chance (.5)."
        }
      ),
      theme = theme_upsidedown()
    ) &
    theme(legend.position = "bottom")
  fname <- paste0("sheet_et_", sheet, "_", phase, if (points) "_points" else "", ".png")
  ggsave(
    file.path(et_figure_dir, "sheets", if (points) "points" else "bars", fname),
    sheet_fig,
    width = 13,
    height = 9.5,
    dpi = 300
  )
})

# Behavioural figures ----
behav_data <- readRDS(here::here(cfg$paths$behav_data))
behav_data <- set_factor_levels(behav_data, cfg)

behav_figure_dir <- file.path(here::here(cfg$paths$out_figures), "behavioural")
dir.create(behav_figure_dir, recursive = TRUE, showWarnings = FALSE)

# Accuracy, RT, IES: collapsed and by object
behav_figure_combos <- expand_grid(
  dv = map_chr(cfg$dvs_behav, "name"),
  by_object = c(FALSE, TRUE)
)
pwalk(behav_figure_combos, function(dv, by_object) {
  fig <- make_behav_figure(behav_data, dv, by_object = by_object, cfg = cfg)
  fname <- paste0(
    "fig_behav_",
    dv,
    if (by_object) "_by_object" else "",
    ".png"
  )
  ggsave(
    file.path(behav_figure_dir, fname),
    fig,
    width = if (by_object) 8.5 else 5.5,
    height = 5,
    dpi = 300
  )
})

# Delta IES (inverted - upright)
delta_data <- behav_data |>
  select(participant_id, object, position, change_type, ies) |>
  pivot_wider(names_from = position, values_from = ies) |>
  mutate(delta_ies = inverted - upright)

# Within-subject SE computed separately for each change type (one per panel)
delta_bar_data <- map(levels(behav_data$change_type), function(change_type_name) {
  delta_data |>
    filter(change_type == change_type_name) |>
    cousineau_morey_se("delta_ies", object) |>
    mutate(change_type = change_type_name, .before = 1)
}) |>
  list_rbind()

delta_fig <- ggplot(delta_bar_data, aes(x = object, y = mean)) +
  geom_col(width = 0.66, colour = "black", linewidth = 0.4, fill = "grey85") +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  geom_errorbar(
    aes(ymin = mean - se_within, ymax = mean + se_within),
    width = 0.18,
    linewidth = 0.5
  ) +
  facet_wrap(
    vars(change_type),
    labeller = as_labeller(str_to_title)
  ) +
  scale_x_discrete(labels = str_to_title) +
  labs(
    x = "Object",
    y = "Delta IES (inverted - upright, ms)",
    caption = paste(
      "Positive = worse inverted. Error bars = within-subject",
      "(Cousineau-Morey) SE, per change type."
    )
  ) +
  theme_upsidedown()
fname <- "fig_behav_delta_ies.png"
ggsave(
  file.path(behav_figure_dir, fname),
  delta_fig,
  width = 8,
  height = 4.6,
  dpi = 300
)

# Run log ----
write(
  c(
    "",
    "# 03_figures.R",
    paste0("date: ", format(Sys.time(), "%Y-%m-%d %H:%M"))
  ),
  here::here(cfg$paths$run_log),
  append = TRUE
)
