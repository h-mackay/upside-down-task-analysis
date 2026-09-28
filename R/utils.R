# Shared helpers: config, media codes, number formatting, factor levels

# Config ----
load_config <- function() {
  config::get(file = here::here("config.yml"))
}

# Media codes ----
normalize_media_code <- function(x) {
  x |>
    stringr::str_trim() |>
    stringr::str_to_upper() |>
    stringr::str_remove("\\.(BMP|PNG|JPG|JPEG)$")
}

# 1 = "U" (upside down = inverted), 2 = object, 3 = change type (O/S), 4 = variant
media_re <- "^(U?)([BCHM])([OS])([0-9]+)$"

# NA if code doesn't match
extract_object <- function(x) {
  obj <- stringr::str_match(normalize_media_code(x), media_re)[, 3]
  unname(c(B = "banana", C = "chess", H = "hammer", M = "mug")[obj])
}

# Formatting ----
# APA: "< .001"; else 3 decimals (no leading zero)
fmt_p <- function(p) {
  dplyr::case_when(
    is.na(p) ~ NA_character_,
    p < .001 ~ "< .001",
    .default = sub("^0", "", sprintf("%.3f", p))
  )
}

# Degrees of freedom: whole-number = integers; corrected (fractional) = 2 decimals
fmt_dof <- function(x) {
  dplyr::if_else(x %% 1 == 0, sprintf("%d", as.integer(x)), sprintf("%.2f", x))
}

# APA: 2 decimals; no leading zero
fmt_r <- function(x) {
  dplyr::case_when(
    is.na(x) ~ NA_character_,
    .default = sub("^(-?)0", "\\1", sprintf("%.2f", x))
  )
}

# Means ----
# NA when all values missing (instead of NaN)
mean_valid <- function(x) {
  if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
}

# Factors ----
set_factor_levels <- function(df, cfg) {
  f <- purrr::map(cfg$factors, ~ unlist(.x, use.names = FALSE))
  df |>
    dplyr::mutate(
      dplyr::across(dplyr::any_of("object"), ~ factor(.x, levels = f$object)),
      dplyr::across(dplyr::any_of("position"), ~ factor(.x, levels = f$position)),
      dplyr::across(dplyr::any_of("change_type"), ~ factor(.x, levels = f$change_type)),
      dplyr::across(dplyr::any_of("phase"), ~ factor(.x, levels = f$phase)),
      dplyr::across(dplyr::any_of("aoi"), ~ factor(.x, levels = f$aoi))
    )
}
