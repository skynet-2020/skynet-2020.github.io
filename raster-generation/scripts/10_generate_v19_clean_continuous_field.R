#!/usr/bin/env Rscript

# V19: clean reconstruction from V8's underlying spatial components.
# No protected core, expansion mask, piecewise remapping, or post-hoc blending.
# One analytic distance function controls the center and every outward transition.

required_packages <- c("viridisLite", "png")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop("Install required packages: ", paste(missing_packages, collapse = ", "))
}

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1) stop("Run this file with Rscript.")
script_path <- normalizePath(sub("^--file=", "", script_arg))
raster_root <- dirname(dirname(script_path))
output_dir <- file.path(raster_root, "outputs", "landscape")
study_dir <- file.path(output_dir, "raster-studies")
dir.create(study_dir, recursive = TRUE, showWarnings = FALSE)

input_path <- file.path(output_dir, "03h_radial-organic-expanded-core-v8.rds")
if (!file.exists(input_path)) stop("The preserved V8 field is missing: ", input_path)

source <- readRDS(input_path)
distance <- source$warped_distance
corner <- source$corner_component
broad <- source$broad_component
extensions <- source$dark_extensions
height <- nrow(distance)
width <- ncol(distance)

# Normalize the organic distance only at the extreme exterior. The exponential
# rise starts with a shallow, nonzero derivative and accelerates continuously.
# Consequently, violet occupies a long physical interval without a flat core.
distance_limit <- stats::quantile(distance, 0.995, names = FALSE)
d <- distance / distance_limit
acceleration <- 2.85
radial <- expm1(acceleration * d) / expm1(acceleration)

# Heterogeneity amplitude grows continuously with distance. The center retains
# only broad, quiet variation; the exterior progressively recovers more texture.
broad_normalized <- (broad - mean(broad)) / stats::sd(as.vector(broad))
texture_growth <- 0.010 + 0.060 * (1 - exp(-(distance / 0.98)^3))
texture <- broad_normalized * texture_growth

# Corner lift is continuous everywhere and subordinate through the interior.
# It allows the final registers to accelerate at the margins without modifying
# or blending the central radial derivative.
corner_lift <- 0.105 * corner^1.7

# Restrained dark irregularity increases outward with the same smooth gate.
extension_gate <- 1 - exp(-(distance / 1.05)^3)
extension_term <- 0.022 * extensions * extension_gate

raw_field <- radial + texture + corner_lift - extension_term

# Fixed robust limits are used only for display. They do not reshape the field.
display_limits <- stats::quantile(raw_field, c(0.005, 0.995), names = FALSE)
field_v19 <- raw_field
palette_begin <- source$parameters$palette_begin
palette_end <- source$parameters$palette_end

render_fixed <- function(values, path) {
  display <- (values - display_limits[1]) / diff(display_limits)
  display[display < 0] <- 0
  display[display > 1] <- 1
  palette <- viridisLite::viridis(
    256, option = "D", begin = palette_begin, end = palette_end
  )
  index <- pmax(1L, pmin(256L, floor(display * 255) + 1L))
  lookup <- grDevices::col2rgb(palette) / 255
  image <- array(0, c(height, width, 3L))
  for (channel in seq_len(3L)) {
    image[, , channel] <- matrix(lookup[channel, index], height, width)
  }
  png::writePNG(image, path)
}

aggregate_grid <- function(values, target_rows, target_cols) {
  row_group <- pmin(target_rows, ceiling(seq_len(nrow(values)) * target_rows / nrow(values)))
  col_group <- pmin(target_cols, ceiling(seq_len(ncol(values)) * target_cols / ncol(values)))
  sums <- rowsum(values, row_group, reorder = FALSE)
  coarse <- t(rowsum(t(sums), col_group, reorder = FALSE))
  coarse <- coarse / outer(
    tabulate(row_group, nbins = target_rows),
    tabulate(col_group, nbins = target_cols)
  )
  coarse[row_group, col_group, drop = FALSE]
}

render_fixed(field_v19, file.path(output_dir, "12_clean-continuous-field-v19.png"))
render_fixed(
  aggregate_grid(field_v19, 90L, 160L),
  file.path(study_dir, "12_v19_small-cells_160x90.png")
)
render_fixed(
  aggregate_grid(field_v19, 45L, 80L),
  file.path(study_dir, "12_v19_medium-cells_80x45.png")
)

saveRDS(
  list(
    field = field_v19,
    warped_distance = distance,
    radial_component = radial,
    texture_component = texture,
    corner_component = corner_lift,
    seed = source$seed,
    width = width,
    height = height,
    parameters = list(
      source_geometry = "v8 underlying components",
      distance_limit = distance_limit,
      acceleration = acceleration,
      texture_center_amplitude = 0.010,
      texture_outer_amplitude = 0.070,
      texture_growth_radius = 0.98,
      corner_weight = 0.105,
      extension_weight = 0.022,
      reference_display_limits = display_limits,
      palette_begin = palette_begin,
      palette_end = palette_end
    )
  ),
  file.path(output_dir, "12_clean-continuous-field-v19.rds")
)

cat("V19 clean continuous field complete\n")
cat("No protected core or post-hoc transition remapping\n")
