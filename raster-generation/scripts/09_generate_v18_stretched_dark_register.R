#!/usr/bin/env Rscript

# V18: preserve V17's center, but dramatically lengthen the surrounding dark
# register. A single analytic rational curve avoids shelves and knots. Spatial
# C2 blends keep the central piece unchanged and return to V17 near the corners.

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

input_path <- file.path(output_dir, "10_v8-spatial-transition-v17.rds")
if (!file.exists(input_path)) stop("The preserved V17 field is missing: ", input_path)

source <- readRDS(input_path)
field_v17 <- source$field
distance <- source$warped_distance
limits <- source$parameters$reference_display_limits

x <- (field_v17 - limits[1]) / diff(limits)
x[x < 0] <- 0
x[x > 1] <- 1

smootherstep <- function(value) {
  value[value < 0] <- 0
  value[value > 1] <- 1
  value^3 * (value * (value * 6 - 15) + 10)
}

# This curve has identical endpoints and changes its rate continuously. It
# strongly compresses the violet/blue value progression, then catches up rapidly
# through the upper register. Unlike a power curve, its slope at zero is 1, so
# it does not create a larger flat minimum.
stretch_strength <- 3.0
stretched <- x / (1 + stretch_strength * x * (1 - x))

# V17 remains exactly unchanged through the central piece. The stretched field
# takes over gradually only in the surrounding transition zone.
inner_blend <- smootherstep((distance - 0.34) / (0.64 - 0.34))

# Return continuously to V17 toward the extreme exterior, retaining the fast,
# heterogeneous corners and their brightest-yellow footprint.
outer_return <- smootherstep((distance - 1.36) / (1.76 - 1.36))
effect_weight <- inner_blend * (1 - outer_return)
x_v18 <- x * (1 - effect_weight) + stretched * effect_weight
field_v18 <- limits[1] + x_v18 * diff(limits)

height <- nrow(field_v18)
width <- ncol(field_v18)
palette_begin <- source$parameters$palette_begin
palette_end <- source$parameters$palette_end

render_fixed <- function(values, path) {
  display <- (values - limits[1]) / diff(limits)
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

render_fixed(field_v18, file.path(output_dir, "11_v17-stretched-dark-register-v18.png"))
render_fixed(
  aggregate_grid(field_v18, 90L, 160L),
  file.path(study_dir, "11_v18_small-cells_160x90.png")
)
render_fixed(
  aggregate_grid(field_v18, 45L, 80L),
  file.path(study_dir, "11_v18_medium-cells_80x45.png")
)

source$field <- field_v18
source$parameters$source_version <- "v17"
source$parameters$dark_register_stretch_strength <- stretch_strength
source$parameters$unchanged_center_distance <- 0.34
source$parameters$full_stretch_distance <- 0.64
source$parameters$outer_return_distance <- c(1.36, 1.76)
saveRDS(source, file.path(output_dir, "11_v17-stretched-dark-register-v18.rds"))

cat("V18 stretched-dark-register study complete\n")
cat("Source: V17; center preserved through distance 0.34\n")
