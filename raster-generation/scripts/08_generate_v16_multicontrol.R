#!/usr/bin/env Rscript

# V16: three independently controlled changes from preserved V8:
#   1. broaden the dark core and its dark surroundings;
#   2. smooth the center more strongly than the exterior;
#   3. reduce the spatial footprint of the brightest yellow while preserving it.
# All masks and value curves are continuously differentiable.

required_packages <- c("terra", "viridisLite", "png")
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
field_v8 <- source$field
distance <- source$warped_distance
limits <- source$parameters$reference_display_limits
make_v17 <- "--v17" %in% commandArgs(trailingOnly = TRUE)

x <- (field_v8 - limits[1]) / diff(limits)
x[x < 0] <- 0
x[x > 1] <- 1

# Gaussian smoothing is blended spatially: strong near the text area and fading
# continuously to zero toward the exterior. Exterior heterogeneity stays crisp.
sigma <- if (make_v17) 42 else 10
if (make_v17) {
  # Coarse averaging followed by bilinear interpolation creates a genuinely
  # broad smoothing scale without an impractically large convolution kernel.
  source_raster <- terra::rast(x)
  coarse_raster <- terra::aggregate(source_raster, fact = 40, fun = mean)
  smoothed <- terra::as.matrix(
    terra::resample(coarse_raster, source_raster, method = "bilinear"),
    wide = TRUE
  )
} else {
  axis <- seq(-4 * sigma, 4 * sigma)
  gaussian <- outer(axis, axis, function(a, b) exp(-(a^2 + b^2) / (2 * sigma^2)))
  gaussian <- gaussian / sum(gaussian)
  smoothed <- terra::as.matrix(
    terra::focal(terra::rast(x), w = gaussian, na.policy = "omit", fillvalue = NA),
    wide = TRUE
  )
}
smoothed[!is.finite(smoothed)] <- x[!is.finite(smoothed)]
center_smoothing_weight <- if (make_v17) {
  0.94 * exp(-(distance / 0.95)^6)
} else {
  0.78 * exp(-(distance / 0.72)^4)
}
x_center <- x * (1 - center_smoothing_weight) +
  smoothed * center_smoothing_weight

# A broad multiplicative low-value envelope deepens the core AND its surrounding
# violet/blue field. Multiplication preserves an internal gradient rather than
# replacing the center with a single flat target value.
if (make_v17) {
  # Stronger and much broader than V16. The exponent itself changes smoothly
  # over space: low-value progress is slow through most of the interior, then
  # accelerates continuously as it returns to V8 toward the corners.
  local_exponent <- 1 + 1.05 * exp(-(distance / 1.22)^4)
  x_dark <- x_center^local_exponent
} else {
  dark_envelope <- 0.22 * exp(-(distance / 1.02)^3)
  x_dark <- x_center * (1 - dark_envelope)
}

# Smoothly compress only the upper register. The power curve reaches exactly 1,
# so the brightest viridis yellow survives, but fewer pixels reach that range.
# Quintic smootherstep introduces the compression with no visible threshold.
yellow_start <- 0.60
yellow_full <- 0.82
if (make_v17) {
  # Leave the exterior register uncompressed. Combined with the spatial return
  # above, this makes the corner transitions faster and less smoothed.
  x_v16 <- x_dark
} else {
  t <- (x_dark - yellow_start) / (yellow_full - yellow_start)
  t[t < 0] <- 0
  t[t > 1] <- 1
  upper_weight <- t^3 * (t * (t * 6 - 15) + 10)
  yellow_compressed <- x_dark^1.18
  x_v16 <- x_dark * (1 - upper_weight) + yellow_compressed * upper_weight
}
field_v16 <- limits[1] + x_v16 * diff(limits)

height <- nrow(field_v16)
width <- ncol(field_v16)
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
  row_counts <- tabulate(row_group, nbins = target_rows)
  col_counts <- tabulate(col_group, nbins = target_cols)
  coarse <- coarse / outer(row_counts, col_counts)
  coarse[row_group, col_group, drop = FALSE]
}

version <- if (make_v17) "v17" else "v16"
prefix <- if (make_v17) "10" else "09"
render_fixed(
  field_v16,
  file.path(output_dir, sprintf("%s_v8-spatial-transition-%s.png", prefix, version))
)
render_fixed(
  aggregate_grid(field_v16, 90L, 160L),
  file.path(study_dir, sprintf("%s_%s_small-cells_160x90.png", prefix, version))
)
render_fixed(
  aggregate_grid(field_v16, 45L, 80L),
  file.path(study_dir, sprintf("%s_%s_medium-cells_80x45.png", prefix, version))
)

source$field <- field_v16
source$parameters$source_version <- "v8"
source$parameters$center_smoothing_sigma <- sigma
source$parameters$center_smoothing_max_weight <- if (make_v17) 0.94 else 0.78
source$parameters$dark_envelope_strength <- if (make_v17) NA_real_ else 0.22
source$parameters$dark_envelope_radius <- if (make_v17) NA_real_ else 1.02
source$parameters$spatial_exponent_strength <- if (make_v17) 1.05 else NA_real_
source$parameters$spatial_exponent_radius <- if (make_v17) 1.22 else NA_real_
source$parameters$yellow_compression_exponent <- if (make_v17) 1 else 1.18
source$parameters$yellow_compression_range <- c(yellow_start, yellow_full)
saveRDS(
  source,
  file.path(output_dir, sprintf("%s_v8-spatial-transition-%s.rds", prefix, version))
)

cat(toupper(version), " spatial-transition study complete\n", sep = "")
cat("Source: V8\n")
