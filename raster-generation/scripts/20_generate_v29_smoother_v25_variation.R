#!/usr/bin/env Rscript

# V29: preserve V25's inner ellipse, radial transition, and viridis mapping while
# modestly increasing the autocorrelation of its V19-derived exterior variation.

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
make_v30 <- "--v30" %in% commandArgs(trailingOnly = TRUE)

v19_path <- file.path(output_dir, "12_clean-continuous-field-v19.rds")
v21_path <- file.path(output_dir, "14_radial-dominant-field-v21.rds")
v25_path <- file.path(output_dir, "18_v21-core-v19-contours-v25.rds")
if (!all(file.exists(c(v19_path, v21_path, v25_path)))) {
  stop("The preserved V19, V21, and V25 fields are required.")
}

v19 <- readRDS(v19_path)
v21 <- readRDS(v21_path)
v25 <- readRDS(v25_path)
field_v21 <- v21$field
distance_v21 <- v21$warped_distance
distance_v19 <- v19$warped_distance
height <- v21$height
width <- v21$width
display_limits <- v21$parameters$reference_display_limits

# A coarse mean followed by bilinear reconstruction supplies the low-pass field.
# Only 25% of each V19-derived source is replaced, making this a restrained
# correlation-length increase rather than a redesign.
low_pass <- function(values, factor = 24L) {
  source_raster <- terra::rast(values)
  coarse <- terra::aggregate(source_raster, fact = factor, fun = mean)
  smoothed <- terra::as.matrix(
    terra::resample(coarse, source_raster, method = "bilinear"),
    wide = TRUE
  )
  smoothed[!is.finite(smoothed)] <- values[!is.finite(smoothed)]
  smoothed
}

smoothing_share <- if (make_v30) 0.85 else 0.25
smoothing_factor <- if (make_v30) 64L else 24L

# Smooth the logarithmic deformation, not the distance itself. This retains
# positive nested contours and V19's broad organic geometry.
safe_v21 <- pmax(distance_v21, 1e-8)
safe_v19 <- pmax(distance_v19, 1e-8)
log_deformation <- log(safe_v19 / safe_v21)
smooth_log_deformation <- low_pass(log_deformation, factor = smoothing_factor)
log_deformation_v29 <-
  (1 - smoothing_share) * log_deformation +
  smoothing_share * smooth_log_deformation
distance_v19_smooth <- distance_v21 * exp(log_deformation_v29)

# Apply the same restrained smoothing to V19's broad residual variation.
residual_v19 <- v19$field - v19$radial_component
smooth_residual_v19 <- low_pass(residual_v19, factor = smoothing_factor)
residual_v19_v29 <-
  (1 - smoothing_share) * residual_v19 +
  smoothing_share * smooth_residual_v19

# Reuse V25's exact continuous outward envelope and residual strength.
inner_ellipse <- v25$parameters$preserved_inner_ellipse
envelope_curvature <- v25$parameters$outward_envelope_curvature
v19_residual_gain <- v25$parameters$v19_residual_gain
outer_distance <- stats::quantile(distance_v21, 0.995, names = FALSE)
u <- (distance_v21 - inner_ellipse) / (outer_distance - inner_ellipse)
u[u < 0] <- 0
u[u > 1] <- 1
outward_envelope <-
  (1 - exp(-envelope_curvature * u^2)) /
  (1 - exp(-envelope_curvature))

distance_v29 <- exp(
  (1 - outward_envelope) * log(safe_v21) +
    outward_envelope * log(pmax(distance_v19_smooth, 1e-8))
)
distance_v29[distance_v21 == 0] <- 0

acceleration <- v21$parameters$acceleration
d <- distance_v29 / v21$parameters$distance_limit
radial_v29 <- expm1(acceleration * d) / expm1(acceleration)

residual_v21 <- field_v21 - v21$radial_component
residual_v29 <-
  (1 - outward_envelope) * residual_v21 +
  outward_envelope * v19_residual_gain * residual_v19_v29

field_v29 <- radial_v29 + residual_v29
inside <- distance_v21 <= inner_ellipse
field_v29[inside] <- field_v21[inside]
distance_v29[inside] <- distance_v21[inside]

palette_begin <- v21$parameters$palette_begin
palette_end <- v21$parameters$palette_end

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

version <- if (make_v30) "v30" else "v29"
prefix <- if (make_v30) "23" else "22"
render_fixed(
  field_v29,
  file.path(output_dir, sprintf("%s_v25-smoother-variation-%s.png", prefix, version))
)
render_fixed(
  aggregate_grid(field_v29, 90L, 160L),
  file.path(study_dir, sprintf("%s_%s_small-cells_160x90.png", prefix, version))
)
render_fixed(
  aggregate_grid(field_v29, 45L, 80L),
  file.path(study_dir, sprintf("%s_%s_medium-cells_80x45.png", prefix, version))
)

v25$field <- field_v29
v25$warped_distance <- distance_v29
v25$radial_component <- radial_v29
v25$parameters$source_version <- "v25"
v25$parameters$variation_smoothing_share <- smoothing_share
v25$parameters$variation_smoothing_factor <- smoothing_factor
saveRDS(
  v25,
  file.path(output_dir, sprintf("%s_v25-smoother-variation-%s.rds", prefix, version))
)

cat(toupper(version), " smoother V25 variation complete\n", sep = "")
cat("Low-pass share:", smoothing_share, "\n")
cat("Low-pass factor:", smoothing_factor, "\n")
cat("Maximum absolute change inside:", max(abs(field_v29[inside] - field_v21[inside])), "\n")
