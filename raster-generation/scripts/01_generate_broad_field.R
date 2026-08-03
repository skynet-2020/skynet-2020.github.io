#!/usr/bin/env Rscript

# Layer 1: broad isotropic background field
#
# This script intentionally creates only the lowest-frequency spatial layer.
# It does not add the central basin, local minima, boundary warping, or fine
# texture. Those are separate stages so their effects remain inspectable.

required_packages <- c("ambient", "viridisLite", "png")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Install the required packages before running this script: ",
    paste(missing_packages, collapse = ", ")
  )
}

# Resolve paths from the script location so the script works from any directory.
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1) {
  stop("Run this file with Rscript so its output path can be resolved.")
}

script_path <- normalizePath(sub("^--file=", "", script_arg))
raster_root <- dirname(dirname(script_path))
output_dir <- file.path(raster_root, "outputs", "landscape")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Work at a moderate 16:9 resolution. The final field can later be regenerated
# at production resolution using the same seed and parameters.
width <- 1600L
height <- 900L
seed <- 20260803L

set.seed(seed)

# Two low-frequency simplex octaves establish broad regional structure without
# introducing the medium- and fine-scale layers reserved for later phases.
broad_field <- ambient::noise_simplex(
  dim = c(height, width),
  frequency = 0.0018,
  fractal = "fbm",
  octaves = 2,
  lacunarity = 1.8,
  gain = 0.42,
  pertubation = "none"
)

# Quantile clipping prevents isolated extremes from controlling the diagnostic
# color stretch. This changes display values only, not spatial geometry.
display_limits <- stats::quantile(
  broad_field,
  probs = c(0.01, 0.99),
  names = FALSE
)

display_field <- pmin(pmax(broad_field, display_limits[1]), display_limits[2])
display_field <- ambient::normalise(display_field)

# Use full standard viridis for this structural diagnostic. Palette truncation
# will be tuned only after the complete numeric surface exists.
palette <- viridisLite::viridis(256, option = "D", begin = 0, end = 1)
color_index <- pmax(1L, pmin(256L, floor(display_field * 255) + 1L))
rgb_lookup <- grDevices::col2rgb(palette) / 255

rgb_image <- array(0, dim = c(height, width, 3L))
for (channel in seq_len(3L)) {
  rgb_image[, , channel] <- matrix(
    rgb_lookup[channel, color_index],
    nrow = height,
    ncol = width
  )
}

field_path <- file.path(output_dir, "01_broad-field.rds")
image_path <- file.path(output_dir, "01_broad-field.png")

saveRDS(
  list(
    field = broad_field,
    seed = seed,
    width = width,
    height = height,
    parameters = list(
      generator = "simplex",
      frequency = 0.0018,
      fractal = "fbm",
      octaves = 2,
      lacunarity = 1.8,
      gain = 0.42
    )
  ),
  field_path
)

png::writePNG(rgb_image, image_path)

cat("Layer 1 complete\n")
cat("Seed:", seed, "\n")
cat("Dimensions:", width, "x", height, "\n")
cat("Field range:", paste(signif(range(broad_field), 5), collapse = " to "), "\n")
cat("Image:", image_path, "\n")
cat("Numeric field:", field_path, "\n")
