#!/usr/bin/env Rscript

# Layer 2 (continuous model): long-tailed central attractor
#
# The text region is not a masked object. It is the deepest, calmest portion of
# a continuous potential that influences the entire canvas and has no explicit
# outer boundary.

required_packages <- c("ambient", "viridisLite", "png")
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
input_path <- file.path(output_dir, "01_broad-field.rds")
if (!file.exists(input_path)) stop("Run 01_generate_broad_field.R first.")

layer_1 <- readRDS(input_path)
broad_field <- layer_1$field
height <- nrow(broad_field)
width <- ncol(broad_field)

x <- seq(-1, 1, length.out = width)
y <- seq(-1, 1, length.out = height)
x_grid <- matrix(rep(x, each = height), nrow = height, ncol = width)
y_grid <- matrix(rep(y, times = width), nrow = height, ncol = width)

# Broad aspect-aware distance, subsequently warped so equal-value contours do
# not form concentric ellipses.
radius_x <- 0.66
radius_y <- 0.76
base_distance <- sqrt((x_grid / radius_x)^2 + (y_grid / radius_y)^2)

set.seed(layer_1$seed + 201L)
warp_large <- ambient::noise_simplex(
  c(height, width), frequency = 0.00115, fractal = "none"
)
set.seed(layer_1$seed + 202L)
warp_support <- ambient::noise_simplex(
  c(height, width), frequency = 0.00205, fractal = "none"
)
warp <- 0.8 * warp_large + 0.2 * warp_support
warp_limits <- stats::quantile(warp, c(0.02, 0.98), names = FALSE)
warp <- pmin(pmax(warp, warp_limits[1]), warp_limits[2])
warp <- ambient::normalise(warp, to = c(-1, 1))

# Exponential scaling keeps distance positive and creates broad unequal reach.
warped_distance <- base_distance * exp(0.24 * warp)

# This attractor never switches off: its influence decays gradually to every
# edge of the canvas. The exponent below 2 creates a longer tail than a Gaussian.
attractor_strength <- 2.5
attractor_decay <- 1.2
attractor_exponent <- 1.25
central_attractor <- -attractor_strength * exp(
  -(warped_distance / attractor_decay)^attractor_exponent
)

# Calm the center by smoothly reducing broad-field amplitude. This is not a
# bounded mask: the Gaussian attenuation has an infinite, gradual tail.
calm_radius <- 0.7
calm_weight <- exp(-(warped_distance / calm_radius)^2)
background_amplitude <- 1 - 0.88 * calm_weight

combined_field <- broad_field * background_amplitude + central_attractor

render_field <- function(field, path) {
  limits <- stats::quantile(field, c(0.01, 0.99), names = FALSE)
  display <- ambient::normalise(pmin(pmax(field, limits[1]), limits[2]))
  palette <- viridisLite::viridis(256, option = "D")
  index <- pmax(1L, pmin(256L, floor(display * 255) + 1L))
  lookup <- grDevices::col2rgb(palette) / 255
  image <- array(0, c(height, width, 3L))
  for (channel in seq_len(3L)) {
    image[, , channel] <- matrix(lookup[channel, index], height, width)
  }
  png::writePNG(image, path)
}

field_path <- file.path(output_dir, "02_continuous-attractor.rds")
image_path <- file.path(output_dir, "02_continuous-attractor.png")
saveRDS(
  list(
    field = combined_field,
    broad_field = broad_field,
    central_attractor = central_attractor,
    warped_distance = warped_distance,
    calm_weight = calm_weight,
    seed = layer_1$seed,
    width = width,
    height = height,
    parameters = list(
      radius_x = radius_x,
      radius_y = radius_y,
      attractor_strength = attractor_strength,
      attractor_decay = attractor_decay,
      attractor_exponent = attractor_exponent,
      calm_radius = calm_radius,
      central_background_suppression = 0.88
    )
  ),
  field_path
)
render_field(combined_field, image_path)

cat("Continuous Layer 2 complete\n")
cat("Image:", image_path, "\n")
cat("Numeric field:", field_path, "\n")
