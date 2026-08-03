#!/usr/bin/env Rscript

# Layer 3 (continuous model): canvas-scale organizing field
#
# The organizer is added to the long-tailed attractor without compositing a
# separate core. Heterogeneity is smoothly attenuated near the text area, making
# the center calm while preserving one continuous field across the canvas.

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
input_path <- file.path(output_dir, "02_continuous-attractor.rds")
if (!file.exists(input_path)) stop("Run 02_add_central_basin.R first.")

layer_2 <- readRDS(input_path)
broad_field <- layer_2$broad_field
base_field <- layer_2$field
warped_distance <- layer_2$warped_distance
height <- nrow(base_field)
width <- ncol(base_field)

# Broad fields produce the canvas-scale arms, lobes, corridors, and gaps.
set.seed(layer_2$seed + 311L)
organizer_large <- ambient::noise_simplex(
  c(height, width), frequency = 0.001, fractal = "none"
)
set.seed(layer_2$seed + 312L)
organizer_support <- ambient::noise_simplex(
  c(height, width), frequency = 0.0021, fractal = "none"
)
organizer_noise <- 0.76 * organizer_large + 0.24 * organizer_support
noise_limits <- stats::quantile(organizer_noise, c(0.02, 0.98), names = FALSE)
organizer_noise <- pmin(pmax(organizer_noise, noise_limits[1]), noise_limits[2])
organizer_noise <- ambient::normalise(organizer_noise, to = c(-1, 1))

# A smooth nonlinear mapping makes negative influence widespread but unequal.
# It has no hard threshold and therefore creates no patch perimeter.
organizer_influence <- stats::plogis(2.8 * (organizer_noise + 0.12))

# Organizer variation fades smoothly—not abruptly—toward the text region. This
# leaves the central attractor calm while allowing its surrounding arms to carry
# the same field outward across the canvas.
organizer_calm_radius <- 0.7
organizer_gate <- 1 - 0.9 * exp(
  -(warped_distance / organizer_calm_radius)^2
)

organizer_strength <- 1.25
organizer_effect <- -organizer_strength * organizer_influence * organizer_gate
combined_field <- base_field + organizer_effect

render_viridis <- function(field, path) {
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

# Gradient magnitude exposes any remaining enclosing seam around the text area.
dx <- cbind(combined_field[, -1] - combined_field[, -width], 0)
dy <- rbind(combined_field[-1, ] - combined_field[-height, ], 0)
gradient <- sqrt(dx^2 + dy^2)
gradient_display <- ambient::normalise(
  pmin(gradient, stats::quantile(gradient, 0.995, names = FALSE))
)
gradient_image <- array(
  rep(gradient_display, 3L),
  dim = c(height, width, 3L)
)

field_path <- file.path(output_dir, "03_continuous-organizer.rds")
image_path <- file.path(output_dir, "03_continuous-organizer.png")
gradient_path <- file.path(output_dir, "03_gradient-diagnostic.png")

saveRDS(
  list(
    field = combined_field,
    broad_field = broad_field,
    central_attractor = layer_2$central_attractor,
    organizer_noise = organizer_noise,
    organizer_influence = organizer_influence,
    organizer_gate = organizer_gate,
    warped_distance = warped_distance,
    seed = layer_2$seed,
    width = width,
    height = height,
    parameters = c(
      layer_2$parameters,
      list(
        organizer_strength = organizer_strength,
        organizer_calm_radius = organizer_calm_radius
      )
    )
  ),
  field_path
)
render_viridis(combined_field, image_path)
png::writePNG(gradient_image, gradient_path)

cat("Continuous Layer 3 complete\n")
cat("Image:", image_path, "\n")
cat("Gradient diagnostic:", gradient_path, "\n")
cat("Numeric field:", field_path, "\n")
