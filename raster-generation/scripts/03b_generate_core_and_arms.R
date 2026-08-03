#!/usr/bin/env Rscript

# Alternative Layer 3: explicit scale hierarchy
#
# Candidate v2 uses three deliberately separated scales:
#   1. a broad, calm central plateau for text;
#   2. long, wide, smoothly connected arms organizing the whole canvas;
#   3. broad-field variation that can sculpt the arms outside the plateau.
# Short-scale outer heterogeneity is intentionally deferred to the next layer.

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

# Smooth coordinate warping prevents the arms from reading as geometric rays.
set.seed(layer_1$seed + 331L)
warp_x <- ambient::noise_simplex(
  c(height, width), frequency = 0.00125, fractal = "none"
)
set.seed(layer_1$seed + 332L)
warp_y <- ambient::noise_simplex(
  c(height, width), frequency = 0.00135, fractal = "none"
)
warp_x <- ambient::normalise(warp_x, to = c(-1, 1))
warp_y <- ambient::normalise(warp_y, to = c(-1, 1))
x_warped <- x_grid + 0.11 * warp_x
y_warped <- y_grid + 0.11 * warp_y

# A super-Gaussian plateau gives the text region a broad, nearly uniform floor.
# Its transition remains smooth and is obscured by overlapping connected arms.
core_radius_x <- 0.5
core_radius_y <- 0.58
core_distance <- sqrt(
  (x_warped / core_radius_x)^2 +
    (y_warped / core_radius_y)^2
)
core_transition_exponent <- 2
core_plateau <- exp(-(core_distance^core_transition_exponent))

# Long, wide arm definitions. Unequal angles, lengths, widths, and strengths
# avoid symmetry and collectively have no preferred global orientation.
arm_specs <- data.frame(
  angle = c(12, 76, 139, 207, 278, 326) * pi / 180,
  length = c(1.55, 1.35, 1.5, 1.42, 1.48, 1.36),
  width = c(0.24, 0.2, 0.23, 0.22, 0.25, 0.19),
  strength = c(0.92, 0.74, 0.86, 0.77, 0.9, 0.7)
)

arm_union <- matrix(0, nrow = height, ncol = width)
for (i in seq_len(nrow(arm_specs))) {
  spec <- arm_specs[i, ]
  along <- x_warped * cos(spec$angle) + y_warped * sin(spec$angle)
  across <- -x_warped * sin(spec$angle) + y_warped * cos(spec$angle)

  # Each broad ridge begins inside the plateau and extends toward a margin.
  arm <- exp(
    -0.5 * ((along - spec$length / 2) / (spec$length * 0.52))^2 -
      0.5 * (across / spec$width)^2
  )
  arm <- arm * stats::plogis(9 * (along + 0.2)) * spec$strength

  # Smooth probabilistic union avoids hard pmax() seams where arms overlap.
  arm_union <- 1 - (1 - arm_union) * (1 - arm)
}

# The plateau and arms are one connected phenomenon. This smooth union creates
# no pasted central object and guarantees all arms originate within it.
organizing_structure <- 1 - (1 - core_plateau) * (1 - arm_union)

# Retain very little broad-field variation in the plateau, then restore it
# gradually outside. This allows limited broad incursions without sacrificing
# the white-text region.
calm_transition_exponent <- 1.8
calm_weight <- exp(
  -(core_distance / 1.22)^calm_transition_exponent
)
heterogeneity_amplitude <- 1 - 0.94 * calm_weight

structure_depth <- 2.25
background_strength <- 0.78
combined_field <-
  -structure_depth * organizing_structure +
  background_strength * broad_field * heterogeneity_amplitude

limits <- stats::quantile(combined_field, c(0.01, 0.99), names = FALSE)
display <- ambient::normalise(pmin(pmax(combined_field, limits[1]), limits[2]))
palette_begin <- 0.02
palette_end <- 0.58
palette <- viridisLite::viridis(
  256,
  option = "D",
  begin = palette_begin,
  end = palette_end
)
index <- pmax(1L, pmin(256L, floor(display * 255) + 1L))
lookup <- grDevices::col2rgb(palette) / 255
image <- array(0, c(height, width, 3L))
for (channel in seq_len(3L)) {
  image[, , channel] <- matrix(lookup[channel, index], height, width)
}

field_path <- file.path(output_dir, "03d_low-contrast-core-and-arms-v4.rds")
image_path <- file.path(output_dir, "03d_low-contrast-core-and-arms-v4.png")
saveRDS(
  list(
    field = combined_field,
    broad_field = broad_field,
    core_plateau = core_plateau,
    arm_union = arm_union,
    organizing_structure = organizing_structure,
    core_distance = core_distance,
    seed = layer_1$seed,
    width = width,
    height = height,
    parameters = list(
      core_radius_x = core_radius_x,
      core_radius_y = core_radius_y,
      core_transition_exponent = core_transition_exponent,
      calm_transition_exponent = calm_transition_exponent,
      structure_depth = structure_depth,
      background_strength = background_strength,
      palette_begin = palette_begin,
      palette_end = palette_end,
      arm_specs = arm_specs
    )
  ),
  field_path
)
png::writePNG(image, image_path)

cat("Alternative Layer 3 v4 complete\n")
cat("Image:", image_path, "\n")
cat("Numeric field:", field_path, "\n")
