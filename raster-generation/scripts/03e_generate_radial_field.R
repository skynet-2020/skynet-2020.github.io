#!/usr/bin/env Rscript

# Alternative Layer 3 v5: predominantly radial composition
#
# Radial organization controls the composition. Broad spatial warping, residual
# low-frequency heterogeneity, and restrained outward dark extensions prevent a
# perfect bullseye. Bright values are explicitly biased toward the corners.

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

# Very broad noise controls unequal transition distances around the center.
set.seed(layer_1$seed + 351L)
radial_warp_large <- ambient::noise_simplex(
  c(height, width), frequency = 0.00105, fractal = "none"
)
set.seed(layer_1$seed + 352L)
radial_warp_support <- ambient::noise_simplex(
  c(height, width), frequency = 0.0019, fractal = "none"
)
radial_warp <- 0.8 * radial_warp_large + 0.2 * radial_warp_support
radial_warp <- ambient::normalise(radial_warp, to = c(-1, 1))

# The aspect-aware distance makes a wide central text region. Warping perturbs
# the contours without allowing heterogeneity to replace radial organization.
radius_x <- 0.68
radius_y <- 0.82
base_distance <- sqrt(
  (x_grid / radius_x)^2 +
    (y_grid / radius_y)^2
)
warped_distance <- base_distance * exp(0.13 * radial_warp)

# Smooth radial rise: broad and dark near the center, progressively higher
# outward, with no hard core boundary.
radial_component <- warped_distance^1.32
radial_component <- ambient::normalise(
  pmin(radial_component, stats::quantile(radial_component, 0.995, names = FALSE))
)

# This term is high only when both |x| and |y| are large. It concentrates the
# brightest values in corners rather than creating bright horizontal/vertical
# edge bands.
corner_component <- (abs(x_grid) * abs(y_grid))^1.45
corner_component <- ambient::normalise(corner_component)

# A small amount of the original broad field restores non-monotonic spatial
# variation. Its influence fades almost completely inside the text region.
calm_radius <- 0.62
calm_weight <- exp(-(warped_distance / calm_radius)^2.4)
heterogeneity_gate <- 1 - 0.94 * calm_weight
broad_component <- ambient::normalise(broad_field, to = c(-1, 1))

# Restrained broad dark extensions preserve some of 03b's outward-reaching
# character without appearing as six discrete arms.
set.seed(layer_1$seed + 353L)
extension_noise <- ambient::noise_simplex(
  c(height, width),
  frequency = 0.00155,
  fractal = "fbm",
  octaves = 2,
  lacunarity = 1.8,
  gain = 0.4
)
extension_noise <- ambient::normalise(extension_noise, to = c(-1, 1))
dark_extensions <- pmax(extension_noise, 0) * heterogeneity_gate

# Approximate compositional balance: 76% radial/corner organization and 24%
# broad irregularity. Corner emphasis is part of the radial organization.
field <-
  0.62 * radial_component +
  0.14 * corner_component +
  0.18 * broad_component * heterogeneity_gate -
  0.06 * dark_extensions

# Further compress residual variation at the center without creating a mask:
# the same infinite-tailed calm weight simply reduces local amplitude.
central_floor <- stats::quantile(field[warped_distance < 0.18], 0.2, names = FALSE)
central_compression <- 0.94
central_darkening <- 0.11
field <- central_floor +
  (field - central_floor) * (1 - central_compression * calm_weight) -
  central_darkening * calm_weight

# Preserve v7's display scale before modifying the center. Using these fixed
# limits for v8 guarantees that unchanged outer values retain identical colors.
reference_display_limits <- stats::quantile(
  field,
  c(0.01, 0.99),
  names = FALSE
)

# Establish a genuinely broad dark interior, then compress the color transition
# into a more peripheral band. At the nominal radii below, expansion_inner spans
# about 60% of canvas width and 71% of canvas height. The adjustment becomes
# exactly zero before the outer field and corners.
expansion_inner <- 0.66
expansion_outer <- 1.16
expansion_strength <- 0.72
transition_position <- pmin(
  1,
  pmax(
    0,
    (warped_distance - expansion_inner) /
      (expansion_outer - expansion_inner)
  )
)
smoothstep <- transition_position^2 * (3 - 2 * transition_position)
expansion_weight <- expansion_strength * (1 - smoothstep)
expanded_dark_target <- stats::quantile(
  field[warped_distance < 0.3],
  0.25,
  names = FALSE
)
field <- field * (1 - expansion_weight) +
  expanded_dark_target * expansion_weight

# Four broad directional arms extend the low-value interior beyond the visible
# canvas. They have no in-frame endpoints: once established, each ridge remains
# active until it is cropped by an edge. Unequal angles, widths, curvature, and
# strength prevent a symmetric spoke or starburst pattern.
arm_x <- x_grid + 0.08 * radial_warp
arm_y <- y_grid + 0.08 * extension_noise
arm_specs <- data.frame(
  angle = c(7, 84, 187, 276) * pi / 180,
  width = c(0.25, 0.22, 0.27, 0.23),
  curvature = c(0.07, -0.09, 0.06, -0.08),
  strength = c(0.72, 0.62, 0.68, 0.64)
)

arm_union <- matrix(0, nrow = height, ncol = width)

# Robust color scaling retains nearly full viridis, including controlled yellow
# at the corner maxima.
limits <- reference_display_limits
display <- ambient::normalise(pmin(pmax(field, limits[1]), limits[2]))
palette_begin <- 0
palette_end <- 0.96
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

field_path <- file.path(output_dir, "03k_radial-organic-midpoint-v11.rds")
image_path <- file.path(output_dir, "03k_radial-organic-midpoint-v11.png")
saveRDS(
  list(
    field = field,
    radial_component = radial_component,
    corner_component = corner_component,
    broad_component = broad_component,
    dark_extensions = dark_extensions,
    warped_distance = warped_distance,
    seed = layer_1$seed,
    width = width,
    height = height,
    parameters = list(
      radius_x = radius_x,
      radius_y = radius_y,
      calm_radius = calm_radius,
      central_compression = central_compression,
      central_darkening = central_darkening,
      expansion_inner = expansion_inner,
      expansion_outer = expansion_outer,
      expansion_strength = expansion_strength,
      arm_specs = arm_specs,
      reference_display_limits = reference_display_limits,
      radial_warp_strength = 0.13,
      radial_exponent = 1.32,
      radial_weight = 0.62,
      corner_weight = 0.14,
      broad_weight = 0.18,
      extension_weight = 0.06,
      palette_begin = palette_begin,
      palette_end = palette_end
    )
  ),
  field_path
)
png::writePNG(image, image_path)

cat("Alternative Layer 3 radial v11 complete\n")
cat("Image:", image_path, "\n")
cat("Numeric field:", field_path, "\n")
