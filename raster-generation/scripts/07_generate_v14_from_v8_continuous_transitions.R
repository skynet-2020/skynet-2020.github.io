#!/usr/bin/env Rscript

# V14: rebuild from preserved V8 with continuously changing transition lengths.
# There are no piecewise color shelves. A smooth accelerating warp lengthens the
# low and middle-low transitions; a C2-smooth return to identity preserves V8's
# lightest-yellow colors and footprint.

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
field_v8 <- source$field
limits <- source$parameters$reference_display_limits
if (is.null(limits)) stop("V8 does not contain fixed reference display limits.")

x <- (field_v8 - limits[1]) / diff(limits)
x[x < 0] <- 0
x[x > 1] <- 1

# Pass --v15 for a restrained increase in the same continuous effect. Both
# versions always start from V8, so V15 does not compound V14's transformation.
make_v15 <- "--v15" %in% commandArgs(trailingOnly = TRUE)

# Smooth convex mapping. The identity component guarantees a nonzero gradient
# through the darkest region; the exponential component makes transition speed
# accelerate continuously rather than changing at bands or thresholds.
curve_strength <- if (make_v15) 2.05 else 1.70
identity_share <- if (make_v15) 0.24 else 0.38
accelerated <- identity_share * x +
  (1 - identity_share) * expm1(curve_strength * x) / expm1(curve_strength)

# Preserve the exact upper register. quintic_smootherstep has zero first and
# second derivatives at both endpoints, so the return to V8 is visually smooth.
return_start <- 0.70
return_end <- 0.91
t <- (x - return_start) / (return_end - return_start)
t[t < 0] <- 0
t[t > 1] <- 1
upper_blend <- t^3 * (t * (t * 6 - 15) + 10)
remapped <- accelerated * (1 - upper_blend) + x * upper_blend
field_v14 <- limits[1] + remapped * diff(limits)

height <- nrow(field_v14)
width <- ncol(field_v14)
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

version <- if (make_v15) "v15" else "v14"
prefix <- if (make_v15) "08" else "07"
smooth_name <- sprintf("%s_v8-continuous-transitions-%s.png", prefix, version)
rds_name <- sprintf("%s_v8-continuous-transitions-%s.rds", prefix, version)

render_fixed(field_v14, file.path(output_dir, smooth_name))
render_fixed(
  aggregate_grid(field_v14, 90L, 160L),
  file.path(study_dir, sprintf("%s_%s_small-cells_160x90.png", prefix, version))
)
render_fixed(
  aggregate_grid(field_v14, 45L, 80L),
  file.path(study_dir, sprintf("%s_%s_medium-cells_80x45.png", prefix, version))
)

source$field <- field_v14
source$parameters$source_version <- "v8"
source$parameters$transition_method <- "smooth exponential acceleration with C2 upper return"
source$parameters$curve_strength <- curve_strength
source$parameters$identity_share <- identity_share
source$parameters$upper_return <- c(return_start, return_end)
saveRDS(source, file.path(output_dir, rds_name))

cat(toupper(version), " continuous-transition study complete\n", sep = "")
cat("Source: V8\n")
cat("Upper register identical to V8 from:", return_end, "to 1\n")
