#!/usr/bin/env Rscript

# V13: targeted transition redistribution based on V11.
# Preserve the darkest core, spatially stretch violet -> blue-purple and
# blue-purple -> teal, then accelerate the remaining green/yellow register.

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

input_path <- file.path(output_dir, "03k_radial-organic-midpoint-v11.rds")
if (!file.exists(input_path)) stop("The preserved V11 field is missing: ", input_path)

source <- readRDS(input_path)
field_v11 <- source$field
limits <- source$parameters$reference_display_limits
if (is.null(limits)) stop("V11 does not contain fixed reference display limits.")

normalized <- (field_v11 - limits[1]) / diff(limits)
normalized[normalized < 0] <- 0
normalized[normalized > 1] <- 1

# Input and output anchors in normalized viridis space. The initial segment is
# unchanged (slope 1), preventing expansion of the flattest dark core. The two
# middle-low segments have shallow slopes and therefore occupy much more of the
# canvas. The upper register catches up rapidly and still reaches full yellow.
input_knots  <- c(0.00, 0.15, 0.55, 0.80, 1.00)
output_knots <- c(0.00, 0.15, 0.33, 0.52, 1.00)
remapped <- matrix(
  stats::approx(input_knots, output_knots, xout = as.vector(normalized), ties = "ordered")$y,
  nrow = nrow(normalized),
  ncol = ncol(normalized)
)
field_v13 <- limits[1] + remapped * diff(limits)

height <- nrow(field_v13)
width <- ncol(field_v13)
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

render_fixed(field_v13, file.path(output_dir, "06_piecewise-transitions-v13.png"))
render_fixed(
  aggregate_grid(field_v13, 90L, 160L),
  file.path(study_dir, "06_v13_small-cells_160x90.png")
)
render_fixed(
  aggregate_grid(field_v13, 45L, 80L),
  file.path(study_dir, "06_v13_medium-cells_80x45.png")
)

source$field <- field_v13
source$parameters$source_version <- "v11"
source$parameters$transition_method <- "anchored piecewise-linear transform"
source$parameters$transition_input_knots <- input_knots
source$parameters$transition_output_knots <- output_knots
saveRDS(source, file.path(output_dir, "06_piecewise-transitions-v13.rds"))

cat("V13 piecewise-transition study complete\n")
cat("Input knots: ", paste(input_knots, collapse = ", "), "\n", sep = "")
cat("Output knots: ", paste(output_knots, collapse = ", "), "\n", sep = "")
