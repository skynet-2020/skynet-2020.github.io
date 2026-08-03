#!/usr/bin/env Rscript

# V28: preserve V25's numeric field and replace viridis with a continuous,
# perceptually interpolated indigo-to-cherry palette.

required_packages <- c("png")
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

input_path <- file.path(output_dir, "18_v21-core-v19-contours-v25.rds")
if (!file.exists(input_path)) stop("The preserved V25 field is required: ", input_path)

v25 <- readRDS(input_path)
field_v28 <- v25$field
height <- v25$height
width <- v25$width
display_limits <- v25$parameters$reference_display_limits
center_color <- "#3B2F68"
outer_color <- "#A63446"

# CIE Lab interpolation produces a continuous perceptual transition rather than
# independently mixing display RGB channels.
palette <- grDevices::colorRampPalette(
  c(center_color, outer_color),
  space = "Lab"
)(256)

render_fixed <- function(values, path) {
  display <- (values - display_limits[1]) / diff(display_limits)
  display[display < 0] <- 0
  display[display > 1] <- 1
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

render_fixed(field_v28, file.path(output_dir, "21_v25-indigo-cherry-v28.png"))
render_fixed(
  aggregate_grid(field_v28, 90L, 160L),
  file.path(study_dir, "21_v28_small-cells_160x90.png")
)
render_fixed(
  aggregate_grid(field_v28, 45L, 80L),
  file.path(study_dir, "21_v28_medium-cells_80x45.png")
)

v25$parameters$source_version <- "v25"
v25$parameters$palette_type <- "two-endpoint CIE Lab interpolation"
v25$parameters$center_color <- center_color
v25$parameters$outer_color <- outer_color
v25$parameters$recommended_text_color <- "#E8E6DF"
saveRDS(v25, file.path(output_dir, "21_v25-indigo-cherry-v28.rds"))

cat("V28 indigo-to-cherry study complete\n")
cat("Center:", center_color, "\n")
cat("Exterior:", outer_color, "\n")
