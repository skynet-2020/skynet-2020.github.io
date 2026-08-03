#!/usr/bin/env Rscript

# V35: preserve V25's numeric field and raster structure exactly while
# replacing viridis with a continuous warm-white-to-indigo palette.

if (!requireNamespace("png", quietly = TRUE)) {
  stop("Install required package: png")
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
field <- v25$field
height <- v25$height
width <- v25$width
display_limits <- v25$parameters$reference_display_limits
center_color <- "#E8E6DF"
outer_color <- "#3B2F68"
text_color <- "#A63446"

palette <- grDevices::colorRampPalette(
  c(center_color, outer_color),
  space = "Lab"
)(256)

render_fixed <- function(values, path) {
  display <- (values - display_limits[1]) / diff(display_limits)
  display <- pmax(0, pmin(1, display))
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

render_fixed(field, file.path(output_dir, "28_v25-warm-white-indigo-v35.png"))
render_fixed(
  aggregate_grid(field, 90L, 160L),
  file.path(study_dir, "28_v35_small-cells_160x90.png")
)
render_fixed(
  aggregate_grid(field, 45L, 80L),
  file.path(study_dir, "28_v35_medium-cells_80x45.png")
)

v25$parameters$source_version <- "v25"
v25$parameters$palette_type <- "two-endpoint CIE Lab interpolation"
v25$parameters$center_color <- center_color
v25$parameters$outer_color <- outer_color
v25$parameters$recommended_text_color <- text_color
saveRDS(v25, file.path(output_dir, "28_v25-warm-white-indigo-v35.rds"))

cat("V35 V25 warm-white-to-indigo conversion complete\n")
cat("Center:", center_color, "\n")
cat("Exterior:", outer_color, "\n")
cat("Text:", text_color, "\n")
