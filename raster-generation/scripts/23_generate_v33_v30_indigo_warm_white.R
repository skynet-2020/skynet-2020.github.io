#!/usr/bin/env Rscript

# V33: preserve V30's numeric field and raster structure exactly while
# replacing viridis with a continuous indigo-to-warm-white palette.

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

input_path <- file.path(output_dir, "23_v25-smoother-variation-v30.rds")
if (!file.exists(input_path)) stop("The preserved V30 field is required: ", input_path)

v30 <- readRDS(input_path)
field <- v30$field
height <- v30$height
width <- v30$width
display_limits <- v30$parameters$reference_display_limits
center_color <- "#3B2F68"
outer_color <- "#E8E6DF"
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

render_fixed(field, file.path(output_dir, "26_v30-indigo-warm-white-v33.png"))
render_fixed(
  aggregate_grid(field, 90L, 160L),
  file.path(study_dir, "26_v33_small-cells_160x90.png")
)
render_fixed(
  aggregate_grid(field, 45L, 80L),
  file.path(study_dir, "26_v33_medium-cells_80x45.png")
)

v30$parameters$source_version <- "v30"
v30$parameters$palette_type <- "two-endpoint CIE Lab interpolation"
v30$parameters$center_color <- center_color
v30$parameters$outer_color <- outer_color
v30$parameters$recommended_text_color <- text_color
saveRDS(v30, file.path(output_dir, "26_v30-indigo-warm-white-v33.rds"))

cat("V33 V30 indigo-to-warm-white conversion complete\n")
cat("Center:", center_color, "\n")
cat("Exterior:", outer_color, "\n")
cat("Text:", text_color, "\n")
