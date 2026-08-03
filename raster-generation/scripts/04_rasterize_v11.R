#!/usr/bin/env Rscript

# Raster-cell studies derived from the preserved v11 numeric field.
#
# Outputs:
#   - three uniformly aggregated grids;
#   - adaptive quadtree with fine center and coarse margins;
#   - adaptive quadtree with coarse center and fine margins.
#
# The source values, viridis palette, and v11 display limits remain fixed.

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
output_dir <- file.path(raster_root, "outputs", "landscape", "raster-studies")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

input_path <- file.path(
  raster_root,
  "outputs",
  "landscape",
  "03k_radial-organic-midpoint-v11.rds"
)
if (!file.exists(input_path)) stop("The preserved v11 RDS is missing: ", input_path)

source <- readRDS(input_path)
field <- source$field
distance <- source$warped_distance
height <- nrow(field)
width <- ncol(field)

display_limits <- source$parameters$reference_display_limits
palette_begin <- source$parameters$palette_begin
palette_end <- source$parameters$palette_end

if (is.null(display_limits)) {
  stop("v11 does not contain its fixed reference display limits.")
}

render_fixed <- function(values, path) {
  display <- ambient::normalise(
    pmin(pmax(values, display_limits[1]), display_limits[2]),
    from = display_limits,
    to = c(0, 1)
  )
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
  png::writePNG(image, path)
}

aggregate_grid <- function(values, target_rows, target_cols) {
  row_breaks <- floor(seq(0, nrow(values), length.out = target_rows + 1L))
  col_breaks <- floor(seq(0, ncol(values), length.out = target_cols + 1L))
  coarse <- matrix(NA_real_, target_rows, target_cols)

  for (row in seq_len(target_rows)) {
    row_ids <- seq.int(row_breaks[row] + 1L, row_breaks[row + 1L])
    for (col in seq_len(target_cols)) {
      col_ids <- seq.int(col_breaks[col] + 1L, col_breaks[col + 1L])
      coarse[row, col] <- mean(values[row_ids, col_ids])
    }
  }

  output_rows <- pmin(
    target_rows,
    pmax(1L, ceiling(seq_len(nrow(values)) * target_rows / nrow(values)))
  )
  output_cols <- pmin(
    target_cols,
    pmax(1L, ceiling(seq_len(ncol(values)) * target_cols / ncol(values)))
  )
  coarse[output_rows, output_cols, drop = FALSE]
}

uniform_specs <- list(
  subtle = c(rows = 90L, cols = 160L),
  medium = c(rows = 45L, cols = 80L),
  bold = c(rows = 23L, cols = 40L)
)

for (name in names(uniform_specs)) {
  spec <- uniform_specs[[name]]
  aggregated <- aggregate_grid(field, spec[["rows"]], spec[["cols"]])
  path <- file.path(
    output_dir,
    sprintf(
      "04_uniform-%s_%dx%d.png",
      name,
      spec[["cols"]],
      spec[["rows"]]
    )
  )
  render_fixed(aggregated, path)
}

# Recursively subdivide aligned square blocks. Leaf sizes are powers of two, so
# transitions remain structured and read as a true multiresolution raster.
adaptive_quadtree <- function(values, distance_field, direction) {
  output <- matrix(NA_real_, nrow(values), ncol(values))

  desired_size <- function(d) {
    if (direction == "fine-center") {
      if (d < 0.55) return(8L)
      if (d < 0.95) return(16L)
      if (d < 1.30) return(32L)
      return(64L)
    }

    if (d < 0.55) return(64L)
    if (d < 0.95) return(32L)
    if (d < 1.30) return(16L)
    8L
  }

  fill_block <- function(row_start, col_start, size) {
    row_end <- min(row_start + size - 1L, nrow(values))
    col_end <- min(col_start + size - 1L, ncol(values))
    center_row <- floor((row_start + row_end) / 2)
    center_col <- floor((col_start + col_end) / 2)
    target <- desired_size(distance_field[center_row, center_col])

    if (size > target && size > 8L) {
      half <- as.integer(size / 2L)
      fill_block(row_start, col_start, half)
      if (col_start + half <= ncol(values)) {
        fill_block(row_start, col_start + half, half)
      }
      if (row_start + half <= nrow(values)) {
        fill_block(row_start + half, col_start, half)
      }
      if (
        row_start + half <= nrow(values) &&
          col_start + half <= ncol(values)
      ) {
        fill_block(row_start + half, col_start + half, half)
      }
      return(invisible(NULL))
    }

    block_mean <- mean(values[row_start:row_end, col_start:col_end])
    output[row_start:row_end, col_start:col_end] <<- block_mean
    invisible(NULL)
  }

  for (row_start in seq(1L, nrow(values), by = 64L)) {
    for (col_start in seq(1L, ncol(values), by = 64L)) {
      fill_block(row_start, col_start, 64L)
    }
  }

  output
}

adaptive_fine_center <- adaptive_quadtree(field, distance, "fine-center")
adaptive_coarse_center <- adaptive_quadtree(field, distance, "coarse-center")

render_fixed(
  adaptive_fine_center,
  file.path(output_dir, "04_adaptive-fine-center_coarse-margins.png")
)
render_fixed(
  adaptive_coarse_center,
  file.path(output_dir, "04_adaptive-coarse-center_fine-margins.png")
)

cat("Raster studies complete\n")
cat("Source: v11\n")
cat("Output directory:", output_dir, "\n")
