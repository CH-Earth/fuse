#!/usr/bin/env Rscript
#
# make_elevation_bands.R
#
# Create a FUSE elevation-band NetCDF file from:
#   1. a DEM GeoTIFF;
#   2. a catchment polygon; and
#   3. a specified elevation-band width.
#
# Usage:
#
#   Rscript make_elevation_bands.R \
#       <dem.tif> \
#       <catchment.shp> \
#       <band_width_m> \
#       <output.nc>
#
# Example:
#
#   Rscript make_elevation_bands.R \
#       CAN_05BB001_dem.tif \
#       CAN_05BB001_lumped.shp \
#       100 \
#       CAN_05BB001_elev_bands.nc
#

library(terra)
library(ncdf4)

# ----------------------------------------------------------------------
# Read and validate command-line arguments
# ----------------------------------------------------------------------

debug <- TRUE

if (debug) {

  dem_file    <- "../CAN_07BB001/geospatial/tif/CAN_05BB001_merit_hydro_elv.tif"
  basin_file  <- "../CAN_07BB001/geospatial/shp/CAN_05BB001_lumped.shp"
  band_width  <- 100
  output_file <- "../CAN_07BB001/geospatial/bands/CAN_05BB001_elev_bands.nc"

} else {

  args <- commandArgs(trailingOnly = TRUE)
  
  if (length(args) != 4) {
    stop(
      paste(
        "Usage:",
        "Rscript make_elevation_bands.R",
        "<dem.tif> <catchment.shp> <band_width_m> <output.nc>"
      )
    )
  }
  
  dem_file <- args[1]
  basin_file <- args[2]
  band_width <- as.numeric(args[3])
  output_file <- args[4]
  
  if (!file.exists(dem_file)) {
    stop("DEM file not found: ", dem_file)
  }
  
  if (!file.exists(basin_file)) {
    stop("Catchment shapefile not found: ", basin_file)
  }
  
  if (!is.finite(band_width) || band_width <= 0) {
    stop("The elevation-band width must be a positive number.")
  }

} # debug switch

# ----------------------------------------------------------------------
# Read the DEM and catchment polygon
# ----------------------------------------------------------------------

dem <- rast(dem_file)
basin <- vect(basin_file)

if (nlyr(dem) != 1) {
  stop("The DEM must contain exactly one raster layer.")
}

# Transform the basin boundary to the DEM coordinate reference system.
if (!same.crs(dem, basin)) {
  basin <- project(basin, crs(dem))
}

# ----------------------------------------------------------------------
# Clip and mask the DEM to the catchment
# ----------------------------------------------------------------------

dem_basin <- crop(dem, basin)
dem_basin <- mask(dem_basin, basin)

valid_cells <- global(!is.na(dem_basin), "sum", na.rm = TRUE)[1, 1]

if (!is.finite(valid_cells) || valid_cells == 0) {
  stop("The catchment does not overlap any valid DEM cells.")
}

# ----------------------------------------------------------------------
# Define fixed-width elevation bands
# ----------------------------------------------------------------------

elev_min <- global(dem_basin, "min", na.rm = TRUE)[1, 1]
elev_max <- global(dem_basin, "max", na.rm = TRUE)[1, 1]

lower_limit <- floor(elev_min / band_width) * band_width
upper_limit <- ceiling(elev_max / band_width) * band_width

# Ensure that the maximum DEM value falls within the final interval.
if (upper_limit <= elev_max) {
  upper_limit <- upper_limit + band_width
}

breaks <- seq(
  from = lower_limit,
  to = upper_limit,
  by = band_width
)

number_of_bands <- length(breaks) - 1

# terra classification matrix:
# lower bound, upper bound, resulting band number.
classification <- cbind(
  breaks[-length(breaks)],
  breaks[-1],
  seq_len(number_of_bands)
)

band_raster <- classify(
  dem_basin,
  classification,
  include.lowest = TRUE,
  right = FALSE,
  others = NA
)

# Explicitly place values equal to the final upper boundary in the last
# band, although this should rarely be necessary.
band_raster <- ifel(
  dem_basin == upper_limit,
  number_of_bands,
  band_raster
)

# ----------------------------------------------------------------------
# Calculate cell areas
# ----------------------------------------------------------------------
#
# cellSize() gives appropriate cell areas for both projected and
# longitude-latitude rasters.

cell_area <- cellSize(dem_basin, unit = "m")

band_results <- vector("list", number_of_bands)

for (band_id in seq_len(number_of_bands)) {

  in_band <- band_raster == band_id

  area_m2 <- global(
    ifel(in_band, cell_area, NA),
    "sum",
    na.rm = TRUE
  )[1, 1]

  mean_elevation <- global(
    ifel(in_band, dem_basin, NA),
    "mean",
    na.rm = TRUE
  )[1, 1]

  band_results[[band_id]] <- data.frame(
    elevation_band = band_id,
    lower_elev = breaks[band_id],
    upper_elev = breaks[band_id + 1],
    area_m2 = area_m2,
    mean_elev = mean_elevation
  )
}

band_table <- do.call(rbind, band_results)

# Remove elevation intervals containing no DEM cells.
band_table <- band_table[
  is.finite(band_table$area_m2) &
    band_table$area_m2 > 0,
  ,
  drop = FALSE
]

if (nrow(band_table) == 0) {
  stop("No non-empty elevation bands were produced.")
}

# Renumber after removing empty bands.
band_table$elevation_band <- seq_len(nrow(band_table))

band_table$area_frac <- (
  band_table$area_m2 / sum(band_table$area_m2)
)

# The existing FUSE input format uses precipitation fractions equal to
# the corresponding elevation-band area fractions.
band_table$prec_frac <- band_table$area_frac

# Check that fractions sum to one within numerical precision.
if (abs(sum(band_table$area_frac) - 1) > 1e-10) {
  stop("Elevation-band area fractions do not sum to one.")
}

# Force exact closure (put roundoff error in the last band)
band_table$area_frac[nrow(band_table)] <-
    1 - sum(band_table$area_frac[-nrow(band_table)])

# ----------------------------------------------------------------------
# Determine representative catchment longitude and latitude
# ----------------------------------------------------------------------
#
# Use the catchment centroid, transformed to geographic coordinates.

basin_lonlat <- project(basin, "EPSG:4326")
basin_centroid <- centroids(aggregate(basin_lonlat))

centroid_coordinates <- crds(basin_centroid)

longitude <- centroid_coordinates[1, 1]
latitude <- centroid_coordinates[1, 2]

# ----------------------------------------------------------------------
# Define NetCDF dimensions
# ----------------------------------------------------------------------

longitude_dim <- ncdim_def(
  name = "longitude",
  units = "degreesE",
  vals = longitude,
  longname = "longitude",
  create_dimvar = TRUE
)

latitude_dim <- ncdim_def(
  name = "latitude",
  units = "degreesN",
  vals = latitude,
  longname = "latitude",
  create_dimvar = TRUE
)

elevation_band_dim <- ncdim_def(
  name = "elevation_band",
  units = "-",
  vals = band_table$elevation_band,
  longname = "elevation_band",
  create_dimvar = TRUE
)

# ----------------------------------------------------------------------
# Define NetCDF variables
# ----------------------------------------------------------------------
#
# Define the variable dimensions.
dims <- list(
  longitude_dim,
  latitude_dim,
  elevation_band_dim
)

area_frac_var <- ncvar_def(
  name = "area_frac",
  units = "-",
  dim = dims,
  missval = NA_real_,
  longname = "Fraction of the catchment covered by each elevation band",
  prec = "double"
)

mean_elev_var <- ncvar_def(
  name = "mean_elev",
  units = "m asl",
  dim = dims,
  missval = NA_real_,
  longname = "Mean elevation of each elevation band",
  prec = "double"
)

prec_frac_var <- ncvar_def(
  name = "prec_frac",
  units = "-",
  dim = dims,
  missval = NA_real_,
  longname = paste(
    "Fraction of catchment precipitation that falls on each",
    "elevation band - same as area_frac"
  ),
  prec = "double"
)

# ----------------------------------------------------------------------
# Create and populate the NetCDF file
# ----------------------------------------------------------------------

nc <- nc_create(
  filename = output_file,
  vars = list(
    area_frac_var,
    mean_elev_var,
    prec_frac_var
  ),
  force_v4 = FALSE
)

number_of_output_bands <- nrow(band_table)

area_array <- array(
  band_table$area_frac,
  dim = c(number_of_output_bands, 1, 1)
)

mean_array <- array(
  band_table$mean_elev,
  dim = c(number_of_output_bands, 1, 1)
)

prec_array <- array(
  band_table$prec_frac,
  dim = c(number_of_output_bands, 1, 1)
)

ncvar_put(nc, "area_frac", area_array)
ncvar_put(nc, "mean_elev", mean_array)
ncvar_put(nc, "prec_frac", prec_array)

# ----------------------------------------------------------------------
# Add global attributes
# ----------------------------------------------------------------------

creation_time <- format(
  Sys.time(),
  format = "Created %Y/%m/%d %H:%M:%S"
)

ncatt_put(
  nc,
  varid = 0,
  attname = "author",
  attval = "Martyn Clark"
)

ncatt_put(
  nc,
  varid = 0,
  attname = "date",
  attval = creation_time
)

ncatt_put(
  nc,
  varid = 0,
  attname = "institution",
  attval = "University of Calgary"
)

ncatt_put(
  nc,
  varid = 0,
  attname = "source_dem",
  attval = basename(dem_file)
)

ncatt_put(
  nc,
  varid = 0,
  attname = "catchment_boundary",
  attval = basename(basin_file)
)

ncatt_put(
  nc,
  varid = 0,
  attname = "elevation_band_width_m",
  attval = band_width
)

nc_close(nc)

# Prevent on.exit() from trying to close an already closed object.
on.exit(NULL, add = FALSE)

# ----------------------------------------------------------------------
# Report summary
# ----------------------------------------------------------------------

cat("Created:", output_file, "\n")
cat("Catchment longitude:", longitude, "\n")
cat("Catchment latitude: ", latitude, "\n")
cat("DEM elevation range:", elev_min, "to", elev_max, "m\n")
cat("Elevation-band width:", band_width, "m\n")
cat("Number of non-empty bands:", number_of_output_bands, "\n")
cat(
  "Represented catchment area:",
  sum(band_table$area_m2) / 1e6,
  "km2\n"
)
cat(
  "Sum of area fractions:",
  format(sum(band_table$area_frac), digits = 16),
  "\n"
)
