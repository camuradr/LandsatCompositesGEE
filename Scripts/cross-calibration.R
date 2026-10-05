
### Landsat sensors cross-calibration


library(pacman)
library(reticulate)
library(ggplot2)

# Load packages for data handling etc.
library(sf)
library(dplyr)
library(tidyr)
library(purrr)
library(data.table)
library(stringr)
library(rgee)

# Load LandsatTS package
library(LandsatTS)
library(jsonlite)  # For parsing the coordinate string
library(terra)
library(raster)

## 1. Set up google earth engine environment

rgee_environment_dir = "C:\\Users\\camr1n24\\AppData\\Local\\ESRI\\conda\\envs\\rgee_py\\"
Sys.setenv(RETICULATE_PYTHON = rgee_environment_dir)
Sys.setenv(EARTHENGINE_PYTHON = rgee_environment_dir)

# Load the Earth Engine Python API
#ee$Authenticate(auth_mode='notebook')
#ee$Initialize(project='ee-leafborealforests')  # <-- EDIT THIS FOR YOUR PROJECT
#ee$String('Hello from the Earth Engine servers!')$getInfo()

# Intialize the Earth Engine with rgee
ee_Initialize(email = "cesaramuradr1n24@gmail.com", drive = TRUE)
#ee_Initialize(email = "cesaramuradr@gmail.com", drive = TRUE)

ee_Authenticate()

# clean environment
rm(list = ls()); gc()

## 2. Read and load sample polygons COORDINATES

#polygon_coordinates <- read.csv("C:/Users/camr1n24/Documents/LandsatTS/Sample_Polygon_Coordinates.csv", header = TRUE, stringsAsFactors = FALSE)

samples_path <- "G:/My Drive/GEE_exports/PolygonSamples/StableForestSamples"
samples_files <- list.files(samples_path,
                         pattern = "\\.csv",
                         full.names = TRUE)

samples_csv_list <- lapply(samples_files, read.csv)

polygon_coordinates <- bind_rows(samples_csv_list)

polygon_coordinates$coords_list <- lapply(polygon_coordinates$coordinates, fromJSON)

# Flatten the coordinates
flat_coords <- polygon_coordinates %>%
  mutate(poly_id = row_number()) %>%
  rowwise() %>%
  mutate(coords = list(matrix(unlist(coords_list), ncol = 10, byrow = TRUE))) %>%
  ungroup()

# Expand to long format
long_coords <- flat_coords %>%
  dplyr::select(poly_id, coords) %>%
  rowwise() %>%
  mutate(coords_df = list(as.data.frame(coords))) %>%
  unnest(cols = c(coords_df)) %>%
  rename(x = V1, y = V2)

#Rename columns
colnames(long_coords) <- c("ID","Coordinates","X1", "X2","X3","X4","X5","Y1","Y2","Y3","Y4","Y5") 

# Initialize an empty list to store polygons
polygons_list <- list()

# Loop over each row (polygon) in long_coords
for (i in 1:nrow(long_coords)) {
  row <- long_coords[i, ]
  
  # Get all columns with names starting with X or Y
  x_cols <- grep("^X\\d+$", names(row), value = TRUE)
  y_cols <- grep("^Y\\d+$", names(row), value = TRUE)
  
  # Sort them to ensure the order is correct (X1, X2,...)
  x_vals <- as.numeric(unlist(row[ , x_cols[order(x_cols)]]))
  y_vals <- as.numeric(unlist(row[ , y_cols[order(y_cols)]]))
  
  # Combine into coordinate matrix
  coords <- matrix(c(x_vals, y_vals), ncol = 2)
  
  # Create polygon
  poly <- st_polygon(list(coords))
  
  # Add to list
  polygons_list[[i]] <- poly
}

# Convert to an sf object
sf_polygons <- st_sf(id = 1:length(polygons_list), geometry = st_sfc(polygons_list), crs = 4326)
#sf_polygons_subset <- sf_polygons[sf_polygons$id %in% c(1, 2, 3, 4, 5), ] # Select first 5 polygons

# Plot to verify
shapefile_path <- "C:/Users/camr1n24/Documents/SpatialAutocorrelation/Shapefiles/IntactForest_Quebec_ShieldTaiga.shp"
ShieldTaiga <- vect(shapefile_path)
plot(ShieldTaiga, border = "red", lwd = 1)
plot(sf_polygons, add = TRUE, col = "blue", border = "blue", lwd = 0.5)

# Make sure sf_polygons is of class sf with valid CRS
# sf_polygons <- st_read(...) or already created

# Create a list to store results
pixel_list_all <- list()

# Loop through each polygon
for (i in 1:nrow(sf_polygons)) {
  message("Processing polygon ", i)
  
  # Extract single polygon as sf object
  single_poly <- sf_polygons[i, ]
  
  # Get pixel centers within that polygon
  pixel_centers <- lsat_get_pixel_centers(single_poly, lsat_WRS2_scene_bounds = "C:/Users/camr1n24/Documents/LandsatTS/WRS-2_bound_world_0.kml")
  
  # Store result
  pixel_list_all[[i]] <- pixel_centers
}

# Save the list
saveRDS(pixel_list_all, file = "pixel_list.rds")

# Later you can load it back:
loaded_pixel_list <- readRDS("pixel_list.rds")



# 3. Export time-series using lsat_export_ts()

split_list_in_parts <- function(input_list, split_size) {
  n <- length(input_list)
  indices <- seq(1, n, by = split_size)
  split_lists <- lapply(indices, function(start_idx) {
    end_idx <- min(start_idx + split_size - 1, n)
    input_list[start_idx:end_idx]
  })
  return(split_lists)
}

pixel_sublists <- split_list_in_parts(pixel_list_all, 90)

#Select input list to process

index <- 5  # Change this to select a different sublist
input_pixel_list <- pixel_sublists[[index]]

# Loop through each pixel list
for (i in seq_along(input_pixel_list)) {
  message("Exporting time series for polygon ", i)
  
  # Extract single polygon as sf object
  single_list <- pixel_list_all[[i]]
  
  # Export time-series using lsat_export_ts()
  task_list <- lsat_export_ts(single_list, start_doy = 152, end_doy = 243, start_date = "1984-01-01", end_date = "2024-12-31", buffer_dist = 30,
                              drive_export_dir = "lsatTS_boreal", file_prefix = "lsatTS_export")
  
}

# Load csv file

pixel_samples_path <- "D:/GEE/PixelDataset(SecondRunFromPolygons)/CompletePixelSamples"
pixel_samples <- list.files(pixel_samples_path, pattern = "\\.csv$", full.names = TRUE)

#Read each csv file, count unique sample_id and add to total count (Check total number f)
total_pixels <- 0
for (i in seq_along(pixel_samples)) {
  cat("Processing file", i, "\n")
  # Read one file
  dt <- fread(pixel_samples[[i]])
  
  # Get unique sample_ids and how many there are
  unique_ids <- unique(dt$sample_id)
  n_pixels <- length(unique_ids)
  
  # Print the number of unique sample_ids
  total_pixels <- total_pixels + n_pixels
}

## Update sample_id values so that each pixel across all files gets a unique global ID. All rows corresponding to the same pixel (within a file) keep the same ID.Then merge all files into one final data frame (lsat.dt).
### Requirements:
### Rename files (ctrl+a to selelect all and then rename) and place them on a single folder

pixel_offset <- 0 # Initialize global counter
result_list <- vector("list", length(pixel_samples))

for (i in seq_along(pixel_samples)) {
  # Read one file
  dt <- fread(pixel_samples[[i]])
  
  # Get unique sample_ids and how many there are
  unique_ids <- unique(dt$sample_id)
  n_pixels <- length(unique_ids)
  
  # Create new global pixel IDs
  new_ids <- paste0("pixel_", seq(pixel_offset + 1, pixel_offset + n_pixels))
  
  # Create a mapping: old_id -> new_id
  id_map <- setNames(new_ids, unique_ids)
  
  # Replace sample_id with new global IDs
  dt$sample_id <- id_map[dt$sample_id]
  
  # Save and update pixel_offset
  result_list[[i]] <- dt
  pixel_offset <- pixel_offset + n_pixels
}

#### Combine all data frames into one

lsat.dt <- rbindlist(result_list)
# Save the RDS file
dt_folder_path <- "D:/GEE/PixelDataset(SecondRunFromPolygons)/CompletePixelSamples"
saveRDS(lsat.dt, file = file.path(dt_folder_path, "lsatdt.rds"))
lsat.dt <- readRDS(file.path(dt_folder_path, "lsatdt.rds"))

n_total <- nrow(lsat.dt) # Check number of rows in the data table
n_subsets <- 5 #number of subsets
subset_indices <- split(1:n_total, ceiling(seq_along(1:n_total) / ceiling(n_total / n_subsets)))

output_dir <- "D:/GEE/PixelDataset(SecondRunFromPolygons)"
#lsat.dt <- lsat_format_data(lsat.dt)
for (i in seq_len(n_subsets)) {
  message("Processing subset ", i, " of ", n_subsets)
  
  rows <- subset_indices[[i]]
  dtsubset <- lsat.dt[rows, ]  # This subset will never overlap or miss rows
  
  # Format and save the subset
  formatted_subset <- lsat_format_data(dtsubset)
  saveRDS(formatted_subset, file = file.path(output_dir, paste0("formatted_subset_", i, ".rds")))
  
  # Clean up memory
  rm(dtsubset, formatted_subset); gc()
  
  # Remove R's temporary files to free disk space
  temp_files <- list.files(tempdir(), full.names = TRUE, recursive = TRUE)
  unlink(list.files(tempdir(), full.names = TRUE, recursive = TRUE), recursive = TRUE, force = TRUE)
}

formatted_files <- list.files(output_dir, pattern = "\\.rds$", full.names = TRUE)

lsat.dt1 <- readRDS(file.path(output_dir, "formatted_subset_1.rds"))
lsat.dt2 <- readRDS(file.path(output_dir, "formatted_subset_2.rds"))
lsat.dt3 <- readRDS(file.path(output_dir, "formatted_subset_3.rds"))
lsat.dt4 <- readRDS(file.path(output_dir, "formatted_subset_4.rds"))
lsat.dt5 <- readRDS(file.path(output_dir, "formatted_subset_5.rds"))

formatted_data <- do.call(rbind, lapply(formatted_files, readRDS))
saveRDS(formatted_data, file = file.path(output_dir, "formattedlsatdt.rds"))
formatted_data <- readRDS(file.path(output_dir, "formattedlsatdt.rds"))
pixels<- unique(formatted_data$sample.id)
cleanlsat.dt <- lsat_clean_data(formatted_data, geom.max = 15, cloud.max = 80, sza.max = 60, filter.cfmask.snow = T, filter.cfmask.water = T, filter.jrc.water = T)
saveRDS(cleanlsat.dt, file = file.path(output_dir, "cleanformattedlsatdt.rds"))

# #######Subset tests#######
# testlsat.dt1 <- lsat_format_data(testlsat.dt1)
# testlsat.dt1 <- lsat_clean_data(testlsat.dt1, geom.max = 15, cloud.max = 80, sza.max = 60, filter.cfmask.snow = T, filter.cfmask.water = T, filter.jrc.water = T)
# 
lsat.dt1 <- lsat_format_data(lsat.dt1)
clsat.dt1 <- lsat_clean_data(lsat.dt1, geom.max = 15, cloud.max = 80, sza.max = 60, filter.cfmask.snow = T, filter.cfmask.water = T, filter.jrc.water = T)
saveRDS(clsat.dt1, file = file.path(output_dir, "clsatdt1.rds"))
clsat.dt1 <- readRDS(file.path(output_dir, "clsatdt1.rds"))

lsat.dt2 <- lsat_format_data(lsat.dt2)
clsat.dt2 <- lsat_clean_data(lsat.dt2, geom.max = 15, cloud.max = 80, sza.max = 60, filter.cfmask.snow = T, filter.cfmask.water = T, filter.jrc.water = T)
saveRDS(clsat.dt2, file = file.path(output_dir, "clsatdt2.rds"))
clsat.dt2 <- readRDS(file.path(output_dir, "clsatdt2.rds"))

lsat.dt3 <- lsat_format_data(lsat.dt3)
clsat.dt3 <- lsat_clean_data(lsat.dt3, geom.max = 15, cloud.max = 80, sza.max = 60, filter.cfmask.snow = T, filter.cfmask.water = T, filter.jrc.water = T)
saveRDS(clsat.dt3, file = file.path(output_dir, "clsatdt3.rds"))
clsat.dt3 <- readRDS(file.path(output_dir, "clsatdt3.rds"))

lsat.dt4 <- lsat_format_data(lsat.dt4)
clsat.dt4 <- lsat_clean_data(lsat.dt4, geom.max = 15, cloud.max = 80, sza.max = 60, filter.cfmask.snow = T, filter.cfmask.water = T, filter.jrc.water = T)
saveRDS(clsat.dt4, file = file.path(output_dir, "clsatdt4.rds"))
clsat.dt4 <- readRDS(file.path(output_dir, "clsatdt4.rds"))

lsat.dt5 <- lsat_format_data(lsat.dt5)
clsat.dt5 <- lsat_clean_data(lsat.dt5, geom.max = 15, cloud.max = 80, sza.max = 60, filter.cfmask.snow = T, filter.cfmask.water = T, filter.jrc.water = T)
saveRDS(clsat.dt5, file = file.path(output_dir, "clsatdt5.rds"))
clsat.dt5 <- readRDS(file.path(output_dir, "clsatdt5.rds"))

#Combine all data tables into one
cleanlsat.dt <- rbind(clsat.dt1, clsat.dt2, clsat.dt3, clsat.dt4, clsat.dt5) #, clsat.dt5
#####################################

data.summary.dt <- lsat_summarize_data(cleanlsat.dt)
data.summary.dt


# Store all data frames in a list
lsat_list <- list(lsat.dt1, lsat.dt2, lsat.dt3, lsat.dt4, lsat.dt5) #, lsat.dt5

# Combine all cleaned data frames
nlsat.dt <- do.call(rbind, lsat_list)

# Compute NDVI or other vegetation index
cleanlsat.dt <- lsat_calc_spectral_index(cleanlsat.dt, si = 'ndvi')

# Cross-calibrate NDVI among sensors using an approach based on Random Forest machine learning
lsat.dt <- lsat_calibrate_rf(lsat.dt, 
                             band.or.si = 'evi', 
                             doy.rng = 151:242, #01 June to 31 August
                             min.obs = 5, 
                             frac.train = 0.75,
                             train.with.highlat.data = F,
                             overwrite.col = F, 
                             write.output = F)

cleanlsat.dt <- lsat_calibrate_poly(cleanlsat.dt, 
                               band.or.si = 'nir', 
                               doy.rng = 151:242, 
                               min.obs = 5, #Minimum number of paired, seasonally-matched observations from Landsat 7 and Landsat 5/8 required to include a sampling sample.
                               frac.train = 0.75,# Fraction of the data used for training the model. The rest is used for validation.
                               train.with.highlat.data = F,
                               overwrite.col = F, 
                               write.output = F)

saveRDS(cleanlsat.dt, file = file.path(output_dir, "cleanformattedlsatdtxcalv2.rds"))
cleanlsat.dt <- readRDS(file.path("D:/GEE/PixelDataset", "cleanformattedlsatdt.rds"))

########### Select 10,000 random sample locations #############

# Remove non-masked years
fcleanlsat.dt <- cleanlsat.dt %>%
  filter(year >= 1985 & year <= 2020)

# Get max and min values from column
max_val <- max(cleanlsat.dt$nir, na.rm = TRUE)
min_val <- min(cleanlsat.dt$nir, na.rm = TRUE)

n_samples <- 10000  
unique_ids <- unique(fcleanlsat.dt$sample.id) # Get unique sample IDs

set.seed(555)  # Randomly sample n IDs
sampled_ids <- sample(unique_ids, n_samples)

sampled_data <- fcleanlsat.dt[fcleanlsat.dt$sample.id %in% sampled_ids, ]
unique_pixels <- unique(sampled_data$sample.id) #Get a List of Unique Pixels

sampled_data <- lsat_calibrate_poly(sampled_data, 
                                    band.or.si = 'nir', 
                                    doy.rng = 151:242, 
                                    min.obs = 5, #Minimum number of paired, seasonally-matched observations from Landsat 7 and Landsat 5/8 required to include a sampling sample.
                                    frac.train = 0.75,# Fraction of the data used for training the model. The rest is used for validation.
                                    train.with.highlat.data = F,
                                    trim = F,
                                    overwrite.col = F, 
                                    write.output = F)

#Remove columns
sampled_data <- sampled_data %>%
  select(-c("red.xcal.1"))

nirvalues <- cleanlsat.dt[, c("sample.id", "nir", "satellite", "ndvi")]
#boxplot of nir values
ggplot(nirvalues, aes(x = satellite, y = nir)) +
  geom_boxplot() +
  labs(title = "NIR Values by Satellite", x = "Satellite", y = "NIR Value") +
  theme_minimal()

# Calculate yearly mean for both NDVI and calibrated NDVI
ndvi_summary <- sampled_data %>%
  group_by(year) %>%
  summarise(
    mean_ndvi = mean(ndvi, na.rm = TRUE),
    mean_ndvi_xcal = mean(ndvi.xcal, na.rm = TRUE)
  )

#Export as CSV
write.csv(ndvi_summary, "ndvi_summary.csv", row.names = FALSE)

ggplot(polylsat.dt, aes(x = year)) +
  geom_smooth(aes(y = ndvi, color = "NDVI"), method = "loess") +
  geom_smooth(aes(y = ndvi.xcal, color = "NDVI (XCalibrated)"), method = "loess") +
  labs(x = "Year", y = "NDVI", title = "Smoothed NDVI Trends", color = "Legend") +
  theme_minimal()

ggplot(ndvi_summary, aes(x = year)) +
  geom_line(aes(y = mean_ndvi, color = "NDVI"), size = 1) + 
  geom_line(aes(y = mean_ndvi_xcal, color = "NDVI (XCalibrated)"), size = 1) +
  labs(
    title = "Mean NDVI vs. Calibrated NDVI Over Time",
    x = "Year",
    y = "Mean NDVI",
    color = "Legend"
  ) +
  theme_minimal()

#################################################################################
################################ Masked AOI #####################################
#################################################################################
# Load the raster for the AOI

BorealDomainAfile <- "C:/Users/camr1n24/Documents/LandsatTS/BorealForestDomain_EasternCanadianShieldTaigaA.tif"
BorealDomainBfile <- "C:/Users/camr1n24/Documents/LandsatTS/BorealForestDomain_EasternCanadianShieldTaigaB.tif"
BorealDomainA <- rast(BorealDomainAfile)
BorealDomainB <- rast(BorealDomainBfile)

BorealDomain <- merge(BorealDomainA, BorealDomainB)
#make zero values NA
BorealDomain[BorealDomain == 0] <- NA
plot(BorealDomain)

#Pixel Count
# Total number of pixels (including NA)
total_pixels <- ncell(BorealDomain)
print(total_pixels)
non_na_pixels <- global(!is.na(BorealDomain), "sum", na.rm = TRUE) #Terra
non_na_pixels <- cellStats(!is.na(r), sum) #Raster
print(non_na_pixels)

aoi.r <- raster(BorealDomain)
plot(aoi.r)

n.samples <- 10000
pts <- sampleRandom(aoi.r, size = n.samples*5, xy = T, sp = T, na.rm=T)
pts <- pts[sample(1:nrow(pts), n.samples, replace = F),]

pts$sample_id <- paste0('pixel_', 1:nrow(pts)) #Column has to be named "sample_id"
plot(pts, add=T)

#NDVI_Timeseries_2yrs <- raster("C:/Users/camr1n24/Documents/LAI_Timeseries/Composites/Two-year_1984-2022/NDVI_1984-2022_2yrs.tif")
agg.r <- aggregate(BorealDomain, fact = 10000, fun = "max", na.rm = TRUE)
aoi.300km.r <- raster(agg.r) #::aggregate(aoi.r, fact = 10000, fun = 'max')
vals <- values(aoi.300km.r)
vals <- vals[is.na(vals)==F]
aoi.300km.r[aoi.300km.r==1] <- 1:length(vals)
plot(aoi.300km.r)
plot(agg.r)

pts$cluster <- paste0('cluster_', raster::extract(aoi.300km.r, pts))
length(unique(pts$cluster))

pts.sf <- st_as_sf(pts)
plot(pts.sf, add = T, col = 'red', pch = 16, cex = 0.3) #pch for solid circle, cex for size

# EXTRACT LANDSAT DATA USING GEE ====================================================================================

lsat_export_ts(pts.sf, start_doy = 152, end_doy = 243, start_date = "1984-01-01", end_date = "2024-12-31", buffer_dist = 30,
               chunks_from = 'cluster', drive_export_dir = "lsatTS_boreal", file_prefix = "lsatTS_export")

## DID NOT DELETED RECENT FIRES (2023) so exclude samples of the years were disturbances were not masked

exported_clusters_path <- "D:/GEE/SamplesBufferFullMask/cluster_samples"

# Vector of folder paths
exported_clusters <- list(exported_clusters_path)

# Use lapply to get CSV files from each folder
temp_files <- lapply(exported_clusters, list.files, pattern = "\\.csv$", full.names = TRUE)

lsat.dt <- do.call("rbind", lapply(temp_files[[1]], fread))
saveRDS(lsat.dt, file = file.path("D:/GEE/SamplesBufferFullMask", "lsatdt.rds"))
total_pixels <- unique(lsat.dt$sample_id) # Check the years present in the data

n_total <- nrow(lsat.dt) # Check number of rows in the data table
n_subsets <- 3 #number of subsets
subset_indices <- split(1:n_total, ceiling(seq_along(1:n_total) / ceiling(n_total / n_subsets)))

output_dir <- "D:/GEE/SamplesBufferFullMask/formatted_data"

for (i in seq_len(n_subsets)) {
  message("Processing subset ", i, " of ", n_subsets)
  
  rows <- subset_indices[[i]]
  dtsubset <- lsat.dt[rows, ]  # This subset will never overlap or miss rows
  
  # Format and save the subset
  formatted_subset <- lsat_format_data(dtsubset)
  saveRDS(formatted_subset, file = file.path(output_dir, paste0("formatted_subset_", i, ".rds")))
  
  # Clean up memory
  rm(dtsubset, formatted_subset); gc()
  
  # Remove R's temporary files to free disk space
  temp_files <- list.files(tempdir(), full.names = TRUE, recursive = TRUE)
  unlink(list.files(tempdir(), full.names = TRUE, recursive = TRUE), recursive = TRUE, force = TRUE)
}

formatted_files <- list.files(output_dir, pattern = "\\.rds$", full.names = TRUE)

formatted_data <- do.call(rbind, lapply(formatted_files, readRDS))
saveRDS(formatted_data, file = file.path(output_dir, "formattedlsatdt.rds"))

cleanlsat.dt <- lsat_clean_data(formatted_data, geom.max = 15, cloud.max = 80, sza.max = 60, filter.cfmask.snow = T, filter.cfmask.water = T, filter.jrc.water = T)
saveRDS(cleanlsat.dt, file = file.path(output_dir, "cleanformattedlsatdt.rds"))
bufcleanlsat.dt <- readRDS(file.path("D:/GEE/SamplesBufferedPartiallyMasked", "cleanformattedlsatdt.rds"))

### Loop to avoid crash

n_subsets <- 3
clean_list <- vector("list", n_subsets)

for (i in seq_len(n_subsets)) {
  
  # Read original formatted subset
  lsat_dt <- readRDS(file.path(output_dir, paste0("formatted_subset_", i, ".rds")))
  
  # Clean the data
  clsat_dt <- lsat_clean_data(
    lsat_dt,
    geom.max = 15,
    cloud.max = 80,
    sza.max = 60,
    filter.cfmask.snow = TRUE,
    filter.cfmask.water = TRUE,
    filter.jrc.water = TRUE
  )
  
  # Save cleaned file
  cleaned_file <- file.path(output_dir, paste0("clsatdt", i, ".rds"))
  saveRDS(clsat_dt, cleaned_file)
  
  # Read it back in (if needed)
  clean_list[[i]] <- readRDS(cleaned_file)
}

# Merge all cleaned subsets into one dataframe
cleanlsat.dt <- do.call(rbind, clean_list)
saveRDS(cleanlsat.dt, file = file.path(output_dir, "cleanformattedlsatdt.rds"))

## Because it used buffered
#neighborslat.dt <- lsat_neighborhood_mean(cleanlsat.dt) 

total_years <- unique(cleanlsat.dt$year) # Check the years present in the data
total_pixels <- unique(cleanlsat.dt$sample.id) # Check the years present in the data

# Remove non-masked years
# bufcleanlsat.dt <- bufcleanlsat.dt %>%
#   filter(year >= 1984 & year <= 2020) # Adjust the years as needed

# Get max and min values from column
max_val <- max(cleanlsat.dt$nir, na.rm = TRUE)
min_val <- min(cleanlsat.dt$nir, na.rm = TRUE)

n_samples <- 10000  
unique_ids <- unique(cleanlsat.dt$sample.id) # Get unique sample IDs

sampled_ids <- sample(unique_ids, n_samples)
sampled_data <- cleanlsat.dt[cleanlsat.dt$sample.id %in% sampled_ids, ]
sample_pixels <- unique(sampled_data$sample.id)

calibratedlsat.dt <- lsat_calibrate_poly(cleanlsat.dt, 
                                    band.or.si = 'nir', 
                                    doy.rng = 151:242, 
                                    min.obs = 5, #Minimum number of paired, seasonally-matched observations from Landsat 7 and Landsat 5/8 required to include a sampling sample.
                                    frac.train = 0.75,# Fraction of the data used for training the model. The rest is used for validation.
                                    train.with.highlat.data = F,
                                    trim = T,
                                    overwrite.col = F, 
                                    write.output = F)
