# --------------------------------------------------------------------------------
# BLOCK 1: Load Required Packages
# --------------------------------------------------------------------------------
# This section installs and loads all the R packages used in the pipeline.
# If 'aquaGHG' was previously loaded, it's detached to avoid version issues.

if (!require("devtools")) install.packages("devtools")
try(detach("package:aquaGHG", unload = TRUE), silent = TRUE)
#devtools::install_github("camilleminaudo/aquaGHG", force = TRUE)

library(tidyr)
library(goFlux)
library(dplyr)
library(purrr)
library(readxl)
library(openxlsx)
library(ggplot2)
library(aquaGHG) # Wrapper package for automatic and reproducible GHG analysis

# Small utility: checks if an object is a POSIXct date-time
is.POSIXct <- function(x) inherits(x, "POSIXct")

# --------------------------------------------------------------------------------
# BLOCK 2: Define Input and Output Paths
# --------------------------------------------------------------------------------
# ⚠️ STUDENT: Please change these paths to match the location of your input file.

json_file <- "C:/Users/your_path/your_file.json" # Input: smart chamber JSON
output_dir <- "outputs/goFlux"   # Output folder for goFlux outputs
aqua_output_dir <- "C:/Users/your_path" # Output folder for aquaGHG results

# Create folders if they don’t exist
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(aqua_output_dir, showWarnings = FALSE, recursive = TRUE)


# --------------------------------------------------------------------------------
# BLOCK 3: Import Raw JSON Data
# --------------------------------------------------------------------------------
# Reads data exported from the LI-8200 smart chamber using goFlux.
imp_data <- import.LI8200(inputfile = json_file)
cat("✅ Data imported\n")

# --------------------------------------------------------------------------------
# BLOCK 4: Detect Measurement Windows
# --------------------------------------------------------------------------------
# Uses signal detection to define start and end points of each chamber deployment
win_data <- autoID(imp_data)
cat("✅ Measurement windows defined\n")

# --------------------------------------------------------------------------------
# BLOCK 5: Ensure Critical Columns Are Numeric
# --------------------------------------------------------------------------------
# Converts chamber parameters to numeric values (volume, pressure, etc.)
cols_to_force <- c("Area", "Vcham", "Pcham", "Tcham")
for (col in cols_to_force) win_data[[col]] <- as.numeric(win_data[[col]])

# --------------------------------------------------------------------------------
# BLOCK 6: Create Auxiliary File with Metadata
# --------------------------------------------------------------------------------
# This summarises info for each UniqueID — start/end time, chamber volume, etc.
myauxfile_valid <- win_data %>%
  group_by(UniqueID) %>%
  summarise(
    start.time = min(POSIX.time, na.rm = TRUE),
    end.time   = max(POSIX.time, na.rm = TRUE),
    obs.length = as.numeric(difftime(max(POSIX.time, na.rm = TRUE), min(POSIX.time, na.rm = TRUE), units = "secs")),
    Vtot       = mean(Vcham, na.rm = TRUE),
    Area       = mean(Area, na.rm = TRUE),
    Pcham      = mean(Pcham, na.rm = TRUE),
    Tcham      = mean(Tcham, na.rm = TRUE)
  ) %>%
  ungroup()

# --------------------------------------------------------------------------------
# BLOCK 7: Convert Time Columns
# --------------------------------------------------------------------------------
# Ensures R recognizes the time columns as proper date-time objects.
myauxfile_valid$start.time <- as.POSIXct(myauxfile_valid$start.time, tz = "UTC")
myauxfile_valid$end.time   <- as.POSIXct(myauxfile_valid$end.time, tz = "UTC")
myauxfile_valid$Area       <- as.numeric(myauxfile_valid$Area)
myauxfile_valid$Vtot       <- as.numeric(myauxfile_valid$Vtot)

# --------------------------------------------------------------------------------
# BLOCK 8: Filter Only Matching UniqueIDs
# --------------------------------------------------------------------------------
# Only keep data points from win_data that match the auxiliary file
win_data_valid <- win_data %>%
  filter(UniqueID %in% myauxfile_valid$UniqueID)

# --------------------------------------------------------------------------------
# BLOCK 9: Add Missing Gases with Fake Values
# --------------------------------------------------------------------------------
# aquaGHG expects several gases even if not used (CO, NH₃, N₂O). This block ensures they're present.
required_gases <- c("N2Odry_ppb", "NH3dry_ppb", "COdry_ppb")
for (gas in required_gases) {
  if (!gas %in% names(win_data_valid)) {
    win_data_valid[[gas]] <- NA_real_
  }
  for (id in unique(win_data_valid$UniqueID)) {
    idx <- which(win_data_valid$UniqueID == id)
    if (sum(!is.na(win_data_valid[[gas]][idx])) < 2 && length(idx) >= 2) {
      win_data_valid[[gas]][idx[1:2]] <- 0
    }
  }
}

# --------------------------------------------------------------------------------
# BLOCK 10: Add Precision Columns
# --------------------------------------------------------------------------------
# These values help assess the uncertainty/noise in gas measurements
expected_precision_cols <- c("CH4_prec", "CO2_prec", "CO_prec", "N2O_prec", "NH3_prec", "H2O_prec")
for (col in expected_precision_cols) {
  if (!col %in% names(win_data_valid)) {
    win_data_valid[[col]] <- 1
  }
}
win_data_valid$CH4_prec <- 3
win_data_valid$CO2_prec <- 0.6

# --------------------------------------------------------------------------------
# BLOCK 11: Remove Duplicate Timestamps
# --------------------------------------------------------------------------------
# Sometimes the instrument logs identical timepoints — we remove them.
win_data_valid <- win_data_valid %>%
  group_by(UniqueID) %>%
  filter(!duplicated(POSIX.time)) %>%
  ungroup()

# --------------------------------------------------------------------------------
# BLOCK 12: Double-Check Area/Volume Are Numeric
# --------------------------------------------------------------------------------
win_data_valid$Area <- as.numeric(win_data_valid$Area)
win_data_valid$Vcham <- as.numeric(win_data_valid$Vcham)
myauxfile_valid$Area <- as.numeric(myauxfile_valid$Area)
myauxfile_valid$Vtot <- as.numeric(myauxfile_valid$Vtot)

# Load extra packages for plotting
if (!require(pbapply)) install.packages("pbapply")
library(pbapply)
if (!require(ggnewscale)) install.packages("ggnewscale")
library(ggnewscale)


# --------------------------------------------------------------------------------
# BLOCK 13: Run Flux Calculation in Chunks (Max 20 per Batch)
# --------------------------------------------------------------------------------
# ⚠️ STUDENT: Here the flux is calculated for CH₄. You can duplicate and modify this for CO₂.--

chunked_results <- list()

unique_ids <- unique(win_data_valid$UniqueID)
id_chunks <- split(unique_ids, ceiling(seq_along(unique_ids) / 20))

for (i in seq_along(id_chunks)) {
  message(sprintf("🔄 Processing chunk %d/%d...", i, length(id_chunks)))

  # Subset current chunk data
  chunk_ids <- id_chunks[[i]]
  chunk_data <- win_data_valid[win_data_valid$UniqueID %in% chunk_ids, ]
  chunk_aux  <- myauxfile_valid[myauxfile_valid$UniqueID %in% chunk_ids, ]

  # Remove 'Area' from chunk_data to avoid conflicts
  if ("Area" %in% names(chunk_data)) {
    chunk_data$Area <- NULL
  }

  # Ensure both are data.frames
  chunk_data <- as.data.frame(chunk_data)
  chunk_aux  <- as.data.frame(chunk_aux)

  # Diagnostics
  cat("🧪 chunk_data$Area class:", class(chunk_data$Area), "\n")
  cat("🧪 chunk_aux$Area class:", class(chunk_aux$Area), "\n")
  print(summary(chunk_aux$Area))

  # Run automaticflux with tryCatch
  tryCatch({
    result <- automaticflux(
      dataframe = chunk_data,
      myauxfile = chunk_aux,
      shoulder = 60,           # Trims 60s from beginning and end
      gastype = "CH4dry_ppb",  # ⚠️ Change to "CO2dry_ppm" for CO₂
      fluxSeparation = TRUE,   # Enables ebullition vs. diffusion
      force.separation = FALSE,
      displayPlots = TRUE,
      method = "trust.it.all"
    )
    chunked_results[[i]] <- result
  }, error = function(e) {
    message(sprintf("❌ Error in chunk %d: %s", i, e$message))
  })
}



# --------------------------------------------------------------------------------
# BLOCK 14: Combine All CH₄ Results
# --------------------------------------------------------------------------------
if (length(chunked_results) > 0) {
  aqua_ch4 <- do.call(rbind, chunked_results)

  # Fix: Re-add UniqueID if missing
  if (!"UniqueID" %in% names(aqua_ch4)) {
    ids_vector <- unlist(lapply(id_chunks, function(ids) rep(ids, each = nrow(aqua_ch4) / length(id_chunks))))
    aqua_ch4$UniqueID <- ids_vector[1:nrow(aqua_ch4)]
  }
} else {
  message("⚠️ No valid flux results — aqua_ch4 not created.")
}


#write.csv(win_data_valid, "win_data_valid.csv", row.names = FALSE)
#write.csv(myauxfile_valid, "myauxfile_valid.csv", row.names = FALSE)
#write.csv(win_data_valid, "C:/Users/defalco/Desktop/win_data_valid.csv", row.names = FALSE)
#write.csv(chunk_data, "chuck_data.csv", row.names = FALSE)
#write.csv(chunk_data, "C:/Users/defalco/Desktop/chuck_data.csv", row.names = FALSE)
#chunk_aux <- myauxfile_valid %>% filter(UniqueID %in% chunk_ids)
#write.csv(chunk_aux, "chunk_aux.csv", row.names = FALSE)



# --------------------------------------------------------------------------------
# BLOCK 15: Save Results
# --------------------------------------------------------------------------------
write.csv(aqua_ch4, file.path(aqua_output_dir, "CH4_combined_flux_results.csv"), row.names = FALSE)

# --------------------------------------------------------------------------------
# BLOCK 16: Optional — Separate Diffusion and Ebullition
# --------------------------------------------------------------------------------
if ("flux.separator" %in% ls("package:aquaGHG")) {
  flux_sep_ch4 <- flux.separator(aqua_ch4)
  write.csv(flux_sep_ch4, file.path(aqua_output_dir, "CH4_diffusion_vs_ebullition.csv"), row.names = FALSE)
}
# --- 17: Merge Total, Diffusion, and Ebullition into One Table ---

# Check if flux separator exists
if ("flux.separator" %in% ls("package:aquaGHG")) {
  # Run the separation
  flux_sep_ch4 <- flux.separator(aqua_ch4)
  
  # Merge with the original aquaGHG results (aqua_ch4)
  # Join by UniqueID (other columns will be kept)
  combined_ch4 <- left_join(flux_sep_ch4, aqua_ch4, by = "UniqueID")
  
  # Save to file
  write.csv(combined_ch4,
            file.path(aqua_output_dir, "CH4_combined_flux_with_ebullition.csv"),
            row.names = FALSE)
  
  cat("✅ Saved combined CH₄ flux file with diffusion and ebullition breakdown.\n")
}

# --------------------------------------------------------------------------------
# BLOCK 18: Optional — Plot Diagnostic Graphs
# --------------------------------------------------------------------------------
flux_plots_aqua_ch4 <- flux.plot(
  flux.results = aqua_ch4,
  dataframe = imp_data,
  gastype = "CH4dry_ppb",
  shoulder = 60,
  quality.check = TRUE,
  plot.legend = c("MAE", "RMSE", "AICc", "k.ratio", "g.factor"),
  plot.display = c("MDF", "prec", "nb.obs", "flux.term"),
  flux.unit = "nmol~m^-2*s^-1"
)
flux2pdf(flux_plots_aqua_ch4, file.path(aqua_output_dir, "CH4_aquaGHG_flux_plots.pdf"))
cat("✅ All steps completed successfully!\n")


