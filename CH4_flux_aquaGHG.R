# --- 1. Load Required Packages ---
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
library(aquaGHG)
is.POSIXct <- function(x) inherits(x, "POSIXct")

# --- 2. Define Paths ---
json_file <- "C:/Users/defalco/OneDrive - BGU/python_projects/floating_chamber_ch4/smart_chamber/data/Yeruham_Smart_Chamber_12032025.json"
output_dir <- "outputs/goFlux"
aqua_output_dir <- "C:/Users/defalco/OneDrive - BGU/python_projects/floating_chamber_ch4/smart_chamber/GoFlux/aquafiles"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(aqua_output_dir, showWarnings = FALSE, recursive = TRUE)

# --- 3. Import JSON Data ---
imp_data <- import.LI8200(inputfile = json_file)
cat("✅ Data imported\n")

# --- 4. Detect Measurement Windows ---
win_data <- autoID(imp_data)
cat("✅ Measurement windows defined\n")

# --- 5. Ensure critical columns are numeric ---
cols_to_force <- c("Area", "Vcham", "Pcham", "Tcham")
for (col in cols_to_force) win_data[[col]] <- as.numeric(win_data[[col]])

# --- 6. Create Auxiliary File Automatically ---
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

# --- 7. Clean types ---
myauxfile_valid$start.time <- as.POSIXct(myauxfile_valid$start.time, tz = "UTC")
myauxfile_valid$end.time   <- as.POSIXct(myauxfile_valid$end.time, tz = "UTC")
myauxfile_valid$Area       <- as.numeric(myauxfile_valid$Area)
myauxfile_valid$Vtot       <- as.numeric(myauxfile_valid$Vtot)

# --- 8. Clean and filter win_data_valid based on myauxfile_valid UniqueIDs ---
win_data_valid <- win_data %>%
  filter(UniqueID %in% myauxfile_valid$UniqueID)

# --- 9. Patch required gases ---
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

# --- 10. Add precision columns ---
expected_precision_cols <- c("CH4_prec", "CO2_prec", "CO_prec", "N2O_prec", "NH3_prec", "H2O_prec")
for (col in expected_precision_cols) {
  if (!col %in% names(win_data_valid)) {
    win_data_valid[[col]] <- 1
  }
}
win_data_valid$CH4_prec <- 3
win_data_valid$CO2_prec <- 0.6

# --- 11. Remove duplicate timestamps ---
win_data_valid <- win_data_valid %>%
  group_by(UniqueID) %>%
  filter(!duplicated(POSIX.time)) %>%
  ungroup()

# --- 12. Ensure numeric consistency again ---
win_data_valid$Area <- as.numeric(win_data_valid$Area)
win_data_valid$Vcham <- as.numeric(win_data_valid$Vcham)
myauxfile_valid$Area <- as.numeric(myauxfile_valid$Area)
myauxfile_valid$Vtot <- as.numeric(myauxfile_valid$Vtot)

if (!require(pbapply)) install.packages("pbapply")
library(pbapply)
if (!require(ggnewscale)) install.packages("ggnewscale")
library(ggnewscale)


# --- 13. Run CH₄ Flux in chunks of ≤20 UniqueIDs ---

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
      shoulder = 60,
      gastype = "CH4dry_ppb",
      fluxSeparation = TRUE,
      force.separation = FALSE,
      displayPlots = TRUE,
      method = "trust.it.all"
    )
    chunked_results[[i]] <- result
  }, error = function(e) {
    message(sprintf("❌ Error in chunk %d: %s", i, e$message))
  })
}



# Combine all CH₄ results only if results exist
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



# --- 14. Save CH₄ results ---
write.csv(aqua_ch4, file.path(aqua_output_dir, "CH4_combined_flux_results.csv"), row.names = FALSE)

# --- 15. (Optional) Diffusion vs. Ebullition separation ---
if ("flux.separator" %in% ls("package:aquaGHG")) {
  flux_sep_ch4 <- flux.separator(aqua_ch4)
  write.csv(flux_sep_ch4, file.path(aqua_output_dir, "CH4_diffusion_vs_ebullition.csv"), row.names = FALSE)
}

# --- 16. Plot CH₄ fluxes ---
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


