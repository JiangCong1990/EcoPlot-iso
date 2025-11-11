### Conceptual plot-scale isotope SoilWaterBalance model (isoEcoPlot) ###
## Christian Birkel, OACG-UCR | Revised by Cong Jiang, 2024.11 ###

# Clean workspace
graphics.off()  # Clear all plots
rm(list = ls(all = TRUE))  # Remove all objects
cat("\014")  # Clear console

# Load required libraries
library(parallel)
library(data.table)
library(FME)         # Calibration and uncertainty estimation
library(tidyverse)   # Data manipulation
library(ncdf4)

# ---------------------------------------
# 1. Load Input Data
# ---------------------------------------
# Main input file
inp <- read.csv("Demnitz_inpmod_Fo.csv", na.strings = "", sep = ",")
names(inp) <- c("Date", "P", "NR", "AT", "RH", "LAI", "PET", "SW_10", "SW_30", 
                "SW_100", "TF", "Tobs", "P_D", "SW_10_D", "SW_30_D", "SW_100_D",
                "TempS_10", "TempS_30", "TempS_100")

# LAI input file
inp_lai <- read.csv("Demnitz_inpmod_Fo_LAI.csv", na.strings = "", sep = ",")
names(inp_lai) <- c("Date", "LAI_1", "LAI_2", "LAI_3")

# Parameter file
para <- read.csv("isoEcoPlot_LHres_sel_100sims.csv", na.strings = "")[1:20]
names(para) <- c("rE", "alpha", "Smax", "Ic", "ks1", "ks2", "k3", "GWmax", "Lmax", 
                 "g1", "g2", "g3", "PF_Scale", "INTp", "stoSp", "gwSp", "lowSp", 
                 "k", "x", "beta")

# ---------------------------------------
# 2. Simulation Function
# ---------------------------------------
source("out_SWBiso_Forest_Scen.R")  # Load SWMsim function

Q_sim <- function(pa, inp) {
  model_output <- SWMsim(rE = pa[1], alpha = pa[2], Smax = pa[3], Ic = pa[4],
                         ks1 = pa[5], ks2 = pa[6], ks3 = pa[7], GWmax = pa[8],
                         Lmax = pa[9], g1 = pa[10], g2 = pa[11], g3 = pa[12],
                         PF_Scale = pa[13], INTp = pa[14], stoSp = pa[15], 
                         gwSp = pa[16], lowSp = pa[17], k = pa[18], x = pa[19], 
                         beta = pa[20], inp = inp)
  return(as.matrix(model_output))
}

# >>> Add here <<<
library(compiler)
enableJIT(3)
SWMsim <- cmpfun(SWMsim)
Q_sim  <- cmpfun(Q_sim)
# <<< Done >>>

# ---------------------------------------
# 3. Parallel Simulation Function
# ---------------------------------------
## Run parallel simulations
#cl <- makeCluster(num_cores)
#clusterExport(cl, varlist = c("Q_sim", "inp", "SWMsim"), envir = environment())

sensRange_parallel <- function(func, parInput, inp, cl) {
  # Define dimensions
  time_steps <- nrow(inp) - 366
  num_sims <- nrow(parInput)
  num_vars <- 27  # Number of output variables

  
  # Initialize smaller 3D results array for the current simulation
  results_array <- array(NA, dim = c(num_sims, time_steps, num_vars))

  # Run parallel simulations
 # cl <- makeCluster(num_cores)
 # clusterExport(cl, varlist = c("Q_sim", "inp", "SWMsim"), envir = environment())

  results <- parLapply(cl, seq_len(num_sims), function(i) {
    pa <- parInput[i, ]
    return(func(pa, inp))  # Ensure function call works correctly
  })

#  stopCluster(cl)

  # Populate results_array
  for (i in seq_len(num_sims)) {
    results_array[i, , ] <- results[[i]]
  }

  return(results_array)
}

# ---------------------------------------
# 4. Helper Functions
# ---------------------------------------
# Initialize the NetCDF file for results storage

time_steps <- as.Date(inp_lai$Date[367:nrow(inp_lai)])  # Ensure Date format
reference_date <- min(time_steps)  # Reference date
time_days <- as.numeric(difftime(time_steps, reference_date, units = "days"))  # Convert to days

initialize_results_netcdf <- function(filename, num_sims, inp_lai, num_vars, forest_types, scaling_factors, forest_ages) {
  
  dim_sim <- ncdim_def("Simulations", "index", seq_len(num_sims))

  time_steps <- as.Date(inp_lai$Date[367:nrow(inp_lai)])  # Ensure Date format
  reference_date <- min(time_steps)  # Reference date
  time_days <- as.numeric(difftime(time_steps, reference_date, units = "days"))  # Convert to days
   
  dim_time <- ncdim_def("Time",  paste("days since", reference_date), time_days, longname = "Simulation Time Steps")
  dim_vars <- ncdim_def("Variables", "names", seq_len(num_vars))
  dim_forest <- ncdim_def("Forest_Types", "category", seq_along(forest_types))
  dim_scale <- ncdim_def("Scaling_Factors", "multiplier", scaling_factors)
  dim_age <- ncdim_def("Forest_Ages", "years", forest_ages)

  results_var <- ncvar_def("results_array", units = "units",
                           dim = list(dim_sim, dim_time, dim_vars, dim_forest, dim_scale, dim_age),
                           prec = "float", longname = "Simulation Results Array")

  nc <- nc_create(filename, list(results_var), force_v4 = TRUE)  # Use NetCDF4 format to allow compression
  nc_close(nc)
}

# Save a slice to NetCDF
save_results_slice <- function(filename, results_slice, forest_index, scale_index, age_index) {
  nc <- nc_open(filename, write = TRUE)
  ncvar_put(nc, "results_array", results_slice,
            start = c(1, 1, 1, forest_index, scale_index, age_index),
            count = c(dim(results_slice)[1], dim(results_slice)[2], dim(results_slice)[3], 1, 1, 1))
  nc_close(nc)
}

# Compute mean results from NetCDF slices
compute_mean_results <- function(filename, num_sims, num_rows, num_vars, forest_types, scaling_factors, forest_ages) {
  mean_results <- array(NA, dim = c(num_rows, num_vars, length(forest_types), length(scaling_factors), length(forest_ages)))
  nc <- nc_open(filename)

  for (i in seq_along(forest_types)) {
    for (j in seq_along(scaling_factors)) {
      for (k in seq_along(forest_ages)) {

        results_slice <- ncvar_get(nc, "results_array",
                                   start = c(1, 1, 1, i, j, k),
                                   count = c(num_sims, num_rows, num_vars, 1, 1, 1))

        # Compute slice mean and update mean_results

        slice_mean <- apply(results_slice, c(2, 3), mean, na.rm = TRUE)

        mean_results[, , i, j, k] <- slice_mean
      }
    }
  }
  nc_close(nc)

  return(mean_results)
}

# ---------------------------------------
# 5. Scenario Loop and Results Storage
# ---------------------------------------
forest_types <- c("LAI_1", "LAI_2", "LAI_3")
scaling_factors <- seq(0.2, 1.8, by = 0.2)
forest_ages <- c(0.0, 0.5, 1.0, 2.0)

# ---------------------------------------
# Load parameter files and compute means for LAI_1, LAI_2 and LAI_3
# ---------------------------------------

# Read and calculate mean values for LAI_1 (Agroforestry)
para_agro <- read.csv("/data/scratch/jiangcong/Ecolpt_root_agroforestry_longterm_beta/isoEcoPlot_LHres_sel_100sims.csv")
means_agro <- colMeans(para_agro[, c("rE", "alpha", "INTp")], na.rm = TRUE)
medians_agro <- apply(para_agro[, c("rE", "alpha", "INTp")], 2, median, na.rm = TRUE)

# Read and calculate mean values for LAI_2 (Broadleaf forest)
para_brofo <- read.csv("/data/scratch/jiangcong/Ecolpt_root_forest_longterm_beta/isoEcoPlot_LHres_sel_100sims.csv")
means_brofo <- colMeans(para_brofo[, c("rE", "alpha", "INTp")], na.rm = TRUE)
medians_brofo <- apply(para_brofo[, c("rE", "alpha", "INTp")], 2, median, na.rm = TRUE)

# Read and calculate mean values for LAI_3 (Confier Forest)
para_confo <- read.csv("/data/scratch/jiangcong/Ecolpt_root_forestB_longterm_beta/isoEcoPlot_LHres_sel_100sims.csv")
means_confo <- colMeans(para_confo[, c("rE", "alpha", "INTp")], na.rm = TRUE)
medians_confo <- apply(para_confo[, c("rE", "alpha", "INTp")], 2, median, na.rm = TRUE)

# Store the constants in a list for use in the loop
param_overrides <- list(
  LAI_1 = list(rE = medians_agro["rE"], alpha = medians_agro["alpha"], INTp = medians_agro["INTp"]),
  LAI_2 = list(rE = medians_brofo["rE"], alpha = medians_brofo["alpha"], INTp = medians_brofo["INTp"]),
  LAI_3 = list(rE = medians_confo["rE"], alpha = medians_confo["alpha"], INTp = medians_confo["INTp"])
)





# Function: select random rows within 5–95% range for given parameters
filter_and_sample <- function(df, vars, n = 100) {
  # Compute 5th and 95th percentiles for each variable
  q_low  <- apply(df[, vars], 2, quantile, probs = 0.05, na.rm = TRUE)
  q_high <- apply(df[, vars], 2, quantile, probs = 0.95, na.rm = TRUE)
  
  # Filter rows where all parameters fall within 5–95% range
  in_range <- apply(df[, vars], 1, function(row) all(row >= q_low & row <= q_high))
  df_filtered <- df[in_range, vars]
  
  # Randomly sample rows (with replacement if fewer than n available)
  df_sampled <- df_filtered[sample(1:nrow(df_filtered), n, replace = TRUE), ]
  return(df_sampled)
}

# Apply for each forest type
vars <- c("rE", "alpha", "INTp")

para_agro  <- read.csv("/data/scratch/jiangcong/Ecolpt_root_agroforestry_longterm_beta/isoEcoPlot_LHres_sel_100sims.csv")
para_brofo <- read.csv("/data/scratch/jiangcong/Ecolpt_root_forest_longterm_beta/isoEcoPlot_LHres_sel_100sims.csv")
para_confo <- read.csv("/data/scratch/jiangcong/Ecolpt_root_forestB_longterm_beta/isoEcoPlot_LHres_sel_100sims.csv")

# Randomly sample within 5–95% percentile range
sample_agro  <- filter_and_sample(para_agro, vars, n = 1)
sample_brofo <- filter_and_sample(para_brofo, vars, n = 1)
sample_confo <- filter_and_sample(para_confo, vars, n = 1)

# Store in list for scenario loop
param_overrides <- list(
  LAI_1 = sample_agro,
  LAI_2 = para_brofo[, c("rE", "alpha", "INTp")],  # Broadleaf: use original calibrated sets
  LAI_3 = sample_confo
)




num_rows <- nrow(inp_lai) - 366
num_sims <- nrow(para)
num_cores <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK", unset = 8))
num_vars <- 27

# Initialize results NetCDF file
results_file <- "Modelsim_results_array.nc"
initialize_results_netcdf(results_file, num_sims, inp_lai, num_vars, forest_types, scaling_factors, forest_ages)

# Run parallel simulations
cl <- makeCluster(num_cores)
clusterExport(cl, varlist = c("Q_sim", "inp", "SWMsim"), envir = environment())


#for (i in seq_along(forest_types)) {

#  current_lai <- forest_types[i]

#  for (j in seq_along(scaling_factors)) {
#    for (k in seq_along(forest_ages)) {

#      inp$LAI <- inp_lai[[forest_types[i]]] * scaling_factors[j]
#      para$beta <- forest_ages[k]

      # Conditionally override parameters if current LAI is in the override list
#      if (current_lai %in% names(param_overrides)) {
#        override_params <- param_overrides[[current_lai]]
#        para$rE <- override_params$rE
#        para$alpha <- override_params$alpha
#        para$INTp <- override_params$INTp
#      }

      # Define the log file path
#      log_file <- "Ecoplot.log"

      # Check if the log file exists and delete it
#      if (file.exists(log_file)) {
#          file.remove(log_file)  # Delete the previous log file
#      }

      # Run parallel simulations
#      Modelsim <- sensRange_parallel(Q_sim, as.matrix(para), inp, num_cores)

      # Save results slice incrementally
#      save_results_slice(results_file, Modelsim, i, j, k)
#    }
#  }
#}

      # Define the log file path
      log_file <- "Ecoplot.log"

      # Check if the log file exists and delete it
      if (file.exists(log_file)) {
          file.remove(log_file)  # Delete the previous log file
      }

# ============================================
# SLURM multi-task setup (for --ntasks=10)
# ============================================

# Detect task ID and total number of tasks
task_id <- as.numeric(Sys.getenv("SLURM_PROCID", unset = 0)) + 1
num_tasks <- as.numeric(Sys.getenv("SLURM_NTASKS", unset = 1))

message("Running SLURM task ", task_id, " of ", num_tasks)

# Define full scenario grid

scenarios <- expand.grid(
  forest_type = forest_types,
  scale = scaling_factors,
  age = forest_ages
)

scenarios$forest_type <- as.character(scenarios$forest_type)


# Divide total 108 scenarios among tasks
chunk_size <- ceiling(nrow(scenarios) / num_tasks)
start_row  <- (task_id - 1) * chunk_size + 1
end_row    <- min(task_id * chunk_size, nrow(scenarios))
my_scenarios <- scenarios[start_row:end_row, ]

message("This task handles ", nrow(my_scenarios),
        " scenarios (rows ", start_row, "–", end_row, ")")

# Each task writes its own output file
results_file <- sprintf("Modelsim_results_task%02d.nc", task_id)
initialize_results_netcdf(results_file, num_sims, inp_lai, num_vars,
                          forest_types, scaling_factors, forest_ages)



for (s in seq_len(nrow(my_scenarios))) {

  current_lai <- my_scenarios$forest_type[s]
  j_factor    <- my_scenarios$scale[s]
  k_age       <- my_scenarios$age[s]

  message("Task ", task_id, " | ", current_lai,
          " | LAI scale=", j_factor, " | age=", k_age)

  inp$LAI <- inp_lai[[current_lai]] * j_factor
  para$beta <- k_age

  if (current_lai %in% names(param_overrides)) {
    veg_param_set <- param_overrides[[current_lai]]

    ensemble_sum <- 0
    count <- 0
    for (m in 1:nrow(veg_param_set)) {
      para_mod <- para
      para_mod$rE    <- veg_param_set$rE[m]
      para_mod$alpha <- veg_param_set$alpha[m]
      para_mod$INTp  <- veg_param_set$INTp[m]

      Modelsim <- sensRange_parallel(Q_sim, as.matrix(para_mod), inp, cl)
      ensemble_sum <- ensemble_sum + Modelsim
      count <- count + 1
    }

    ensemble_mean <- ensemble_sum / count
    save_results_slice(results_file, ensemble_mean,
                       which(forest_types == current_lai),
                       which(scaling_factors == j_factor),
                       which(forest_ages == k_age))
  }
}

'''

for (i in seq_along(forest_types)) {

  current_lai <- forest_types[i]

  message("Processing forest type: ", current_lai)

  # Temporary storage for all scale–age combinations for this forest
  forest_storage <- list()


  for (j in seq_along(scaling_factors)) {
    for (k in seq_along(forest_ages)) {

      # Prepare LAI and root distribution parameter
      inp$LAI <- inp_lai[[current_lai]] * scaling_factors[j]
      para$beta <- forest_ages[k]

      # Initialize a container for ensemble results (100 vegetation parameter samples)
#      ensemble_results <- array(NA, dim = c(100,  num_sims, nrow(inp) - 366, num_vars))

      # Conditionally override vegetation-type parameters
      if (current_lai %in% names(param_overrides)) {

        veg_param_set <- param_overrides[[current_lai]]  # 100×3 matrix

        ensemble_sum <- 0
        count <- 0
        # Loop over 100 vegetation parameter samples
        for (m in 1:nrow(veg_param_set)) {

          # Copy parameter table and override vegetation-specific parameters
          para_mod <- para
          para_mod$rE    <- veg_param_set$rE[m]
          para_mod$alpha <- veg_param_set$alpha[m]
          para_mod$INTp  <- veg_param_set$INTp[m]

          # Run parallel simulations
          Modelsim <- sensRange_parallel(Q_sim, as.matrix(para_mod), inp, cl)
           ensemble_sum <- ensemble_sum + Modelsim
           count <- count + 1
           
          # Compute mean across parameter ensemble dimension (per simulation)
          #ensemble_results[m, , ,] <- Modelsim
        }
        ensemble_mean <- ensemble_sum / count

        forest_storage[[paste0("s", j, "_a", k)]] <- ensemble_mean

        # Compute ensemble mean across 100 vegetation parameter samples
        #ensemble_mean <- apply(ensemble_results, c(2, 3, 4), mean, na.rm = TRUE)

        # Save ensemble mean slice
      #  save_results_slice(results_file, ensemble_mean, i, j, k)
      }
    }
  }

  # Write all results for this forest type at once
  message("Writing results for forest type: ", current_lai)
  for (j in seq_along(scaling_factors)) {
    for (k in seq_along(forest_ages)) {
      key <- paste0("s", j, "_a", k)
      if (!is.null(forest_storage[[key]])) {
        save_results_slice(results_file, forest_storage[[key]], i, j, k)
      }
    }
  }

  rm(forest_storage)
  gc()


}

'''

stopCluster(cl)


# Compute mean results from saved slices
mean_results <- compute_mean_results(results_file, num_sims, num_rows, num_vars, forest_types, scaling_factors, forest_ages)

# Save mean results to NetCDF
mean_file <- "Modelsim_results_mean.nc"
dim_time <- ncdim_def("Time",  paste("days since", reference_date), time_days, longname = "Simulation Time Steps")
dim_vars <- ncdim_def("Variables", "names", seq_len(num_vars))
dim_forest <- ncdim_def("Forest_Types", "category", seq_along(forest_types))
dim_scale <- ncdim_def("Scaling_Factors", "multiplier", scaling_factors)
dim_age <- ncdim_def("Forest_Ages", "years", forest_ages)

mean_var <- ncvar_def("mean_results", units = "units",
                      dim = list(dim_time, dim_vars, dim_forest, dim_scale, dim_age),
                      prec = "float", longname = "Mean Results over Simulations")

ncfile_mean <- nc_create(mean_file, list(mean_var), force_v4 = TRUE)  # Use NetCDF4 format to allow compression
ncvar_put(ncfile_mean, mean_var, mean_results)
nc_close(ncfile_mean)
cat("Results saved to", results_file, "and", mean_file, "\n")
