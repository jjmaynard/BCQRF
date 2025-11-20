library(targets)
library(dplyr)
library(here)
library(yaml)

tar_source()

# Set global options if needed, e.g.,
tar_option_set(packages = c("here", "dplyr", "magrittr"))

# Load configuration parameters
config <- yaml::read_yaml(here("config.yml"))


list(
  # Target 1: loads global_cor_matrices from RDS
  tar_target(
    global_cor_matrices,
    readRDS(here::here("raw-data/global_cor_matrices.rds")),
    format = "rds"  # recommended for storing R objects
  ),

  # Target 2: loads global_cor_txt_matrices from RDS
  tar_target(
    global_cor_txt_matrices,
    readRDS(here::here("raw-data/global_cor_texture_matrices.rds")),
    format = "rds"
  ),

  # Target 3: Download SSURGO Data for soil properties within AOI
  tar_target(
    ssurgo_data,
    {
      download_ssurgo_tabular(aoi_wkt = config$aoi_wkt, properties = config$properties)
    }
  ),

  # Target 4: Simulate component composition per mukey
  tar_target(
    ssurgo_comppct_sim,
    {
      print(names(ssurgo_data))

      ssurgo_data %>%
        dplyr::group_by(mukey) %>%
        dplyr::group_modify(~ {
          # Attach the grouping key to .x so that sim_component_comp() can use it.
          .x$mukey <- .y$mukey
          out <- sim_component_comp(.x, n_simulations = 1000)
          # Remove the grouping variable from the output if it is present.
          if ("mukey" %in% colnames(out)) {
            out <- dplyr::select(out, -mukey)
          }
          out
        }) %>%
        dplyr::ungroup()
    }
  ),

  # Target 5: Join simulation results back to ssurgo_data
  tar_target(
    ssurgo_data_sim,
    {
      ssurgo_data %>%
        left_join(ssurgo_comppct_sim %>% select(cokey, sim_comppct), by = "cokey")
    }
  ),

  # Target 6: Remove organic layers if present
  tar_target(
    ssurgo_data_clean,
    {
      if (any(grepl("O", ssurgo_data_sim$hzname))) {
        remove_organic_layer(ssurgo_data_sim)
      } else {
        ssurgo_data_sim
      }
    }
  ),

  # Target 7: Create a SoilProfileCollection from ssurgo_data
  tar_target(
    ssurgo_profile_collection,
    {
      sp_df <- ssurgo_data_clean %>%
        # Convert horizon boundary columns to numeric, if they aren't already
        dplyr::mutate(
          hzdept_r = as.numeric(hzdept_r),
          hzdepb_r = as.numeric(hzdepb_r)
        )

      # Ensure the profile identifier exists. If not, create one.
      # For example, if 'id' is missing, you might create it from 'cokey':
      if (!"id" %in% colnames(sp_df)) {
        sp_df <- sp_df %>% dplyr::mutate(id = cokey)
      }

      # Set depths; this converts the data.frame into a SoilProfileCollection
      aqp::depths(sp_df) <- id ~ hzdept_r + hzdepb_r
      aqp::hzdesgnname(sp_df) <- "hzname"

      # Return the SoilProfileCollection
      sp_df
    }
  ),

  # Target 8: Simulate soil profile horizon depth variability in parallel
  tar_target(
    ssurgo_depth_simulation,
    simulate_profile_depths_by_collection_parallel(
      soil_collection = ssurgo_profile_collection,
      seed = 123,
      n_cores = 6
    )
  ),

  # Target 9: Simulate soil properties by Cokey
  tar_target(
    ssurgo_soil_sim_df,
    simulate_soil_properties(ssurgo_profile_collection, global_cor_matrices, global_cor_txt_matrices)
  ),

  # Target 10: Reshape the data into depth-wise matrices
  tar_target(
    final_processed_data,
    {
      adjust_soil_data(ssurgo_soil_sim_df, properties=config$properties)
    }
  ),

  # Target 11: Combine simulated depths w/ simulated properties
  tar_target(
    ssurgo_depth_property_sim,
    {
      # Extract simulated depth data
      ssurgo_depth_simulation_data <- ssurgo_depth_simulation@horizons

      # Extract key columns from ssurgo_data
      ssurgo_df <- ssurgo_data@horizons %>%
        dplyr::select(cokey, hzname, hzdept_r, hzdepb_r)

      # Join simulated properties with ssurgo_df based on cokey and depth columns
      final_processed_data <- final_processed_data %>%
        dplyr::left_join(
          ssurgo_df,
          by = c("cokey" = "cokey", "hzdept_r" = "hzdept_r", "hzdepb_r" = "hzdepb_r")
        )

      # Combine the simulated depth data with the processed properties
      ssurgo_depth_property_sim <- ssurgo_depth_simulation_data %>%
        dplyr::left_join(
          final_processed_data %>%
            dplyr::select(unique_id, hzname, bulk_density_third_bar, water_retention_third_bar, water_retention_15_bar, rfv, sand_total, silt_total, clay_total),
          by = c("id" = "unique_id", "hzname" = "hzname")
        )

      # Optionally view the result in an interactive session:
      # View(ssurgo_depth_property_sim)

      # Set depth information for the soil profile object
      depths(ssurgo_depth_property_sim) <- id ~ hzdept_r + hzdepb_r

      # Return the final combined data
      ssurgo_depth_property_sim
    }
  )

)

