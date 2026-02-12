library(devtools)
install_git("https://github.com/JonJup/pulseR.git")
library(pulseR)


library(sfarrow)
library(arrow)
library(dplyr)
library(sf)
library(mclust)

x <- st_read_parquet("../paper/data/rawData/parquet/angerman.parquet")
y <- read_parquet("../paper/data/rawData/catchments/Angerman_w_variables.parquet")


x2 <- filter(x, ID %in% y$ID)
z <- left_join(x2, y, by = "ID")
z$saturated_soil_water_content <- round(z$saturated_soil_water_content, 3)

# Are there NAs n Z?
if (any(is.na(z))){
        
        vars_with_na <- names(z)[colSums(is.na(z)) > 0]
        # Print status
        cat("Imputing variables:", paste(vars_with_na, collapse = ", "), "\n")
        # 4. Loop through each variable containing NAs
        for (var in vars_with_na) {
                
                # Identify which rows (indices) have missing values for this variable
                na_indices <- which(is.na(z[[var]]))
                
                # Iterate through each missing row
                for (i in na_indices) {
                        nb <- st_nearest_feature(x = z[i, ], z)
                        z[[var]][[i]] <- z[[var]][nb]
                        
                }
        }
        
}

## Regions 
g <- pulseR::dev_polygon_to_network(z)
g2 <- pulseR::add_edge_weight(g, id_col = "ID")
r  <- pulseR::skater_con(g2, final_regions = 3, n.rst = 100)
pulseR::plot_fuzzy_area(g$polygons,fuzzy_memberships = r, plot = "dominant")

## Rivers
ca <- pulseR::get_core_regions(x = cbind(g$polygons, r), membership_cols = paste0("X", 1:3), cutoff = 0.8)
rt <- pulseR::river_types(core_regions = ca, membership_cols = paste0("X", 1:3), non_value_cols = "ID", n_river_types = 1:20, membership_id = r, all_data = z)
pulseR::plot_fuzzy_area(g$polygons,fuzzy_memberships = rt$rivertypes, plot = "dominant")
pulseR::plot_river_centroids(models = rt$models)
## Landscapes 
