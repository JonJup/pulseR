library(devtools)
install_git("https://github.com/JonJup/pulseR.git")
library(pulseR)


library(sfarrow)
library(arrow)
library(dplyr)
library(sf)


x <- st_read_parquet("../paper/data/rawData/parquet/Skjern.parquet")
y <- read_parquet("../paper/data/rawData/catchments/Skjern_w_variables.parquet")


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


g2 <- dev_polygon_to_network(z)
g2 <- add_edge_weight(g2, id_col = "ID")
r  <- skater_con(g2, final_regions = 3, n.rst = 100)
