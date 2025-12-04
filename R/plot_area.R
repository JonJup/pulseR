plot_fuzzy_area <- function(polygons, fuzzy_memberships, plot = NULL){
        
        n_types = ncol(fuzzy_memberships)
        colnames(fuzzy_memberships) <- paste0("class", 1:n_types)
        if (is.null(plot)){
                stop("Please provide one of the following options to the plot argument:\nRGB\nternary1")
        }
        
        pd <- cbind(polygons, fuzzy_memberships)
        
        if(plot == "RGB"){
                if (n_types != 3) stop("RGB only works with 3 types")
                pd$color <- rgb(pd$class1, pd$class2, pd$class3)
                ggplot(pd) +
                        geom_sf(aes(fill = color), color = NA) +
                        theme_void() +
                        guides(fill = "none")
                        
        } else if(plot == "ternary1"){
                if (n_types != 3) stop("ternary only works with 3 types")
                data_df <- pd %>% st_drop_geometry()
                
                ggtern::ggtern(data_df, aes(x = class1, y = class2, z = class3)) +
                        geom_point(alpha = 0.6, color = "darkblue", size = 2) +
                        ggtern::theme_rgbw() + # A clean white theme
                        labs(
                                title = "Fuzzy Classification in Ternary Space",
                                x = "Class A",
                                y = "Class B",
                                z = "Class C"
                        )
        } else if(plot == "ternary2"){
                if (n_types != 3) stop("ternary only works with 3 types")
                tric_obj <- tricolore::Tricolore(pd, 
                                                 p1 = "class1", 
                                                 p2 = "class2", 
                                                 p3 = "class3")
                
                # 2. Add the generated colors back to your SF object
                pd$tric_color <- tric_obj$rgb
                
                # 3. Plot the map
                map_plot <- ggplot(pd) +
                        geom_sf(aes(fill = tric_color), color = NA) +
                        scale_fill_identity() +
                        theme_void() +
                        labs(title = "Map with Ternary Legend")
                
                # 4. Add the legend (tricolore generates a separate legend plot)
                # You can use cowplot or patchwork to combine them
                cowplot::plot_grid(map_plot, tric_obj$key, ncol = 2, rel_widths = c(3, 1))
        } else if (plot == "small multiples"){
                
                
                data_long <- pd %>%
                        pivot_longer(cols = starts_with("class"), 
                                     names_to = "Class", 
                                     values_to = "Probability") %>%  
                        select(Class, Probability)
                
                ggplot(data_long) +
                        geom_sf(aes(fill = Probability), color = NA) +
                        scale_fill_viridis_c(option = "magma") +
                        facet_wrap(~Class, ncol = 5) +
                        theme_void() +
                        labs(title = "Fuzzy Membership Probability by Class")
        } else if (plot == "dominant"){
                data_summary <- pd |>
                        #dplyr::select(tidyselect::starts_with("class")) |>
                        dplyr::rowwise() %>%
                        dplyr::mutate(
                                Max_Prob = max(c_across(starts_with("class"))),
                                Dominant_Class = names(.)[which.max(dplyr::c_across(tidyselect::starts_with("class")))]
                        ) %>%
                        dplyr::ungroup()
                
                ggplot(data_summary) +
                        geom_sf(aes(fill = Dominant_Class, alpha = Max_Prob)) +
                        scale_alpha_continuous(range = c(0.3, 1)) + # Ensure low probs are still slightly visible
                        theme_minimal() +
                        labs(title = "Dominant Class Modulated by Certainty")
        }
        
}
