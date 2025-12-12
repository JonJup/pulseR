#' Visualize Fuzzy Classification of Geographic Areas
#'
#' This function visualizes fuzzy classification memberships over a set of geographic polygons.
#' It supports multiple visualization methods, including RGB mapping, ternary plots,
#' small multiples, and dominant class maps modulated by uncertainty.
#'
#' @param polygons An \code{sf} object representing the geographic areas (e.g., river segments, catchments).
#'   Must contain a geometry column.
#' @param fuzzy_memberships A matrix or data frame containing the membership probabilities for each class.
#'   Rows must correspond to the rows in \code{polygons}.
#' @param plot Character string. Specifies the type of plot to generate. Options are:
#'   \describe{
#'     \item{\code{"RGB"}}{Maps the first three classes to Red, Green, and Blue channels. Requires exactly 3 classes.}
#'     \item{\code{"ternary1"}}{Produces a standard ternary plot of the data points (non-spatial).}
#'     \item{\code{"ternary2"}}{Produces a choropleth map with a ternary color scale legend (requires \code{tricolore}).}
#'     \item{\code{"small multiples"}}{Facets the map by class, showing probability intensity for each class separately.}
#'     \item{\code{"dominant"}}{Maps the class with the highest probability, with transparency scaled by the certainty (max probability).}
#'   }
#' @param auto_rename Logical. If \code{TRUE} (default), renames the columns of \code{fuzzy_memberships} to
#'   "class1", "class2", ..., "classN". If \code{FALSE}, existing column names are preserved.
#'
#' @return A \code{ggplot2} object (or a cowplot object for "ternary2").
#'
#' @details
#' \strong{Plot Types:}
#' \itemize{
#'   \item \strong{RGB:} Useful for visualizing the transition zones between exactly three distinct types.
#'   \item \strong{Ternary:} "ternary1" plots the data distribution in barycentric coordinates. "ternary2" maps these coordinates spatially.
#'   \item \strong{Small Multiples:} Best for inspecting the spatial distribution of every single class individually.
#'   \item \strong{Dominant:} Provides a "hard" classification map but retains information about uncertainty via transparency (alpha).
#' }
#'
#' @importFrom sf st_drop_geometry
#' @importFrom ggplot2 ggplot aes geom_sf theme_void guides theme_minimal labs scale_fill_identity scale_fill_viridis_c scale_alpha_continuous facet_wrap
#' @importFrom dplyr mutate rowwise c_across select ungroup bind_cols all_of
#' @importFrom tidyr pivot_longer
#' @importFrom ggtern ggtern theme_rgbw
#' @importFrom tricolore Tricolore
#' @importFrom cowplot plot_grid
#' @importFrom grDevices rgb
#'
#' @export
plot_fuzzy_area <- function(polygons, fuzzy_memberships, plot = NULL, auto_rename = TRUE) {
        
        if (is.null(plot)) {
                stop("Please provide one of the following options to the plot argument:\n 'RGB', 'ternary1', 'ternary2', 'small multiples', 'dominant'")
        }
        
        # Ensure inputs are compatible
        if (nrow(polygons) != nrow(fuzzy_memberships)) {
                stop("The number of rows in 'polygons' and 'fuzzy_memberships' must match.")
        }
        
        n_types <- ncol(fuzzy_memberships)
        
        # Rename columns if requested
        if (auto_rename) {
                colnames(fuzzy_memberships) <- paste0("class", 1:n_types)
        }
        
        # Bind geometry and data
        # We use bind_cols carefully to preserve sf class
        pd <- dplyr::bind_cols(polygons, as.data.frame(fuzzy_memberships))
        
        # --- Plotting Logic ---
        
        if (plot == "RGB") {
                if (n_types != 3) stop("RGB plot only works with exactly 3 types.")
                
                # Calculate RGB colors based on membership
                pd$color <- rgb(pd[[colnames(fuzzy_memberships)[1]]], 
                                pd[[colnames(fuzzy_memberships)[2]]], 
                                pd[[colnames(fuzzy_memberships)[3]]])
                
                p <- ggplot2::ggplot(pd) +
                        ggplot2::geom_sf(ggplot2::aes(fill = color), color = NA) +
                        ggplot2::scale_fill_identity() +
                        ggplot2::theme_void() +
                        ggplot2::guides(fill = "none")
                
                return(p)
                
        } else if (plot == "ternary1") {
                if (n_types != 3) stop("Ternary plots only work with exactly 3 types.")
                
                data_df <- sf::st_drop_geometry(pd)
                class_names <- colnames(fuzzy_memberships)
                
                p <- ggtern::ggtern(data_df, ggplot2::aes(x = .data[[class_names[1]]], 
                                                          y = .data[[class_names[2]]], 
                                                          z = .data[[class_names[3]]])) +
                        ggplot2::geom_point(alpha = 0.6, color = "darkblue", size = 2) +
                        ggtern::theme_rgbw() + 
                        ggplot2::labs(
                                title = "Fuzzy Classification in Ternary Space",
                                x = class_names[1],
                                y = class_names[2],
                                z = class_names[3]
                        )
                
                return(p)
                
        } else if (plot == "ternary2") {
                if (n_types != 3) stop("Ternary plots only work with exactly 3 types.")
                if (!requireNamespace("tricolore", quietly = TRUE)) stop("Package 'tricolore' is required for this plot.")
                if (!requireNamespace("cowplot", quietly = TRUE)) stop("Package 'cowplot' is required for this plot.")
                
                # Generate ternary colors
                data_df <- sf::st_drop_geometry(pd)
                class_names <- colnames(fuzzy_memberships)
                
                tric_obj <- tricolore::Tricolore(data_df, 
                                                 p1 = class_names[1], 
                                                 p2 = class_names[2], 
                                                 p3 = class_names[3])
                
                pd$tric_color <- tric_obj$rgb
                
                map_plot <- ggplot2::ggplot(pd) +
                        ggplot2::geom_sf(ggplot2::aes(fill = tric_color), color = NA) +
                        ggplot2::scale_fill_identity() +
                        ggplot2::theme_void() +
                        ggplot2::labs(title = "Map with Ternary Legend")
                
                # Combine map and legend
                p <- cowplot::plot_grid(map_plot, tric_obj$key, ncol = 2, rel_widths = c(3, 1))
                
                return(p)
                
        } else if (plot == "small multiples") {
                
                # Convert to long format for faceting
                class_cols <- colnames(fuzzy_memberships)
                
                data_long <- pd |>
                        tidyr::pivot_longer(cols = dplyr::all_of(class_cols), 
                                            names_to = "Class", 
                                            values_to = "Probability")
                
                p <- ggplot2::ggplot(data_long) +
                        ggplot2::geom_sf(ggplot2::aes(fill = Probability), color = NA) +
                        ggplot2::scale_fill_viridis_c(option = "magma") +
                        ggplot2::facet_wrap(~Class, ncol = 5) +
                        ggplot2::theme_void() +
                        ggplot2::labs(title = "Fuzzy Membership Probability by Class")
                
                return(p)
                
        } else if (plot == "dominant") {
                
                class_cols <- colnames(fuzzy_memberships)
                
                # Calculate dominant class and max probability
                # We apply rowwise logic only to the class columns
                data_summary <- pd |>
                        dplyr::rowwise() |>
                        dplyr::mutate(
                                Max_Prob = max(dplyr::c_across(dplyr::all_of(class_cols))),
                                Dominant_Class = class_cols[which.max(dplyr::c_across(dplyr::all_of(class_cols)))]
                        ) |>
                        dplyr::ungroup()
                
                p <- ggplot2::ggplot(data_summary) +
                        ggplot2::geom_sf(ggplot2::aes(fill = Dominant_Class, alpha = Max_Prob), color = NA) +
                        ggplot2::scale_alpha_continuous(range = c(0.3, 1), name = "Certainty") +
                        ggplot2::theme_minimal() +
                        ggplot2::labs(title = "Dominant Class Modulated by Certainty")
                
                return(p)
                
        } else {
                stop("Invalid plot type selected.")
        }
}