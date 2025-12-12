#' @importFrom utils globalVariables
utils::globalVariables(c(
        # plot_fuzzy_area variables
        "color", 
        "tric_color", 
        "Probability", 
        "Dominant_Class", 
        "Max_Prob", 
        ".data",
        
        # plot_river_centroids variables
        "PC1", 
        "PC2", 
        "Region", 
        "TypeID",
        
        # polygon_to_network variables
        "from", 
        "connected", 
        "to"
))