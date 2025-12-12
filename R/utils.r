# --- Helper Function: Get Internal Edges & Boundary Polygons ---
process_file_internals <- function(file_path, global_offset, temp_dir) {
        
        # 1. Load Data
        polys <- sf::st_read(file_path, quiet = TRUE)
        n <- nrow(polys)
        
        # 2. Assign Global IDs
        # We use these IDs for the final graph
        polys$global_id <- (global_offset + 1):(global_offset + n)
        
        # 3. Calculate Internal Neighbors (Queen)
        nb <- spdep::poly2nb(polys, queen = TRUE, snap = 1e-6)
        
        # 4. Extract Internal Edges (List method - RAM safe)
        edges <- lapply(1:n, function(i) {
                nbrs <- nb[[i]]
                valid <- nbrs[nbrs > i & nbrs != 0] # Keep only one direction
                if(length(valid) > 0) {
                        return(data.frame(from = polys$global_id[i], 
                                          to = polys$global_id[valid]))
                }
                return(NULL)
        })
        edge_df <- do.call(rbind, edges)
        
        # 5. Identify "Rim" Polygons (for stitching later)
        # Heuristic: Polygons that intersect the Bounding Box of the file usually 
        # contain the set of polygons that might touch a neighbor file.
        # For irregular shapes, we just keep geometries.
        
        # To save RAM, we only save the geometry and global_id of polygons 
        # that touch the BBox border.
        bbox <- sf::st_bbox(polys)
        bbox_poly <- sf::st_as_sfc(bbox)
        
        # Convert bbox boundary to lines to catch touching polygons
        bbox_boundary <- sf::st_cast(bbox_poly, "MULTILINESTRING")
        
        # Filter: Keep polygons touching the file's bounding box
        rim_indices <- sf::st_intersects(polys, bbox_boundary, sparse = FALSE)[,1]
        rim_polys <- polys[rim_indices, c("global_id")] # Keep only ID and geometry
        
        # 6. Save outputs
        # Save edges to a CSV (append mode not ideal for CSV, better to list and bind later, 
        # but for massive data, writing to disk is safer)
        return(list(
                edges = edge_df,
                rim = rim_polys,
                bbox = bbox,
                count = n,
                file = file_path
        ))
}
