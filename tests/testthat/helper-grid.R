# A small square-grid sf landscape with two environmental gradients and an
# integer ID column. Used by tests that need the polygon -> graph pipeline.
make_grid_sf <- function(nx = 8L, ny = 8L, seed = 1L) {
        .local_seed(seed)
        bb <- sf::st_bbox(c(xmin = 0, ymin = 0, xmax = nx, ymax = ny),
                      crs = sf::st_crs("+proj=laea +lat_0=52 +lon_0=10 +ellps=GRS80 +units=m +no_defs"))
        g  <- sf::st_make_grid(sf::st_as_sfc(bb), n = c(nx, ny))
        xy <- sf::st_coordinates(sf::st_centroid(g))
        sf::st_sf(ID   = seq_len(nx * ny),
                  temp = 10 * xy[, 1] + stats::rnorm(nx * ny),
                  prec = 0.1 * xy[, 2] + stats::rnorm(nx * ny, sd = 0.05),
                  geometry = g)
}
