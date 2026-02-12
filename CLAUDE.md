# CLAUDE.md - pulseR Development Guide

## Project Overview

pulseR is an R package for **ecotypology creation** — spatially constrained fuzzy clustering to classify geographic areas (river catchments, segments) into ecological types based on environmental variables. It implements a pipeline from spatial polygon networks through fuzzy regionalization to river type classification using Gaussian Mixture Models.

- **Version**: 0.0.0.9000 (early development)
- **Author**: Jonathan Jupke
- **License**: MIT

## Repository Structure

```
pulseR/
├── R/                    # Package source code (11 files)
├── man/                  # Auto-generated roxygen2 documentation (.Rd files)
├── vignettes/            # Vignette templates (minimal content)
├── testing/              # Manual example/test scripts (not unit tests)
├── DESCRIPTION           # Package metadata and dependencies
├── NAMESPACE             # Auto-generated exports/imports (do not edit)
├── README.md             # Minimal readme
├── LICENSE / LICENSE.md  # MIT license
├── .Rbuildignore         # R build exclusions
└── pulseR.Rproj          # RStudio project configuration
```

No compiled code (`src/`), no formal test suite (`tests/testthat/`), no CI/CD configuration.

## Build & Development Commands

```bash
# Install package dependencies
Rscript -e 'install.packages(c("sf", "spdep", "igraph", "cluster", "mclust", "parallelDist", "dplyr", "tidyr", "purrr", "tidyselect", "ggplot2", "ggtern", "cowplot", "tricolore", "knitr", "rmarkdown"))'

# Build and check package
R CMD build .
R CMD check pulseR_0.0.0.9000.tar.gz

# Or via devtools (preferred workflow)
Rscript -e 'devtools::document()'   # Regenerate NAMESPACE and man/ from roxygen2
Rscript -e 'devtools::load_all()'   # Load package for interactive development
Rscript -e 'devtools::check()'      # Full R CMD check
Rscript -e 'devtools::install()'    # Install locally
```

### Key Development Notes

- **NAMESPACE is auto-generated** by roxygen2. Never edit it by hand. Run `devtools::document()` after changing any `@export`, `@importFrom`, or other roxygen tags.
- **man/*.Rd files are auto-generated** by roxygen2. Never edit them by hand.
- The RStudio project file (`pulseR.Rproj`) configures roxygen to auto-generate `rd`, `collate`, and `namespace`.
- No formal test suite exists. The `testing/` directory contains manual example scripts that depend on external data files not included in the repository.

## Architecture & Workflow Pipeline

The package implements a hierarchical spatial ecological classification pipeline:

```
sf Polygons
    │  polygon_to_network() / dev_polygon_to_network()
    ▼
Spatial Network (igraph) with polygon data
    │  add_edge_weight()
    ▼
Weighted Network (edges weighted by environmental dissimilarity)
    │  skater_con()
    ▼
Fuzzy Region Memberships (matrix: nodes × regions)
    ├── plot_fuzzy_area()           → Visualize regions
    └── get_core_regions()          → Extract high-membership contiguous areas
            │  river_types()
            ▼
        Fuzzy River Types (weighted GMM probabilities)
            ├── plot_fuzzy_area()       → Visualize types
            ├── plot_river_centroids()  → PCA of cluster centroids
            └── bhattacharyya_dist()    → Distribution similarity
```

## Public API (10 Exported Functions)

| Function | File | Purpose |
|---|---|---|
| `polygon_to_network()` | `R/polygons_to_network.R` | Create spatial network from sf polygons (Queen contiguity) |
| `dev_polygon_to_network()` | `R/dev_polygons_to_network.R` | Memory-optimized version avoiding `nb2mat` |
| `add_edge_weight()` | `R/add_edge_weight.R` | Weight edges by scaled Euclidean distance of environmental variables |
| `connect_components()` | `R/connect_components.R` | Connect disconnected graph components via minimum spanning tree |
| `skater_con()` | `R/skate.R` | Fuzzy spatially-constrained clustering (SKATER ensemble + PAM + fuzzy memberships) |
| `get_core_regions()` | `R/get_core_regions.R` | Extract contiguous regions exceeding a membership threshold |
| `river_types()` | `R/river_types.R` | Train regional GMMs and classify river types with weighted probabilities |
| `plot_fuzzy_area()` | `R/plot_area.R` | Visualize fuzzy classifications (RGB, ternary, small multiples, dominant) |
| `plot_river_centroids()` | `R/meta_cluster_pca.R` | PCA visualization of cluster centroids |
| `bhattacharyya_dist()` | `R/meta_cluster_pca.R` | Bhattacharyya distance between normal distributions |

## Code Conventions

### Style

- **Indentation**: 8-space tabs (matching RStudio project config)
- **Naming**: `snake_case` for functions and variables (e.g., `skater_con`, `fuzzy_memberships`, `core_region_id`)
- **Assignment**: `<-` operator (not `=`)
- **Pipe**: Native R pipe `|>` in newer code (e.g., `plot_area.R`), `purrr::map2` for functional patterns
- **Explicit namespacing**: Functions from dependencies are called with `package::function()` syntax in function bodies (e.g., `igraph::components()`, `sf::st_drop_geometry()`)
- **Return values**: Explicit `return()` at end of functions

### Roxygen2 Documentation

All exported functions use comprehensive roxygen2 documentation with markdown enabled:

- `#'` comment prefix for all documentation blocks
- `@param` with type info and detailed descriptions
- `@return` with `\describe{}` blocks for structured outputs
- `@details` with `\enumerate{}` or `\itemize{}` for algorithm steps
- `@examples` wrapped in `\dontrun{}` (examples require external data)
- `@importFrom package function` for each specific import (no blanket `@import`)
- `@seealso` cross-references to related functions and dependency functions
- `@references` for academic citations where applicable
- Mathematical notation via `\eqn{}` and `\deqn{}`
- `@export` tag on all public functions

### Data Structures

- **Input polygons**: `sf` objects with environmental variable columns
- **Network output**: List with `$graph` (igraph) and `$polygons` (sf) components
- **Fuzzy memberships**: Numeric matrix (rows = nodes, columns = clusters), rows sum to 1
- **Core regions**: Named list of sf objects, one per class
- **River types output**: List with `$rivertypes` (weighted probability matrix) and `$models` (fitted Mclust objects)

### NSE (Non-Standard Evaluation)

Global variables used in ggplot2/dplyr expressions are declared in `R/globals.R` via `utils::globalVariables()` to avoid R CMD check NOTEs.

## Key Dependencies

| Package | Role |
|---|---|
| `sf` | Spatial data (polygons, geometries, spatial predicates) |
| `spdep` | Spatial weights, neighbor structures, Queen/Rook contiguity |
| `igraph` | Graph construction, spanning trees, component detection, MST |
| `cluster` | PAM (Partitioning Around Medoids) clustering |
| `mclust` | Gaussian Mixture Models for river type classification |
| `parallelDist` | Parallel Hamming distance for coassociation matrices |
| `dplyr` / `tidyr` / `purrr` / `tidyselect` | Data manipulation |
| `ggplot2` / `ggtern` / `cowplot` / `tricolore` | Visualization |

## Typical Usage Pattern

From the testing scripts (`testing/angerman.R`):

```r
# 1. Load and prepare spatial data
z <- st_read_parquet("catchments.parquet")

# 2. Build spatial network
g <- dev_polygon_to_network(z)

# 3. Weight edges by environmental dissimilarity
g2 <- add_edge_weight(g, id_col = "ID")

# 4. Fuzzy regionalization
r <- skater_con(g2, final_regions = 3, n.rst = 100)

# 5. Visualize regions
plot_fuzzy_area(g$polygons, fuzzy_memberships = r, plot = "dominant")

# 6. Extract core regions and classify river types
ca <- get_core_regions(cbind(g$polygons, r), membership_cols = paste0("X", 1:3), cutoff = 0.8)
rt <- river_types(core_regions = ca, membership_id = r, all_data = z,
                  membership_cols = paste0("X", 1:3), non_value_cols = "ID")

# 7. Visualize river types
plot_fuzzy_area(g$polygons, fuzzy_memberships = rt$rivertypes, plot = "dominant")
plot_river_centroids(models = rt$models)
```

## Important Caveats

- `dev_polygon_to_network()` is the preferred version over `polygon_to_network()` — it avoids a memory leak caused by `spdep::nb2mat` on large datasets.
- The `river_types()` function has hardcoded handling for geological compositional variables (`area_sediment`, `area_calcareous`, `area_siliceous`) to avoid covariance singularity.
- Ternary plots (`plot = "ternary1"` / `"ternary2"`) and RGB plots require exactly 3 classes.
- No test data is shipped with the package. Testing scripts reference external parquet files.
