#' Classify local Types using Regional Gaussian Mixture Models
#'
#' Trains a Gaussian Mixture Model (GMM) on each of several user-defined
#' "core regions" and projects those classifications onto a larger target
#' dataset. For every observation the predicted local-type probabilities from
#' each regional model are weighted by that observation's (fuzzy or crisp)
#' membership in the corresponding region and concatenated into a single wide
#' probability matrix.
#'
#' @param graph The output of \code{\link{polygon_to_network}}. Contains a, 
#'   \code{igraph} graph object and an \code{sf} object containing the complete 
#'   dataset to be classified. 
#'   Any \code{sf} geometry column is dropped automatically.
#' @param core_regions A non-empty list of \code{sf} objects or
#'   \code{data.frame}s. Each element holds the training observations for one
#'   region. All elements must contain the same predictor columns.
#' @param cutoff numeric in range \code{0,1}`. Sets the membership degree threshold 
#'   for deriving core regions. Directly passed to \code{\link{get_core_regions}} 
#'   inside the function.
#' @param n_local_types Integer scalar or vector passed to the \code{G} argument
#'   of \code{\link[mclust]{Mclust}}; the number(s) of mixture components to
#'   consider. When a vector is supplied, \code{Mclust} selects the best value
#'   by BIC. Default \code{1:9}.
#' @param regions The output of \code{\link{get_regions}} or \code{\link{tune_regions}}. 
#' @param non_value_cols Character vector of identifier columns to drop before
#'   training (e.g. primary keys). Default \code{"ID"}. The column
#'   \code{"core_region_id"} is always dropped from the training data.
#' @param crisp Logical. When \code{TRUE} both the regional memberships and the
#'   per-region local-type probabilities are hardened to 0/1 (each observation
#'   assigned to its single most probable region and, within a region, its
#'   single most probable local type) before weighting. Default \code{FALSE}.
#' @param verbose Logical. When \code{TRUE} progress over regions is reported
#'   via \code{\link{message}}. Default \code{FALSE}.
#' @param modelNames Character vector indicating the models to be fitted in the 
#'   EM phase of clustering. Is passed directly to mclust::Mclust, see there for
#   for more details
#'
#' @details
#' Pre-processing applied before fitting:
#' \itemize{
#'   \item \code{sf} geometry is removed with \code{\link[sf]{st_drop_geometry}}
#'     (a no-op for plain data frames).
#'   \item Identifier and prior-membership columns are removed.
#'   \item Columns that are constant in \emph{any} region are removed from
#'     \emph{all} regions, since a zero-variance predictor produces a singular
#'     covariance matrix in \code{Mclust}.
#' }
#'
#' For region \eqn{r} and local type \eqn{k} the returned weighted probability
#' is \eqn{P(\mathrm{Type}_k \mid \mathrm{Region}_r)\, W_r}, where \eqn{W_r} is
#' the membership weight in column \eqn{r} of \code{membership_id}.
#'
#' @return A named list with:
#' \describe{
#'   \item{\code{localtypes}}{Numeric matrix of membership-weighted local-type
#'     probabilities, one row per observation in \code{graph$polygons}. Columns are
#'     named \code{"Region{i}localType{k}"}.}
#'   \item{\code{localtypes_raw}}{List (one element per region, named
#'     \code{"Region{i}"}) of the \emph{unweighted} predicted probability
#'     matrices.}
#'   \item{\code{models}}{List of the fitted \code{Mclust} objects, named by
#'     region.}
#' }
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' make_region <- function(n, shift) {
#'   data.frame(
#'     ID             = seq_len(n),
#'     core_region_id = 1,
#'     X1 = stats::runif(n), X2 = stats::runif(n), # prior-membership cols (dropped)
#'     temp     = stats::rnorm(n, shift),
#'     slope    = stats::rnorm(n, -shift),
#'     rainfall = stats::rnorm(n)
#'   )
#' }
#' core_regions  <- list(make_region(100, 0), make_region(100, 3))
#' graph <- list()
#' regions <- list()
#' graph$polygons <- rbind(make_region(50, 0), make_region(50, 3))
#' 
#' memb <- matrix(stats::runif(nrow(graph$polygons) * 2), ncol = 2)
#' memb <- memb / rowSums(memb)                    # fuzzy weights summing to 1
#' regions$memberships <- memb
#' 
#' res <- get_local_types(graph = graph,
#'                        core_regions = core_regions,
#'                        regions = regions,
#'                        n_local_types = 1:4)
#' head(res$localtypes)
#' }
#'
#' @importFrom sf st_drop_geometry
#' @importFrom mclust Mclust mclustBIC
#' @importFrom stats predict complete.cases
#'
#' @export
get_local_types <- function(graph,
                        core_regions = NULL,
                        cutoff = NULL,
                        n_local_types = 1:9,
                        regions,
                        non_value_cols  = "ID",
                        crisp           = FALSE,
                        verbose         = FALSE,
                        modelNames = NULL) {
        
        ## ---- 0. Input validation -------------------------------------------*
        if (is.null(core_regions) & is.null(cutoff)){
                stop("Either core_regions or cutoff must be specified.")
        } else if (!is.null(core_regions) & !is.null(cutoff)){
                stop("You should only specify core_regions (a list of sf objects created with get_core_regions()) or cutoff. The latter option creates core_regions on the fly in the function.")
        }
        if (!is.data.frame(graph$polygons)) {
                stop("`graph$polygons` must be a data.frame or sf object.", call. = FALSE)
        }
        if (!all(vapply(core_regions, is.data.frame, logical(1)))) {
                stop("Every element of `core_regions` must be a data.frame or sf object.",
                     call. = FALSE)
        }
        if (!is.logical(crisp) || length(crisp) != 1L || is.na(crisp)) {
                stop("`crisp` must be a single TRUE/FALSE value.", call. = FALSE)
        }
        if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
                stop("`verbose` must be a single TRUE/FALSE value.", call. = FALSE)
        }
        if (!is.numeric(n_local_types) || length(n_local_types) == 0L ||
            any(n_local_types < 1L) || any(n_local_types != as.integer(n_local_types))) {
                stop("`n_local_types` must be one or more positive integers.", call. = FALSE)
        }
        if (!is.null(non_value_cols) && !is.character(non_value_cols)) {
                stop("`non_value_cols` must be a character vector or NULL.", call. = FALSE)
        }
        # if (!is.null(membership_cols) && !is.character(membership_cols)) {
        #         stop("`membership_cols` must be a character vector or NULL.", call. = FALSE)
        # }

        n_regions <- ncol(regions$memberships)
        ## Coerce membership weights to a numeric matrix for predictable indexing,
        ## then check that its dimensions are consistent with the other inputs.
        membership_id <- as.matrix(regions$memberships)
        if (!is.numeric(membership_id)) {
                stop("`membership_id` must be numeric.", call. = FALSE)
        }
        if (nrow(membership_id) != nrow(graph$polygons)) {
                stop(sprintf(
                        "`membership_id` has %d rows but `graph$polygons` has %d; they must match.",
                        nrow(membership_id), nrow(graph$polygons)), call. = FALSE)
        }
        
        ## Local helper: harden each row of a numeric matrix to a 0/1 indicator of
        ## its single largest value (ties -> first column; an all-NA row -> all 0).
        ## Vectorised and shape-stable for any number of columns (incl. one).
        harden_rows <- function(m) {
                m   <- as.matrix(m)
                out <- matrix(0, nrow = nrow(m), ncol = ncol(m))
                j   <- max.col(replace(m, is.na(m), -Inf), ties.method = "first")
                out[cbind(seq_len(nrow(m)), j)] <- 1
                all_na <- rowSums(!is.na(m)) == 0L
                if (any(all_na)) out[all_na, ] <- 0
                out
        }

        ## ---- 1. Pre-processing --------------------------------------------- *
        ## Finding Core regions 
        if (is.null(core_regions)){
                # No core regions where provided so we will compute them here. 
                core_regions <- get_core_regions(
                        graph = graph, 
                        regions = regions,
                        cutoff = cutoff
                        )
                
        }
        
        ## Drop sf geometry (no-op on a plain data.frame) and standardise to data.frame.
        train_list   <- lapply(core_regions, function(x) as.data.frame(sf::st_drop_geometry(x)))
        predict_data <- as.data.frame(sf::st_drop_geometry(graph$polygons))
        
        ## Columns that must never enter a model: the region key, the prior-membership
        ## columns and any user-supplied identifier columns.
        drop_cols  <- unique(c("core_region_id", non_value_cols, paste0("region",1:n_regions)))
        train_list <- lapply(train_list, function(x) x[, !names(x) %in% drop_cols, drop = FALSE])
        
        ## Prediction data only needs the identifier columns removed; the predictors
        ## actually used are selected by name per region further below.
        predict_data <- predict_data[, !names(predict_data) %in% non_value_cols, drop = FALSE]
        
        ## All regions must share the same predictor set (the global constant-column
        ## logic below assumes this).
        ref_cols <- names(train_list[[1L]])
        if (!all(vapply(train_list, function(x) setequal(names(x), ref_cols), logical(1)))) {
                stop("All core regions must contain the same predictor columns.", call. = FALSE)
        }
        
        ## ---- 1b. Optional crisp (hard) regional membership ------------------------
        if (crisp) membership_id <- harden_rows(membership_id)
        
        ## ---- 2. Drop constant & collinear compositional predictors ----------------
        
        ## A column constant within ANY region is dropped from ALL regions: a
        ## zero-variance predictor yields a singular covariance matrix in Mclust.
        ## NB: a column that is entirely NA also has a single unique value and is
        ## therefore also flagged here.
        const_names <- character(0)
        for (i in seq_len(n_regions)) {
                n_unique    <- vapply(train_list[[i]], function(x) length(unique(x)), integer(1))
                const_names <- union(const_names, names(n_unique)[n_unique <= 1L])
        }
        
        ## Remove constant variables by name
        if (length(const_names)) {
                train_list <- lapply(train_list, function(x) x[, !names(x) %in% const_names, drop = FALSE])
        }
        
        ## Predictors must be numeric and complete; otherwise Mclust fails cryptically.
        for (i in seq_len(n_regions)) {
                if (ncol(train_list[[i]]) == 0L) {
                        stop(sprintf("Region %d has no predictor columns left after pre-processing.", i),
                             call. = FALSE)
                }
                non_num <- !vapply(train_list[[i]], is.numeric, logical(1))
                if (any(non_num)) {
                        stop(sprintf("Region %d has non-numeric predictor column(s): %s.",
                                     i, paste(names(train_list[[i]])[non_num], collapse = ", ")),
                             call. = FALSE)
                }
                if (anyNA(train_list[[i]])) {
                        stop(sprintf("Region %d training data contains missing values; Mclust requires complete cases.", i),
                             call. = FALSE)
                }
        }
        
        ## ---- 3. Fit one Gaussian Mixture Model per region -------------------------
        ## (No scaling: Mclust models the full covariance structure, so standardising
        ## the predictors is unnecessary. The variables used to train each region are
        ## recorded explicitly rather than recovered from the fitted object, which is
        ## robust to the univariate case where `parameters$mean` has no rownames.)
        models     <- vector("list", n_regions)
        train_vars <- vector("list", n_regions)
        
        for (i in seq_len(n_regions)) {
                if (verbose) message(sprintf("Fitting GMM for region %d/%d ...", i, n_regions))
                
                train_vars[[i]] <- names(train_list[[i]])
                
                ## IMPORTANT: Mclust() rebuilds its arguments into a call to mclustBIC() and
                ## evaluates that call in the calling frame (this function's frame), not in
                ## mclust's own namespace. So mclustBIC must be reachable from here!
                fit <- mclust::Mclust(train_list[[i]], G = n_local_types, modelNames = modelNames)
                if (is.null(fit)) {
                        stop(sprintf("Mclust could not fit a model for region %d over G = %s.",
                                     i, paste(range(n_local_types), collapse = ":")), call. = FALSE)
                }
                models[[i]] <- fit
        }
        names(models) <- paste0("Region", seq_len(n_regions))
        
        ## ---- 4. Predict on `graph$polygons` and weight by regional membership -----------
        raw      <- vector("list", n_regions)
        weighted <- vector("list", n_regions)
        
        for (i in seq_len(n_regions)) {
                vars_i <- train_vars[[i]]
                
                ## Every training variable must be present in the target data.
                missing_vars <- setdiff(vars_i, names(predict_data))
                if (length(missing_vars)) {
                        stop(sprintf("`graph$polygons` is missing variable(s) required by region %d: %s.",
                                     i, paste(missing_vars, collapse = ", ")), call. = FALSE)
                }
                
                model_data <- predict_data[, vars_i, drop = FALSE]
                
                non_num <- !vapply(model_data, is.numeric, logical(1))
                if (any(non_num)) {
                        stop(sprintf("`all_data` column(s) used by region %d are not numeric: %s.",
                                     i, paste(vars_i[non_num], collapse = ", ")), call. = FALSE)
                }
                
                ## Mclust prediction needs complete, finite rows. Flag offenders in a
                ## vectorised way (avoids a slow row-wise apply on large datasets) and
                ## report how many and where.
                mat <- as.matrix(model_data)
                bad <- !stats::complete.cases(mat) | rowSums(!is.finite(mat)) > 0L
                if (any(bad)) {
                        offenders <- which(bad)
                        stop(sprintf(
                                "%d row(s) in `graph$polygons` have missing/non-finite values for region %d predictors (%s). First offending row(s): %s.",
                                length(offenders), i, paste(vars_i, collapse = ", "),
                                paste(offenders[seq_len(min(6L, length(offenders)))], collapse = ", ")),
                             call. = FALSE)
                }
                
                ## predict.Mclust is dispatched through stats::predict. $z holds
                ##  the posterior cluster probabilities as an n-by-G matrix.
                pred <- stats::predict(models[[i]], newdata = model_data)$z
                ## guarantee matrix form (e.g. when G = 1)
                pred <- as.matrix(pred)   
                
                if (crisp) {
                        pred <- harden_rows(pred)
                }
                
                raw[[i]] <- pred
                
                ## Weight every local-type column by this region's membership weight. The
                ## length-n vector recycles down the columns, scaling each row by its weight.
                w <- pred * membership_id[, i]
                colnames(w) <- paste0("R", i, "\u00b7","T", seq_len(ncol(w)))
                weighted[[i]] <- w
        }
        names(raw) <- paste0("Region", seq_len(n_regions))
        
        ## ---- 5. Assemble output ---------------------------------------------------
        out <- do.call(cbind, weighted)
        
        list(
                localtypes     = out,
                localtypes_raw = raw,
                models         = models
        )
}


#' Determine Core Regions from Fuzzy Membership Degrees
#'
#' @description
#' Identifies *core regions* for each region from per-polygon membership
#' degrees. A core region is a set of spatially contiguous polygons whose
#' membership degree for that class meets or exceeds a threshold (`cutoff`).
#'
#' @param graph The output of \code{\link{polygon_to_network}}.
#' @param regions The result of \code{\link{get_regions}} or \code{.$best_result} of the result 
#'   \code{\link{tune_regions}}. 
#' @param cutoff Single numeric in `[0, 1]`. Minimum membership degree for a
#'   polygon to be included in a class's core region. Default `0.5`.
#' @param queen Logical scalar. If `TRUE` (default) neighbours are defined by
#'   queen contiguity (polygons sharing an edge *or* a single vertex). If
#'   `FALSE`, rook contiguity is used (only polygons sharing an edge of
#'   positive length).
#'
#' @return
#' A named `list` with one element per membership column. Each element is an
#' `sf` object containing the polygons that belong to that class's core
#' region(s), with an additional **first** column `core_region_id` (integer)
#' labelling distinct contiguous regions. Classes with no qualifying polygons
#' yield a zero-row `sf` object with the same columns. The `core_region_id`
#' labels are arbitrary integers and are only comparable within a single class.
#'
#' @details
#' For each membership column the function:
#' \enumerate{
#'   \item selects polygons whose membership degree is `>= cutoff`
#'     (missing values are treated as *not* meeting the threshold);
#'   \item builds a sparse spatial neighbour list using the requested
#'     contiguity rule;
#'   \item derives an undirected graph and extracts its connected components,
#'     each of which is one contiguous core region;
#'   \item returns the selected polygons annotated with their region id.
#' }
#'
#' Contiguity is evaluated with DE-9IM relate patterns: `"F***T****"` for queen
#' (boundaries meet at a point *or* a line) and `"F***1****"` for rook
#' (boundaries meet along a line only). The neighbor list is kept sparse so the
#' routine scales to large polygon sets without materializing an
#' \eqn{n \times n} adjacency matrix. Contiguity is only meaningful for
#' polygonal geometries.
#'
#' @examples
#' \dontrun{
#' # 'fuzzy_data' is an sf object with one membership column per class.
#' core_regions <- get_core_regions(
#'   graph = g,
#'   regions = regions,
#'   cutoff = 0.6,
#'   queen = TRUE
#' )
#'
#' # Core regions for class_A (sf object, 'core_region_id' is the first column):
#' core_regions$class_A
#' }
#'
#' @importFrom sf st_drop_geometry st_relate
#' @importFrom igraph graph_from_adj_list components simplify
#' @export
get_core_regions <- function(graph,
                             regions,
                             cutoff = 0.5,
                             queen = TRUE,
                             typicality = FALSE) {
        
        
        x <- graph$polygons
        U <- regions$memberships
        Typi <- regions$typicality
        ## ---- Argument validation -----------------------------------*
        
        
        # `x` must be an sf object; we rely on its geometry and sf subsetting rules.
        if (!inherits(x, "sf")) {
                stop("'x' must be an sf object.", call. = FALSE)
        }
        
        # `cutoff` must be a single, non-missing number in the unit interval, because
        # membership degrees are fuzzy weights in [0, 1]. (An NA here would otherwise
        # crash the range test below, which combines NAs with `||`.)
        if (!is.numeric(cutoff) || length(cutoff) != 1L || is.na(cutoff)) {
                stop("'cutoff' must be a single non-missing numeric value.", call. = FALSE)
        }
        if (cutoff < 0 || cutoff > 1) {
                stop("'cutoff' must be between 0 and 1.", call. = FALSE)
        }
        
        # `queen` selects the contiguity rule and must be an unambiguous flag.
        if (!is.logical(queen) || length(queen) != 1L || is.na(queen)) {
                stop("'queen' must be a single logical value (TRUE or FALSE).",
                     call. = FALSE)
        }
        
        # Warn before silently overwriting a pre-existing column of the same name.
        if ("core_region_id" %in% names(x)) {
                warning("'x' already contains a 'core_region_id' column; ",
                        "it will be overwritten in the output.", call. = FALSE)
        }
        
        # Resolve the active geometry column by name
        geom_col <- attr(x, "sf_column")
        #x <- x[,geom_col]
        colnames(U) <- paste0("region",1:ncol(U))
        x <- cbind(U,x,Typi)
        rangeBound <- all(U <= 1) & all(U>=0)
        if (rangeBound==FALSE){
                warning("Some membership values fall outside [0, 1]; check that ",
                        "'membership_cols' hold fuzzy membership degrees.", call. = FALSE)
        }
        
        ## ---- Local helpers -------------------------------------------*
        
        # DE-9IM relate pattern for the chosen contiguity rule.
        #   queen  "F***T****": interiors disjoint, boundaries meet at a point OR line
        #   rook   "F***1****": interiors disjoint, boundaries meet along a line only
        # (The leading "F" excludes self-matches, since a polygon's interior
        #  intersects itself, so no manual diagonal removal is required.)
        contiguity_pattern <- if (queen) "F***T****" else "F***1****"
        
        # Move `core_region_id` to the front, preserve all other attribute columns,
        # and keep the (possibly non-"geometry") geometry column last. Used for both
        # populated and empty results so every list element shares one column order.
        place_id_first <- function(obj) {
                middle <- setdiff(names(obj), c("core_region_id", geom_col))
                obj[, c("core_region_id", middle, geom_col)]
        }
        
        ## ---- Per-class core region detection --------------------------*
        
        core_regions_list <- vector("list", ncol(U))
        names(core_regions_list) <- colnames(U)#paste0("Region",1:ncol(U))
        
        for (class_col in seq_along(colnames(U))) {
                
                # Polygons meeting the threshold. Coercing NA comparisons to FALSE excludes
                # missing memberships; otherwise NA logical indices would inject all-NA
                # "phantom" rows when used to subset `x`.
                
                x2 <- st_drop_geometry(x)
                itername <- colnames(U)[class_col]
                x3 <- x2[, "Typi"]
                x2 <- x2[,itername]
                 
                if (!typicality) {
                        core_mask <- !is.na(x2) & x2 >= cutoff
                } else {
                        core_mask <- !is.na(x2) & x2 >= cutoff & x3 > cutoff
                }
                
                        
                        
                if (!any(core_mask)) {
                        message("No polygons meet the cutoff for class '", class_col, "'.")
                        empty <- x[0, , drop = FALSE]
                        empty$core_region_id <- integer(0)
                        core_regions_list[[class_col]] <- place_id_first(empty)
                        next
                }
                
                core_polys <- x[core_mask, , drop = FALSE]
                
                if (nrow(core_polys) == 1L) {
                        # A lone qualifying polygon is trivially its own region; skip the (for
                        # n == 1, unnecessary) graph construction.
                        core_polys$core_region_id <- 1L
                } else {
                        nb <- sf::st_relate(core_polys,
                                            pattern = contiguity_pattern,
                                            sparse = TRUE)
                        
                        adj_graph <- igraph::simplify(
                                igraph::graph_from_adj_list(unclass(nb), mode = "all")
                        )
                        
                        comp <- igraph::components(adj_graph)
                        membership <- as.integer(comp$membership)
                        
                        # Keep ONLY the largest contiguous region. 
                        areas <- as.numeric(sf::st_area(core_polys))
                        largest_id <- as.integer(names(which.max(tapply(areas, membership, sum))))
                        core_polys  <- core_polys[membership == largest_id, , drop = FALSE]
                        core_polys$core_region_id <- class_col
                }
                
                core_regions_list[[class_col]] <- place_id_first(core_polys)
        }
        
        core_regions_list
}

