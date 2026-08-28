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
#'   for more details
#' @param typicality Logical scalar. If `FALSE` (default) the membership matrix 
#'   \code{.$memberships} is used to weight results of mixuter models If `TRUE` the 
#'   typicality matrix \code{.$typicality} is used instead. 
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
                        modelNames = NULL,
                        typicality = FALSE) {
        
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
        
        
        membership_id <- if (!typicality){
                as.matrix(regions$memberships)
        } else {
                as.matrix(regions$typicality)         
                }
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
                        graph   = graph, 
                        regions = regions,
                        cutoff  = cutoff
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
        
        ## ---- 1b. Optional crisp (hard) regional membership ------------------------*
        if (crisp) membership_id <- harden_rows(membership_id)
        
        ## ---- 2. Drop constant & collinear compositional predictors ----------------*
        
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
        
        ## ---- 3. Fit one Gaussian Mixture Model per region -------------------------*
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
        
        ## ---- 4. Predict on `graph$polygons` and weight by regional membership -----------*
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
                
                ## Weight every local-type column by this region's membership 
                ## weight or its typicality. The length-n vector recycles down 
                ## the columns, scaling each row by its weight.
                 
                w <- pred * membership_id[, i]
                colnames(w) <- paste0("R", i, "\u00b7","T", seq_len(ncol(w)))
                weighted[[i]] <- w
        }
        names(raw) <- paste0("Region", seq_len(n_regions))
        
        ## ---- 5. Assemble output ---------------------------------------------------*
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
#' @param typicality Logical scalar. If `FALSE` (default) the membership matrix 
#'   is used to compute core regions, i.e., the threshold in `cutoff` is used to 
#'   select rows in the membership matrix \code{.$memberships}. If `TRUE` the 
#'   typicality matrix \code{.$typicality} is used instead. 
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
        
        if (isTRUE(typicality)){
                U <- regions$typicality
        } else {
                U <- regions$memberships 
        }
        x <- graph$polygons
        # Typi <- regions$typicality
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
        x <- cbind(U,x)
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
               # x3 <- x2[, "Typi"]
                x2 <- x2[,itername]
                 
                # if (!typicality) {
                        core_mask <- !is.na(x2) & x2 >= cutoff
                # } else {
                        # core_mask <- !is.na(x2) & x2 >= cutoff & x3 > cutoff
                # }
                
                        
                        
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


#' Collapse local types across regions into meta-types
#'
#' Cuts the cross-region meta-clustering (see [plot_local_similarity()]) and
#' merges local types that fall in the same meta-cluster by summing their
#' membership columns. Meta-types are labelled sequentially (`meta.1`, `meta.2`,
#' ...) in increasing order of their meta-cluster id.
#'
#' Run this after inspecting [plot_local_similarity()] to commit to a particular
#' number of meta-types (`k`) or a dendrogram cut height (`h`).
#'
#' @param local A fitted local-type object (as returned by `get_local_types()`).
#'   Only used to (re)compute `similarity` when it is not supplied.
#' @param memberships The wide local-type membership matrix (catchments x types).
#'   Column names must resolve to the internal meta labels (`R<region>\u00b7T<type>`);
#'   columns named `Region<n>localType<g>` are relabelled to that form
#'   automatically.
#' @param k,h Cut the meta-clustering into `k` meta-types or at dendrogram height
#'   `h`. Supply exactly one of the two; supplying neither is an error.
#' @param similarity Optional precomputed `localtype_similarity` (from
#'   `.localtype_meta()`); recomputed from `local` when `NULL`.
#' @param distance Distance used only when `similarity` is recomputed: `"mean"`
#'   (component means), `"hellinger"` (mean and shape), or `"wasserstein"`.
#' @param sep Separator assumed between region and type in the membership column
#'   labels; used when building the mismatch error message.
#' @param ... Currently unused; reserved for future arguments.
#'
#' @return A list with components:
#'   \describe{
#'     \item{memberships}{Catchments x meta-types matrix of summed memberships.
#'       Rows are *not* renormalised, so they need not sum to 1.}
#'     \item{hard}{Hard meta-type label per catchment (argmax of the normalised
#'       memberships).}
#'     \item{entropy}{Shannon entropy of each catchment's normalised membership
#'       vector; higher values flag more ambiguous assignments.}
#'     \item{row_sums}{Row sums of `memberships` before normalisation.}
#'     \item{mapping}{Data frame linking each original type (`region`, `type`,
#'       `label`) to its meta-cluster (`meta`) and label (`meta_label`).}
#'     \item{n_before,n_after}{Counts of local types before, and meta-types
#'       after, collapsing.}
#'     \item{k,h}{The cut arguments that were used.}
#'   }
#' @seealso [plot_local_similarity()], [plot_local_embedding()]
#' @export
collapse_local <- function(local, memberships, k = NULL, h = NULL, similarity = NULL,
                           distance = c("mean", "hellinger", "wasserstein"),
                           sep = ".", ...) {
        distance <- match.arg(distance)
        if (is.null(k) && is.null(h)) stop("supply k or h.")
        # Build the cross-region similarity / meta-clustering unless one was passed in
        if (is.null(similarity)) similarity <- .localtype_meta(local, distance = distance, sep = sep)
        
        # Cut the meta-dendrogram: grp is a named integer vector (names = type labels,
        # value = meta-cluster id)
        grp  <-
                if (!is.null(k)){
                        stats::cutree(similarity$hclust, k = k)
                } else {
                        stats::cutree(similarity$hclust, h = h)
                }
        info <- similarity$info
        info$meta <- grp[info$label]                   # attach each label's meta-cluster id
        
        # name each meta-type sequentially: meta.1, meta.2, ... in ascending id order
        ids <- sort(unique(info$meta))
        med <- paste("meta", seq_along(ids), sep = ".")
        names(med) <- ids
        info$meta_label <- med[as.character(info$meta)]
        
        # Align membership column names to the internal labels: e.g. "Region3localType2"
        # becomes "R3\u00b7T2". Columns already in that form are left untouched by sub().
        colnames(memberships) <-
                sub(pattern = "^Region([0-9]+)localType([0-9]+)$", replacement = "R\\1\u00b7T\\2", x = colnames(memberships))
        
        # NOTE: the labels use a middle dot ("\u00b7"); `sep` only formats the message below.
        if (!all(info$label %in% colnames(memberships)))
                stop("membership columns must match the meta labels (region", sep, "type). ",
                     "Confirm the get_local_types() column naming and pass a matching matrix or relabel.")
        
        # Restrict and order the membership matrix to the known labels, then map each
        # column to the meta-type it belongs to
        M   <- as.matrix(memberships[, info$label, drop = FALSE])
        key <- info$meta_label[match(colnames(M), info$label)]
        
        # Sum the member columns within each meta-type -> catchments x meta-types
        collapsed <- vapply(split(seq_len(ncol(M)), key),
                            function(j) rowSums(M[, j, drop = FALSE]), numeric(nrow(M)))
        
        collapsed <- collapsed[, unique(key), drop = FALSE]  # restore first-appearance order
        # rows are not guaranteed to sum to 1
        rs   <- rowSums(collapsed)
        # normalise only for hard label and entropy (the returned memberships stay raw)
        prob <- collapsed / ifelse(rs > 0, rs, 1)
        hard <- colnames(collapsed)[max.col(prob, ties.method = "first")]  # argmax meta-type
        ent  <- -rowSums(ifelse(prob > 0, prob * log(prob), 0))            # Shannon entropy
        
        list(memberships = collapsed, hard = hard, entropy = ent, row_sums = rs,
             mapping = info[order(info$meta), c("region", "type", "label", "meta", "meta_label")],
             n_before = nrow(info), n_after = ncol(collapsed), k = k, h = h)
}

#' Plot cross-region similarity of local types
#'
#' Compares the GMM components of every region's local types in a common,
#' standardised predictor space and draws a meta-clustered distance heatmap.
#' Off-diagonal blocks that mix regions flag local types that are candidates for
#' collapsing with [collapse_local()].
#'
#' @param local A fitted local-type object (uses `local$models`).
#' @param distance Component distance: `"mean"` (standardised component means,
#'   matches the fingerprint plot), `"hellinger"` (mean and shape), or
#'   `"wasserstein"`.
#' @param k,h Optional meta-clustering cut used *only* to colour the axis labels:
#'   `k` meta-types or dendrogram height `h`.
#' @param groups Optional grouping that overrides `k`/`h` for label colouring.
#'   Either a `collapse_local()` result (its `$mapping` is used) or a named
#'   vector mapping each label to a group.
#' @param triangle If `TRUE` (default) draw only the upper triangle; otherwise
#'   the full off-diagonal matrix. The diagonal is always dropped.
#' @param ... Passed to the internal `.localtype_meta()` (e.g. `scale`, `ridge`).
#'
#' @return A `ggplot` object: the cross-region similarity heatmap.
#' @seealso [collapse_local()], [plot_local_embedding()]
#' @export
plot_local_similarity <- function(local, distance = c("mean", "hellinger", "wasserstein"),
                                  k = NULL, h = NULL, groups = NULL, triangle = TRUE, ...) {
        .region_require("ggplot2")
        distance <- match.arg(distance)
        meta <- .localtype_meta(local, distance = distance, ...)
        ord  <- meta$hclust$order                      # dendrogram order, so blocks read along the diagonal
        labs <- meta$info$label[ord]
        n    <- length(labs)
        M    <- meta$matrix[ord, ord]
        
        ## reshape to long form; keep upper triangle, diagonal always dropped
        gi   <- expand.grid(i = seq_len(n), j = seq_len(n))
        gi$d <- M[cbind(gi$i, gi$j)]
        gi   <- if (triangle) gi[gi$j > gi$i, ] else gi[gi$i != gi$j, ]
        gi$a <- factor(labs[gi$i], levels = labs)
        gi$b <- factor(labs[gi$j], levels = labs)
        
        ## groups: a collapse_local() result, a named vector, or an internal cut
        grp <- if (is.list(groups) && !is.null(groups$mapping))
                stats::setNames(groups$mapping$meta_label, groups$mapping$label)[labs]
        else if (!is.null(groups)) groups[labs]
        else if (!is.null(k))      stats::cutree(meta$hclust, k = k)[labs]
        else if (!is.null(h))      stats::cutree(meta$hclust, h = h)[labs]
        else                       NULL
        
        ## colour each label by its group; merged types share a colour, singletons stay grey
        lab_col <- rep("grey35", n)
        if (!is.null(grp)) {
                tab   <- table(grp)
                multi <- names(tab)[tab >= 2]          # only groups with >= 2 members get a colour
                if (length(multi)) {
                        pal <- stats::setNames(grDevices::hcl.colors(length(multi), "Dark 3"), multi)
                        hit <- !is.na(grp) & grp %in% multi
                        lab_col[hit] <- pal[as.character(grp[hit])]
                }
        }
        
        ggplot2::ggplot(gi, ggplot2::aes(.data$a, .data$b, fill = .data$d)) +
                ggplot2::geom_tile() +
                ggplot2::scale_fill_gradient(low = "#08306B", high = "grey92", name = "distance") +
                ggplot2::scale_x_discrete(drop = FALSE) +
                ggplot2::scale_y_discrete(drop = FALSE) +
                ggplot2::coord_equal() +
                ggplot2::labs(title = "Cross-region local-type similarity",
                              subtitle = sprintf("%s distance between GMM components; shared label colour marks a collapse group", distance),
                              x = NULL, y = NULL) +
                theme_eco() +
                ggplot2::theme(panel.grid  = ggplot2::element_blank(),
                               axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5,
                                                                   colour = lab_col),
                               axis.text.y = ggplot2::element_text(colour = lab_col))
}

#' Per-component covariance matrices from an mclust fit
#'
#' Normalises the several shapes `mclust` uses to store variances (per-component
#' `sigma`, shared `Sigma`, or univariate `sigmasq`) into a list of `G` labelled
#' covariance matrices. Any `NULL` element in the result is the caller's signal
#' to fall back to a mean-only distance.
#'
#' @param fit A fitted `mclust` model for one region's local types.
#' @param vn Character vector of variable names used to label the matrices.
#' @return List of length `G` of d x d covariance matrices (entries may be NULL).
#' @noRd
.gmm_sigma <- function(fit, vn) {
        v <- fit$parameters$variance
        G <- if (is.null(dim(fit$parameters$mean))) {  # univariate: mean is a length-G vector
                length(fit$parameters$mean)
        }  else {
                ncol(fit$parameters$mean)
        }
        out <- vector("list", G)
        if (!is.null(v$sigma)) {                       # d x d x G (most multivariate models)
                for (g in seq_len(G)) { S <- v$sigma[, , g]; dimnames(S) <- list(vn, vn); out[[g]] <- S }
        } else if (!is.null(v$Sigma)) {                # shared covariance across components
                S <- v$Sigma; dimnames(S) <- list(vn, vn); for (g in seq_len(G)) out[[g]] <- S
        } else if (!is.null(v$sigmasq)) {              # univariate variance(s)
                s2 <- v$sigmasq; if (length(s2) == 1L) s2 <- rep(s2, G)
                for (g in seq_len(G)) out[[g]] <- matrix(s2[g], 1, 1, dimnames = list(vn, vn))
        }
        out                                            # any NULL element triggers the mean-only fallback
}

#' Symmetric matrix square root via eigendecomposition
#'
#' Returns the symmetric PSD square root of `S`, clamping tiny negative
#' eigenvalues to 0 for numerical safety.
#'
#' @param S A symmetric (covariance) matrix.
#' @return The symmetric square-root matrix of `S`.
#' @noRd
.mat_sqrt <- function(S) {
        e <- eigen(S, symmetric = TRUE); e$vectors %*% diag(sqrt(pmax(e$values, 0)), nrow(S)) %*% t(e$vectors)
}

#' Distance between two GMM components
#'
#' Computes the `"mean"` (Euclidean between standardised means), `"wasserstein"`
#' (closed-form Gaussian W2), or `"hellinger"` (via the Bhattacharyya
#' coefficient) distance. Falls back to the mean distance whenever either
#' covariance is missing. A ridge proportional to the mean variance stabilises
#' the matrix inverse / square root.
#'
#' @param a,b Component descriptors with standardised mean `$m` and optional
#'   standardised covariance `$S`.
#' @param distance One of `"mean"`, `"hellinger"`, `"wasserstein"`.
#' @param ridge Relative ridge added to the covariance diagonals.
#' @return A scalar distance.
#' @noRd
.pair_dist <- function(a, b, distance, ridge) {
        dm <- a$m - b$m
        # mean distance, or fall back to it when either covariance is unavailable
        if (distance == "mean" || is.null(a$S) || is.null(b$S)) return(sqrt(sum(dm^2)))
        Sa <- a$S; Sb <- b$S
        rg <- ridge * mean(c(diag(Sa), diag(Sb)))      # ridge scaled to the typical variance
        Sa <- Sa + diag(rg, nrow(Sa)); Sb <- Sb + diag(rg, nrow(Sb))
        if (distance == "wasserstein") {               # closed-form 2-Wasserstein between Gaussians
                ra    <- .mat_sqrt(Sa)
                inner <- .mat_sqrt(ra %*% Sb %*% ra)
                return(sqrt(max(0, sum(dm^2) + sum(diag(Sa + Sb - 2 * inner)))))
        }
        Sbar <- (Sa + Sb) / 2                          # Hellinger via Bhattacharyya
        quad <- as.numeric(t(dm) %*% solve(Sbar, dm)) / 8                  # mean-separation term
        ld   <- as.numeric(determinant(Sbar, TRUE)$modulus) -             # covariance-mismatch term
                0.5 * (as.numeric(determinant(Sa, TRUE)$modulus) + as.numeric(determinant(Sb, TRUE)$modulus))
        sqrt(max(0, 1 - exp(-(quad + 0.5 * ld))))
}

#' Cross-region meta-clustering of local types
#'
#' Pools every region's GMM components into one standardised predictor space,
#' builds the pairwise component-distance matrix, and average-links them into a
#' meta-clustering. This is the engine behind both [plot_local_similarity()] and
#' [collapse_local()].
#'
#' @param lt A fitted local-type object (uses `lt$models`).
#' @param distance Component distance; see [plot_local_similarity()].
#' @param scale If `TRUE`, standardise predictors by their pooled SD so that no
#'   single variable dominates the distance.
#' @param sep Reserved separator for region/type labels (labels currently use a
#'   middle dot regardless of this value).
#' @param ridge Relative ridge forwarded to `.pair_dist()`.
#' @return A `localtype_similarity` object: the `dist` and full `matrix`, the
#'   `hclust` tree, per-label `info`, the common `vars`, and the `distance` used.
#' @noRd
.localtype_meta <- function(lt, distance = c("mean", "hellinger", "wasserstein"),
                            scale = TRUE, sep = ".", ridge = 1e-6) {
        distance <- match.arg(distance)
        
        comps  <- list()                               # one entry per (region, type) component
        pooled <- list()                               # per-region predictor tables, for SD scaling
        # Loop over Regions
        for (nm in names(lt$models)) {
                # Extract regional model
                fit <- lt$models[[nm]]
                # Means for each river type
                mu  <- fit$parameters$mean
                if (is.null(dim(mu))) {                # univariate fit: promote to a 1-row matrix
                        vn <- colnames(fit$data); if (is.null(vn)) vn <- "var"
                        mu <- matrix(mu, nrow = 1L, dimnames = list(vn[1], NULL))
                }
                # Variable names
                vn <- rownames(mu)
                # Add replacement V1...Vn if no rownames are given
                if (is.null(vn)) {
                        vn <- paste0("V", seq_len(nrow(mu)))
                        rownames(mu) <- vn
                }
                # run local helper which extracts covariance matrix
                sig <- .gmm_sigma(fit, vn)
                # register each component under a "R<region>\u00b7T<type>" label
                for (g in seq_len(ncol(mu))){
                        comps[[paste0(sub("Region","R", nm),"\u00b7", "T", g)]] <-
                                list(region = nm, type = g, mean = mu[, g], sigma = sig[[g]], vars = vn)
                }
                
                if (!is.null(fit$data)) pooled[[nm]] <- as.data.frame(fit$data)
        }
        if (length(comps) < 2L) stop("need at least two local types to compare.")
        
        vars <- Reduce(intersect, lapply(comps, `[[`, "vars"))   # common predictor space
        if (!length(vars)) stop("regional models share no common predictors.")
        
        # Per-variable scaling factor: pooled SD when scaling, else 1 (no scaling)
        sd_v <- stats::setNames(rep(1, length(vars)), vars)
        if (scale && length(pooled)) {
                pooled  <- do.call(rbind, lapply(pooled, function(d) d[, vars, drop = FALSE]))
                sd_v[]  <- vapply(vars, function(v) stats::sd(pooled[[v]], na.rm = TRUE), numeric(1))
                sd_v[!is.finite(sd_v) | sd_v <= 0] <- 1  # guard against zero / NA SDs
        }
        D <- diag(1 / sd_v, length(vars))              # diagonal standardising transform
        
        # Standardise each component's mean (and covariance) into the common space
        for (lab in names(comps)) {
                comps[[lab]]$m <- as.numeric(D %*% comps[[lab]]$mean[vars])
                S <- comps[[lab]]$sigma
                comps[[lab]]$S <- if (is.null(S)) NULL else D %*% S[vars, vars, drop = FALSE] %*% D
        }
        
        # Fill the symmetric pairwise distance matrix
        labs <- names(comps)
        P <- length(labs)
        M <- matrix(0, P, P, dimnames = list(labs, labs))
        for (a in seq_len(P - 1L)) for (b in (a + 1L):P){
                M[a, b] <- M[b, a] <- .pair_dist(comps[[a]],
                                                 comps[[b]],
                                                 distance,
                                                 ridge)
        }
        
        
        info <- data.frame(label  = labs,
                           region = vapply(comps, `[[`, "", "region"),
                           type   = vapply(comps, `[[`, 0L, "type"),
                           row.names = NULL, stringsAsFactors = FALSE)
        structure(list(dist = stats::as.dist(M), matrix = M,
                       hclust = stats::hclust(stats::as.dist(M), method = "average"),
                       info = info, vars = vars, distance = distance),
                  class = "localtype_similarity")
}

#' Plot local types in a reduced environmental space
#'
#' PCA-embeds all samples, then overlays original-type centroids and, optionally,
#' the collapsed-type centroids from [collapse_local()]. Useful for checking
#' whether a proposed collapse merges types that genuinely overlap in predictor
#' space.
#'
#' @param local Fitted local-type object; supplies the predictor set via
#'   `local$models`.
#' @param data All-sample predictor table, row-aligned to `memberships`.
#' @param memberships Wide local-type membership matrix; the argmax of each row
#'   gives that sample's original type.
#' @param collapse Optional `collapse_local()` result; colours samples by merged
#'   type and adds the collapsed centroids (diamonds).
#' @param dims 2 (static ggplot2) or 3 (interactive plotly).
#' @param n_max Maximum number of samples drawn per original type in the scatter.
#' @param scale If `TRUE`, scale predictors to unit variance before the PCA.
#' @param seed Optional RNG seed for the per-type sample thinning.
#' @param samples If `TRUE`, draw the sample points; if `FALSE` they are made
#'   fully transparent so only the centroids show (2D only).
#' @param ... Currently unused; reserved for future arguments.
#'
#' @return For `dims = 2` a `ggplot` object; for `dims = 3` a `plotly` object.
#' @seealso [collapse_local()], [plot_local_similarity()]
#' @export
plot_local_embedding <- function(local, data, memberships, collapse = NULL,
                                 dims = 2, n_max = 100, scale = TRUE, seed = NULL, samples = TRUE, ...) {
        .region_require("ggplot2")
        if (!dims %in% c(2L, 3L)) stop("dims must be 2 or 3.")
        if (missing(data) || missing(memberships))
                stop("supply `data` (all-sample predictors) and `memberships` (wide river-type matrix).")
        memberships <- as.matrix(memberships)
        data <- as.data.frame(data)
        if (nrow(memberships) != nrow(data)) stop("`data` and `memberships` need the same rows.")
        
        # Predictor set: the model's variables, intersected with what `data` provides
        mvars <- rownames(local$models[[1]]$parameters$mean)
        if (is.null(mvars)) mvars <- names(Filter(is.numeric, data))
        mvars <- intersect(mvars, colnames(data))
        X <- as.matrix(data[, mvars, drop = FALSE])
        X <- X[, apply(X, 2, stats::sd) > 0, drop = FALSE]   # drop constant columns (PCA can't use them)
        
        # Original (pre-collapse) type per sample = argmax membership column
        orig <- colnames(memberships)[max.col(memberships, ties.method = "first")]
        coll <- NULL
        if (!is.null(collapse)) {
                # Map each original type to its collapsed (meta) type
                map <- stats::setNames(collapse$mapping$meta_label, collapse$mapping$label)
                if (!all(orig %in% names(map)))
                        warning("some membership columns are not in the collapse mapping; check the column naming.")
                coll <- unname(map[orig])
        }
        grp <- if (is.null(coll)) orig else coll       # colour samples by collapsed type if available
        
        ## region for each original-type label: from the mapping, else from per-region block sizes
        col_region <- if (!is.null(collapse)) {
                stats::setNames(collapse$mapping$region, collapse$mapping$label)
        } else {
                # rebuild region-per-column from the number of components in each region's model
                gper <- vapply(local$models, function(f) { mu <- f$parameters$mean
                if (is.null(dim(mu))) length(mu) else ncol(mu) }, integer(1))
                rcv  <- rep(names(local$models), gper)
                if (length(rcv) == ncol(memberships)) stats::setNames(rcv, colnames(memberships)) else NULL
        }
        have_region <- !is.null(col_region)
        if (!have_region) warning("could not resolve region per type; crosses drawn in one colour.")
        
        # compute PCA on the (optionally scaled) predictor matrix
        pc <- stats::prcomp(X, center = TRUE, scale. = scale)
        sc <- as.data.frame(pc$x[, seq_len(dims), drop = FALSE]); names(sc) <- paste0("PC", seq_len(dims))
        
        # helper: per-label centroid in PC space (ignores NA labels)
        cen <- function(lab) {
                ok <- !is.na(lab); stats::aggregate(sc[ok, , drop = FALSE], list(label = lab[ok]), mean)
        }
        cen_o <- cen(orig)                             # original-type centroids
        if (have_region) cen_o$region <- unname(col_region[cen_o$label])
        cen_c <- if (!is.null(coll)) cen(coll) else NULL   # collapsed-type centroids (optional)
        
        # Thin to at most n_max samples per original type so dense types don't dominate
        if (!is.null(seed)) set.seed(seed)
        pick <- unlist(lapply(split(seq_len(nrow(sc)), orig),
                              function(ix) if (length(ix) > n_max) sample(ix, n_max) else ix),
                       use.names = FALSE)
        pts <- cbind(sc[pick, , drop = FALSE], group = grp[pick])
        
        if (dims == 2L) {
                # `samples = FALSE` keeps the points but makes them invisible (centroids only)
                if (samples == TRUE){
                        pnt.alpha = 0.5
                } else {
                        pnt.alpha = 0
                }
                p <- ggplot2::ggplot(pts, ggplot2::aes(.data$PC1, .data$PC2)) +
                        ggplot2::geom_point(ggplot2::aes(fill = .data$group), shape = 21, stroke = 0,
                                            size = 1.6, alpha = pnt.alpha)
                # crosses = original centroids, coloured by region when known
                p <- p + if (have_region)
                        ggplot2::geom_point(data = cen_o, ggplot2::aes(colour = .data$region),
                                            shape = 4, size = 2.6, stroke = 1)
                else
                        ggplot2::geom_point(data = cen_o, shape = 4, size = 2.6, stroke = 1, colour = "grey20")
                p <- p +
                        ggplot2::geom_text(data = cen_o, ggplot2::aes(label = .data$label),
                                           size = 2.6, colour = "grey20", vjust = -0.8, check_overlap = TRUE) +
                        ggplot2::labs(title = "Local types in reduced environmental space",
                                      subtitle = sprintf("PCA of %d predictors; crosses are original centroids%s",
                                                         ncol(X),
                                                         if (!is.null(cen_c)) ", diamonds are collapsed centroids" else ""),
                                      fill = if (is.null(coll)) "type" else "collapsed type", colour = "region") +
                        theme_eco()
                # diamonds = collapsed centroids, overlaid on top
                if (!is.null(cen_c))
                        p <- p +
                        ggplot2::geom_point(data = cen_c, shape = 18, size = 4, colour = "black") +
                        ggplot2::geom_text(data = cen_c, ggplot2::aes(label = .data$label),
                                           size = 3, fontface = "bold", vjust = 1.6, check_overlap = TRUE)
                return(p)
        }
        
        ## 3D (interactive)
        .region_require("plotly")
        fig <- plotly::add_markers(plotly::plot_ly(), data = pts, x = ~PC1, y = ~PC2, z = ~PC3,
                                   color = ~group, marker = list(size = 3, opacity = 0.5), name = "samples")
        fig <- if (have_region)
                plotly::add_markers(fig, data = cen_o, x = ~PC1, y = ~PC2, z = ~PC3, color = ~region,
                                    text = ~label, marker = list(size = 5, symbol = "x"), name = "centroids")
        else
                plotly::add_markers(fig, data = cen_o, x = ~PC1, y = ~PC2, z = ~PC3, text = ~label,
                                    marker = list(size = 5, symbol = "x", color = "grey30"), name = "centroids")
        if (!is.null(cen_c))
                fig <- plotly::add_markers(fig, data = cen_c, x = ~PC1, y = ~PC2, z = ~PC3, text = ~label,
                                           marker = list(size = 8, symbol = "diamond", color = "black"),
                                           name = "collapsed centroids")
        plotly::layout(fig, title = "Local types in reduced environmental space",
                       scene = list(xaxis = list(title = "PC1"), yaxis = list(title = "PC2"),
                                    zaxis = list(title = "PC3")))
}

