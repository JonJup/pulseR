#' Store the environment matrix on a graph with its own edge weights
#'
#' The default Level 1 settings, \code{cut_rule = "ssd"} and
#' \code{refine_boundary = TRUE}, cut trees and move boundary units on the
#' within-region sum of squares of the environment matrix stored on the graph
#' (graph attribute \code{env}). \code{\link{add_edge_weight}} stores this
#' matrix. \code{set_env()} stores it on a graph whose edge weights come from
#' elsewhere, matching rows to vertices by ID.
#'
#' @param graph An \code{igraph} object, one vertex per unit, usually with an
#'   edge attribute \code{weight}.
#' @param env Data frame, matrix or \pkg{sf} object, one row per unit. Extra
#'   rows are dropped when rows are matched by ID.
#' @param env_id Name of the ID column in \code{env}. \code{NULL} (default)
#'   uses a column named \code{ID} or \code{id} if there is one. Rows are
#'   matched to the vertex attribute \code{polygon_id} (set by
#'   \code{\link{build_network}}). Without an ID column or without
#'   \code{polygon_id}, row \eqn{i} is taken to be vertex \eqn{i}, with a
#'   warning.
#' @param variables Character vector of columns to use. \code{NULL} (default)
#'   uses every numeric column except the ID column and columns named
#'   \code{true_*}.
#' @param scale Logical. Standardise each column to mean 0 and SD 1, dropping
#'   constant columns. Default \code{FALSE}: the matrix is stored as given, so
#'   transformed predictors keep their fitted weighting. Set \code{TRUE} for
#'   raw predictors on different scales.
#' @param verbose Logical. Report the size of the stored matrix and how well
#'   Euclidean distances on it reproduce the edge weights. Default \code{TRUE}.
#'
#' @return \code{graph} with graph attribute \code{env}: a numeric matrix, one
#'   row per vertex in vertex order, with no missing values. An \code{env}
#'   already on \code{graph} is replaced.
#'
#' @section Bring your own weights:
#' \code{\link{get_regions}} and \code{\link{tune_regions}} take any graph
#' with a numeric \code{weight} edge attribute. The edge weights decide which
#' spanning trees are sampled. Under the defaults, the environment matrix
#' decides where each tree is cut and which boundary units move. Both
#' defaults stop if \code{graph} has no environment matrix. There are two
#' ways forward:
#' \itemize{
#'   \item Store a matrix with \code{set_env()}. For weights from a
#'     generalised dissimilarity model (GDM), use the I-spline-transformed
#'     predictors with \code{scale = FALSE}. GDM distance is the Manhattan
#'     distance between transformed predictors passed through a link
#'     function, whereas the cuts and the refinement use squared Euclidean
#'     distance on the same predictors. The two orderings are close but not
#'     identical.
#'   \item Use only the edge weights: \code{cut_rule = "propensity"} and
#'     \code{refine_boundary = FALSE}.
#' }
#' The \code{env} argument of \code{\link{get_regions}} and
#' \code{\link{tune_regions}} does not replace the stored matrix. It sets only
#' the space of the validity indices \code{CH_FIT} and \code{R2_FIT}, which
#' default to the stored matrix.
#'
#' \code{\link{compute_null_baseline}} permutes environment rows and recomputes
#' the weights as Euclidean distances, so it needs weights that this matrix
#' reproduces (correlation above 0.999). With \code{verbose = TRUE},
#' \code{set_env()} reports this correlation. Weights from another model, such
#' as GDM, usually do not pass this check.
#'
#' @seealso \code{\link{add_edge_weight}}, \code{\link{get_regions}}.
#'
#' @examples
#' g   <- build_network(pulse_grid, verbose = FALSE)$graph
#' env <- sf::st_drop_geometry(pulse_grid)
#' ## Edge weights from another source; here a stand-in on one variable.
#' ep  <- igraph::ends(g, igraph::E(g), names = FALSE)
#' igraph::E(g)$weight <- abs(env$temp[ep[, 1]] - env$temp[ep[, 2]])
#' g   <- set_env(g, env, variables = c("temp", "prec"), scale = TRUE)
#' dim(igraph::graph_attr(g, "env"))
#'
#' @export
set_env <- function(graph, env, env_id = NULL, variables = NULL,
                    scale = FALSE, verbose = TRUE) {
        if (!inherits(graph, "igraph"))
                stop("`graph` must be an igraph object.", call. = FALSE)
        .assert_flag(scale, "scale")
        .assert_flag(verbose, "verbose")
        if (inherits(env, "sf")) env <- sf::st_drop_geometry(env)
        if (!is.data.frame(env) && !is.matrix(env))
                stop("`env` must be a data frame, matrix or sf object.", call. = FALSE)
        env <- as.data.frame(env, stringsAsFactors = FALSE)
        n   <- igraph::vcount(graph)

        # ---- IDs ------------------------------------------------------------
        if (!is.null(env_id)) {
                if (!is.character(env_id) || length(env_id) != 1L || !env_id %in% names(env))
                        stop("`env_id` must name one column of `env`.", call. = FALSE)
        } else {
                env_id <- .id_col(env)
                if (is.na(env_id)) env_id <- NULL
        }
        env_ids <- if (is.null(env_id)) NULL else env[[env_id]]
        v_ids   <- igraph::vertex_attr(graph, "polygon_id")
        if (is.null(env_ids) || is.null(v_ids))
                warning(sprintf(paste0(
                        "No %s to match on, so row i of `env` is taken to be vertex i. ",
                        "Pass `env_id` and a graph from build_network() to match by ID."),
                        if (is.null(env_ids)) "ID column in `env`" else "`polygon_id` on `graph`"),
                        call. = FALSE)

        # ---- Columns --------------------------------------------------------
        if (!is.null(variables)) {
                miss <- setdiff(variables, names(env))
                if (length(miss))
                        stop(sprintf("Column(s) not found in `env`: %s.",
                                     paste(miss, collapse = ", ")), call. = FALSE)
                bad <- variables[!vapply(env[variables], is.numeric, logical(1))]
                if (length(bad))
                        stop(sprintf("Column(s) in `variables` are not numeric: %s.",
                                     paste(bad, collapse = ", ")), call. = FALSE)
                cols <- variables
        } else {
                cols <- .drop_truth_cols(setdiff(names(env), env_id))
                num  <- vapply(env[cols], is.numeric, logical(1))
                if (verbose && any(!num))
                        message(sprintf("Ignoring non-numeric column(s) %s.",
                                        paste(cols[!num], collapse = ", ")))
                cols <- cols[num]
        }
        if (!length(cols))
                stop("`env` has no numeric columns to store.", call. = FALSE)
        X <- as.matrix(env[, cols, drop = FALSE])
        storage.mode(X) <- "double"

        # ---- Rows -----------------------------------------------------------
        X <- .align_env_rows(X, env_ids, v_ids, n, verbose = verbose)
        if (anyNA(X)) {
                n_na <- colSums(is.na(X))
                stop(sprintf(paste0(
                        "`env` has %d missing value(s) (column(s) %s). The \"ssd\" cut ",
                        "rule and the boundary refinement need complete rows: impute ",
                        "or drop the columns first."), sum(n_na),
                        paste(colnames(X)[n_na > 0], collapse = ", ")), call. = FALSE)
        }
        .warn_id_like(X, "`env_id` or `variables`")

        if (scale) {
                sds   <- apply(X, 2L, stats::sd)
                const <- !(sds > 0)
                if (any(const)) {
                        warning(sprintf("Dropping constant column(s) %s, which cannot be scaled.",
                                        paste(colnames(X)[const], collapse = ", ")), call. = FALSE)
                        X <- X[, !const, drop = FALSE]
                        if (!ncol(X)) stop("Every column of `env` is constant.", call. = FALSE)
                }
                X <- base::scale(X)          # keeps scaled:center/scale, as add_edge_weight()
        }
        dimnames(X) <- list(igraph::vertex_attr(graph, "name"), colnames(X))

        if (verbose) {
                if (!is.null(igraph::graph_attr(graph, "env")))
                        message("Replacing the environment matrix already stored on `graph`.")
                message(sprintf("Stored a %d x %d environment matrix (%s).",
                                nrow(X), ncol(X), paste(colnames(X), collapse = ", ")))
                w <- igraph::edge_attr(graph, "weight")
                if (!is.null(w) && length(w) > 1L) {
                        r <- suppressWarnings(stats::cor(.edge_weights_from_env(graph, X), w))
                        message(sprintf(paste0(
                                "Correlation between Euclidean distances on it and E(graph)$weight: ",
                                "%.3f.%s"), r,
                                if (is.finite(r) && r > 0.999) "" else
                                        " compute_null_baseline() needs > 0.999 and will stop."))
                }
        }
        igraph::set_graph_attr(graph, "env", X)
}

#' Check the environment matrix the Level 1 defaults need
#'
#' \code{cut_rule = "ssd"} and \code{refine_boundary = TRUE} read
#' \code{graph_attr(graph, "env")}. Called at the top of every entry point so a
#' missing or broken matrix stops before any sampling, and again where the
#' matrix is read.
#' @param graph An \code{igraph} object.
#' @param cut_rule Character scalar.
#' @param refine_boundary Logical scalar.
#' @return The matrix invisibly, or \code{NULL} if neither setting needs it.
#' @keywords internal
#' @noRd
.require_env <- function(graph, cut_rule = "propensity", refine_boundary = FALSE) {
        need <- c(if (identical(cut_rule, "ssd")) "cut_rule = \"ssd\"",
                  if (isTRUE(refine_boundary)) "refine_boundary = TRUE")
        if (!length(need)) return(invisible(NULL))
        what <- paste(need, collapse = " and ")
        X    <- igraph::graph_attr(graph, "env")
        if (is.null(X))
                stop(sprintf(paste0(
                        "%s need%s the environment matrix graph_attr(graph, \"env\"), ",
                        "and `graph` has none. Either build `graph` with add_edge_weight(), ",
                        "store a matrix with set_env() (for edge weights from elsewhere, ",
                        "e.g. GDM), or use only the edge weights with ",
                        "cut_rule = \"propensity\", refine_boundary = FALSE. The `env` ",
                        "argument sets only the space of CH_FIT and R2_FIT. See ?set_env."),
                        what, if (length(need) == 1L) "s" else ""), call. = FALSE)
        if (!is.matrix(X) || !is.numeric(X))
                stop(sprintf(paste0(
                        "%s needs graph_attr(graph, \"env\") to be a numeric matrix, ",
                        "not %s. Store it with set_env()."), what, class(X)[1L]), call. = FALSE)
        if (nrow(X) != igraph::vcount(graph))
                stop(sprintf(paste0(
                        "%s needs graph_attr(graph, \"env\") with one row per vertex: it ",
                        "has %d rows, `graph` has %d vertices. Store it with set_env(), ",
                        "which matches rows by ID."), what, nrow(X), igraph::vcount(graph)),
                     call. = FALSE)
        if (anyNA(X))
                stop(sprintf(paste0(
                        "%s needs complete graph_attr(graph, \"env\"); it has %d missing ",
                        "value(s). Impute or drop them first."), what, sum(is.na(X))),
                     call. = FALSE)
        invisible(X)
}
