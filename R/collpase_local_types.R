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