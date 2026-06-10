#' Collapse local types across regions
#'
#' Cuts the cross-region meta-clustering and merges local types that fall in the
#' same meta-cluster, summing their membership columns. Meta-types are named
#' after their medoid component.
#'
#' @param lt A fitted local-type object.
#' @param memberships The wide local-type membership matrix (catchments x types),
#'   columns named to match the meta labels (region.type).
#' @param k,h Cut the meta-clustering into `k` meta-types or at height `h`.
#' @param similarity Optional precomputed `localtype_similarity`; recomputed if NULL.
collapse_localtypes <- function(lt, memberships, k = NULL, h = NULL, similarity = NULL,
                                distance = c("mean", "hellinger", "wasserstein"),
                                sep = ".", ...) {
        distance <- match.arg(distance)
        if (is.null(k) && is.null(h)) stop("supply k or h.")
        if (is.null(similarity)) similarity <- .localtype_meta(lt, distance = distance, sep = sep)
        
        grp  <- 
                if (!is.null(k)){
                        stats::cutree(similarity$hclust, k = k)      
                } else {
                        stats::cutree(similarity$hclust, h = h)
                }
        info <- similarity$info
        info$meta <- grp[info$label]
        
        # name each meta-type sequentially: meta.1, meta.2, ...
        
        ids <- sort(unique(info$meta))
        med <- paste("meta", seq_along(ids), sep = ".")
        names(med) <- ids
        info$meta_label <- med[as.character(info$meta)]
        
        ### CRUD ##### 
        # Mmat <- similarity$matrix
        # med  <- tapply(info$label, info$meta, function(ls)
        #         if (length(ls) == 1L)
        #                 ls
        #         else
        #                 ls[which.min(rowSums(Mmat[ls, ls, drop = FALSE]))])
        # med <- paste("meta", 1:length(unique(info$meta)),sep=".") 
        # 
        # info$meta_label <- med[as.character(info$meta)]
        
        if (!all(info$label %in% colnames(memberships)))
                stop("membership columns must match the meta labels (region", sep, "type). ",
                     "Confirm the river_types() column naming and pass a matching matrix or relabel.")
        
        M   <- as.matrix(memberships[, info$label, drop = FALSE])
        key <- info$meta_label[match(colnames(M), info$label)]
        
        collapsed <- vapply(split(seq_len(ncol(M)), key),
                            function(j) rowSums(M[, j, drop = FALSE]), numeric(nrow(M)))
        
        collapsed <- collapsed[, unique(key), drop = FALSE]
        # rows are not guaranteed to sum to 1
        rs   <- rowSums(collapsed)
        # normalise only for hard label and entropy
        prob <- collapsed / ifelse(rs > 0, rs, 1)
        hard <- colnames(collapsed)[max.col(prob, ties.method = "first")]
        ent  <- -rowSums(ifelse(prob > 0, prob * log(prob), 0))
        
        list(memberships = collapsed, hard = hard, entropy = ent, row_sums = rs,
             mapping = info[order(info$meta), c("region", "type", "label", "meta", "meta_label")],
             n_before = nrow(info), n_after = ncol(collapsed), k = k, h = h)
}

#' Cross-region similarity of local types
#'
#' Compares the GMM components across all ecoregions in a common standardised
#' predictor space and draws a meta-clustered distance heatmap. Off-diagonal
#' blocks that mix regions are candidates for collapsing.
#'
#' @param lt A fitted local-type object (uses `lt$models`).
#' @param distance "mean" (standardised component means, matches the fingerprint
#'   plot), "hellinger" (mean and shape), or "wasserstein".
#' @param k Optional number of meta-types; draws the cut as block separators.
#' @return Invisibly, a `localtype_similarity` object (distance, hclust, info).
plot_localtype_similarity <- function(rt, distance = c("mean", "hellinger", "wasserstein"),
                                      k = NULL, h = NULL, groups = NULL, triangle = TRUE, ...) {
        .eco_require("ggplot2")
        distance <- match.arg(distance)
        meta <- .localtype_meta(rt, distance = distance, ...)
        ord  <- meta$hclust$order
        labs <- meta$info$label[ord]
        n    <- length(labs)
        M    <- meta$matrix[ord, ord]
        
        ## upper triangle, diagonal always dropped
        gi   <- expand.grid(i = seq_len(n), j = seq_len(n))
        gi$d <- M[cbind(gi$i, gi$j)]
        gi   <- if (triangle) gi[gi$j > gi$i, ] else gi[gi$i != gi$j, ]
        gi$a <- factor(labs[gi$i], levels = labs)
        gi$b <- factor(labs[gi$j], levels = labs)
        
        ## groups: a collapse_localtypes() result, a named vector, or an internal cut
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
                multi <- names(tab)[tab >= 2]
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

.gmm_sigma <- function(fit, vn) {
        v <- fit$parameters$variance
        G <- if (is.null(dim(fit$parameters$mean))) {
                length(fit$parameters$mean)
        }  else {
                ncol(fit$parameters$mean)
        }
        out <- vector("list", G)
        if (!is.null(v$sigma)) {                       # d x d x G (most multivariate models)
                for (g in seq_len(G)) { S <- v$sigma[, , g]; dimnames(S) <- list(vn, vn); out[[g]] <- S }
        } else if (!is.null(v$Sigma)) {                # shared covariance
                S <- v$Sigma; dimnames(S) <- list(vn, vn); for (g in seq_len(G)) out[[g]] <- S
        } else if (!is.null(v$sigmasq)) {              # univariate
                s2 <- v$sigmasq; if (length(s2) == 1L) s2 <- rep(s2, G)
                for (g in seq_len(G)) out[[g]] <- matrix(s2[g], 1, 1, dimnames = list(vn, vn))
        }
        out                                            # any NULL element triggers the mean-only fallback
}

.mat_sqrt <- function(S) {
        e <- eigen(S, symmetric = TRUE); e$vectors %*% diag(sqrt(pmax(e$values, 0)), nrow(S)) %*% t(e$vectors)
}

.pair_dist <- function(a, b, distance, ridge) {
        dm <- a$m - b$m
        if (distance == "mean" || is.null(a$S) || is.null(b$S)) return(sqrt(sum(dm^2)))
        Sa <- a$S; Sb <- b$S
        rg <- ridge * mean(c(diag(Sa), diag(Sb)))
        Sa <- Sa + diag(rg, nrow(Sa)); Sb <- Sb + diag(rg, nrow(Sb))
        if (distance == "wasserstein") {
                ra    <- .mat_sqrt(Sa)
                inner <- .mat_sqrt(ra %*% Sb %*% ra)
                return(sqrt(max(0, sum(dm^2) + sum(diag(Sa + Sb - 2 * inner)))))
        }
        Sbar <- (Sa + Sb) / 2                          # Hellinger via Bhattacharyya
        quad <- as.numeric(t(dm) %*% solve(Sbar, dm)) / 8
        ld   <- as.numeric(determinant(Sbar, TRUE)$modulus) -
                0.5 * (as.numeric(determinant(Sa, TRUE)$modulus) + as.numeric(determinant(Sb, TRUE)$modulus))
        sqrt(max(0, 1 - exp(-(quad + 0.5 * ld))))
}

.localtype_meta <- function(lt, distance = c("mean", "hellinger", "wasserstein"),
                            scale = TRUE, sep = ".", ridge = 1e-6) {
        distance <- match.arg(distance)
        
        comps  <- list()
        pooled <- list()
        # Loop over Regions 
        for (nm in names(lt$models)) {
                # Extract regional model
                fit <- lt$models[[nm]]
                # Means for each river type 
                mu  <- fit$parameters$mean
                if (is.null(dim(mu))) {
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
                for (g in seq_len(ncol(mu))){
                        comps[[paste0(nm,"localType", g)]] <-
                                list(region = nm, type = g, mean = mu[, g], sigma = sig[[g]], vars = vn)
                }
                        
                if (!is.null(fit$data)) pooled[[nm]] <- as.data.frame(fit$data)
        }
        if (length(comps) < 2L) stop("need at least two local types to compare.")
        
        vars <- Reduce(intersect, lapply(comps, `[[`, "vars"))   # common predictor space
        if (!length(vars)) stop("regional models share no common predictors.")
        
        sd_v <- stats::setNames(rep(1, length(vars)), vars)
        if (scale && length(pooled)) {
                pooled  <- do.call(rbind, lapply(pooled, function(d) d[, vars, drop = FALSE]))
                sd_v[]  <- vapply(vars, function(v) stats::sd(pooled[[v]], na.rm = TRUE), numeric(1))
                sd_v[!is.finite(sd_v) | sd_v <= 0] <- 1
        }
        D <- diag(1 / sd_v, length(vars))
        
        for (lab in names(comps)) {
                comps[[lab]]$m <- as.numeric(D %*% comps[[lab]]$mean[vars])
                S <- comps[[lab]]$sigma
                comps[[lab]]$S <- if (is.null(S)) NULL else D %*% S[vars, vars, drop = FALSE] %*% D
        }
        
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
#' Local types in a reduced environmental space
#'
#' PCA-embeds all samples, then overlays original-type centroids and, optionally,
#' collapsed-type centroids from `collapse_localtypes()`. Useful for checking
#' whether a proposed collapse merges types that actually overlap.
#'
#' @param rt Fitted river-type object (used for the predictor set).
#' @param data All-sample predictor table, row-aligned to `memberships`.
#' @param memberships Wide river-type membership matrix; argmax gives each sample's type.
#' @param collapse Optional `collapse_localtypes()` result; colours samples by merged
#'   type and adds collapsed centroids.
#' @param dims 2 or 3.
#' @param n_max Maximum samples drawn per original type for the scatter.
plot_localtype_embedding <- function(rt, data, memberships, collapse = NULL,
                                     dims = 2, n_max = 100, scale = TRUE, seed = NULL, samples = TRUE, ...) {
        .eco_require("ggplot2")
        if (!dims %in% c(2L, 3L)) stop("dims must be 2 or 3.")
        if (missing(data) || missing(memberships))
                stop("supply `data` (all-sample predictors) and `memberships` (wide river-type matrix).")
        memberships <- as.matrix(memberships)
        data <- as.data.frame(data)
        if (nrow(memberships) != nrow(data)) stop("`data` and `memberships` need the same rows.")
        
        mvars <- rownames(rt$models[[1]]$parameters$mean)
        if (is.null(mvars)) mvars <- names(Filter(is.numeric, data))
        mvars <- intersect(mvars, colnames(data))
        X <- as.matrix(data[, mvars, drop = FALSE])
        X <- X[, apply(X, 2, stats::sd) > 0, drop = FALSE]
        
        orig <- colnames(memberships)[max.col(memberships, ties.method = "first")]
        coll <- NULL
        if (!is.null(collapse)) {
                map <- stats::setNames(collapse$mapping$meta_label, collapse$mapping$label)
                if (!all(orig %in% names(map)))
                        warning("some membership columns are not in the collapse mapping; check the column naming.")
                coll <- unname(map[orig])
        }
        grp <- if (is.null(coll)) orig else coll
        
        ## region for each original-type label: from the mapping, else from per-region block sizes
        col_region <- if (!is.null(collapse)) {
                stats::setNames(collapse$mapping$region, collapse$mapping$label)
        } else {
                gper <- vapply(rt$models, function(f) { mu <- f$parameters$mean
                if (is.null(dim(mu))) length(mu) else ncol(mu) }, integer(1))
                rcv  <- rep(names(rt$models), gper)
                if (length(rcv) == ncol(memberships)) stats::setNames(rcv, colnames(memberships)) else NULL
        }
        have_region <- !is.null(col_region)
        if (!have_region) warning("could not resolve region per type; crosses drawn in one colour.")
        
        # compute PCA
        pc <- stats::prcomp(X, center = TRUE, scale. = scale)
        sc <- as.data.frame(pc$x[, seq_len(dims), drop = FALSE]); names(sc) <- paste0("PC", seq_len(dims))
        
        cen <- function(lab) {
                ok <- !is.na(lab); stats::aggregate(sc[ok, , drop = FALSE], list(label = lab[ok]), mean)
        }
        cen_o <- cen(orig)
        if (have_region) cen_o$region <- unname(col_region[cen_o$label])
        cen_c <- if (!is.null(coll)) cen(coll) else NULL
        
        if (!is.null(seed)) set.seed(seed)
        pick <- unlist(lapply(split(seq_len(nrow(sc)), orig),
                              function(ix) if (length(ix) > n_max) sample(ix, n_max) else ix),
                       use.names = FALSE)
        pts <- cbind(sc[pick, , drop = FALSE], group = grp[pick])
        
        if (dims == 2L) {
                if (samples == TRUE){
                        pnt.alpha = 0.5
                } else {
                        pnt.alpha = 0
                }
                p <- ggplot2::ggplot(pts, ggplot2::aes(.data$PC1, .data$PC2)) +
                        ggplot2::geom_point(ggplot2::aes(fill = .data$group), shape = 21, stroke = 0,
                                            size = 1.6, alpha = pnt.alpha)
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
                if (!is.null(cen_c))
                        p <- p +
                        ggplot2::geom_point(data = cen_c, shape = 18, size = 4, colour = "black") +
                        ggplot2::geom_text(data = cen_c, ggplot2::aes(label = .data$label),
                                           size = 3, fontface = "bold", vjust = 1.6, check_overlap = TRUE)
                return(p)
        }
        
        ## 3D
        .eco_require("plotly")
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
