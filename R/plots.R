# =============================================================================
# pulseR — visualisation for the two-level typology
#
#   Level 1  pulse_regions      N units x k regions
#   Level 2  pulse_local_types  N units x sum(G_r) local types, nested in regions
#
# Design contract
# ---------------
# 1. ONE palette per fit. pulse_palette() is exported and both levels consume
#    it, so a region map and a local-type map from the same fit are readable
#    side by side. Local-type colours are shades of their parent region's hue.
# 2. Colours are ORDERED BY SIMILARITY, not by label index. Region 3 and
#    region 17 get neighbouring hues if the ensemble keeps confusing them.
#    Arbitrary hues across 20+ European regions are unreadable.
# 3. Uncertainty is baked into the fill, not applied as ggplot `alpha`. Alpha
#    composites against whatever is underneath and breaks the moment a second
#    layer or a basemap appears. Every map here precomputes a hex fill and uses
#    scale_*_identity(), which also keeps large-N plots cheap.
# 4. Geometry lives outside the fit objects (they only ever saw a graph), so it
#    is attached once with attach_geometry() and validated for row alignment.
#
# Depends: ggplot2, sf, grDevices, stats. Optional: terra (rasterisation).
# =============================================================================


# -----------------------------------------------------------------------------
# Geometry
# -----------------------------------------------------------------------------

#' Attach polygon geometry to a fitted typology
#'
#' The fitting functions see only a graph, so the polygons have to be supplied
#' before anything can be mapped. This validates row alignment once and stores
#' the geometry as an attribute, so every later plotting call is a one-liner.
#'
#' @param x A \code{"pulse_regions"} or \code{"pulse_local_types"} object.
#' @param geom An \code{sf} object with one row per unit, in graph vertex order.
#' @param simplify Numeric or \code{NULL}. If numeric, \code{sf::st_simplify}
#'   tolerance (map units) applied once and cached. Continental coverages carry
#'   far more vertices than any figure resolves; simplifying here rather than at
#'   draw time is the difference between seconds and minutes per panel.
#' @return \code{x}, with geometry attached.
#' @export
attach_geometry <- function(x, geom, simplify = NULL) {
        if (!inherits(geom, "sf")) stop("`geom` must be an sf object.", call. = FALSE)
        n <- .n_units(x)
        if (nrow(geom) != n)
                stop(sprintf("`geom` has %d rows but the fit has %d units. Rows must be in graph vertex order.",
                             nrow(geom), n), call. = FALSE)
        geom <- geom[, attr(geom, "sf_column"), drop = FALSE]
        if (!is.null(simplify))
                geom <- sf::st_simplify(geom, dTolerance = simplify, preserveTopology = TRUE)
        attr(x, "pulse_geometry") <- geom
        x
}

#' @keywords internal
#' @noRd
.get_geom <- function(x, geom = NULL) {
        g <- if (!is.null(geom)) geom else attr(x, "pulse_geometry")
        if (is.null(g))
                stop("No geometry. Call attach_geometry(x, polygons) first, or pass `geom =`.",
                     call. = FALSE)
        if (nrow(g) != .n_units(x))
                stop("Geometry rows do not match the number of units.", call. = FALSE)
        g
}


# -----------------------------------------------------------------------------
# Palettes
# -----------------------------------------------------------------------------

#' Similarity-ordered colours for the regions of a fit
#'
#' Regions are embedded in one or two dimensions from a between-region
#' dissimilarity and the embedding is mapped into HCL space, so perceptually
#' close colours mean ecologically close regions. Compute this once and pass it
#' to every plotting call for the fit; the default palette in each function is
#' this same call, so figures agree unless you deliberately override.
#'
#' @param x A \code{"pulse_regions"} object.
#' @param method \code{"auto"} (default), \code{"ramp1d"}, \code{"hue2d"} or
#'   \code{"qualitative"}. \code{"auto"} uses \code{"qualitative"} for
#'   \eqn{k \le 8} — at that size distinguishability beats similarity ordering —
#'   and \code{"hue2d"} above it.
#' @param dissimilarity Optional k x k matrix. Defaults to
#'   \eqn{1 - \mathrm{region\_coassoc}} from the fit (how often the ensemble
#'   put two regions' members together), falling back to Euclidean distance
#'   between environmental centroids if \code{env} is supplied, and to no
#'   ordering if neither is available.
#' @param env Optional predictor matrix, one row per unit, used for centroids
#'   when the fit carries no co-association matrix.
#' @param chroma,luminance Numeric. HCL chroma at maximum distinctness and
#'   luminance of the palette. Defaults 65 and 65 print and project acceptably.
#' @param hue_range Numeric length 2, degrees, for \code{"ramp1d"}.
#' @return Character vector of \code{k} hex colours, named by region index.
#' @export
pulse_palette <- function(x,
                          method        = c("auto", "ramp1d", "hue2d", "qualitative"),
                          dissimilarity = NULL,
                          env           = NULL,
                          chroma        = 65,
                          luminance     = 65,
                          hue_range     = c(10, 340)) {
        method <- match.arg(method)
        k <- if (inherits(x, "pulse_regions")) x$k else length(x$n_local_types)
        
        D <- .region_dissimilarity(x, dissimilarity, env, k)
        if (method == "auto") method <- if (k <= 8 || is.null(D)) "qualitative" else "hue2d"
        if (is.null(D) && method != "qualitative") {
                warning("No region dissimilarity available; falling back to a qualitative palette.",
                        call. = FALSE)
                method <- "qualitative"
        }
        
        cols <- switch(
                method,
                qualitative = grDevices::hcl(h = seq(15, 375, length.out = k + 1)[seq_len(k)],
                                             c = chroma + 30, l = luminance),
                ramp1d = {
                        s <- stats::cmdscale(stats::as.dist(D), k = 1)[, 1]
                        h <- hue_range[1] + diff(hue_range) * .rescale01(s)
                        grDevices::hcl(h = h, c = chroma + 20, l = luminance)
                },
                hue2d = {
                        cfg <- stats::cmdscale(stats::as.dist(D), k = 2)
                        ## Angle -> hue, radius -> chroma. Regions at the centre of the
                        ## ordination are intermediate and correctly come out muted;
                        ## extreme regions get saturated, unmistakable colours.
                        h <- (atan2(cfg[, 2], cfg[, 1]) * 180 / pi) %% 360
                        r <- sqrt(rowSums(cfg^2))
                        cc <- 20 + (chroma + 20 - 20) * .rescale01(r)
                        grDevices::hcl(h = h, c = cc, l = luminance)
                })
        stats::setNames(cols, seq_len(k))
}

#' @keywords internal
#' @noRd
.region_dissimilarity <- function(x, dissimilarity, env, k) {
        if (!is.null(dissimilarity)) return(as.matrix(dissimilarity))
        B <- x$fit$region_coassoc
        if (!is.null(B) && all(dim(B) == c(k, k))) {
                D <- 1 - B
                D[!is.finite(D)] <- 1
                diag(D) <- 0
                return((D + t(D)) / 2)
        }
        if (!is.null(env) && inherits(x, "pulse_regions")) {
                X <- scale(as.matrix(env))
                G <- rowsum(X, x$hard_labels) / as.vector(tabulate(x$hard_labels, k))
                return(as.matrix(stats::dist(G)))
        }
        NULL
}

#' Nested colours for local types
#'
#' Every local type keeps its parent region's hue and varies in luminance, so a
#' local-type map still reads as the region map at a glance. Without this the
#' two levels look like unrelated results.
#'
#' @param x A \code{"pulse_local_types"} object.
#' @param region_colours Character vector of region colours, normally from
#'   \code{\link{pulse_palette}} on the parent \code{"pulse_regions"} object.
#' @param l_range Numeric length 2. Luminance span within a region.
#' @return Character vector of hex colours, one per column of \code{x$memberships}.
#' @export
local_palette <- function(x, region_colours, l_range = c(80, 38)) {
        stopifnot(inherits(x, "pulse_local_types"))
        key <- x$key
        hcl <- t(grDevices::col2rgb(region_colours))
        hcl <- grDevices::convertColor(hcl / 255, "sRGB", "Luv")
        h   <- (atan2(hcl[, 3], hcl[, 2]) * 180 / pi) %% 360
        cc  <- sqrt(hcl[, 2]^2 + hcl[, 3]^2)
        
        out <- character(nrow(key))
        for (r in unique(key$region)) {
                i  <- which(key$region == r)
                G  <- length(i)
                ll <- if (G == 1L) mean(l_range) else seq(l_range[1], l_range[2], length.out = G)
                out[i] <- grDevices::hcl(h = h[r], c = cc[r], l = ll)
        }
        stats::setNames(out, key$name)
}

#' Fade a colour toward neutral in proportion to confidence
#'
#' Hue carries the class, chroma carries the confidence. A unit the ensemble
#' never settled on comes out near-grey, which is the honest rendering for a
#' probabilistic typology and makes transition zones visible without a second
#' panel.
#'
#' @param cols Character vector of hex colours.
#' @param w Numeric confidence in \code{[0, 1]}, same length.
#' @param floor Numeric. Chroma retained at \code{w = 0}, as a fraction.
#' @param to_l Numeric. Luminance the fully-uncertain end moves toward.
#' @return Character vector of hex colours.
#' @export
fade_by_confidence <- function(cols, w, floor = 0.12, to_l = 88) {
        w <- pmin(pmax(as.numeric(w), 0), 1)
        w[!is.finite(w)] <- 0
        luv <- grDevices::convertColor(t(grDevices::col2rgb(cols)) / 255, "sRGB", "Luv")
        h   <- (atan2(luv[, 3], luv[, 2]) * 180 / pi) %% 360
        cc  <- sqrt(luv[, 2]^2 + luv[, 3]^2)
        f   <- floor + (1 - floor) * w
        grDevices::hcl(h = h, c = cc * f, l = luv[, 1] + (to_l - luv[, 1]) * (1 - w))
}

#' @keywords internal
#' @noRd
.rescale01 <- function(v) {
        r <- range(v, finite = TRUE)
        if (diff(r) == 0) return(rep(0.5, length(v)))
        (v - r[1]) / diff(r)
}


# -----------------------------------------------------------------------------
# Shared map builder
# -----------------------------------------------------------------------------

#' @keywords internal
#' @noRd
.map <- function(geom, fill, title = NULL, subtitle = NULL, border = NA,
                 legend = NULL, raster = c("auto", "never", "always"),
                 raster_px = 2000) {
        raster <- match.arg(raster)
        use_raster <- raster == "always" ||
                (raster == "auto" && nrow(geom) > 5e4 && requireNamespace("terra", quietly = TRUE))
        
        if (use_raster) {
                if (!requireNamespace("terra", quietly = TRUE))
                        stop("Rasterisation needs the 'terra' package.", call. = FALSE)
                v  <- terra::vect(geom)
                bb <- as.vector(sf::st_bbox(geom))
                ar <- (bb[4] - bb[2]) / (bb[3] - bb[1])
                rr <- terra::rast(v, ncols = raster_px, nrows = round(raster_px * ar))
                idx <- terra::rasterize(v, rr, field = seq_len(nrow(geom)))
                df  <- as.data.frame(idx, xy = TRUE)
                names(df)[3] <- "i"
                df$fill <- fill[df$i]
                p <- ggplot2::ggplot(df, ggplot2::aes(.data$x, .data$y, fill = .data$fill)) +
                        ggplot2::geom_raster() + ggplot2::coord_equal()
        } else {
                g <- geom
                g$fill <- fill
                p <- ggplot2::ggplot(g) +
                        ggplot2::geom_sf(ggplot2::aes(fill = .data$fill),
                                         colour = border, linewidth = 0.05)
        }
        p +
                ggplot2::scale_fill_identity(guide = legend) +
                ggplot2::labs(title = title, subtitle = subtitle, x = NULL, y = NULL) +
                ggplot2::theme_void(base_size = 11) +
                ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"),
                               legend.position = "right")
}

#' @keywords internal
#' @noRd
.legend_df <- function(labels, cols) {
        data.frame(label = factor(labels, levels = labels), colour = unname(cols),
                   stringsAsFactors = FALSE)
}


# -----------------------------------------------------------------------------
# Level 1: regions
# -----------------------------------------------------------------------------

#' Plot a fitted regionalisation
#'
#' @param x A \code{"pulse_regions"} object.
#' @param what One of:
#'   \describe{
#'     \item{\code{"map"}}{Hard regions, faded by confidence. The headline figure.}
#'     \item{\code{"crisp"}}{Hard regions at full saturation, no fading.}
#'     \item{\code{"confidence"}}{The \code{confidence} quantity as a continuous surface.}
#'     \item{\code{"entropy"}}{Normalised membership entropy.}
#'     \item{\code{"transition"}}{Units whose best and second-best regions are
#'       within \code{margin} — the fuzzy boundaries themselves.}
#'     \item{\code{"membership"}}{One faceted panel per region (see \code{which}).}
#'     \item{\code{"coassociation"}}{k x k between-region co-association heatmap.}
#'     \item{\code{"ordination"}}{Regions in co-association space.}
#'     \item{\code{"tuning"}}{Score against each searched axis.}
#'   }
#' @param geom An \code{sf} object, if not already attached.
#' @param palette Region colours; defaults to \code{\link{pulse_palette}(x)}.
#' @param confidence Which quantity drives the fade and the confidence surface.
#'   \code{"certainty"} (default): own-region composition membership
#'   \eqn{\Pi_{i,\mathrm{hard}_i}}. \code{"margin"}: best minus second-best
#'   \eqn{\Pi}, which isolates boundaries. \code{"cohesion"}: own-region
#'   cohesion. Cohesion shrinks with region size, so large homogeneous regions
#'   look uncertain; use it only to compare units within one region.
#'   \code{"membership"} is accepted as an alias for \code{"certainty"}.
#' @param which Integer vector of regions, for \code{what = "membership"}.
#' @param margin Numeric, for \code{what = "transition"}. Default 0.1.
#' @param raster,raster_px Rasterisation control, passed to the map builder.
#' @param ... Unused.
#' @return A \code{ggplot} object.
#' @export
plot_regions <- function(x,
                         what       = c("map", "crisp", "confidence", "entropy",
                                        "transition", "membership", "coassociation",
                                        "ordination", "tuning"),
                         geom       = NULL,
                         palette    = NULL,
                         confidence = c("certainty", "margin", "cohesion"),
                         which      = NULL,
                         margin     = 0.1,
                         raster     = "auto",
                         raster_px  = 2000,
                         ...) {
        .need("ggplot2")
        stopifnot(inherits(x, "pulse_regions"))
        what       <- match.arg(what)
        if (identical(confidence, "membership")) confidence <- "certainty"
        confidence <- match.arg(confidence)
        if (is.null(palette)) palette <- pulse_palette(x)
        
        needs_geom <- what %in% c("map", "crisp", "confidence", "entropy",
                                  "transition", "membership")
        g <- if (needs_geom) .get_geom(x, geom) else NULL
        conf <- .confidence_vec(x, confidence)
        
        switch(
                what,
                
                map = .map(g,
                           fade_by_confidence(palette[x$hard_labels], conf),
                           title    = sprintf("River ecotypes (k = %d)", x$k),
                           subtitle = "hue = region, saturation = confidence",
                           raster = raster, raster_px = raster_px) +
                        .region_key(palette),
                
                crisp = .map(g, unname(palette[x$hard_labels]),
                             title = sprintf("River ecotypes (k = %d)", x$k),
                             raster = raster, raster_px = raster_px) +
                        .region_key(palette),
                
                confidence = .continuous_map(g, conf,
                                             sprintf("Assignment confidence (%s)", confidence),
                                             limits = c(0, 1), raster, raster_px),
                
                entropy = .continuous_map(g, x$entropy, "Membership entropy",
                                          limits = c(0, 1), raster, raster_px,
                                          direction = -1),
                
                transition = {
                        M    <- .conf_matrix(x, confidence)
                        srt  <- t(apply(M, 1, function(r) sort(r, decreasing = TRUE)[1:2]))
                        gap  <- srt[, 1] - srt[, 2]
                        fill <- ifelse(gap <= margin, "#B2182B", "#EEEEEE")
                        .map(g, fill,
                             title    = "Transition zones",
                             subtitle = sprintf("best minus second-best <= %.2f (%.1f%% of units)",
                                                margin, 100 * mean(gap <= margin)),
                             raster = raster, raster_px = raster_px)
                },
                
                membership = {
                        M  <- .conf_matrix(x, confidence)
                        wh <- if (is.null(which)) seq_len(ncol(M)) else which
                        if (length(wh) > 12)
                                warning("Faceting ", length(wh),
                                        " regions; pass `which =` for a legible figure.",
                                        call. = FALSE)
                        df <- do.call(rbind, lapply(wh, function(r) {
                                d <- g; d$value <- M[, r]
                                d$region <- factor(sprintf("Region %d", r),
                                                   levels = sprintf("Region %d", wh))
                                d
                        }))
                        ggplot2::ggplot(df) +
                                ggplot2::geom_sf(ggplot2::aes(fill = .data$value),
                                                 colour = NA) +
                                ggplot2::facet_wrap(~ region) +
                                ggplot2::scale_fill_viridis_c(limits = c(0, 1), option = "mako",
                                                              direction = -1, name = confidence) +
                                ggplot2::labs(title = "Per-region membership surfaces") +
                                ggplot2::theme_void(base_size = 10)
                },
                
                coassociation = {
                        B <- x$fit$region_coassoc
                        if (is.null(B)) stop("Fit carries no `region_coassoc`.", call. = FALSE)
                        ord <- .seriate(1 - B)
                        df  <- expand.grid(i = seq_len(nrow(B)), j = seq_len(ncol(B)))
                        df$value <- as.vector(B[ord, ord])
                        df$i <- factor(ord[df$i], levels = ord)
                        df$j <- factor(ord[df$j], levels = ord)
                        ggplot2::ggplot(df, ggplot2::aes(.data$i, .data$j, fill = .data$value)) +
                                ggplot2::geom_tile() +
                                ggplot2::scale_fill_viridis_c(name = "co-association",
                                                              option = "rocket", direction = -1) +
                                ggplot2::coord_equal() +
                                ggplot2::labs(title = "Between-region co-association",
                                              subtitle = "ordered by similarity; bright off-diagonal = regions the ensemble confuses",
                                              x = "region", y = "region") +
                                ggplot2::theme_minimal(base_size = 10)
                },
                
                ordination = {
                        D <- .region_dissimilarity(x, NULL, NULL, x$k)
                        if (is.null(D)) stop("No region dissimilarity available.", call. = FALSE)
                        cfg <- as.data.frame(stats::cmdscale(stats::as.dist(D), k = 2))
                        names(cfg) <- c("Axis1", "Axis2")
                        cfg$region <- factor(seq_len(x$k))
                        cfg$size   <- as.vector(tabulate(x$hard_labels, x$k))
                        ggplot2::ggplot(cfg, ggplot2::aes(.data$Axis1, .data$Axis2)) +
                                ggplot2::geom_point(ggplot2::aes(colour = .data$region,
                                                                 size = .data$size)) +
                                ggplot2::geom_text(ggplot2::aes(label = .data$region),
                                                   size = 3, vjust = -1.2) +
                                ggplot2::scale_colour_manual(values = palette, guide = "none") +
                                ggplot2::scale_size_area(name = "units", max_size = 10) +
                                ggplot2::coord_equal() +
                                ggplot2::labs(title = "Regions in co-association space") +
                                ggplot2::theme_minimal(base_size = 10)
                },
                
                tuning = {
                        tl <- x$tuning_log
                        if (is.null(tl)) stop("This fit was not tuned.", call. = FALSE)
                        axes <- .TUNE_AXES[vapply(.TUNE_AXES,
                                                  function(a) length(unique(tl[[a]])) > 1L,
                                                  logical(1))]
                        if (!length(axes)) stop("No axis varied during tuning.", call. = FALSE)
                        df <- do.call(rbind, lapply(axes, function(a)
                                data.frame(axis = a, value = tl[[a]], score = tl$score)))
                        best <- data.frame(axis = axes,
                                           value = unlist(x$config[axes]),
                                           score = x$score)
                        ggplot2::ggplot(df, ggplot2::aes(.data$value, .data$score)) +
                                ggplot2::geom_point(alpha = 0.5) +
                                ggplot2::geom_vline(data = best,
                                                    ggplot2::aes(xintercept = .data$value),
                                                    linetype = 2, colour = "#B2182B") +
                                ggplot2::facet_wrap(~ axis, scales = "free_x") +
                                ggplot2::labs(title = sprintf("Tuning (%s, %s)", x$tuning_metric, x$tuning_strategy),
                                              subtitle = "each point is one configuration; dashed line = selected",
                                              x = NULL, y = "score (higher is better)") +
                                ggplot2::theme_bw(base_size = 10)
                })
}

#' @keywords internal
#' @noRd
.conf_matrix <- function(x, confidence) {
        if (confidence != "cohesion") return(x$memberships)
        M <- x$fit$cohesion
        if (is.null(M)) {
                warning("No cohesion matrix on this fit; using memberships.", call. = FALSE)
                return(x$memberships)
        }
        M
}

#' Per-unit confidence of a pulse_regions fit
#'
#' \code{"certainty"}: own-region \eqn{\Pi}. \code{"margin"}: best minus
#' second-best \eqn{\Pi}. \code{"cohesion"}: own-region cohesion. Hard
#' labels can be \code{NA} (uncovered units under an exact prior); those units
#' and \code{NA} rows get \code{NA}.
#' @keywords internal
#' @noRd
.confidence_vec <- function(x, confidence = c("certainty", "margin", "cohesion")) {
        confidence <- match.arg(confidence)
        M <- .conf_matrix(x, confidence)
        if (confidence == "margin") {
                if (ncol(M) < 2L) return(rep(1, nrow(M)))
                return(apply(M, 1, function(r) {
                        if (anyNA(r)) return(NA_real_)
                        s <- sort(r, decreasing = TRUE); s[1] - s[2]
                }))
        }
        h <- x$hard_labels
        out <- rep(NA_real_, nrow(M))
        ok <- !is.na(h)
        out[ok] <- M[cbind(which(ok), h[ok])]
        out
}

#' @keywords internal
#' @noRd
.continuous_map <- function(g, value, title, limits, raster, raster_px, direction = 1) {
        g$value <- value
        ggplot2::ggplot(g) +
                ggplot2::geom_sf(ggplot2::aes(fill = .data$value), colour = NA) +
                ggplot2::scale_fill_viridis_c(limits = limits, option = "mako",
                                              direction = direction, name = NULL) +
                ggplot2::labs(title = title) +
                ggplot2::theme_void(base_size = 11)
}

#' @keywords internal
#' @noRd
.region_key <- function(palette) {
        ## A discrete legend on top of scale_fill_identity: draw invisible points
        ## carrying the factor, so the identity fill stays cheap.
        df <- .legend_df(names(palette), palette)
        list(ggplot2::geom_point(data = df,
                                 ggplot2::aes(x = -Inf, y = -Inf, colour = .data$label),
                                 alpha = 0, show.legend = TRUE, inherit.aes = FALSE),
             ggplot2::scale_colour_manual(values = stats::setNames(unname(palette),
                                                                   names(palette)),
                                          name = "region"),
             ggplot2::guides(colour = ggplot2::guide_legend(
                     override.aes = list(alpha = 1, size = 3, shape = 15))))
}

#' @keywords internal
#' @noRd
.seriate <- function(D) {
        s <- stats::cmdscale(stats::as.dist(D), k = 1)[, 1]
        order(s)
}


# -----------------------------------------------------------------------------
# Level 2: local types
# -----------------------------------------------------------------------------

#' Plot fitted local types
#'
#' @param x A \code{"pulse_local_types"} object.
#' @param what One of:
#'   \describe{
#'     \item{\code{"map"}}{All local types, shaded within region hues.}
#'     \item{\code{"region"}}{One region's local types, everything else greyed
#'       out (see \code{which}). The readable way to show nesting when
#'       \eqn{\sum G_r} runs past about 30 classes.}
#'     \item{\code{"facet"}}{One panel per region, that region's types only.}
#'     \item{\code{"confidence"}}{Local-type confidence surface.}
#'     \item{\code{"bic"}}{Weighted BIC against G, per region — the model
#'       selection audit. A curve still rising at the largest candidate means
#'       \code{n_local_types} was capped too low for that region.}
#'     \item{\code{"profiles"}}{Membership-weighted predictor means per local
#'       type, z-scored per predictor. Requires \code{env}.}
#'     \item{\code{"sizes"}}{Effective size of each local type.}
#'   }
#' @param geom An \code{sf} object, if not already attached.
#' @param regions The parent \code{"pulse_regions"} object. Supplying it keeps
#'   the hues identical to the level-1 figures.
#' @param palette Local-type colors; defaults to \code{\link{local_palette}}
#'   built on \code{\link{pulse_palette}(regions)}.
#' @param confidence \code{"own"} (default) uses \eqn{z_{irg}}, the 
#'   within-region probability of the assigned local type, ignoring how sure 
#'   Level 1 is about the region; \code{"joint"} uses \eqn{w_{ir} z_{irg}} of 
#'   the assigned local type: Level 1
#'   region membership × within-region type probability. Use \code{"joint"} for 
#'   any figure that stands alone.
#' @param which Integer region, for \code{what = "region"}.
#' @param env Predictor matrix, for \code{what = "profiles"}.
#' @param raster,raster_px Rasterization control.
#' @param ... Unused.
#' @return A \code{ggplot} object.
#' @export
plot_local_types <- function(x,
                             what       = c("map", "region", "facet", "confidence",
                                            "bic", "profiles", "sizes"),
                             geom       = NULL,
                             regions    = NULL,
                             palette    = NULL,
                             confidence = c("own", "joint"),
                             which      = NULL,
                             env        = NULL,
                             raster     = "auto",
                             raster_px  = 2000,
                             ...) {
        .need("ggplot2")
        stopifnot(inherits(x, "pulse_local_types"))
        what       <- match.arg(what)
        confidence <- match.arg(confidence)
        
        if (is.null(palette)) {
                rc <- if (!is.null(regions)) pulse_palette(regions) else
                        stats::setNames(grDevices::hcl(
                                h = seq(15, 375, length.out = length(x$n_local_types) + 1)[
                                        seq_len(length(x$n_local_types))],
                                c = 95, l = 65), seq_along(x$n_local_types))
                palette <- local_palette(x, rc)
        }
        conf <- .local_confidence(x, regions, confidence)
        needs_geom <- what %in% c("map", "region", "facet", "confidence")
        g <- if (needs_geom) .get_geom(x, geom) else NULL
        
        switch(
                what,
                
                map = .map(g, fade_by_confidence(palette[x$hard_labels], conf),
                           title    = sprintf("Local types (%d across %d regions)",
                                              ncol(x$memberships), length(x$n_local_types)),
                           subtitle = sprintf("hue = region, shade = local type, saturation = %s confidence",
                                              confidence),
                           raster = raster, raster_px = raster_px),
                
                region = {
                        if (is.null(which)) stop("`what = \"region\"` needs `which = <region>`.",
                                                 call. = FALSE)
                        cols <- x$key$column[x$key$region == which]
                        inr  <- x$hard_labels %in% cols
                        fill <- rep("#F0F0F0", length(inr))
                        fill[inr] <- fade_by_confidence(palette[x$hard_labels[inr]], conf[inr])
                        .map(g, fill,
                             title    = sprintf("Region %d: %d local types", which, length(cols)),
                             subtitle = "other regions greyed out",
                             raster = raster, raster_px = raster_px)
                },
                
                facet = {
                        df <- do.call(rbind, lapply(sort(unique(x$key$region)), function(r) {
                                cols <- x$key$column[x$key$region == r]
                                inr  <- x$hard_labels %in% cols
                                d <- g
                                d$fill <- ifelse(inr,
                                                 fade_by_confidence(palette[x$hard_labels], conf),
                                                 "#F2F2F2")
                                d$region <- factor(sprintf("Region %d", r))
                                d
                        }))
                        ggplot2::ggplot(df) +
                                ggplot2::geom_sf(ggplot2::aes(fill = .data$fill), colour = NA) +
                                ggplot2::facet_wrap(~ region) +
                                ggplot2::scale_fill_identity() +
                                ggplot2::labs(title = "Local types by region") +
                                ggplot2::theme_void(base_size = 10)
                },
                
                confidence = .continuous_map(g, conf,
                                             sprintf("Local-type confidence (%s)", confidence),
                                             limits = c(0, 1), raster, raster_px),
                
                bic = {
                        B  <- x$BIC
                        df <- expand.grid(G = as.integer(rownames(B)), region = colnames(B),
                                          stringsAsFactors = FALSE)
                        df$BIC <- as.vector(B)
                        sel <- data.frame(region = colnames(B), G = as.integer(x$n_local_types))
                        ggplot2::ggplot(df, ggplot2::aes(.data$G, .data$BIC)) +
                                ggplot2::geom_line(na.rm = TRUE) +
                                ggplot2::geom_point(na.rm = TRUE, size = 1) +
                                ggplot2::geom_vline(data = sel,
                                                    ggplot2::aes(xintercept = .data$G),
                                                    linetype = 2, colour = "#B2182B") +
                                ggplot2::facet_wrap(~ region, scales = "free_y") +
                                ggplot2::labs(title = "Weighted BIC by number of local types",
                                              subtitle = "dashed line = selected; a curve still rising at the right edge means n_local_types was capped too low",
                                              x = "G", y = "BIC") +
                                ggplot2::theme_bw(base_size = 9)
                },
                
                profiles = {
                        if (is.null(env)) stop("`what = \"profiles\"` needs `env`.", call. = FALSE)
                        X <- as.matrix(env[, x$variables, drop = FALSE])
                        U <- x$memberships
                        ## Membership-weighted centroid: the fuzzy analogue of a group
                        ## mean, so partially-assigned units contribute partially.
                        cent <- (t(U) %*% X) / pmax(colSums(U), .Machine$double.eps)
                        cent <- scale(cent)
                        df <- expand.grid(type = rownames(cent) %||% x$key$name,
                                          variable = colnames(cent), stringsAsFactors = FALSE)
                        df$z <- as.vector(cent)
                        df$type <- factor(df$type, levels = x$key$name)
                        ggplot2::ggplot(df, ggplot2::aes(.data$variable, .data$type, fill = .data$z)) +
                                ggplot2::geom_tile() +
                                ggplot2::scale_fill_gradient2(name = "z", low = "#2166AC",
                                                              mid = "white", high = "#B2182B") +
                                ggplot2::labs(title = "Local-type environmental profiles",
                                              subtitle = "membership-weighted means, z-scored per predictor",
                                              x = NULL, y = NULL) +
                                ggplot2::theme_minimal(base_size = 9) +
                                ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
                },
                
                sizes = {
                        df <- x$key
                        df$n_eff <- colSums(x$memberships)
                        df$name  <- factor(df$name, levels = rev(df$name))
                        ggplot2::ggplot(df, ggplot2::aes(.data$n_eff, .data$name,
                                                         fill = .data$name)) +
                                ggplot2::geom_col() +
                                ggplot2::scale_fill_manual(values = palette, guide = "none") +
                                ggplot2::labs(title = "Effective size of each local type",
                                              x = "sum of memberships", y = NULL) +
                                ggplot2::theme_minimal(base_size = 9)
                })
}

#' @keywords internal
#' @noRd
.local_confidence <- function(x, regions, confidence) {
        if (confidence == "joint") {
                ## certainty_raw is w_ir * z_irg of the assigned local type. Level 1
                ## membership is already in it, so it is not multiplied again.
                return(.first_non_null(x$certainty_raw, x$certainty))
        }
        ## "own": z_irg of the assigned type, conditional on its region, i.e.
        ## Zw[i, g] / w_ir. Since z rows sum to 1 within a region, w_ir is the
        ## sum of Zw over that region's columns. No refit or extra slot needed.
        Zw <- x$weighted
        if (is.null(Zw) || is.null(x$key)) return(x$certainty)
        reg  <- x$key$region
        hard <- x$hard_labels
        Wr   <- t(rowsum(t(Zw), reg))                      # N x R region weights
        w_i  <- Wr[cbind(seq_len(nrow(Zw)), match(as.character(reg[hard]), colnames(Wr)))]
        z    <- Zw[cbind(seq_len(nrow(Zw)), hard)] / w_i
        z[!is.finite(z)] <- NA_real_                       # all-zero rows; fade maps NA -> 0
        z
}

#' @keywords internal
#' @noRd
`%||%` <- function(a, b) if (is.null(a)) b else a


# -----------------------------------------------------------------------------
# Both levels together
# -----------------------------------------------------------------------------

#' Region and local-type maps side by side
#'
#' Uses one palette for both panels, so the nesting is visible: the right panel
#' is the left panel subdivided into shades of the same hues.
#'
#' @param regions A \code{"pulse_regions"} object.
#' @param local A \code{"pulse_local_types"} object fitted on it.
#' @param geom An \code{sf} object, if not already attached.
#' @param confidence Passed to both panels; \code{"joint"} on the right.
#' @param ... Passed to both plotting calls.
#' @return A \code{patchwork} object if \pkg{patchwork} is installed, otherwise
#'   a list of two \code{ggplot} objects.
#' @export
plot_region_local <- function(regions, local, geom = NULL, confidence = "joint", ...) {
        .need("ggplot2")
        pal_r <- pulse_palette(regions)
        pal_l <- local_palette(local, pal_r)
        g <- .get_geom(regions, geom)
        p1 <- plot_regions(regions, "map", geom = g, palette = pal_r, ...)
        p2 <- plot_local_types(local, "map", geom = g, regions = regions,
                               palette = pal_l, confidence = confidence, ...)
        if (requireNamespace("patchwork", quietly = TRUE)) p1 + p2 else list(p1, p2)
}


# =============================================================================
# pulseR — visualisation, levels 2b and 3
#
#   Level 2b  pulse_meta_types   N units x M cross-region meta-types
#   Level 3   pulse_mosaic_types       clustered on the COMPLETE-CASE rows only
#
# Follows the contract in the level-1/2 plotting file: one similarity-ordered
# palette per fit, confidence baked into the hex fill rather than applied as
# alpha, geometry attached once with attach_geometry().
#
# Two things differ from the earlier levels and drive the design here:
#
#   * Meta-types deliberately have no parent region, so they CANNOT inherit
#     region hues. Giving them their own palette family is the point, not an
#     oversight: if a meta-type map looked like the region map, the merge had
#     failed.
#   * A mosaic fit covers only the polygons with a complete signature. Isolated
#     units have no mosaic type at all and must render as absent, not as a
#     silent extra class. .scatter_valid() handles that everywhere.
#
# Replaces plot_meta_types() from meta_types.R and the older
# plot_mosaic_signatures() / plot_mosaic_ordination().
# =============================================================================


# -----------------------------------------------------------------------------
# 1. Shared palette machinery
# -----------------------------------------------------------------------------

#' Similarity-ordered HCL palette from a dissimilarity matrix
#'
#' Factored out of \code{pulse_palette()} so every level uses one
#' implementation. Angle in the ordination becomes hue and radius becomes
#' chroma, so classes at the centre of the configuration are intermediate and
#' come out muted while extreme classes get saturated colours.
#'
#' @param D A \eqn{k \times k} dissimilarity, or \code{NULL} for a qualitative
#'   ramp.
#' @param k Number of classes.
#' @param method \code{"auto"}, \code{"ramp1d"}, \code{"hue2d"} or
#'   \code{"qualitative"}.
#' @param chroma,luminance,hue_range HCL controls.
#' @param hue_offset Degrees added to every hue. Use it to give a level its own
#'   colour family so readers do not confuse two maps of the same coverage.
#' @return Character vector of \code{k} hex colours.
#' @keywords internal
#' @noRd
.similarity_palette <- function(D, k,
                                method    = c("auto", "ramp1d", "hue2d", "qualitative"),
                                chroma    = 65,
                                luminance = 65,
                                hue_range = c(10, 340),
                                hue_offset = 0) {
        method <- match.arg(method)
        if (method == "auto") method <- if (k <= 8 || is.null(D)) "qualitative" else "hue2d"
        if (is.null(D) && method != "qualitative") {
                warning("No dissimilarity available; falling back to a qualitative palette.",
                        call. = FALSE)
                method <- "qualitative"
        }
        cols <- switch(
                method,
                qualitative = grDevices::hcl(
                        h = (seq(15, 375, length.out = k + 1)[seq_len(k)] + hue_offset) %% 360,
                        c = chroma + 30, l = luminance),
                ramp1d = {
                        s <- stats::cmdscale(stats::as.dist(D), k = 1)[, 1]
                        h <- hue_range[1] + diff(hue_range) * .rescale01(s)
                        grDevices::hcl(h = (h + hue_offset) %% 360, c = chroma + 20, l = luminance)
                },
                hue2d = {
                        cfg <- stats::cmdscale(stats::as.dist(D), k = 2)
                        h <- (atan2(cfg[, 2], cfg[, 1]) * 180 / pi + hue_offset) %% 360
                        r <- sqrt(rowSums(cfg^2))
                        grDevices::hcl(h = h, c = 20 + chroma * .rescale01(r), l = luminance)
                })
        stats::setNames(cols, seq_len(k))
}

#' Colours for meta-types
#'
#' Ordered by average between-meta-type distance in the space the merge was
#' computed in, so environmentally similar meta-types get neighbouring hues.
#' Offset in hue from the region palette by default, because a meta-type map
#' that resembled the region map would mean the merge had not decoupled class
#' from location.
#'
#' @param x A \code{"pulse_meta_types"} object.
#' @param method,chroma,luminance,hue_offset Passed to the palette builder.
#' @return Named character vector of \code{n_meta_types} hex colours.
#' @export
meta_palette <- function(x, method = "auto", chroma = 65, luminance = 65,
                         hue_offset = 40) {
        stopifnot(inherits(x, "pulse_meta_types"))
        D <- .meta_dissimilarity(x)
        .similarity_palette(D, x$n_meta_types, method, chroma, luminance,
                            hue_offset = hue_offset)
}

#' @keywords internal
#' @noRd
.meta_dissimilarity <- function(x) {
        Dl <- as.matrix(x$distance)
        cl <- x$key$meta
        M  <- x$n_meta_types
        if (M < 2L) return(NULL)
        out <- matrix(0, M, M)
        for (i in seq_len(M)) for (j in seq_len(M)) {
                out[i, j] <- mean(Dl[cl == i, cl == j, drop = FALSE])
        }
        diag(out) <- 0
        (out + t(out)) / 2
}

#' Colours for mosaic types
#'
#' Ordered by the distance between the medoid signatures, in the same metric the
#' fit used.
#'
#' @param x A \code{"pulse_mosaic_types"} object.
#' @param method,chroma,luminance,hue_offset Passed to the palette builder.
#' @return Named character vector of \code{n_mosaic_types} hex colours.
#' @export
mosaic_palette <- function(x, method = "auto", chroma = 60, luminance = 60,
                           hue_offset = 200) {
        stopifnot(inherits(x, "pulse_mosaic_types"))
        D <- tryCatch(
                as.matrix(compute_signature_distance(x$medoid_signatures, x$distance)),
                error = function(e) NULL)
        .similarity_palette(D, x$n_mosaic_types, method, chroma, luminance,
                            hue_offset = hue_offset)
}


# -----------------------------------------------------------------------------
# 2. Helpers these levels need
# -----------------------------------------------------------------------------

#' Unit count of a fitted object
#'
#' A mosaic fit is sized by its \code{valid} vector, not by its membership matrix, because the clustering
#' saw only the complete-case rows.
#' @keywords internal
#' @noRd
.n_units <- function(x) {
        if (inherits(x, "pulse_mosaic_types"))      return(length(x$valid))
        if (inherits(x, "pulse_meta_types"))  return(nrow(x$memberships))
        if (inherits(x, "pulse_local_types")) return(nrow(x$memberships))
        if (inherits(x, "pulse_regions"))     return(nrow(x$memberships))
        stop("Unsupported object.", call. = FALSE)
}

#' Scatter values computed on the clustered subset back onto every polygon
#'
#' Isolated units have no signature and therefore no mosaic type. They get
#' \code{fill_na}, never a class of their own.
#' @keywords internal
#' @noRd
.scatter_valid <- function(values, valid, fill_na = NA) {
        out <- rep(fill_na, length(valid))
        out[valid] <- values
        out
}

#' Rescale a row-maximum membership onto \\[0, 1\\]
#'
#' A row of \eqn{k} memberships summing to one has maximum \eqn{1/k} when the
#' unit is maximally undecided, not 0. Fading on the raw maximum therefore makes
#' every map look uniformly uncertain as \eqn{k} grows. This maps the useful
#' range onto the full colour range.
#' @keywords internal
#' @noRd
.rescale_certainty <- function(own, k) {
        if (k < 2L) return(rep(1, length(own)))
        pmin(pmax((own - 1 / k) / (1 - 1 / k), 0), 1)
}

#' @keywords internal
#' @noRd
.split_pairs <- function(cn) {
        sp <- strsplit(cn, ":", fixed = TRUE)
        list(a = vapply(sp, `[`, character(1), 1L),
             b = vapply(sp, `[`, character(1), 2L))
}


# -----------------------------------------------------------------------------
# 3. Meta-types
# -----------------------------------------------------------------------------

#' Plot merged meta-types
#'
#' @param x A \code{"pulse_meta_types"} object.
#' @param what One of:
#'   \describe{
#'     \item{\code{"map"}}{Meta-types in space, faded by confidence.}
#'     \item{\code{"composition"}}{The diagnostic figure. One bar per meta-type,
#'       segments coloured by the parent region of each local type it absorbed.
#'       Multi-coloured bars mean types recur across regions and the merge has
#'       decoupled class from location; single-coloured bars mean it has not,
#'       and the mosaic level will keep reproducing the region map.}
#'     \item{\code{"overlap"}}{Local-type distance matrix ordered by the merge
#'       tree, with a region annotation strip.}
#'     \item{\code{"selection"}}{Silhouette and ARI-with-region against
#'       \code{n_meta_types}. Choose on the second: cluster tightness is not the
#'       criterion here.}
#'     \item{\code{"confidence"}, \code{"entropy"}}{Uncertainty surfaces.}
#'     \item{\code{"profiles"}}{Membership-weighted predictor means per
#'       meta-type. Needs \code{env} and \code{variables}.}
#'     \item{\code{"sizes"}}{Effective size of each meta-type.}
#'   }
#' @param geom An \code{sf} object, if not attached.
#' @param palette Meta-type colours; defaults to \code{\link{meta_palette}(x)}.
#' @param region_colours Region colours for the annotation in
#'   \code{"composition"} and \code{"overlap"}, normally
#'   \code{\link{pulse_palette}} on the level-1 fit.
#' @param env,variables Predictor matrix and column names for
#'   \code{what = "profiles"}.
#' @param rescale_confidence Map the row maximum from \eqn{[1/M, 1]} onto
#'   \eqn{[0, 1]} before fading. Default \code{TRUE}.
#' @param raster,raster_px Rasterisation control.
#' @param ... Unused.
#' @return A \code{ggplot} object.
#' @export
plot_meta_types <- function(x,
                            what    = c("map", "composition", "overlap", "selection",
                                        "confidence", "entropy", "profiles", "sizes"),
                            geom    = NULL,
                            palette = NULL,
                            region_colours = NULL,
                            env     = NULL,
                            variables = NULL,
                            rescale_confidence = TRUE,
                            raster  = "auto",
                            raster_px = 2000,
                            ...) {
        .need("ggplot2")
        stopifnot(inherits(x, "pulse_meta_types"))
        what <- match.arg(what)
        if (is.null(palette)) palette <- meta_palette(x)
        
        conf <- if (rescale_confidence) .rescale_certainty(x$certainty, x$n_meta_types) else x$certainty
        needs_geom <- what %in% c("map", "confidence", "entropy")
        g <- if (needs_geom) .get_geom(x, geom) else NULL
        
        switch(
                what,
                
                map = .map(g, fade_by_confidence(palette[x$hard_labels], conf),
                           title    = sprintf("Meta-types (%d, merged from %d local types)",
                                              x$n_meta_types, nrow(x$key)),
                           subtitle = sprintf("ARI with region = %.2f; low means these are not region proxies",
                                              x$crossing$ari_with_region),
                           raster = raster, raster_px = raster_px) +
                        .class_key(palette, "meta-type"),
                
                composition = {
                        df <- x$key
                        df$meta_f   <- factor(df$meta, levels = seq_len(x$n_meta_types))
                        df$region_f <- factor(df$region)
                        p <- ggplot2::ggplot(df, ggplot2::aes(y = .data$meta_f,
                                                              fill = .data$region_f)) +
                                ggplot2::geom_bar() +
                                ggplot2::labs(
                                        title = "What each meta-type absorbed",
                                        subtitle = sprintf(
                                                "%.0f%% of meta-types draw on a single region | ARI with region = %.2f",
                                                100 * x$crossing$frac_single_region,
                                                x$crossing$ari_with_region),
                                        x = "local types", y = "meta-type", fill = "region") +
                                ggplot2::theme_minimal(base_size = 10)
                        if (!is.null(region_colours))
                                p <- p + ggplot2::scale_fill_manual(values = region_colours)
                        p
                },
                
                overlap = {
                        D   <- as.matrix(x$distance)
                        ord <- x$hclust$order
                        lab <- rownames(D)[ord]
                        df  <- expand.grid(i = seq_along(ord), j = seq_along(ord))
                        df$value <- as.vector(D[ord, ord])
                        df$i <- factor(lab[df$i], levels = lab)
                        df$j <- factor(lab[df$j], levels = rev(lab))
                        ann <- data.frame(type   = factor(lab, levels = lab),
                                          region = factor(x$key$region[ord]))
                        p <- ggplot2::ggplot(df, ggplot2::aes(.data$i, .data$j,
                                                              fill = .data$value)) +
                                ggplot2::geom_tile() +
                                ggplot2::scale_fill_viridis_c(name = "distance", option = "mako",
                                                              direction = -1) +
                                ggplot2::labs(
                                        title = "Local-type overlap, ordered by the merge tree",
                                        subtitle = "blocks that mix region colours below the axis are types recurring across regions",
                                        x = NULL, y = NULL) +
                                ggplot2::theme_minimal(base_size = 8) +
                                ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5),
                                               panel.grid = ggplot2::element_blank())
                        p <- p + ggplot2::geom_point(
                                data = ann,
                                ggplot2::aes(x = .data$type, y = 0.4, colour = .data$region),
                                inherit.aes = FALSE, size = 1.5, shape = 15)
                        if (is.null(region_colours))
                                p + ggplot2::scale_colour_discrete(name = "region")
                        else p + ggplot2::scale_colour_manual(values = region_colours, name = "region")
                },
                
                selection = {
                        sel <- x$selection
                        df <- rbind(
                                data.frame(n_meta_types = sel$n_meta_types, value = sel$silhouette,
                                           panel = "silhouette (tightness)"),
                                data.frame(n_meta_types = sel$n_meta_types, value = sel$ari_with_region,
                                           panel = "ARI with region (lower is better)"),
                                data.frame(n_meta_types = sel$n_meta_types, value = sel$mean_regions_per_meta,
                                           panel = "mean regions per meta-type"))
                        ggplot2::ggplot(df, ggplot2::aes(.data$n_meta_types, .data$value)) +
                                ggplot2::geom_line() + ggplot2::geom_point(size = 1) +
                                ggplot2::geom_vline(xintercept = x$n_meta_types, linetype = 2,
                                                    colour = "#B2182B") +
                                ggplot2::facet_wrap(~ panel, scales = "free_y", ncol = 1) +
                                ggplot2::labs(title = "Choosing the number of meta-types",
                                              subtitle = "dashed line = selected; the middle panel is the one that matters",
                                              x = "n_meta_types", y = NULL) +
                                ggplot2::theme_bw(base_size = 9)
                },
                
                confidence = .continuous_map(g, conf,
                                             "Meta-type confidence", c(0, 1), raster, raster_px),
                
                entropy = .continuous_map(g, x$entropy, "Meta-type membership entropy",
                                          c(0, 1), raster, raster_px, direction = -1),
                
                profiles = {
                        if (is.null(env)) stop("`what = \"profiles\"` needs `env`.", call. = FALSE)
                        vs <- if (is.null(variables)) colnames(env) else variables
                        X <- as.matrix(env[, vs, drop = FALSE])
                        U <- x$memberships
                        cent <- scale((t(U) %*% X) / pmax(colSums(U), .Machine$double.eps))
                        df <- expand.grid(type = colnames(U), variable = colnames(cent),
                                          stringsAsFactors = FALSE)
                        df$z <- as.vector(cent)
                        df$type <- factor(df$type, levels = colnames(U))
                        ggplot2::ggplot(df, ggplot2::aes(.data$variable, .data$type,
                                                         fill = .data$z)) +
                                ggplot2::geom_tile() +
                                ggplot2::scale_fill_gradient2(name = "z", low = "#2166AC",
                                                              mid = "white", high = "#B2182B") +
                                ggplot2::labs(title = "Meta-type environmental profiles",
                                              subtitle = "membership-weighted means, z-scored per predictor",
                                              x = NULL, y = NULL) +
                                ggplot2::theme_minimal(base_size = 9) +
                                ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
                },
                
                sizes = {
                        df <- data.frame(type = factor(colnames(x$memberships),
                                                       levels = rev(colnames(x$memberships))),
                                         n_eff = colSums(x$memberships),
                                         n_local = as.vector(table(factor(x$key$meta,
                                                                          levels = seq_len(x$n_meta_types)))))
                        ggplot2::ggplot(df, ggplot2::aes(.data$n_eff, .data$type,
                                                         fill = .data$type)) +
                                ggplot2::geom_col() +
                                ggplot2::geom_text(ggplot2::aes(label = sprintf("%d local", .data$n_local)),
                                                   hjust = -0.1, size = 2.5) +
                                ggplot2::scale_fill_manual(values = stats::setNames(
                                        palette, colnames(x$memberships)), guide = "none") +
                                ggplot2::labs(title = "Effective size of each meta-type",
                                              x = "sum of memberships", y = NULL) +
                                ggplot2::theme_minimal(base_size = 9)
                })
}


# -----------------------------------------------------------------------------
# 4. Mosaic types
# -----------------------------------------------------------------------------

#' Plot fitted mosaic types
#'
#' @param x A \code{"pulse_mosaic_types"} object.
#' @param what One of:
#'   \describe{
#'     \item{\code{"map"}}{Mosaic types in space, faded by confidence. Polygons
#'       with no signature render in \code{na_colour}, not as a class.}
#'     \item{\code{"confidence"}, \code{"entropy"}}{Uncertainty surfaces.}
#'     \item{\code{"signatures"}}{The fingerprint figure: mean class-by-class
#'       co-occurrence per mosaic type, lower triangle, never-occurring
#'       combinations in white. Needs \code{signatures}.}
#'     \item{\code{"ordination"}}{Signatures in two dimensions, coloured by
#'       mosaic type. Needs \code{signatures}.}
#'     \item{\code{"validity"}}{Every validity index against \code{k}, with the
#'       selection marked.}
#'     \item{\code{"saturation"}}{Distribution of pairwise signature distances
#'       against the metric ceiling. Read this before any of the others: a spike
#'       at the ceiling means the partition is close to arbitrary. Needs
#'       \code{signatures}.}
#'     \item{\code{"sizes"}}{Effective size of each mosaic type.}
#'   }
#' @param signatures The signature matrix the fit was built on, row-aligned to
#'   \code{x$valid}.
#' @param geom An \code{sf} object, if not attached.
#' @param palette Mosaic colours; defaults to \code{\link{mosaic_palette}(x)}.
#' @param class_labels Optional shorter names for the classes appearing in the
#'   signature column names, for \code{what = "signatures"}.
#' @param na_colour Fill for polygons with no signature. Default light grey.
#' @param rescale_confidence Map the row maximum from \eqn{[1/k, 1]} onto
#'   \eqn{[0, 1]} before fading. Default \code{TRUE}.
#' @param n_sample Signatures sampled for \code{"ordination"} and
#'   \code{"saturation"}. Default 3000.
#' @param raster,raster_px Rasterisation control.
#' @param ... Unused.
#' @return A \code{ggplot} object.
#' @export
plot_mosaic_types <- function(x,
                              what = c("map", "confidence", "entropy", "signatures",
                                       "ordination", "validity", "saturation", "sizes"),
                              signatures = NULL,
                              geom = NULL,
                              palette = NULL,
                              class_labels = NULL,
                              na_colour = "#EDEDED",
                              rescale_confidence = TRUE,
                              n_sample = 3000L,
                              raster = "auto",
                              raster_px = 2000,
                              ...) {
        .need("ggplot2")
        stopifnot(inherits(x, "pulse_mosaic_types"))
        what <- match.arg(what)
        if (is.null(palette)) palette <- mosaic_palette(x)
        
        k    <- x$n_mosaic_types
        own  <- x$memberships[cbind(seq_len(nrow(x$memberships)), x$hard_labels)]
        conf <- if (rescale_confidence) .rescale_certainty(own, k) else own
        
        needs_geom <- what %in% c("map", "confidence", "entropy")
        g <- if (needs_geom) .get_geom(x, geom) else NULL
        needs_sigs <- what %in% c("signatures", "ordination", "saturation")
        if (needs_sigs && is.null(signatures))
                stop(sprintf("`what = \"%s\"` needs the signature matrix; pass `signatures =`.", what),
                     call. = FALSE)
        S <- if (needs_sigs) as.matrix(signatures)[x$valid, , drop = FALSE] else NULL
        
        switch(
                what,
                
                map = .map(g,
                           .scatter_valid(fade_by_confidence(palette[x$hard_labels], conf),
                                          x$valid, na_colour),
                           title    = sprintf("Mosaic types (k = %d)", k),
                           subtitle = sprintf("configuration of classes in the neighbourhood | %s distance%s",
                                              x$distance,
                                              if (is.finite(x$saturation))
                                                      sprintf(", saturation %.2f", x$saturation) else ""),
                           raster = raster, raster_px = raster_px) +
                        .class_key(palette, "mosaic type"),
                
                confidence = .continuous_map(g, .scatter_valid(conf, x$valid),
                                             "Mosaic-type confidence", c(0, 1),
                                             raster, raster_px),
                
                entropy = .continuous_map(g, .scatter_valid(.row_entropy(x$memberships), x$valid),
                                          "Mosaic-type membership entropy", c(0, 1),
                                          raster, raster_px, direction = -1),
                
                signatures = {
                        pr  <- .split_pairs(colnames(S))
                        cls <- if (is.null(class_labels)) sort(unique(c(pr$a, pr$b))) else class_labels
                        lab <- x$hard_labels
                        rows <- lapply(sort(unique(lab)), function(gg) {
                                m <- colMeans(S[lab == gg, , drop = FALSE], na.rm = TRUE)
                                data.frame(cluster = factor(paste0("Mosaic ", gg)),
                                           a = pr$a, b = pr$b, value = m,
                                           stringsAsFactors = FALSE)
                        })
                        df <- do.call(rbind, rows)
                        ia <- match(df$a, cls); ib <- match(df$b, cls)
                        ## Keep the lower triangle plus the diagonal. The signature already
                        ## holds each unordered pair once, so no mirroring is needed; the
                        ## old version built both halves and then discarded one.
                        swap <- ib < ia
                        df$a[swap] <- pr$b[swap]; df$b[swap] <- pr$a[swap]
                        df$a <- factor(df$a, levels = cls)
                        df$b <- factor(df$b, levels = rev(cls))
                        ## Zero means "never co-occurs"; on a sequential scale that is the
                        ## dark end, which reads as a strong signal. White it out.
                        df$value[df$value <= 0] <- NA_real_
                        ggplot2::ggplot(df, ggplot2::aes(.data$a, .data$b, fill = .data$value)) +
                                ggplot2::geom_tile(colour = "white", linewidth = 0.2) +
                                ggplot2::facet_wrap(~ cluster) +
                                ggplot2::scale_fill_viridis_c(option = "G", direction = -1,
                                                              name = "mean\nco-occur.",
                                                              na.value = "white") +
                                ggplot2::coord_equal() +
                                ggplot2::labs(title = "Mosaic-type co-occurrence fingerprints",
                                              subtitle = "mean class-by-class adjacency within each mosaic type; white = never observed",
                                              x = NULL, y = NULL) +
                                ggplot2::theme_minimal(base_size = 8) +
                                ggplot2::theme(panel.grid = ggplot2::element_blank(),
                                               axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5))
                },
                
                ordination = {
                        i <- if (nrow(S) > n_sample) sort(sample.int(nrow(S), n_sample)) else seq_len(nrow(S))
                        ## Uses the metric the fit used.
                        d   <- compute_signature_distance(S[i, , drop = FALSE], x$distance)
                        mds <- stats::cmdscale(d, k = 2)
                        df  <- data.frame(MDS1 = mds[, 1], MDS2 = mds[, 2],
                                          cluster = factor(x$hard_labels[i]))
                        ggplot2::ggplot(df, ggplot2::aes(.data$MDS1, .data$MDS2,
                                                         colour = .data$cluster)) +
                                ggplot2::geom_point(alpha = 0.5, size = 1.2) +
                                ggplot2::stat_ellipse(level = 0.95, linewidth = 0.6) +
                                ggplot2::scale_colour_manual(values = palette, name = "mosaic") +
                                ggplot2::coord_equal() +
                                ggplot2::labs(title = "Mosaic signature ordination",
                                              subtitle = sprintf("classical MDS of %s distances, %d signatures",
                                                                 x$distance, length(i)),
                                              x = "axis 1", y = "axis 2") +
                                ggplot2::theme_minimal(base_size = 10)
                },
                
                validity = {
                        v <- x$validity
                        cols <- intersect(c("R2_JS", "CH_JS", "GD5", "STAB",
                                            "SILH_fuzzy", "SILH_HARD", "score"), names(v))
                        df <- do.call(rbind, lapply(cols, function(cc)
                                data.frame(k = v$k, value = v[[cc]], index = cc)))
                        df <- df[is.finite(df$value), , drop = FALSE]
                        ggplot2::ggplot(df, ggplot2::aes(.data$k, .data$value)) +
                                ggplot2::geom_line() + ggplot2::geom_point(size = 1) +
                                ggplot2::geom_vline(xintercept = k, linetype = 2, colour = "#B2182B") +
                                ggplot2::facet_wrap(~ index, scales = "free_y") +
                                ggplot2::labs(title = "Mosaic validity indices",
                                              subtitle = sprintf("dashed line = selected k = %d | saturation %.2f",
                                                                 k, x$saturation),
                                              x = "number of mosaic types", y = NULL) +
                                ggplot2::theme_bw(base_size = 9)
                },
                
                saturation = {
                        i  <- sample.int(nrow(S), min(n_sample, nrow(S)))
                        j  <- sample.int(nrow(S), min(n_sample, nrow(S)))
                        ok <- i != j
                        dv <- diag(compute_signature_distance(S[i[ok], , drop = FALSE], x$distance,
                                                      signatures_ref = S[j[ok], , drop = FALSE]))
                        ceil <- switch(x$distance, js = sqrt(log(2)), hellinger = 1,
                                       tv = 1, bray = 1, NA_real_)
                        p <- ggplot2::ggplot(data.frame(d = dv), ggplot2::aes(.data$d)) +
                                ggplot2::geom_histogram(bins = 60, fill = "#4C72B0") +
                                ggplot2::labs(title = "Pairwise signature distances",
                                              subtitle = sprintf("%s metric | %.0f%% of pairs at the ceiling",
                                                                 x$distance, 100 * x$saturation),
                                              x = "distance", y = "pairs") +
                                ggplot2::theme_minimal(base_size = 10)
                        if (is.finite(ceil))
                                p <- p + ggplot2::geom_vline(xintercept = ceil, linetype = 2,
                                                             colour = "#B2182B") +
                                ggplot2::annotate("text", x = ceil, y = Inf, hjust = 1.1,
                                                  vjust = 1.5, size = 3, label = "ceiling")
                        p
                },
                
                sizes = {
                        nm <- paste0("Mosaic ", seq_len(k))
                        df <- data.frame(type = factor(nm, levels = rev(nm)),
                                         n_eff = colSums(x$memberships),
                                         n_hard = as.vector(tabulate(x$hard_labels, k)))
                        ggplot2::ggplot(df, ggplot2::aes(.data$n_eff, .data$type,
                                                         fill = .data$type)) +
                                ggplot2::geom_col() +
                                ggplot2::scale_fill_manual(values = stats::setNames(palette, nm),
                                                           guide = "none") +
                                ggplot2::labs(title = "Effective size of each mosaic type",
                                              x = "sum of memberships", y = NULL) +
                                ggplot2::theme_minimal(base_size = 9)
                })
}


# -----------------------------------------------------------------------------
# 5. Discrete legend on an identity fill
# -----------------------------------------------------------------------------

#' Generalisation of .region_key() to any level
#' @keywords internal
#' @noRd
.class_key <- function(palette, name = "class") {
        df <- data.frame(label = factor(names(palette), levels = names(palette)),
                         stringsAsFactors = FALSE)
        list(ggplot2::geom_point(data = df,
                                 ggplot2::aes(x = -Inf, y = -Inf, colour = .data$label),
                                 alpha = 0, show.legend = TRUE, inherit.aes = FALSE),
             ggplot2::scale_colour_manual(
                     values = stats::setNames(unname(palette), names(palette)), name = name),
             ggplot2::guides(colour = ggplot2::guide_legend(
                     override.aes = list(alpha = 1, size = 3, shape = 15))))
}


# -----------------------------------------------------------------------------
# 6. The whole hierarchy
# -----------------------------------------------------------------------------

#' Every level of the typology side by side
#'
#' Supersedes \code{plot_region_local()}. Pass whichever levels you have; panels are
#' laid out in order. Region and local-type panels share one palette so the
#' nesting is visible. Meta-type and mosaic panels deliberately do not, because
#' they are not nested in the regions and a shared palette would imply they were.
#'
#' @param regions,local,meta,mosaic Fitted objects; any may be \code{NULL}.
#' @param geom An \code{sf} object, if not attached to the objects.
#' @param signatures Signature matrix, only needed if \code{mosaic} is supplied and you
#'   want its subtitle to carry the saturation.
#' @param ... Passed to each plotting call.
#' @return A \code{patchwork} object if \pkg{patchwork} is installed, otherwise a
#'   list of \code{ggplot} objects.
#' @export
plot_all_levels <- function(regions = NULL, local = NULL, meta = NULL, mosaic = NULL,
                           geom = NULL, signatures = NULL, ...) {
        .need("ggplot2")
        ps <- list()
        g  <- geom
        if (is.null(g)) {
                for (o in list(regions, local, meta, mosaic))
                        if (!is.null(o) && !is.null(attr(o, "pulse_geometry"))) {
                                g <- attr(o, "pulse_geometry"); break
                        }
        }
        if (is.null(g)) stop("No geometry found; pass `geom =`.", call. = FALSE)
        
        pal_r <- if (!is.null(regions)) pulse_palette(regions) else NULL
        if (!is.null(regions))
                ps$regions <- plot_regions(regions, "map", geom = g, palette = pal_r, ...)
        if (!is.null(local))
                ps$local <- plot_local_types(local, "map", geom = g, regions = regions,
                                             palette = if (!is.null(pal_r))
                                                     local_palette(local, pal_r) else NULL,
                                             confidence = if (is.null(regions)) "own" else "joint",
                                             ...)
        if (!is.null(meta))
                ps$meta <- plot_meta_types(meta, "map", geom = g, ...)
        if (!is.null(mosaic))
                ps$mosaic <- plot_mosaic_types(mosaic, "map", geom = g, signatures = signatures, ...)
        
        if (!length(ps)) stop("Nothing to plot.", call. = FALSE)
        if (requireNamespace("patchwork", quietly = TRUE)) Reduce(`+`, ps) else ps
}