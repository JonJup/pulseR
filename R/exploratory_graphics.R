# =============================================================================
# ECOTYPOLOGY EXPLORATORY GRAPHICS
# ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
# A cohesive set of exploratory-plot functions for the three nested levels of
# the river ecotypology pipeline:
#
#   Level 1  ECO-REGIONS    -> skater_con()        (WHERE  : spatial regions)
#   Level 2  RIVER TYPES     -> river_types()        (WHAT   : environment within region)
#   Level 3  LANDSCAPE TYPES -> get_mosaic_types() (HOW    : spatial co-occurrence)
#
# Each level gets its own family of plots; a final family visualises the three
# levels *together* (alluvial flow, concordance, side-by-side maps, fuzziness).
#
# Design conventions
#   * One shared minimal theme (theme_eco()) and one scalable, colour-blind
#     -safe qualitative palette (eco_palette()) so all panels read as a set.
#   * Sequential viridis for continuous fields; diverging RdBu for z-scores.
#   * "Uncertainty" is always normalised entropy of the (fuzzy) membership row,
#     so 0 = a crisp assignment and 1 = maximally ambiguous, comparable across
#     levels and across k.
#   * Plot functions RETURN a ggplot/patchwork object (they never print), so the
#     caller can compose, theme, or save them.
#   * Heavy/optional packages are checked at call time with a clear message.
#
# Dependencies: ggplot2 (required by all). sf (maps), patchwork (composites),
#   tidyr + dplyr (reshaping), scales, ggalluvial (alluvial only). The pipeline
#   helpers js_distance_fast() etc. are used when available but never required.
# =============================================================================


# ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
# 0. Shared utilities (theme, palette, label/uncertainty extractors)
# ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾

#' Stop unless every named package is installed
#' @keywords internal
.eco_require <- function(...) {
        pkgs <- c(...)
        missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
        if (length(missing)) {
                stop("This plot needs package(s) not installed: ",
                     paste(missing, collapse = ", "),
                     ".\n  install.packages(c(", paste0('"', missing, '"', collapse = ", "), "))",
                     call. = FALSE)
        }
        invisible(TRUE)
}

#' A clean shared theme for all ecotypology figures
#'
#' @param base_size Base font size.
#' @param map Logical. If TRUE, strips axes/grid for choropleths.
#' @export
theme_eco <- function(base_size = 12, map = FALSE) {
        .eco_require("ggplot2")
        th <- ggplot2::theme_minimal(base_size = base_size) +
                ggplot2::theme(
                        plot.title       = ggplot2::element_text(face = "bold", size = base_size * 1.25),
                        plot.subtitle    = ggplot2::element_text(colour = "grey35"),
                        plot.caption     = ggplot2::element_text(colour = "grey55", size = base_size * 0.8),
                        strip.text       = ggplot2::element_text(face = "bold"),
                        strip.background = ggplot2::element_rect(fill = "grey92", colour = NA),
                        panel.grid.minor = ggplot2::element_blank(),
                        legend.position  = "right",
                        legend.key.height = grid::unit(1, "lines")
                )
        if (map) {
                th <- th + ggplot2::theme(
                        axis.text   = ggplot2::element_blank(),
                        axis.title  = ggplot2::element_blank(),
                        axis.ticks  = ggplot2::element_blank(),
                        panel.grid  = ggplot2::element_blank()
                )
        }
        th
}

#' Scalable colour-blind-safe qualitative palette
#'
#' Uses the 8-colour Okabe-Ito set for k <= 8, then falls back to a perceptually
#' even HCL palette for larger k so arbitrarily many clusters stay distinct.
#'
#' @param n Number of colours.
#' @export
eco_palette <- function(n) {
        okabe <- c("#0072B2", "#E69F00", "#009E73", "#D55E00",
                   "#CC79A7", "#56B4E9", "#F0E442", "#999999")
        if (n <= length(okabe)) return(okabe[seq_len(n)])
        grDevices::hcl.colors(n, palette = "Dark 3")
}

#' Row-wise assignment uncertainty of a membership matrix
#'
#' @param U Numeric membership matrix (rows need not sum to 1; they are
#'   renormalised defensively).
#' @param method "entropy" (normalised Shannon entropy, default), "margin"
#'   (1 - gap between the top two memberships) or "maxprob" (1 - top membership).
#' @return Numeric vector in [0, 1]; higher = more ambiguous. Returns NA for
#'   all-zero / NA rows.
#' @keywords internal
.eco_uncertainty <- function(U, method = c("entropy", "margin", "maxprob")) {
        method <- match.arg(method)
        U <- as.matrix(U)
        k <- ncol(U)
        rs <- rowSums(U)
        ok <- is.finite(rs) & rs > 0
        P <- U
        P[ok, ] <- U[ok, , drop = FALSE] / rs[ok]
        out <- rep(NA_real_, nrow(U))
        if (method == "entropy") {
                H <- -rowSums(P * log(pmax(P, .Machine$double.eps)))
                out[ok] <- (H / log(k))[ok]
        } else {
                srt <- t(apply(P, 1, sort, decreasing = TRUE))
                if (method == "margin") {
                        out[ok] <- (1 - (srt[, 1] - srt[, 2]))[ok]
                } else {
                        out[ok] <- (1 - srt[, 1])[ok]
                }
        }
        pmin(pmax(out, 0), 1)
}

#' Hard label from a membership matrix (argmax)
#' @keywords internal
.eco_hard_from_U <- function(U) max.col(as.matrix(U), ties.method = "first")

#' Parse the wide river_types column names "Region{i}RiverType{k}"
#' @keywords internal
.eco_lt_class_table <- function(rt) {
        cn <- colnames(rt$localtypes)
        m  <- regmatches(cn, regexec("^Region([0-9]+)localType([0-9]+)$", cn))
        reg <- vapply(m, function(x) if (length(x) == 3) x[2] else NA_character_, character(1))
        typ <- vapply(m, function(x) if (length(x) == 3) x[3] else NA_character_, character(1))
        data.frame(col = cn,
                   region    = factor(paste0("Region ", reg), levels = unique(paste0("Region ", reg))),
                   type      = factor(typ),
                   stringsAsFactors = FALSE)
}

#' Dominant (region, local type) label per observation
#' @keywords internal
.eco_lt_hard <- function(lt) {
        M  <- as.matrix(lt$localtypes)
        cn <- colnames(M)
        factor(cn[max.col(M, ties.method = "first")], levels = cn)
}

#' Hard mosaic label, aligned to the clustered (valid) rows
#' @keywords internal
.eco_mosaic_hard <- function(mosaic) {
        cl <- mosaic$clusters
        if (!is.null(cl$pam_result)) return(as.integer(cl$pam_result$clustering)) # fuzzy
        if (!is.null(cl$clustering)) return(as.integer(cl$clustering))            # crisp pam
        if (!is.null(cl$memberships)) return(.eco_hard_from_U(cl$memberships))
        stop("Could not locate hard mosaic labels in `mosaic$clusters`.", call. = FALSE)
}

#' Landscape membership matrix (fuzzy U, or one-hot for crisp), valid rows
#' @keywords internal
.eco_mosaic_U <- function(mosaic) {
        cl <- mosaic$clusters
        if (!is.null(cl$memberships)) return(as.matrix(cl$memberships))
        lab <- .eco_mosaic_hard(mosaic); k <- max(lab)
        U <- matrix(0, length(lab), k); U[cbind(seq_along(lab), lab)] <- 1; U
}

#' Parse "ci:cj" signature column names into class-pair indices
#' @keywords internal
.eco_parse_pairs <- function(cn) {
        sp <- strsplit(cn, ":", fixed = TRUE)
        data.frame(
                a = vapply(sp, `[`, character(1), 1L),
                b = vapply(sp, function(x) if (length(x) >= 2) x[2] else x[1], character(1)),
                stringsAsFactors = FALSE
        )
}

#' Adjusted Rand Index between two label vectors (NA-tolerant)
#' @keywords internal
.eco_ari <- function(a, b) {
        ok <- !is.na(a) & !is.na(b)
        a <- a[ok]; b <- b[ok]
        if (length(a) < 2) return(NA_real_)
        tab <- table(a, b)
        sum_comb <- function(x) sum(choose(x, 2))
        ai <- rowSums(tab); bj <- colSums(tab); n <- sum(tab)
        idx  <- sum_comb(as.vector(tab))
        eidx <- sum_comb(ai) * sum_comb(bj) / choose(n, 2)
        midx <- (sum_comb(ai) + sum_comb(bj)) / 2
        if (midx - eidx == 0) return(0)
        (idx - eidx) / (midx - eidx)
}

#' Attach a vector to an sf layer by row order, with a length check
#' @keywords internal
.eco_attach <- function(polys, values, name) {
        .eco_require("sf")
        if (!inherits(polys, "sf"))
                stop("`polys` must be an sf object whose row order matches the result.", call. = FALSE)
        if (length(values) != nrow(polys))
                stop(sprintf("Length mismatch: %d values vs %d polygons. The sf layer must be ",
                             "row-aligned to the clustering input.", length(values), nrow(polys)),
                     call. = FALSE)
        polys[[name]] <- values
        polys
}

#' Internal: a discrete or continuous choropleth on an sf layer
#' @keywords internal
.eco_choropleth <- function(polys, col, discrete, title, sub, legend, border = "white", alpha = NULL) {
        .eco_require("ggplot2", "sf")
        if (is.null(alpha)){
                g <- ggplot2::ggplot(polys) +
                        ggplot2::geom_sf(ggplot2::aes(fill = .data[[col]]), colour = NA, linewidth = 0.05)
        } else {
                g <- ggplot2::ggplot(polys) +
                        ggplot2::geom_sf(ggplot2::aes(fill = .data[[col]], 
                                                      alpha = .data[[alpha]]),
                                         colour = NA,
                                         linewidth = 0.05) + 
                        labs(alpha = "Membership Degree")
        }

        if (discrete) {
                k <- uniqueN(polys[[col]])
                g <- g + ggplot2::scale_fill_manual(values = eco_palette(k), na.value = "grey85",
                                                    drop = FALSE, name = legend)
        } else {
                g <- g + ggplot2::scale_fill_viridis_c(option = "C", na.value = "grey85", name = legend)
        }
        g + ggplot2::labs(title = title, subtitle = sub) + theme_eco(map = TRUE)
}


# ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
# LEVEL 1: ECO-REGIONS  (output of skater_con())
# ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾

#' Map of fuzzy regions, or their assignment uncertainty
#'
#' The two faces of a fuzzy regionalisation: \code{fill = "region"} shows the
#' hard partition (the regions), \code{fill = "uncertainty"} shows where the
#' partition is *soft* -- the ecotones / transitional catchments that the hard
#' map hides.
#'
#' @param polys An \code{sf} polygon/line layer, row-aligned to the catchments
#'   that were clustered (same order as \code{eco$hard_clusters}).
#' @param eco The list returned by \code{skater_con()}.
#' @param fill "region" (default) or "uncertainty".
#' @param uncertainty_method Passed to the entropy/margin/maxprob helper.
#' @return A ggplot object.
#' @export
plot_region_map <- function(polys, eco, core = NULL, fill = c("region", "uncertainty", "both", "core"),
                               uncertainty_method = "entropy") {
        fill <- match.arg(fill)
        if (fill == "region") {
                lab <- factor(eco$hard_clusters)
                polys <- .eco_attach(polys, lab, ".eco")
                .eco_choropleth(polys, ".eco", discrete = TRUE,
                                title = "Eco-regions",
                                sub = sprintf("%d regions (hard partition) \u00b7 %d catchments",
                                              nlevels(lab), length(lab)),
                                legend = "Eco-region")
        } else if (fill == "uncertainty"){
                u <- .eco_uncertainty(eco$memberships, uncertainty_method)
                polys <- .eco_attach(polys, u, ".unc")
                .eco_choropleth(polys, ".unc", discrete = FALSE,
                                title = "Eco-region assignment uncertainty",
                                sub = "Normalised membership entropy \u00b7 bright = transitional (ecotone)",
                                legend = "Uncertainty")
        } else if (fill == "both"){
                u <- eco$memberships
                pd <- dplyr::bind_cols(polys, as.data.frame(u))
                class_cols <- dplyr::setdiff(names(pd), names(polys))
                
                # Calculate dominant region and membership to dominant region
                data_summary <- dplyr::rowwise(pd)
                data_summary <- dplyr::mutate(data_summary,
                                              Max_Memb = max(dplyr::c_across(dplyr::all_of(class_cols))),
                                              Dominant_Class = class_cols[which.max(dplyr::c_across(dplyr::all_of(class_cols)))])
                data_summary$Dominant_Class <- gsub("V", "", data_summary$Dominant_Class)
                data_summary <- dplyr::ungroup(data_summary)
                #polys <- .eco_attach(polys, u, ".unc")
                .eco_choropleth(data_summary, "Dominant_Class", discrete = TRUE,
                                title = "Continuous Region Membership",
                                alpha = "Max_Memb",
                                sub = "Degree of Membership for highest membership region \n bright = low membership = transitional (ecotone)",
                                legend = "Region")
        } else if (fill == "core"){
                if (is.null(core)){
                        stop("If method core is selected a core object must be provided.")
                }
                u <- eco$memberships
                pd <- dplyr::bind_cols(polys, as.data.frame(u))
                class_cols <- dplyr::setdiff(names(pd), names(polys))
                
                # Calculate dominant region and membership to dominant region
                data_summary <- dplyr::rowwise(pd)
                data_summary <- dplyr::mutate(data_summary,
                                              Max_Memb = max(dplyr::c_across(dplyr::all_of(class_cols))),
                                              Dominant_Class = class_cols[which.max(dplyr::c_across(dplyr::all_of(class_cols)))])
                data_summary <- dplyr::ungroup(data_summary)
                data_summary$Dominant_Class <- gsub("V", "", data_summary$Dominant_Class)
                # Extract IDs of all core regions and subset data_summary
                data_summary <- 
                        dplyr::mutate(
                                data_summary, 
                                Max_Memb = dplyr::if_else(
                                        ID %in% unlist(sapply(core, function(x) x$ID)),
                                        Max_Memb,
                                        0
                                ) 
                        )

                .eco_choropleth(data_summary, "Dominant_Class", discrete = TRUE,
                                title = "Core Regions",
                                alpha = "Max_Memb",
                                sub = "",
                                legend = "Region")
        }
}

#' Eco-region composition and crispness diagnostics
#'
#' Left: region "mass" -- hard catchment count vs soft membership mass (column
#' sums of U); a gap reveals regions that are real but never anyone's argmax.
#' Right: distribution of per-catchment max membership -- how crisp the
#' partition is overall.
#'
#' @param eco Result of \code{skater_con()}.
#' @return A patchwork of two ggplots.
#' @export
plot_region_membership <- function(eco) {
        .eco_require("ggplot2", "patchwork")
        U   <- as.matrix(eco$memberships)
        hard <- eco$hard_clusters
        k   <- ncol(U)
        
        comp <- data.frame(
                region = factor(rep(seq_len(k), 2)),
                kind   = rep(c("Hard count", "Soft mass"), each = k),
                value  = c(as.numeric(table(factor(hard, levels = seq_len(k)))), colSums(U))
        )
        p1 <- ggplot2::ggplot(comp, ggplot2::aes(.data$region, .data$value, fill = .data$kind)) +
                ggplot2::geom_col(position = "dodge", width = 0.7) +
                ggplot2::scale_fill_manual(values = c("Hard count" = "grey55", "Soft mass" = "#0072B2"),
                                           name = NULL) +
                ggplot2::labs(title = "Region size", x = "Eco-region", y = "Catchments") +
                theme_eco() + ggplot2::theme(legend.position = "top")
        
        crisp <- data.frame(maxU = apply(U, 1, max))
        p2 <- ggplot2::ggplot(crisp, ggplot2::aes(.data$maxU)) +
                ggplot2::geom_histogram(bins = 40, fill = "#009E73", colour = "white") +
                ggplot2::geom_vline(xintercept = 1 / k, linetype = 2, colour = "grey40") +
                ggplot2::annotate("text", x = 1 / k, y = Inf, label = "  uniform (1/k)",
                                  hjust = 0, vjust = 1.5, size = 3, colour = "grey40") +
                ggplot2::labs(title = "Partition crispness",
                              subtitle = "Per-catchment maximum membership",
                              x = "max membership", y = "Catchments") +
                theme_eco()
        
        patchwork::wrap_plots(p1, p2, widths = c(1.1, 1)) +
                patchwork::plot_annotation(title = "Eco-region membership structure")
}

#' Hyperparameter tuning surface for the region clustering
#'
#' Heatmap of the selection score over the (skater_regions x final_regions) grid
#' from \code{eco$tuning_log}, with the chosen cell marked. Reveals whether the
#' optimum is a sharp peak or a broad plateau.
#'
#' @param eco Result of a *tuned* \code{skater_con()} run (needs \code{tuning_log}).
#' @param score_col Column to map to fill. Default "score".
#' @return A ggplot object.
#' @export
plot_region_cvi_grid <- function(eco, score_col = "score") {
        .eco_require("ggplot2")
        tl <- eco$tuning_log
        if (is.null(tl)) stop("`eco$tuning_log` is NULL (single-configuration run; nothing to tune).",
                              call. = FALSE)
        need_cols <- c("skater_regions", "final_regions")
        if (!all(need_cols %in% names(tl))) {
                # tolerate alternative column names produced by the tuner
                alt <- c(sr = grep("skater", names(tl), value = TRUE)[1],
                         fr = grep("final",  names(tl), value = TRUE)[1])
                names(tl)[match(alt, names(tl))] <- need_cols
        }
        tl$skater_regions <- factor(tl$skater_regions)
        tl$final_regions  <- factor(tl$final_regions)
        best <- tl[which.max(tl[[score_col]]), , drop = FALSE]
        
        ggplot2::ggplot(tl, ggplot2::aes(.data$final_regions, .data$skater_regions,
                                         fill = .data[[score_col]])) +
                ggplot2::geom_tile(colour = "white", linewidth = 0.4) +
                ggplot2::geom_tile(data = best, fill = NA, colour = "black", linewidth = 1.1) +
                ggplot2::scale_fill_viridis_c(option = "D", name = score_col) +
                ggplot2::labs(title = "Eco-region tuning surface",
                              subtitle = sprintf("Optimum: skater_regions = %s, final_regions = %s",
                                                 best$skater_regions, best$final_regions),
                              x = "final regions (k)", y = "skater regions") +
                theme_eco()
}

#' Ensemble-size (n.rst) stability curve, ggplot version
#'
#' @param eco Result of \code{skater_con()} (uses \code{eco$n.rst_stability}).
#' @return A ggplot object, or NULL with a message if no stability run was done.
#' @export
plot_nrst_stability <- function(eco) {
        .eco_require("ggplot2")
        df <- eco$n.rst_stability
        if (is.null(df)) { message("No n.rst stability data (single n.rst was used)."); return(invisible(NULL)) }
        df <- df[is.finite(df$change), , drop = FALSE]
        knee <- df[df$is_knee, , drop = FALSE]
        ggplot2::ggplot(df, ggplot2::aes(.data$n.rst, .data$change)) +
                ggplot2::geom_line(colour = "grey50") +
                ggplot2::geom_point(size = 2, colour = "grey30") +
                ggplot2::geom_point(data = knee, size = 4, colour = "#D55E00") +
                ggplot2::geom_text(data = knee, ggplot2::aes(label = paste0("knee = ", .data$n.rst)),
                                   vjust = -1, colour = "#D55E00") +
                ggplot2::labs(title = "Ensemble-size stability",
                              subtitle = "Consensus-structure change as trees are added",
                              x = "n.rst (number of spanning trees)", y = "change") +
                theme_eco()
}

#' Combined region tuning dashboard (surface + stability)
#' @param eco Result of \code{skater_con()}.
#' @export
plot_region_tuning <- function(eco) {
        .eco_require("patchwork")
        g  <- tryCatch(plot_region_cvi_grid(eco), error = function(e) NULL)
        s  <- tryCatch(plot_nrst_stability(eco),     error = function(e) NULL)
        ps <- Filter(Negate(is.null), list(g, s))
        if (!length(ps)) stop("Nothing to plot: no tuning_log and no n.rst_stability.", call. = FALSE)
        patchwork::wrap_plots(ps, ncol = length(ps))
}

#' Ternary membership plot for the special case of exactly 3 regions
#'
#' A self-contained simplex projection (no ggtern dependency). Each catchment is
#' a point in the membership triangle; corner = a crisp region, centre = a
#' three-way mixture.
#'
#' @param eco Result of \code{skater_con()} with exactly 3 final regions.
#' @return A ggplot object.
#' @export
plot_membership_ternary <- function(eco) {
        .eco_require("ggplot2")
        U <- as.matrix(eco$memberships)
        if (ncol(U) != 3L) stop("Ternary plot requires exactly 3 regions (k = ", ncol(U), ").",
                                call. = FALSE)
        U <- U / pmax(rowSums(U), .Machine$double.eps)
        x <- U[, 2] + U[, 3] / 2
        y <- U[, 3] * sqrt(3) / 2
        df <- data.frame(x = x, y = y, region = factor(.eco_hard_from_U(U)))
        corners <- data.frame(x = c(0, 1, 0.5), y = c(0, 0, sqrt(3) / 2),
                              lab = c("Region 1", "Region 2", "Region 3"))
        ggplot2::ggplot() +
                ggplot2::geom_polygon(data = corners, ggplot2::aes(.data$x, .data$y),
                                      fill = NA, colour = "grey60") +
                ggplot2::geom_point(data = df, ggplot2::aes(.data$x, .data$y, colour = .data$region),
                                    alpha = 0.5, size = 1) +
                ggplot2::geom_text(data = corners, ggplot2::aes(.data$x, .data$y, label = .data$lab),
                                   vjust = c(2, 2, -1), fontface = "bold") +
                ggplot2::scale_colour_manual(values = eco_palette(3), name = "Argmax region") +
                ggplot2::coord_equal(clip = "off") +
                ggplot2::labs(title = "Eco-region membership simplex") +
                theme_eco(map = TRUE) + ggplot2::theme(legend.position = "right")
}


# ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
# LEVEL 2: RIVER TYPES  (output of get_local_types())
# ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾

#' BIC model-selection curves for every regional GMM
#'
#' For each region, the best BIC at each candidate number of components G
#' (envelope over covariance models), with the selected G highlighted. Shows
#' how many river types each region's mixture "wanted", and how decisive that
#' choice was (sharp peak vs flat ridge).
#'
#' @param rt Result of \code{get_local_types()}.
#' @return A ggplot object (faceted by region).
#' @export
plot_localtype_bic <- function(rt) {
        .eco_require("ggplot2")
        rows <- list()
        for (nm in names(rt$models)) {
                fit <- rt$models[[nm]]
                bic <- fit$BIC
                gvals <- suppressWarnings(as.integer(dimnames(bic)[[1]]))
                mat   <- matrix(as.numeric(bic), nrow = nrow(bic))
                env   <- apply(mat, 1, function(r) suppressWarnings(max(r, na.rm = TRUE)))
                env[!is.finite(env)] <- NA
                rows[[nm]] <- data.frame(region = nm, G = gvals, BIC = env,
                                         selected = gvals == fit$G)
        }
        df  <- do.call(rbind, rows)
        sel <- df[df$selected & is.finite(df$BIC), , drop = FALSE]
        ggplot2::ggplot(df, ggplot2::aes(.data$G, .data$BIC)) +
                ggplot2::geom_line(colour = "grey55") +
                ggplot2::geom_point(colour = "grey40", size = 1.6) +
                ggplot2::geom_point(data = sel, colour = "#D55E00", size = 3.5) +
                ggplot2::facet_wrap(~ region, scales = "free_y") +
                ggplot2::scale_x_continuous(breaks = scales::breaks_width(1)) +
                ggplot2::labs(title = "River-type model selection (per region)",
                              subtitle = "Best BIC at each G \u00b7 orange = selected number of river types",
                              x = "Number of river types (G)", y = "BIC") +
                theme_eco()
}

#' Environmental fingerprint of each river type
#'
#' Heatmap of the GMM component means (z-scored within each predictor across all
#' regions/types) -- rows are predictors, columns are river types, faceted by
#' region. This is the *interpretation* plot: it answers "what does Region 2's
#' river type 3 actually mean ecologically?".
#'
#' @param rt Result of \code{get_local_types()}.
#' @return A ggplot object.
#' @export
# plot_localtype_profiles <- function(rt, top = NULL) {
#         .eco_require("ggplot2")
#         rows <- list()
#         for (nm in names(rt$models)) {
#                 fit <- rt$models[[nm]]
#                 mu  <- fit$parameters$mean
#                 if (is.null(dim(mu))) {
#                         vn <- colnames(fit$data); if (is.null(vn)) vn <- "var"
#                         mu <- matrix(mu, nrow = 1, dimnames = list(vn[1], NULL))
#                 }
#                 vn <- rownames(mu); if (is.null(vn)) vn <- paste0("V", seq_len(nrow(mu)))
#                 for (g in seq_len(ncol(mu))) {
#                         rows[[length(rows) + 1L]] <- data.frame(
#                                 region = nm, type = factor(g), variable = vn, value = mu[, g]
#                         )
#                 }
#         }
#         df <- do.call(rbind, rows)
#         # z-score within each predictor so disparate units are comparable
#         df <- do.call(rbind, lapply(split(df, df$variable), function(d) {
#                 s <- stats::sd(d$value); d$z <- if (is.finite(s) && s > 0) (d$value - mean(d$value)) / s else 0
#                 d
#         }))
#         if (!is.null(top)){
#                 df2 <- dplyr::group_by(df, variable)
#                 df3 <- dplyr::summarize(df2, abs = sum(abs(z)), groups = "variable")
#                 df3 <- dplyr::arrange(df3, abs)
#                 df4 <- tail(df3, top)
#                 df2 <- dplyr::filter(df, variable %in% df4$variable)
#         } else {
#                 df2 <- df
#         }
#        
#         ggplot2::ggplot(df2, ggplot2::aes(.data$type, .data$variable, fill = .data$z)) +
#                 ggplot2::geom_tile(colour = "white", linewidth = 0.3) +
#                 ggplot2::facet_grid(~ region, scales = "free_x", space = "free_x") +
#                 ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "grey95", high = "#B2182B",
#                                               midpoint = 0, name = "z-score") +
#                 ggplot2::labs(title = "River-type environmental fingerprints",
#                               subtitle = "GMM component means, standardised within predictor",
#                               x = "River type", y = NULL) +
#                 theme_eco() + ggplot2::theme(panel.grid = ggplot2::element_blank())
# }
plot_localtype_profiles <- function(rt, top = NULL, cap = 3, order_rows = TRUE) {
        .eco_require("ggplot2")
        
        ## component means in long form, plus core data for a fixed per-predictor scale
        rows <- list(); pooled <- list()
        for (nm in names(rt$models)) {
                fit <- rt$models[[nm]]
                mu  <- fit$parameters$mean
                if (is.null(dim(mu))) {
                        vn <- colnames(fit$data); if (is.null(vn)) vn <- "var"
                        mu <- matrix(mu, nrow = 1L, dimnames = list(vn[1], NULL))
                }
                vn <- rownames(mu); if (is.null(vn)) vn <- paste0("V", seq_len(nrow(mu)))
                for (g in seq_len(ncol(mu))) {
                        rows[[length(rows) + 1L]] <- data.frame(
                                region = nm, type = factor(g), variable = vn, value = mu[, g],
                                stringsAsFactors = FALSE
                        )
                }
                if (!is.null(fit$data)) pooled[[nm]] <- as.data.frame(fit$data)
        }
        df <- do.call(rbind, rows)
        
        ## fixed per-predictor scale (natural variability), kept global so colours stay comparable
        if (length(pooled)) {
                pooled  <- do.call(rbind, pooled)
                scale_v <- vapply(unique(df$variable),
                                  function(v) stats::sd(pooled[[v]], na.rm = TRUE), numeric(1))
                names(scale_v) <- unique(df$variable)
        } else {
                warning("model$data unavailable; falling back to spread of means.")
                scale_v <- vapply(split(df$value, df$variable), stats::sd, numeric(1))
        }
        
        ## centre each predictor WITHIN its region: deviation of each type from its region's average
        s <- scale_v[df$variable]; s[!is.finite(s) | s <= 0] <- NA
        df$z <- (df$value - ave(df$value, df$region, df$variable, FUN = mean)) / s
        df$z[is.na(df$z)] <- 0
        
        ## cap so a stray component cannot blow out the palette (or the ranking)
        if (is.finite(cap)) df$z <- pmax(pmin(df$z, cap), -cap)
        
        ## discriminability = within-region between-type variance, averaged over regions
        rv    <- aggregate(z ~ region + variable, df, stats::var)
        score <- tapply(rv$z, rv$variable, mean, na.rm = TRUE)
        score[!is.finite(score)] <- 0
        
        if (!is.null(top)) {
                keep <- names(sort(score, decreasing = TRUE))[seq_len(min(top, length(score)))]
                df   <- df[df$variable %in% keep, , drop = FALSE]
        }
        if (order_rows) {
                df$variable <- factor(df$variable,
                                      levels = names(sort(score[unique(as.character(df$variable))])))
        }
        
        lim <- if (is.finite(cap)) c(-cap, cap) else NULL
        ggplot2::ggplot(df, ggplot2::aes(.data$type, .data$variable, fill = .data$z)) +
                ggplot2::geom_tile(colour = "white", linewidth = 0.3) +
                ggplot2::facet_grid(~ region, scales = "free_x", space = "free_x") +
                ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "grey95", high = "#B2182B",
                                              midpoint = 0, limits = lim, name = "within-region z") +
                ggplot2::labs(title = "River-type environmental fingerprints",
                              subtitle = "GMM component means, centred within ecoregion, scaled by predictor variability",
                              x = "River type", y = NULL) +
                theme_eco() + ggplot2::theme(panel.grid = ggplot2::element_blank())
}
#' River-type composition by region
#'
#' How much of the dataset each river type captures. \code{weight = "soft"}
#' (default) sums the membership-weighted probabilities (uses the fuzzy
#' information); \code{"hard"} counts argmax assignments.
#'
#' @param rt Result of \code{get_local_types()}.
#' @param weight "soft" or "hard".
#' @return A ggplot object.
#' @export
plot_localtype_composition <- function(rt, weight = c("soft", "hard")) {
        .eco_require("ggplot2")
        weight <- match.arg(weight)
        tab <- .eco_lt_class_table(rt)
        if (weight == "soft") {
                tab$value <- colSums(as.matrix(rt$localtypes))
                ylab <- "Soft membership mass"
        } else {
                hard <- table(.eco_rt_hard(rt))
                tab$value <- as.numeric(hard[tab$col])
                ylab <- "Catchments (argmax)"
        }
        ggplot2::ggplot(tab, ggplot2::aes(.data$region, .data$value, fill = .data$type)) +
                ggplot2::geom_col(width = 0.75, colour = "white", linewidth = 0.2) +
                ggplot2::scale_fill_manual(values = eco_palette(nlevels(tab$type)), name = "River type") +
                ggplot2::labs(title = "River-type composition",
                              subtitle = sprintf("Weighting: %s", weight),
                              x = NULL, y = ylab) +
                theme_eco() + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))
}

#' Map of dominant river type, or river-type assignment uncertainty
#'
#' @param polys An \code{sf} layer row-aligned to \code{rt$localtypes}.
#' @param rt Result of \code{get_local_types()}.
#' @param fill "type" (default) or "uncertainty".
#' @return A ggplot object.
#' @export
plot_localtype_map <- function(polys, lt, fill = c("type", "uncertainty", "both"),
                               uncertainty_method = "entropy") {
        fill <- match.arg(fill)
        if (fill == "type") {
                lab <- .eco_lt_hard(lt)
                lab <- droplevels(lab)
                polys <- .eco_attach(polys, lab, ".lt")
                .eco_choropleth(polys, ".lt", discrete = TRUE,
                                title = "Dominant river type",
                                sub = sprintf("%d realised river-type classes", nlevels(lab)),
                                legend = "River type")
        } else if (fill == "uncertainty"){
                u <- .eco_uncertainty(lt$localtypes, uncertainty_method)
                polys <- .eco_attach(polys, u, ".unc")
                .eco_choropleth(polys, ".unc", discrete = FALSE,
                                title = "River-type assignment uncertainty",
                                sub = "Normalised entropy of weighted river-type probabilities",
                                legend = "Uncertainty")
        } else if (fill == "both"){
                u <- lt$localtypes
                pd <- dplyr::bind_cols(polys, as.data.frame(u))
                class_cols <- dplyr::setdiff(names(pd), names(polys))
                
                # Calculate dominant region and membership to dominant region
                data_summary <- dplyr::rowwise(pd)
                data_summary <- dplyr::mutate(data_summary,
                                              Max_Memb = max(dplyr::c_across(dplyr::all_of(class_cols))),
                                              Dominant_Class = class_cols[which.max(dplyr::c_across(dplyr::all_of(class_cols)))])
                data_summary$Dominant_Class <- gsub("V", "", data_summary$Dominant_Class)
                data_summary <- dplyr::ungroup(data_summary)
                #polys <- .eco_attach(polys, u, ".unc")
                .eco_choropleth(data_summary, "Dominant_Class", discrete = TRUE,
                                title = "Continuous Local Type Membership",
                                alpha = "Max_Memb",
                                sub = "Degree of Membership for highest membership region \n bright = low membership = transitional (ecotone)",
                                legend = "Local type")
        }
}


# ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
# LEVEL 3: MOSAIC TYPES  (output of get_mosaic_types())
# ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾

#' Map of mosaic types (or their uncertainty)
#'
#' Note: \code{get_mosaic_types()} clusters only the complete-case rows, so the
#' \code{mosaic$valid} vector is used to scatter labels back onto the full polygon
#' set; isolated/NA-signature polygons appear in grey.
#'
#' @param polys An \code{sf} layer row-aligned to the *signature* matrix (the
#'   full polygon set; same order as \code{mosaic$valid}).
#' @param mosaic Result of \code{get_mosaic_types()}.
#' @param fill "type" (default) or "uncertainty".
#' @return A ggplot object.
#' @export
plot_mosaic_map <- function(polys, mosaic, fill = c("type", "uncertainty", "both"),
                               uncertainty_method = "entropy") {
        fill  <- match.arg(fill)
        valid <- mosaic$valid
        if (fill == "type") {
                lab <- rep(NA_integer_, length(valid))
                lab[valid] <- .eco_mosaic_hard(mosaic)
                lab <- factor(lab)
                polys <- .eco_attach(polys, lab, ".ls")
                .eco_choropleth(polys, ".ls", discrete = TRUE,
                                title = "Landscape types",
                                sub = sprintf("k = %d co-occurrence clusters \u00b7 grey = isolated polygons",
                                              mosaic$best_k),
                                legend = "Landscape")
        } else if (fill == "uncertainty"){
                u <- rep(NA_real_, length(valid))
                u[valid] <- .eco_uncertainty(.eco_mosaic_U(mosaic), uncertainty_method)
                polys <- .eco_attach(polys, u, ".unc")
                .eco_choropleth(polys, ".unc", discrete = FALSE,
                                title = "Landscape-type uncertainty",
                                sub = "Normalised membership entropy",
                                legend = "Uncertainty")
        } else if (fill == "both"){
                u <- mosaic$clusters$memberships
                pd <- dplyr::bind_cols(polys, as.data.frame(u))
                class_cols <- dplyr::setdiff(names(pd), names(polys))
                
                # Calculate dominant region and membership to dominant region
                data_summary <- dplyr::rowwise(pd)
                data_summary <- dplyr::mutate(data_summary,
                                              Max_Memb = max(dplyr::c_across(dplyr::all_of(class_cols))),
                                              Dominant_Class = class_cols[which.max(dplyr::c_across(dplyr::all_of(class_cols)))])
                data_summary$Dominant_Class <- gsub("V", "", data_summary$Dominant_Class)
                data_summary <- dplyr::ungroup(data_summary)
                #polys <- .eco_attach(polys, u, ".unc")
                .eco_choropleth(data_summary, "Dominant_Class", discrete = TRUE,
                                title = "Continuous Mosaic Type Membership",
                                alpha = "Max_Memb",
                                sub = "Degree of Membership for highest membership mosaic type \n bright = low membership = transitional (ecotone)",
                                legend = "Mosaic type")
        }
}

#' Co-occurrence fingerprints of the mosaic types
#'
#' Each mosiac type is defined by a class co-occurrence distribution. This
#' rebuilds, per cluster, the symmetric class x class mean co-occurrence matrix
#' from the signature columns ("ci:cj") and shows them as small-multiple
#' heatmaps -- the visual "signature" that distinguishes each mosaic type.
#'
#' @param mosaic Result of \code{get_mosaic_types()}.
#' @param sigs The signature matrix passed to \code{get_mosaic_types()}
#'   (row-aligned to \code{mosaic$valid}).
#' @return A ggplot object.
#' @export
plot_mosaic_signatures <- function(mosaic, sigs) {
        .eco_require("ggplot2")
        S   <- as.matrix(sigs)[mosaic$valid, , drop = FALSE]
        lab <- .eco_mosaic_hard(mosaic)
        pr  <- .eco_parse_pairs(colnames(S))
        classes <- sort(unique(c(pr$a, pr$b)))
        ulab <- sort(unique(lab))
        
        rows <- list()
        for (g in ulab) {
                m <- colMeans(S[lab == g, , drop = FALSE], na.rm = TRUE)
                for (j in seq_along(m)) {
                        rows[[length(rows) + 1L]] <- data.frame(
                                cluster = factor(paste0("Landscape ", g)),
                                a = pr$a[j], b = pr$b[j], value = m[j]
                        )
                        if (pr$a[j] != pr$b[j]) {  # mirror to fill the lower triangle
                                rows[[length(rows) + 1L]] <- data.frame(
                                        cluster = factor(paste0("Landscape ", g)),
                                        a = pr$b[j], b = pr$a[j], value = m[j]
                                )
                        }
                }
        }
        df <- do.call(rbind, rows)
        df$a <- factor(df$a, levels = classes)
        df$b <- factor(df$b, levels = rev(classes))
        
        ggplot2::ggplot(df, ggplot2::aes(.data$a, .data$b, fill = .data$value)) +
                ggplot2::geom_tile(colour = "white", linewidth = 0.3) +
                ggplot2::facet_wrap(~ cluster) +
                ggplot2::scale_fill_viridis_c(option = "G", name = "mean\nco-occur.") +
                ggplot2::coord_equal() +
                ggplot2::labs(title = "Landscape-type co-occurrence fingerprints",
                              subtitle = "Mean class \u00d7 class co-occurrence per cluster",
                              x = "class", y = "class") +
                theme_eco() + ggplot2::theme(panel.grid = ggplot2::element_blank(),
                                             axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5))
}

#' Cluster-validity profile for mosaic-type selection
#'
#' The fuzzy/crisp CVIs across candidate k from \code{mosaic$validity}, faceted by
#' index, with the selected k marked. Facet labels flag the optimisation
#' direction (min vs max) so the reader can see *why* best_k won.
#'
#' @param mosaic Result of \code{get_mosaic_types()}.
#' @param metrics Which validity columns to show.
#' @return A ggplot object.
#' @export
plot_mosaic_validity <- function(mosaic,
                                    metrics = c("XB_star", "SIL_F", "STAB",
                                                "SILH_HARD", "MPC")) {
        .eco_require("ggplot2")
        v <- mosaic$validity
        metrics <- intersect(metrics, names(v))
        dir <- c(XB_star = "min", SIL_F = "max", STAB = "max",
                 SILH_HARD = "max", MPC = "max", PE = "min", PC = "max")
        rows <- list()
        for (mt in metrics) {
                if (all(is.na(v[[mt]]))) next
                arrow <- if (!is.na(dir[mt]) && dir[mt] == "min") " (min)" else " (max)"
                rows[[mt]] <- data.frame(k = v$k, value = v[[mt]],
                                         metric = paste0(mt, arrow))
        }
        df <- do.call(rbind, rows)
        ggplot2::ggplot(df, ggplot2::aes(.data$k, .data$value)) +
                ggplot2::geom_vline(xintercept = mosaic$best_k, linetype = 2, colour = "#D55E00") +
                ggplot2::geom_line(colour = "grey55") +
                ggplot2::geom_point(size = 2, colour = "grey25") +
                ggplot2::facet_wrap(~ metric, scales = "free_y") +
                ggplot2::scale_x_continuous(breaks = scales::breaks_width(1)) +
                ggplot2::labs(title = "Mosaic type validity profile",
                              subtitle = sprintf("Selected k = %d (orange)", mosaic$best_k),
                              x = "Number of mosaic types (k)", y = NULL) +
                theme_eco()
}

#' Ordination of mosaic signatures (MDS of Jensen-Shannon distances)
#'
#' Projects the signature space to 2-D so cluster separation is visible. Uses the
#' pipeline's \code{js_distance_fast()} when it is on the search path, otherwise
#' falls back to Euclidean distance on the signatures.
#'
#' @param sigs Signature matrix (row-aligned to \code{mosaic$valid}).
#' @param mosaic Result of \code{get_mosaic_types()}.
#' @param ellipses Draw 95\% normal ellipses per cluster? Default TRUE.
#' @return A ggplot object.
#' @export
plot_mosaic_ordination <- function(sigs, mosaic, ellipses = TRUE) {
        .eco_require("ggplot2")
        S <- as.matrix(sigs)[mosaic$valid, , drop = FALSE]
        d <- if (exists("js_distance_fast", mode = "function")) {
                get("js_distance_fast")(S)
        } else {
                eps <- .Machine$double.eps; P <- (S + eps) / rowSums(S + eps)
                stats::dist(P)
        }
        mds <- stats::cmdscale(d, k = 2)
        df  <- data.frame(MDS1 = mds[, 1], MDS2 = mds[, 2],
                          cluster = factor(.eco_mosaic_hard(mosaic)))
        g <- ggplot2::ggplot(df, ggplot2::aes(.data$MDS1, .data$MDS2, colour = .data$cluster)) +
                ggplot2::geom_point(alpha = 0.6, size = 1.4)
        if (ellipses) g <- g + ggplot2::stat_ellipse(level = 0.95, linewidth = 0.7)
        g + ggplot2::scale_colour_manual(values = eco_palette(nlevels(df$cluster)),
                                         name = "Landscape") +
                ggplot2::labs(title = "Landscape signature ordination",
                              subtitle = "Classical MDS of Jensen-Shannon distances",
                              x = "MDS axis 1", y = "MDS axis 2") +
                theme_eco()
}

#' Joint scale (k-order) x cluster-number sweep heatmap
#'
#' @param sweep_df Output of \code{sweep_k_order()}.
#' @return A ggplot object.
#' @export
plot_mosaic_ksweep <- function(sweep_df) {
        .eco_require("ggplot2")
        df <- sweep_df
        df$k_order    <- factor(df$k_order)
        df$n_clusters <- factor(df$n_clusters)
        best_per <- df[df$n_clusters == df$best, , drop = FALSE]  # selected k per order
        best_per <- df[as.integer(as.character(df$n_clusters)) == df$best, , drop = FALSE]
        ggplot2::ggplot(df, ggplot2::aes(.data$n_clusters, .data$k_order, fill = .data$silhouette)) +
                ggplot2::geom_tile(colour = "white", linewidth = 0.4) +
                ggplot2::geom_point(data = best_per, colour = "white", size = 2) +
                ggplot2::scale_fill_viridis_c(option = "D", name = "silhouette") +
                ggplot2::labs(title = "Neighbourhood-order \u00d7 cluster-number sweep",
                              subtitle = "Hard silhouette width \u00b7 dots = selected k per order",
                              x = "Number of mosaic types (k)", y = "Neighbourhood order") +
                theme_eco()
}


# ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
# CROSS-LEVEL: the three levels together
# ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾

#' Assemble a row-aligned typology table across the three levels
#'
#' All three results must describe the same features in the same row order
#' (region and river-type span the full set; mosaic labels are scattered
#' through \code{mosaic$valid}). Produces one tidy data frame with hard labels and
#' per-level uncertainties -- the substrate for every cross-level plot.
#'
#' @param region,local,mosaic Results of the three pipeline stages (any may be NULL).
#' @return A data.frame with columns \code{region}, \code{river_type},
#'   \code{mosaic}, and \code{u_eco}, \code{u_rt}, \code{u_mosaic}.
#' @export
assemble_typology <- function(region = NULL, local = NULL, mosaic = NULL) {
        lens <- c(region  = if (!is.null(region))  length(region$hard_clusters)      else NA,
                  local   = if (!is.null(local))   nrow(local$localtypes)            else NA,
                  mosaic = if (!is.null(mosaic)) length(mosaic$valid)            else NA)
        N <- unique(stats::na.omit(lens))
        if (length(N) > 1) stop("Inputs are not row-aligned (lengths: ",
                                paste(names(lens), lens, sep = "=", collapse = ", "), ").",
                                call. = FALSE)
        N <- N[1]
        
        out <- data.frame(.row = seq_len(N))
        if (!is.null(region)) {
                out$region <- factor(region$hard_clusters)
                out$u_region     <- .eco_uncertainty(region$memberships)
        }
        if (!is.null(local)) {
                out$local <- droplevels(.eco_rt_hard(local))
                out$u_lt       <- .eco_uncertainty(local$localtypes)
        }
        if (!is.null(mosaic)) {
                ls <- rep(NA_integer_, N); ls[mosaic$valid] <- .eco_mosaic_hard(mosaic)
                out$mosaic <- factor(ls)
                uu <- rep(NA_real_, N); uu[mosaic$valid] <- .eco_uncertainty(.eco_mosaic_U(mosaic))
                out$u_mosaic <- uu
        }
        out$.row <- NULL
        out
}

#' Make alluvial stratum labels unique & readable per axis
#' @keywords internal
.eco_axis_relabel <- function(x, axis) {
        x <- as.character(x)
        if (axis == "local") {
                return(sub(pattern = "^Region([0-9]+)localType([0-9]+)$", replacement = "R\\1\u00b7T\\2", x = x))
        }
        pref <- switch(axis, region = "Region ", mosaic = "Mosaic ",
                       paste0(abbreviate(axis, 4L), ": "))
        paste0(pref, x)
}

#' Alluvial flow across the hierarchy: region -> local type -> mosaic
#'
#' The "money" cross-level figure: stream width = number of catchments, so you
#' can read how each region splits into river types and how those reassemble
#' into mosaic types (i.e. how aligned or crossed the three partitions are).
#'
#' @param typ Output of \code{assemble_typology()}.
#' @param axes Which label columns to chain (in order).
#' @return A ggplot object.
#' @export
plot_hierarchy_alluvial <- function(typ,
                                    axes = c("region", "local", "mosaic")) {
        .eco_require("ggplot2", "ggalluvial")
        axes <- intersect(axes, names(typ))
        if (length(axes) < 2) stop("Need at least two label columns to draw flows.", call. = FALSE)
        d <- typ[stats::complete.cases(typ[, axes, drop = FALSE]), axes, drop = FALSE]
        # Disambiguate identical codes across axes (e.g. region "1" vs mosaic "1")
        # so flows join correctly and stratum labels stay self-describing.
        for (ax in axes) d[[ax]] <- factor(.eco_axis_relabel(d[[ax]], ax))
        d$freq <- 1L
        agg <- stats::aggregate(freq ~ ., data = d, FUN = sum)
        
        mapping <- ggplot2::aes(y = .data$freq)
        for (i in seq_along(axes)) mapping[[paste0("axis", i)]] <- as.name(axes[i])
        
        ggplot2::ggplot(agg, mapping) +
                ggalluvial::geom_alluvium(ggplot2::aes(fill = .data[[axes[1]]]), alpha = 0.6) +
                ggalluvial::geom_stratum(width = 0.28, fill = "grey92", colour = "grey50") +
                ggplot2::geom_text(stat = ggalluvial::StatStratum,
                                   ggplot2::aes(label = ggplot2::after_stat(stratum)), size = 3) +
                ggplot2::scale_x_discrete(limits = axes, expand = c(0.08, 0.08)) +
                ggplot2::scale_fill_manual(values = eco_palette(nlevels(agg[[axes[1]]])),
                                           guide = "none") +
                ggplot2::labs(title = "Typology hierarchy",
                              subtitle = paste(axes, collapse = "  \u2192  "),
                              x = NULL, y = "Catchments") +
                theme_eco()
}

#' Concordance heatmap + ARI between two typology levels
#'
#' Row-normalised contingency between two label sets: each cell is the share of
#' an \code{x}-class that falls in a \code{y}-class. A near-diagonal/blocky map =
#' the partitions agree; a smeared map = they capture orthogonal structure. The
#' Adjusted Rand Index quantifies it in the subtitle.
#'
#' @param typ Output of \code{assemble_typology()}.
#' @param x,y Label columns to compare.
#' @return A ggplot object.
#' @export
plot_hierarchy_concordance <- function(typ, x = "region", y = "mosaic") {
        .eco_require("ggplot2")
        for (cc in c(x, y)) if (!cc %in% names(typ))
                stop("Column '", cc, "' not in the typology table.", call. = FALSE)
        ok <- !is.na(typ[[x]]) & !is.na(typ[[y]])
        tab <- table(typ[[x]][ok], typ[[y]][ok])
        if (x == "local"){
                rownames(tab) <- .eco_axis_relabel(rownames(tab), axis = "local")
        } else if (y == "local"){
                colnames(tab) <- .eco_axis_relabel(colnames(tab), axis = "local")
        }
        prop <- sweep(tab, 1, pmax(rowSums(tab), 1), "/")
        df <- as.data.frame(prop)
        names(df) <- c("xcl", "ycl", "share")
        ari <- .eco_ari(typ[[x]], typ[[y]])
        
        ggplot2::ggplot(df, ggplot2::aes(.data$ycl, .data$xcl, fill = .data$share)) +
                ggplot2::geom_tile(colour = "white", linewidth = 0.3) +
                ggplot2::scale_fill_viridis_c(option = "B", limits = c(0, 1),
                                              name = "row share") +
                ggplot2::labs(title = sprintf("Concordance: %s vs %s", x, y),
                              subtitle = sprintf("Adjusted Rand Index = %.3f", ari),
                              x = y, y = x) +
                theme_eco() + ggplot2::theme(panel.grid = ggplot2::element_blank(),
                                             axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))
}

#' Side-by-side maps of all three levels
#'
#' One geometry, three partitions, shared layout -- the fastest way to *see*
#' whether the levels describe nested or crossing spatial structure.
#'
#' @param polys An \code{sf} layer row-aligned to all three results.
#' @param eco,rt,mosaic The three pipeline results (any may be NULL).
#' @return A patchwork object.
#' @export
plot_hierarchy_maps <- function(polys, eco = NULL, rt = NULL, mosaic = NULL) {
        .eco_require("patchwork")
        ps <- list()
        if (!is.null(eco))  ps$Ecoregion <- plot_region_map(polys, eco) +
                ggplot2::labs(title = "1 \u00b7 Eco-region", subtitle = NULL)
        if (!is.null(rt))   ps$Rivertype <- plot_localtype_map(polys, rt) +
                ggplot2::labs(title = "2 \u00b7 River type", subtitle = NULL)
        if (!is.null(mosaic)) ps$Landscape <- plot_mosaic_map(polys, mosaic) +
                ggplot2::labs(title = "3 \u00b7 Landscape type", subtitle = NULL)
        if (!length(ps)) stop("Supply at least one of `eco`, `rt`, `mosaic`.", call. = FALSE)
        patchwork::wrap_plots(ps, nrow = 1) +
                patchwork::plot_annotation(
                        title = "The three ecotypology levels in space",
                        theme = theme_eco())
}

#' Compare assignment uncertainty across the three levels
#'
#' Which level is crispest? Violins (with quartile boxes) of normalised
#' membership entropy per level -- a higher, fatter violin marks a level whose
#' boundaries are inherently fuzzier.
#'
#' @param typ Output of \code{assemble_typology()}.
#' @return A ggplot object.
#' @export
plot_hierarchy_uncertainty <- function(typ) {
        .eco_require("ggplot2")
        cols <- intersect(c("u_region", "u_lt", "u_mosaic"), names(typ))
        if (!length(cols)) stop("No uncertainty columns found in the typology table.", call. = FALSE)
        labmap <- c(u_region = "Region", u_lt = "Local type", u_mosaic = "Mosaic type")
        rows <- lapply(cols, function(cc)
                data.frame(level = unname(labmap[cc]), u = typ[[cc]]))
        df <- do.call(rbind, rows)
        df <- df[is.finite(df$u), , drop = FALSE]
        df$level <- factor(df$level, levels = labmap[cols])
        
        ggplot2::ggplot(df, ggplot2::aes(.data$level, .data$u, fill = .data$level)) +
                ggplot2::geom_violin(colour = NA, alpha = 0.5, scale = "width") +
                ggplot2::geom_boxplot(width = 0.16, outlier.shape = NA, fill = "white") +
                ggplot2::scale_fill_manual(values = eco_palette(nlevels(df$level)), guide = "none") +
                ggplot2::labs(title = "Assignment uncertainty across levels",
                              subtitle = "Normalised membership entropy (0 = crisp, 1 = ambiguous)",
                              x = NULL, y = "Uncertainty") +
                theme_eco()
}


# ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
# USAGE  (sketch; assumes `polys` is an sf layer row-aligned to each result)
# ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
# eco  <- skater_con(graph, skater_regions = c(20,40), final_regions = 3:6, n.rst = 30)
# rt   <- river_types(all_data, core_regions, membership_id = eco$memberships)
# sigs <- compute_all_signatures(P, graph, k = 2, coords = coords)
# land <- cluster_landscapes(sigs, n_clusters = 2:10)
#
# # Level 1
# plot_region_map(polys, eco)                      # regions
# plot_region_map(polys, eco, fill = "uncertainty")# ecotones
# plot_region_membership(eco)
# plot_region_tuning(eco)
#
# # Level 2
# plot_localtype_bic(rt)
# plot_localtype_profiles(rt)
# plot_localtype_composition(rt)
# plot_localtype_map(polys, rt)
#
# # Level 3
# plot_mosaic_map(polys, land)
# plot_mosaic_signatures(land, sigs)
# plot_mosaic_validity(land)
# plot_mosaic_ordination(sigs, land)
#
# # Together
# typ <- assemble_typology(eco, rt, land)
# plot_hierarchy_alluvial(typ)
# plot_hierarchy_concordance(typ, "region", "landscape")
# plot_hierarchy_maps(polys, eco, rt, land)
# plot_hierarchy_uncertainty(typ)
# =============================================================================