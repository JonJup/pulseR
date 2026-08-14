# ECOTYPOLOGY EXPLORATORY GRAPHICS

# A cohesive set of exploratory-plot functions for the three nested levels of
# the ecotypology pipeline:
#
#   Level 1  REGIONS      -> get_regions()      (WHERE: spatial regions)
#   Level 2  LOCAL TYPES  -> get_local_types()  (WHAT: environment within region)
#   Level 3  MOSAIC TYPES -> get_mosaic_types() (HOW: spatial co-occurrence)
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

# ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
# 0. Shared utilities (theme, palette, label/uncertainty extractors)
# ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾

#' Stop unless every named package is installed
#'
#' @param ... Character package names to check (one or more).
#' @return Invisibly \code{TRUE} if all packages are installed; otherwise stops
#'   with an informative \code{install.packages()} hint.
#' @keywords internal
.region_require <- function(...) {
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
#' @return A \code{ggplot2} theme object to add to a plot.
#' @export
theme_eco <- function(base_size = 12, map = FALSE) {
        .region_require("ggplot2")
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
#' @return A character vector of \code{n} hex colour codes.
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
#' @return Numeric vector in `[0, 1]`; higher = more ambiguous. Returns NA for
#'   all-zero / NA rows.
#' @keywords internal
.region_uncertainty <- function(U, method = c("entropy", "margin", "maxprob")) {
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
#'
#' @param U Numeric membership matrix (observations x classes).
#' @return Integer vector of length \code{nrow(U)} giving the argmax column index
#'   per row (ties broken by first).
#' @keywords internal
.region_hard_from_U <- function(U) max.col(as.matrix(U), ties.method = "first")

#' Parse wide local_types column names
#'
#' @param rt A list with a \code{localtypes} matrix whose column names follow the
#'   \code{R<region>·T<type>} convention.
#' @return A data.frame with columns \code{col} (original name), \code{region}
#'   (factor) and \code{type} (factor), one row per local-type column.
#' @keywords internal
.region_lt_class_table <- function(rt) {
        cn <- colnames(rt$localtypes)
        m  <- regmatches(cn, regexec("^R([0-9]+)\u00b7T([0-9]+)$", cn))
        reg <- vapply(m, function(x) if (length(x) == 3) x[2] else NA_character_, character(1))
        typ <- vapply(m, function(x) if (length(x) == 3) x[3] else NA_character_, character(1))
        data.frame(col = cn,
                   region    = factor(paste0("Region ", reg), levels = unique(paste0("Region ", reg))),
                   type      = factor(typ),
                   stringsAsFactors = FALSE)
}

#' Dominant (region, local type) label per observation
#'
#' @param lt Result of \code{get_local_types()} (uses its \code{localtypes} matrix).
#' @return A factor of length \code{nrow(lt$localtypes)} giving the argmax
#'   (region, local type) label per observation, with levels in column order.
#' @keywords internal
.region_lt_hard <- function(lt) {
        M  <- as.matrix(lt$localtypes)
        cn <- colnames(M)
        cn <- sub(pattern = "^Region([0-9]+)localType([0-9]+)$", replacement = "R\\1\u00b7T\\2", x = cn)
        factor(cn[max.col(M, ties.method = "first")], levels = cn)
}

#' Hard mosaic label, aligned to the clustered (valid) rows
#'
#' @param mosaic Result of \code{get_mosaic_types()}.
#' @return Integer vector of hard cluster labels for the valid (clustered) rows,
#'   sourced from \code{pam_result}, \code{clustering} or \code{memberships} as
#'   available. Stops if none can be found.
#' @keywords internal
.region_mosaic_hard <- function(mosaic) {
        cl <- mosaic$clusters
        if (!is.null(cl$pam_result)) return(as.integer(cl$pam_result$clustering)) # fuzzy
        if (!is.null(cl$clustering)) return(as.integer(cl$clustering))            # crisp pam
        if (!is.null(cl$memberships)) return(.region_hard_from_U(cl$memberships))
        stop("Could not locate hard mosaic labels in `mosaic$clusters`.", call. = FALSE)
}

#' mosaic type membership matrix (fuzzy U, or one-hot for crisp), valid rows
#'
#' @param mosaic Result of \code{get_mosaic_types()}.
#' @return Numeric membership matrix (valid rows x clusters): the fuzzy
#'   \code{memberships} when present, otherwise a one-hot matrix built from the
#'   hard labels.
#' @keywords internal
.region_mosaic_U <- function(mosaic) {
        cl <- mosaic$clusters
        if (!is.null(cl$memberships)) return(as.matrix(cl$memberships))
        lab <- .region_mosaic_hard(mosaic); k <- max(lab)
        U <- matrix(0, length(lab), k); U[cbind(seq_along(lab), lab)] <- 1; U
}

#' Parse "ci:cj" signature column names into class-pair indices
#'
#' @param cn Character vector of signature column names of the form \code{"ci:cj"}
#'   (a bare \code{"ci"} is treated as the self-pair \code{ci:ci}).
#' @return A data.frame with character columns \code{a} and \code{b} holding the
#'   two class codes of each pair.
#' @keywords internal
.region_parse_pairs <- function(cn) {
        sp <- strsplit(cn, ":", fixed = TRUE)
        data.frame(
                a = vapply(sp, `[`, character(1), 1L),
                b = vapply(sp, function(x) if (length(x) >= 2) x[2] else x[1], character(1)),
                stringsAsFactors = FALSE
        )
}

#' Adjusted Rand Index between two label vectors (NA-tolerant)
#'
#' @param a,b Label vectors of equal length; pairs with an \code{NA} in either
#'   vector are dropped before scoring.
#' @return The Adjusted Rand Index as a single numeric, \code{0} for the
#'   degenerate (zero-variance) case, or \code{NA} when fewer than two complete
#'   pairs are available.
#' @keywords internal
.region_ari <- function(a, b) {
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
#'
#' @param polys An \code{sf} object whose row order matches \code{values}.
#' @param values Vector to attach; must have length \code{nrow(polys)}.
#' @param name Character name for the new column.
#' @return The \code{sf} object \code{polys} with column \code{name} added.
#'   Stops if \code{polys} is not \code{sf} or the lengths disagree.
#' @keywords internal
.region_attach <- function(polys, values, name) {
        .region_require("sf")
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
#'
#' @param polys An \code{sf} layer carrying the column to map.
#' @param col Name of the column mapped to fill.
#' @param discrete Logical; \code{TRUE} for a qualitative \code{eco_palette()}
#'   fill, \code{FALSE} for a continuous viridis ("C") scale.
#' @param title,sub Plot title and subtitle.
#' @param legend Legend (fill) title.
#' @param border Polygon border colour. Default \code{"white"}.
#' @param alpha Optional name of a column mapped to the alpha aesthetic
#'   (membership degree); \code{NULL} (default) disables it.
#' @param shuffle Logical; if \code{TRUE}, randomly permute the qualitative
#'   palette (fixed seed) so adjacent classes contrast more. Discrete only.
#' @return A ggplot object.
#' @keywords internal
.region_choropleth <- function(polys, col, discrete, title, sub, legend, 
                               border = "white", alpha = NULL, shuffle = FALSE) {
        .region_require("ggplot2", "sf")
        if (is.null(alpha)){
                g <- ggplot2::ggplot(polys) +
                        ggplot2::geom_sf(ggplot2::aes(fill = .data[[col]]), colour = NA, linewidth = 0.05)
        } else {
                g <- ggplot2::ggplot(polys) +
                        ggplot2::geom_sf(ggplot2::aes(fill = .data[[col]], 
                                                      alpha = .data[[alpha]]),
                                         colour = NA,
                                         linewidth = 0.05) + 
                        ggplot2::labs(alpha = "Membership Degree")
        }
        
        if (discrete) {
                
                k <- data.table::uniqueN(polys[[col]])
                if (shuffle == TRUE){
                        set.seed(1)
                        final_palette <- sample(eco_palette(k))
                } else if (shuffle == FALSE){
                        final_palette <- eco_palette(k)
                }
                g <- g + ggplot2::scale_fill_manual(values = final_palette, na.value = "grey85",
                                                    drop = FALSE, name = legend)
        } else {
                g <- g + ggplot2::scale_fill_viridis_c(option = "C", na.value = "grey85", name = legend)
        }
        g + ggplot2::labs(title = title, subtitle = sub) + theme_eco(map = TRUE)
}


# ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
# LEVEL 1: REGIONS  (output of skater_con())
# ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾

#' Map of fuzzy regions, or their assignment uncertainty
#'
#' The two faces of a fuzzy regionalisation: \code{fill = "region"} shows the
#' hard partition (the regions), \code{fill = "uncertainty"} shows where the
#' partition is *soft* -- the ecotones / transitional catchments that the hard
#' map hides; \code{fill = "both"} shades the dominant region by its membership degree
#' (transitional catchments fade out); \code{fill = "core"} does the same but
#' masks everything outside the supplied core regions.
#'
#' @param polys An \code{sf} polygon/line layer, row-aligned to the catchments
#'   that were clustered (same order as \code{regions$hard_clusters}).
#'   
#' @param regions The list returned by \code{skater_con()}.
#' @param core Optional list of core-region objects (each with an \code{ID}
#'   field). Required only when \code{fill = "core"}; ignored otherwise.
#' @param fill One of "region" (default), "uncertainty", "both", or "core".
#' @param uncertainty_method Passed to the entropy/margin/maxprob helper.
#' @return A ggplot object.
#' @export
plot_region_map <- function(polys, regions, core = NULL, fill = c("region",
                           "uncertainty", "both", "core"), uncertainty_method = 
                           "entropy") {
        fill <- match.arg(fill)
        if (fill == "region") {
                lab <- factor(regions$hard_clusters)
                polys <- .region_attach(polys, lab, ".region")
                .region_choropleth(polys, ".region", discrete = TRUE,
                                   title = "Regions",
                                   sub = sprintf("%d regions (hard partition) \u00b7 %d mapping units",
                                                 nlevels(lab), length(lab)),
                                   legend = "Region")
        } else if (fill == "uncertainty"){
                u <- .region_uncertainty(regions$memberships, uncertainty_method)
                polys <- .region_attach(polys, u, ".unc")
                .region_choropleth(polys, ".unc", discrete = FALSE,
                                   title = "Region assignment uncertainty",
                                   sub = "Normalised membership entropy\u00b7 bright = transitional (ecotone)",
                                   legend = "Uncertainty")
        } else if (fill == "both"){
                u <- regions$memberships
                pd <- dplyr::bind_cols(polys, as.data.frame(u))
                class_cols <- dplyr::setdiff(names(pd), names(polys))
                
                # Calculate dominant region and membership to dominant region
                data_summary <- dplyr::rowwise(pd)
                data_summary <- dplyr::mutate(data_summary,
                                              Max_Memb = max(dplyr::c_across(dplyr::all_of(class_cols))),
                                              Dominant_Class = class_cols[which.max(dplyr::c_across(dplyr::all_of(class_cols)))])
                data_summary$Dominant_Class <- gsub("V", "", data_summary$Dominant_Class)
                data_summary <- dplyr::ungroup(data_summary)
                #polys <- .region_attach(polys, u, ".unc")
                .region_choropleth(data_summary, "Dominant_Class", discrete = TRUE,
                                   title = "Continuous Region Membership",
                                   alpha = "Max_Memb",
                                   sub = "Degree of Membership for highest membership region \n bright = low membership = transitional (ecotone)",
                                   legend = "Region")
        } else if (fill == "core"){
                if (is.null(core)){
                        stop("If method core is selected a core object must be provided.")
                }
                u <- regions$memberships
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
                core_ids <- unlist(sapply(core, function(x) x$ID))
                data_summary$Max_Memb[!data_summary$ID %in% core_ids] <- 0
                
                .region_choropleth(data_summary, "Dominant_Class", discrete = TRUE,
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
#' @param regions Result of \code{skater_con()}.
#' @return A patchwork of two ggplots.
#' @export
plot_region_membership <- function(regions) {
        .region_require("ggplot2", "patchwork")
        U   <- as.matrix(regions$memberships)
        hard <- regions$hard_clusters
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

# ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
# Tuning diagnostics: plotting
# ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾

#' Tuning axes recognised in a `tuning_log`, in canonical order
#' @keywords internal
#' @noRd
.REGION_AXES <- c("intermediate_regions", "weight_exponent",
                  "partitions_per_tree", "n_st", "final_regions")

#' Axis display labels
#' @keywords internal
#' @noRd
.REGION_AXIS_LABELS <- c(
        intermediate_regions = "intermediate regions",
        weight_exponent      = "weight exponent",
        partitions_per_tree  = "partitions per tree",
        n_st                 = "spanning trees",
        final_regions        = "final regions (k)"
)

#' Optimisation direction of each validity index
#'
#' MUST stay in sync with `.derive_score_column()`. Keeping the directions in
#' one place is the point: a silent mismatch between the tuner and the plots
#' would render a voter-agreement panel that is confidently upside down.
#'
#' @keywords internal
#' @noRd
.REGION_CVI_DIRECTION <- c(
        fch             = "max",   # fuzzy Calinski-Harabasz
        Fukuyama_Sugeno = "min",
        gd5             = "max",   # generalised Dunn
        fhv             = "min",   # fuzzy hyper-volume
        STAB            = "max"    # bootstrap ARI
)

#' Extract and normalise a tuning log
#'
#' Accepts either a `"regions_tuning"` object or the bare `tuning_log` data
#' frame, and renames legacy columns (`n.rst` -> `n_st`, `skater_regions` ->
#' `intermediate_regions`) so plots work on logs written by older versions.
#'
#' @keywords internal
#' @noRd
.region_tuning_log <- function(regions) {
        tl <- if (is.data.frame(regions)) regions else regions$tuning_log
        if (is.null(tl)) {
                stop("No `tuning_log` found. This looks like a single-configuration run ",
                     "(no axis had more than one candidate), so there is no tuning ",
                     "surface to plot.", call. = FALSE)
        }
        if (!is.data.frame(tl) || nrow(tl) == 0L) {
                stop("`tuning_log` is empty.", call. = FALSE)
        }
        legacy <- c(n.rst = "n_st", skater_regions = "intermediate_regions")
        for (old in names(legacy)) {
                new <- legacy[[old]]
                if (old %in% names(tl) && !new %in% names(tl)) names(tl)[names(tl) == old] <- new
        }
        if (!any(.REGION_AXES %in% names(tl))) {
                stop("`tuning_log` contains none of the expected axis columns (",
                     paste(.REGION_AXES, collapse = ", "), ").", call. = FALSE)
        }
        tl
}

#' Axes present in a log that actually vary
#' @keywords internal
#' @noRd
.region_active_axes <- function(tl) {
        present <- intersect(.REGION_AXES, names(tl))
        present[vapply(present, function(a) length(unique(tl[[a]])) > 1L, logical(1))]
}

#' Row of the log holding the optimum
#' @keywords internal
#' @noRd
.region_best_row <- function(tl, score_col = "score") {
        if (!score_col %in% names(tl)) {
                stop("Column '", score_col, "' not found in `tuning_log`. Available: ",
                     paste(names(tl), collapse = ", "), ".", call. = FALSE)
        }
        if (all(is.na(tl[[score_col]]))) {
                stop("Column '", score_col, "' is entirely NA.", call. = FALSE)
        }
        tl[which.max(tl[[score_col]]), , drop = FALSE]
}

#' Restrict a log to the slice through the optimum along `free_axes`
#'
#' Every active axis not in `free_axes` is pinned at its value in the best row.
#'
#' @keywords internal
#' @noRd
.region_slice_at_best <- function(tl, free_axes, score_col = "score") {
        best   <- .region_best_row(tl, score_col)
        pinned <- setdiff(.region_active_axes(tl), free_axes)
        keep   <- rep(TRUE, nrow(tl))
        for (a in pinned) keep <- keep & (tl[[a]] == best[[a]])
        list(data = tl[keep, , drop = FALSE], best = best, pinned = pinned)
}

#' Human-readable description of the pinned axes
#' @keywords internal
#' @noRd
.region_pin_label <- function(best, pinned) {
        if (!length(pinned)) return(NULL)
        paste0("held at optimum: ",
               paste(sprintf("%s = %s", .REGION_AXIS_LABELS[pinned],
                             vapply(pinned, function(a) format(best[[a]]), character(1))),
                     collapse = ", "))
}


# ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
# 1. Two-axis tuning surface
# ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾

#' Tuning Surface Over Two Hyper-parameter Axes
#'
#' Heat map of the tuning score over two chosen axes. With more than two axes
#' searched, the remaining ones must be disposed of explicitly, because a tile
#' plot silently overplots duplicated cells and would otherwise show an
#' arbitrary slice.
#'
#' @param regions A `"regions_tuning"` object (or its `tuning_log`).
#' @param x,y Axis names. `NULL` (default) picks the two active axes, with
#'   `final_regions` on x when it varies.
#' @param score_col Column to map to fill. Default `"score"` (higher is better
#'   by construction, including for Borda, which is stored negated).
#' @param others How to handle active axes other than `x` and `y`:
#'   `"slice"` (default) pins them at their value in the best row;
#'   `"facet"` panels over their combinations; `"max"` projects by taking the
#'   best score over them (a profile surface, not a slice).
#' @param max_facets Integer. Refuse to facet beyond this many panels. Default 12.
#'
#' @return A `ggplot` object.
#' @export
plot_region_cvi_grid <- function(regions,
                                 x = NULL, y = NULL,
                                 score_col = "score",
                                 others = c("slice", "facet", "max"),
                                 max_facets = 12L) {
        .region_require("ggplot2")
        others <- match.arg(others)
        
        tl     <- .region_tuning_log(regions)
        active <- .region_active_axes(tl)
        
        if (length(active) < 2L) {
                stop("Only ", length(active), " axis varies in this log (",
                     paste(active, collapse = ", "), "); a surface needs two. ",
                     "Use plot_region_profiles() instead.", call. = FALSE)
        }
        
        # --- Choose the plotted axes ----------------------------------------
        if (is.null(x) && is.null(y)) {
                x <- if ("final_regions" %in% active) "final_regions" else utils::tail(active, 1L)
                y <- setdiff(active, x)[1L]
        } else if (is.null(y)) {
                y <- setdiff(active, x)[1L]
        } else if (is.null(x)) {
                x <- setdiff(active, y)[1L]
        }
        for (a in c(x, y)) {
                if (!a %in% names(tl)) stop("Axis '", a, "' is not a column of `tuning_log`.",
                                            call. = FALSE)
        }
        if (identical(x, y)) stop("`x` and `y` must differ.", call. = FALSE)
        
        rest <- setdiff(active, c(x, y))
        
        # --- Dispose of the remaining axes ----------------------------------
        sub     <- tl
        best    <- .region_best_row(tl, score_col)
        subtitle_extra <- NULL
        facet_by <- NULL
        
        if (length(rest)) {
                if (others == "slice") {
                        sl   <- .region_slice_at_best(tl, free_axes = c(x, y), score_col = score_col)
                        sub  <- sl$data
                        subtitle_extra <- .region_pin_label(sl$best, sl$pinned)
                } else if (others == "facet") {
                        n_panels <- nrow(unique(tl[, rest, drop = FALSE]))
                        if (n_panels > max_facets) {
                                stop(sprintf(paste0(
                                        "Faceting over %s would need %d panels (max_facets = %d). ",
                                        "Use others = \"slice\" or \"max\", or name `x`/`y` differently."),
                                        paste(rest, collapse = " x "), n_panels, max_facets),
                                     call. = FALSE)
                        }
                        sub$.facet <- do.call(paste,
                                              c(lapply(rest, function(a)
                                                      sprintf("%s = %s", a, format(tl[[a]]))),
                                                sep = " | "))
                        facet_by <- ".facet"
                } else {                                   # "max"
                        keys <- do.call(paste, c(sub[, c(x, y), drop = FALSE], sep = "\r"))
                        ord  <- order(keys, -xtfrm(sub[[score_col]]), na.last = TRUE)
                        sub  <- sub[ord, , drop = FALSE]
                        sub  <- sub[!duplicated(keys[ord]), , drop = FALSE]
                        subtitle_extra <- paste0("best over ", paste(rest, collapse = ", "),
                                                 " (projection, not a slice)")
                }
        }
        
        if (!nrow(sub)) {
                stop("The requested slice is empty. This happens when the optimum was ",
                     "never evaluated jointly with the plotted axes, which is normal for ",
                     "coordinate-descent logs; try others = \"max\".", call. = FALSE)
        }
        
        # --- Refuse to draw an overplotted tile map -------------------------
        dup_key <- do.call(paste, c(sub[, c(x, y, facet_by), drop = FALSE], sep = "\r"))
        if (anyDuplicated(dup_key)) {
                stop("Multiple rows share the same (", x, ", ", y, ") cell after handling ",
                     "the other axes; geom_tile() would silently overplot them. ",
                     "Set others = \"slice\" or \"max\".", call. = FALSE)
        }
        
        # Coordinate descent produces cross-shaped, non-factorial logs. Say so
        # rather than letting the user read gaps as failed configurations.
        n_cells <- length(unique(sub[[x]])) * length(unique(sub[[y]]))
        if (nrow(sub) < n_cells) {
                message(sprintf("  %d of %d cells evaluated: the search was not factorial ",
                                "over these axes (expected for sequential / iterative).",
                                nrow(sub), n_cells))
        }
        
        sub$.x <- factor(sub[[x]], levels = sort(unique(sub[[x]])))
        sub$.y <- factor(sub[[y]], levels = sort(unique(sub[[y]])))
        
        # The best row may have been sliced out under others = "facet"/"max";
        # highlight it only if it survived.
        best_key <- do.call(paste, c(best[, c(x, y), drop = FALSE], sep = "\r"))
        sub_key  <- do.call(paste, c(sub[,  c(x, y), drop = FALSE], sep = "\r"))
        best_sub <- sub[sub_key == best_key, , drop = FALSE]
        
        p <- ggplot2::ggplot(sub, ggplot2::aes(.data$.x, .data$.y,
                                               fill = .data[[score_col]])) +
                ggplot2::geom_tile(colour = "white", linewidth = 0.4)
        
        if (nrow(best_sub)) {
                p <- p + ggplot2::geom_tile(data = best_sub, fill = NA,
                                            colour = "black", linewidth = 1.1)
        }
        
        p <- p +
                ggplot2::scale_fill_viridis_c(option = "D", name = score_col) +
                ggplot2::labs(
                        title    = "Region tuning surface",
                        subtitle = paste(c(
                                sprintf("optimum: %s = %s, %s = %s",
                                        .REGION_AXIS_LABELS[x], format(best[[x]]),
                                        .REGION_AXIS_LABELS[y], format(best[[y]])),
                                subtitle_extra), collapse = "\n"),
                        x = unname(.REGION_AXIS_LABELS[x]),
                        y = unname(.REGION_AXIS_LABELS[y])
                ) +
                theme_eco()
        
        if (!is.null(facet_by)) {
                p <- p + ggplot2::facet_wrap(stats::as.formula(paste("~", facet_by)))
        }
        p
}


# ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
# 2. One-dimensional profiles through the optimum
# ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾

#' Score Profiles Along Each Tuning Axis
#'
#' One panel per axis that varies, showing the score along that axis with every
#' other axis pinned at the optimum. This is the general replacement for the
#' two-axis heat map: it works for any number of axes and, unlike a surface,
#' does not require the search to have been factorial, so it reads correctly for
#' `"sequential"` and `"iterative"` runs.
#'
#' A flat profile means the axis does not matter over the range searched. A
#' profile still climbing at an endpoint means the candidate range was too
#' narrow, which is the usual finding for `n_st` and `partitions_per_tree`.
#'
#' @param regions A `"regions_tuning"` object (or its `tuning_log`).
#' @param score_col Column to plot. Default `"score"`.
#' @param axes Optional character vector restricting which axes are shown.
#'
#' @return A `ggplot` object.
#' @export
plot_region_profiles <- function(regions, score_col = "score", axes = NULL) {
        .region_require("ggplot2")
        
        tl     <- .region_tuning_log(regions)
        active <- .region_active_axes(tl)
        if (!length(active)) {
                stop("No axis varies in this log; there is nothing to profile.", call. = FALSE)
        }
        if (!is.null(axes)) {
                unknown <- setdiff(axes, active)
                if (length(unknown)) {
                        warning("Ignoring axes that do not vary: ",
                                paste(unknown, collapse = ", "), call. = FALSE)
                }
                active <- intersect(active, axes)
                if (!length(active)) stop("No requested axis varies.", call. = FALSE)
        }
        
        best <- .region_best_row(tl, score_col)
        
        prof <- do.call(rbind, lapply(active, function(a) {
                sl <- .region_slice_at_best(tl, free_axes = a, score_col = score_col)
                d  <- sl$data
                if (!nrow(d)) return(NULL)
                d <- d[order(d[[a]]), , drop = FALSE]
                data.frame(
                        axis    = a,
                        value   = as.numeric(d[[a]]),
                        score   = d[[score_col]],
                        is_best = d[[a]] == best[[a]],
                        stringsAsFactors = FALSE
                )
        }))
        if (is.null(prof) || !nrow(prof)) {
                stop("No profile rows survived slicing at the optimum.", call. = FALSE)
        }
        
        prof$axis <- factor(prof$axis, levels = intersect(.REGION_AXES, active),
                            labels = unname(.REGION_AXIS_LABELS[intersect(.REGION_AXES, active)]))
        
        ggplot2::ggplot(prof, ggplot2::aes(.data$value, .data$score)) +
                ggplot2::geom_line(linewidth = 0.5, colour = "grey35") +
                ggplot2::geom_point(size = 1.8, colour = "grey20") +
                ggplot2::geom_point(data = prof[prof$is_best, , drop = FALSE],
                                    size = 3.2, shape = 21, stroke = 1.1,
                                    fill = NA, colour = "black") +
                ggplot2::facet_wrap(~ .data$axis, scales = "free_x") +
                ggplot2::labs(
                        title    = "Tuning profiles through the optimum",
                        subtitle = "each panel varies one axis; all others pinned at their best value",
                        x = NULL, y = score_col
                ) +
                theme_eco()
}


# ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
# 3. Voter agreement among the validity indices
# ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾

#' Validity-index Agreement Along One Axis
#'
#' Borda aggregation hides whether the constituent indices agree. This plots
#' each index along a single axis, rescaled to `[0, 1]` and sign-corrected so
#' that higher is better in every case, with the other axes pinned at the
#' optimum. Indices that disagree about the location of the optimum are the
#' signal worth acting on: the Borda winner is only meaningful when the voters
#' are not pulling in opposite directions.
#'
#' @param regions A `"regions_tuning"` object (or its `tuning_log`).
#' @param axis Axis to plot along. `NULL` (default) picks the axis with the most
#'   evaluated candidates.
#' @param score_col Column used to locate the optimum. Default `"score"`.
#'
#' @return A `ggplot` object.
#' @export
plot_region_cvi_votes <- function(regions, axis = NULL, score_col = "score") {
        .region_require("ggplot2")
        
        tl     <- .region_tuning_log(regions)
        active <- .region_active_axes(tl)
        if (!length(active)) stop("No axis varies in this log.", call. = FALSE)
        
        if (is.null(axis)) {
                n_vals <- vapply(active, function(a) length(unique(tl[[a]])), integer(1))
                axis   <- active[which.max(n_vals)]
        }
        if (!axis %in% names(tl)) stop("Axis '", axis, "' is not a column of `tuning_log`.",
                                       call. = FALSE)
        
        cvis <- intersect(names(.REGION_CVI_DIRECTION), names(tl))
        if (!length(cvis)) {
                stop("No validity-index columns found. Expected some of: ",
                     paste(names(.REGION_CVI_DIRECTION), collapse = ", "), ".", call. = FALSE)
        }
        
        sl   <- .region_slice_at_best(tl, free_axes = axis, score_col = score_col)
        d    <- sl$data[order(sl$data[[axis]]), , drop = FALSE]
        best <- sl$best
        if (nrow(d) < 2L) stop("Fewer than two points on this axis after slicing.", call. = FALSE)
        
        rescale01 <- function(v, direction) {
                if (all(is.na(v))) return(rep(NA_real_, length(v)))
                if (direction == "min") v <- -v
                rng <- range(v, na.rm = TRUE)
                if (diff(rng) <= 0) return(rep(0.5, length(v)))
                (v - rng[1L]) / diff(rng)
        }
        
        votes <- do.call(rbind, lapply(cvis, function(cv) {
                data.frame(
                        index = cv,
                        value = as.numeric(d[[axis]]),
                        norm  = rescale01(d[[cv]], .REGION_CVI_DIRECTION[[cv]]),
                        stringsAsFactors = FALSE
                )
        }))
        votes <- votes[!is.na(votes$norm), , drop = FALSE]
        if (!nrow(votes)) stop("All validity indices are NA on this slice.", call. = FALSE)
        
        ggplot2::ggplot(votes, ggplot2::aes(.data$value, .data$norm,
                                            colour = .data$index)) +
                ggplot2::geom_vline(xintercept = as.numeric(best[[axis]]),
                                    linetype = "dashed", colour = "grey40") +
                ggplot2::geom_line(linewidth = 0.6) +
                ggplot2::geom_point(size = 1.6) +
                ggplot2::scale_colour_viridis_d(option = "D", name = NULL, end = 0.9) +
                ggplot2::labs(
                        title    = "Validity-index agreement",
                        subtitle = paste(c(
                                "rescaled to [0, 1], sign-corrected so higher is better",
                                .region_pin_label(best, sl$pinned),
                                sprintf("dashed line: selected %s", .REGION_AXIS_LABELS[axis])),
                                collapse = "\n"),
                        x = unname(.REGION_AXIS_LABELS[axis]),
                        y = "rescaled index"
                ) +
                theme_eco()
}


# ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
# 4. Composite wrapper
# ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾

#' Diagnostic Panel for a Region Tuning Run
#'
#' Assembles whichever panels the tuning log can support. With two axes searched
#' this is a surface plus profiles; with one, profiles alone; with three or more,
#' the surface is a slice through the optimum and the profiles carry the rest.
#'
#' @param regions A `"regions_tuning"` object (or its `tuning_log`).
#' @param panels Character vector; any of `"grid"`, `"profiles"`, `"votes"`.
#'   Panels that the log cannot support are dropped with a message.
#' @param score_col Column to plot. Default `"score"`.
#' @param ncol Columns in the assembled layout. `NULL` lets patchwork decide.
#' @param ... Passed to `plot_region_cvi_grid()`.
#'
#' @return A `patchwork` object, or a single `ggplot` when only one panel applies.
#' @export
plot_region_tuning <- function(regions,
                               panels    = c("grid", "profiles"),
                               score_col = "score",
                               ncol      = NULL,
                               ...) {
        .region_require("patchwork")
        panels <- match.arg(panels, choices = c("grid", "profiles", "votes"),
                            several.ok = TRUE)
        
        tl <- .region_tuning_log(regions)   # errors early with a useful message
        
        builders <- list(
                grid     = function() plot_region_cvi_grid(tl, score_col = score_col, ...),
                profiles = function() plot_region_profiles(tl, score_col = score_col),
                votes    = function() plot_region_cvi_votes(tl, score_col = score_col)
        )
        
        ps   <- list()
        skip <- character(0)
        for (nm in panels) {
                p <- tryCatch(builders[[nm]](), error = function(e) {
                        skip[[length(skip) + 1L]] <<- sprintf("%s (%s)", nm, conditionMessage(e))
                        NULL
                })
                if (!is.null(p)) ps[[nm]] <- p
        }
        
        if (!length(ps)) {
                stop("Nothing to plot. Reasons:\n  - ", paste(skip, collapse = "\n  - "),
                     call. = FALSE)
        }
        if (length(skip)) {
                message("Skipped panel(s):\n  - ", paste(skip, collapse = "\n  - "))
        }
        if (length(ps) == 1L) return(ps[[1L]])
        
        patchwork::wrap_plots(ps, ncol = if (is.null(ncol)) min(2L, length(ps)) else ncol)
}


#' Ternary membership plot for the special case of exactly 3 regions
#'
#' A self-contained simplex projection (no ggtern dependency). Each catchment is
#' a point in the membership triangle; corner = a crisp region, centre = a
#' three-way mixture.
#'
#' @param regions Result of \code{skater_con()} with exactly 3 final regions.
#' @return A ggplot object.
#' @export
plot_membership_ternary <- function(regions) {
        .region_require("ggplot2")
        U <- as.matrix(regions$memberships)
        if (ncol(U) != 3L) stop("Ternary plot requires exactly 3 regions (k = ", ncol(U), ").",
                                call. = FALSE)
        U <- U / pmax(rowSums(U), .Machine$double.eps)
        x <- U[, 2] + U[, 3] / 2
        y <- U[, 3] * sqrt(3) / 2
        df <- data.frame(x = x, y = y, region = factor(.region_hard_from_U(U)))
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
# LEVEL 2: LOCAL TYPES  (output of get_local_types())
# ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾

#' BIC model-selection curves for every regional GMM
#'
#' For each region, the best BIC at each candidate number of components G
#' (envelope over covariance models), with the selected G highlighted. Shows
#' how many local types each region's mixture "wanted", and how decisive that
#' choice was (sharp peak vs flat ridge).
#'
#' @param local Result of \code{get_local_types()}.
#' @return A ggplot object (faceted by region).
#' @export
plot_localtype_bic <- function(local) {
        .region_require("ggplot2")
        rows <- list()
        for (nm in names(local$models)) {
                fit <- local$models[[nm]]
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
                ggplot2::labs(title = "Local type model selection (per region)",
                              subtitle = "Best BIC at each G \u00b7 orange = selected number of local types",
                              x = "Number of local types (G)", y = "BIC") +
                theme_eco()
}

#' Environmental fingerprint of each local type
#'
#' Heatmap of the GMM component means, centred within each ecoregion (deviation
#' of each local type from its region's mean) and scaled by each predictor's
#' natural variability -- rows are predictors, columns are local types, faceted
#' by region. This is the *interpretation* plot: it answers "what does Region 2's
#' local type 3 actually mean ecologically?".
#'
#' @param local Result of \code{get_local_types()}.
#' @param top Optional integer; keep only the \code{top} most discriminating
#'   predictors (ranked by mean within-region between-type variance). Default
#'   \code{NULL} keeps all.
#' @param cap Absolute value at which the within-region z-scores are clipped, so
#'   a single extreme component cannot dominate the colour scale or the ranking.
#'   Default 3.
#' @param order_rows Logical; order predictor rows by discriminability. Default TRUE.
#' @return A ggplot object.
#' @export
plot_local_profiles <- function(local, top = NULL, cap = 3,
                                    order_rows = TRUE) {
        .region_require("ggplot2")
        
        ## component means in long form, plus core data for a fixed per-predictor scale
        rows <- list()
        pooled <- list()
        for (nm in names(local$models)) {
                fit <- local$models[[nm]]
                mu  <- fit$parameters$mean
                if (is.null(dim(mu))) {
                        vn <- colnames(fit$data)
                        if (is.null(vn)) vn <- "var"
                        mu <- matrix(mu, nrow = 1L, dimnames = list(vn[1], NULL))
                }
                vn <- rownames(mu)
                if (is.null(vn)) vn <- paste0("V", seq_len(nrow(mu)))
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
        df$z <- (df$value - stats::ave(df$value, df$region, df$variable, FUN = mean)) / s
        df$z[is.na(df$z)] <- 0
        
        ## cap so a stray component cannot blow out the palette (or the ranking)
        if (is.finite(cap)) df$z <- pmax(pmin(df$z, cap), -cap)
        
        ## discriminability = within-region between-type variance, averaged over regions
        rv    <- stats::aggregate(z ~ region + variable, df, stats::var)
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
                ggplot2::labs(title = "Local type environmental fingerprints",
                              subtitle = "GMM component means, centred within ecoregion, scaled by predictor variability",
                              x = "Local type", y = NULL) +
                theme_eco() + ggplot2::theme(panel.grid = ggplot2::element_blank())
}
#' Local type composition by region
#'
#' How much of the dataset each local type captures. \code{weight = "soft"}
#' (default) sums the membership-weighted probabilities (uses the fuzzy
#' information); \code{"hard"} counts argmax assignments.
#'
#' @param local Result of \code{get_local_types()}.
#' @param weight "soft" or "hard".
#' @return A ggplot object.
#' @export
plot_localtype_composition <- function(local, weight = c("soft", "hard")) {
        .region_require("ggplot2")
        weight <- match.arg(weight)
        tab <- .region_lt_class_table(local)
        if (weight == "soft") {
                tab$value <- colSums(as.matrix(local$localtypes))
                ylab <- "Soft membership mass"
        } else {
                hard <- table(.region_lt_hard(local))
                tab$value <- as.numeric(hard[tab$col])
                ylab <- "Catchments (argmax)"
        }
        ggplot2::ggplot(tab, ggplot2::aes(.data$region, .data$value, fill = .data$type)) +
                ggplot2::geom_col(width = 0.75, colour = "white", linewidth = 0.2) +
                ggplot2::scale_fill_manual(values = eco_palette(nlevels(tab$type)), name = "Local type") +
                ggplot2::labs(title = "Local type composition",
                              subtitle = sprintf("Weighting: %s", weight),
                              x = NULL, y = ylab) +
                theme_eco() + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))
}

#' Map of dominant local type, or local-type assignment uncertainty
#'
#' @param polys An \code{sf} layer row-aligned to \code{local$localtypes}.
#' @param local Result of \code{get_local_types()}.
#' @param fill One of "type" (default), "uncertainty", or "both".
#' @param uncertainty_method Passed to the entropy/margin/maxprob helper.
#' @return A ggplot object.
#' @export
plot_local_map <- function(polys, local, fill = c("type", "uncertainty", "both"),
                           uncertainty_method = "entropy") {
        fill <- match.arg(fill)
        if (fill == "type") {
                lab <- .region_lt_hard(local)
                lab <- droplevels(lab)
                polys <- .region_attach(polys, lab, ".local")
                .region_choropleth(polys, ".local", discrete = TRUE,
                                   title = "Dominant local type",
                                   sub = sprintf("%d realised local type classes", nlevels(lab)),
                                   legend = "Local type", shuffle = TRUE)
        } else if (fill == "uncertainty"){
                u <- .region_uncertainty(local$localtypes, uncertainty_method)
                polys <- .region_attach(polys, u, ".unc")
                .region_choropleth(polys, ".unc", discrete = FALSE,
                                   title = "Local-type assignment uncertainty",
                                   sub = "Normalised entropy of weighted local-type probabilities",
                                   legend = "Uncertainty")
        } else if (fill == "both"){
                u <- local$localtypes
                colnames(u) <- 
                        sub(pattern = "^Region([0-9]+)localType([0-9]+)$", replacement = "R\\1\u00b7T\\2", x = colnames(u))
                pd <- dplyr::bind_cols(polys, as.data.frame(u))
                class_cols <- dplyr::setdiff(names(pd), names(polys))
                
                # Calculate dominant region and membership to dominant region
                data_summary <- dplyr::rowwise(pd)
                data_summary <- dplyr::mutate(data_summary,
                                              Max_Memb = max(dplyr::c_across(dplyr::all_of(class_cols))),
                                              Dominant_Class = class_cols[which.max(dplyr::c_across(dplyr::all_of(class_cols)))])
                data_summary$Dominant_Class <- gsub("V", "", data_summary$Dominant_Class)
                data_summary <- dplyr::ungroup(data_summary)
                #polys <- .region_attach(polys, u, ".unc")
                .region_choropleth(data_summary, "Dominant_Class", discrete = TRUE,
                                   title = "Continuous Local Type Membership",
                                   alpha = "Max_Memb",
                                   sub = "Degree of Membership for highest membership region \n bright = low membership = transitional (ecotone)",
                                   legend = "Local type", shuffle = TRUE)
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
#' @param fill One of "type" (default), "uncertainty", or "both".
#' @param uncertainty_method Passed to the entropy/margin/maxprob helper.
#' @return A ggplot object.
#' @export
plot_mosaic_map <- function(polys, mosaic, fill = c("type", "uncertainty", "both"),
                            uncertainty_method = "entropy") {
        fill  <- match.arg(fill)
        valid <- mosaic$valid
        if (fill == "type") {
                lab <- rep(NA_integer_, length(valid))
                lab[valid] <- .region_mosaic_hard(mosaic)
                lab <- factor(lab)
                polys <- .region_attach(polys, lab, ".ls")
                .region_choropleth(polys, ".ls", discrete = TRUE,
                                   title = "Mosaic types",
                                   sub = sprintf("k = %d co-occurrence clusters \u00b7 grey = isolated polygons",
                                                 mosaic$best_k),
                                   legend = "Mosaic type")
        } else if (fill == "uncertainty"){
                u <- rep(NA_real_, length(valid))
                u[valid] <- .region_uncertainty(.region_mosaic_U(mosaic), uncertainty_method)
                polys <- .region_attach(polys, u, ".unc")
                .region_choropleth(polys, ".unc", discrete = FALSE,
                                   title = "Mosaic-type uncertainty",
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
                #polys <- .region_attach(polys, u, ".unc")
                .region_choropleth(data_summary, "Dominant_Class", discrete = TRUE,
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
        .region_require("ggplot2")
        S   <- as.matrix(sigs)[mosaic$valid, , drop = FALSE]
        lab <- .region_mosaic_hard(mosaic)
        pr  <- .region_parse_pairs(colnames(S))
        classes <- sort(unique(c(pr$a, pr$b)))
        ulab <- sort(unique(lab))
        
        rows <- list()
        for (g in ulab) {
                m <- colMeans(S[lab == g, , drop = FALSE], na.rm = TRUE)
                for (j in seq_along(m)) {
                        rows[[length(rows) + 1L]] <- data.frame(
                                cluster = factor(paste0("Mosaic ", g)),
                                a = pr$a[j], b = pr$b[j], value = m[j]
                        )
                        if (pr$a[j] != pr$b[j]) {  # mirror to fill the lower triangle
                                rows[[length(rows) + 1L]] <- data.frame(
                                        cluster = factor(paste0("Mosaic ", g)),
                                        a = pr$b[j], b = pr$a[j], value = m[j]
                                )
                        }
                }
        }
        df <- do.call(rbind, rows)
        
        ia <- match(as.character(df$a), classes)   # column index in class order
        ib <- match(as.character(df$b), classes)   # row index in class order
        df <- df[ib >= ia, , drop = FALSE]         # lower triangle + diagonal
        
        df$a <- factor(df$a, levels = classes)
        df$b <- factor(df$b, levels = rev(classes))
        
        df$value[df$value == 0] <- NA          # never-occurring combos -> NA
        
        ggplot2::ggplot(df, ggplot2::aes(.data$a, .data$b, fill = .data$value)) +
                ggplot2::geom_tile(colour = "white", linewidth = 0.3) +
                ggplot2::facet_wrap(~ cluster) +
                ggplot2::scale_fill_viridis_c(option = "D", name = "mean\nco-occur.",
                                              na.value = "white") +   # NA -> white
                ggplot2::coord_equal() +
                ggplot2::labs(title = "Mosaic-type co-occurrence fingerprints",
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
        .region_require("ggplot2")
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
#' @param ellipses Draw 95% normal ellipses per cluster? Default TRUE.
#' @return A ggplot object.
#' @export
plot_mosaic_ordination <- function(sigs, mosaic, ellipses = TRUE) {
        .region_require("ggplot2")
        S <- as.matrix(sigs)[mosaic$valid, , drop = FALSE]
        d <- if (exists("js_distance_fast", mode = "function")) {
                get("js_distance_fast")(S)
        } else {
                eps <- .Machine$double.eps; P <- (S + eps) / rowSums(S + eps)
                stats::dist(P)
        }
        mds <- stats::cmdscale(d, k = 2)
        df  <- data.frame(MDS1 = mds[, 1], MDS2 = mds[, 2],
                          cluster = factor(.region_mosaic_hard(mosaic)))
        g <- ggplot2::ggplot(df, ggplot2::aes(.data$MDS1, .data$MDS2, colour = .data$cluster)) +
                ggplot2::geom_point(alpha = 0.6, size = 1.4)
        if (ellipses) g <- g + ggplot2::stat_ellipse(level = 0.95, linewidth = 0.7)
        g + ggplot2::scale_colour_manual(values = eco_palette(nlevels(df$cluster)),
                                         name = "Mosaic") +
                ggplot2::labs(title = "Mosaic signature ordination",
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
        .region_require("ggplot2")
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
#' (region and local-type span the full set; mosaic labels are scattered
#' through \code{mosaic$valid}). Produces one tidy data frame with hard labels and
#' per-level uncertainties -- the substrate for every cross-level plot.
#'
#' @param region,local,mosaic Results of the three pipeline stages (any may be NULL).
#' @return A data.frame with hard-label columns \code{region}, \code{local},
#'   \code{mosaic} and per-level uncertainty columns \code{u_region},
#'   \code{u_lt}, \code{u_mosaic} (only those whose input was supplied).
#' @export
assemble_typology <- function(region = NULL, local = NULL, mosaic = NULL) {
        lens <- c(region  = if (!is.null(region))  length(region$hard_clusters) else NA,
                  local   = if (!is.null(local))   nrow(local$localtypes)       else NA,
                  mosaic  = if (!is.null(mosaic))  length(mosaic$valid)         else NA)
        N <- unique(stats::na.omit(lens))
        if (length(N) > 1) stop("Inputs are not row-aligned (lengths: ",
                                paste(names(lens), lens, sep = "=", collapse = ", "), ").",
                                call. = FALSE)
        N <- N[1]
        
        out <- data.frame(.row = seq_len(N))
        if (!is.null(region)) {
                out$region <- factor(region$hard_clusters)
                out$u_region     <- .region_uncertainty(region$memberships)
        }
        if (!is.null(local)) {
                out$local <- droplevels(.region_lt_hard(local))
                out$u_lt       <- .region_uncertainty(local$localtypes)
        }
        if (!is.null(mosaic)) {
                ls <- rep(NA_integer_, N); ls[mosaic$valid] <- .region_mosaic_hard(mosaic)
                out$mosaic <- factor(ls)
                uu <- rep(NA_real_, N); uu[mosaic$valid] <- .region_uncertainty(.region_mosaic_U(mosaic))
                out$u_mosaic <- uu
        }
        out$.row <- NULL
        out
}

#' Make alluvial stratum labels unique & readable per axis
#'
#' @param x Vector of labels to relabel (coerced to character).
#' @param axis One of \code{"region"}, \code{"local"} or \code{"mosaic"}; any
#'   other value gets an abbreviated prefix.
#' @return A character vector of prefixed, disambiguated labels.
#' @keywords internal
.region_axis_relabel <- function(x, axis) {
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
#' can read how each region splits into local types and how those reassemble
#' into mosaic types (i.e. how aligned or crossed the three partitions are).
#'
#' @param typ Output of \code{assemble_typology()}.
#' @param axes Which label columns to chain (in order).
#' @return A ggplot object.
#' @export
plot_hierarchy_alluvial <- function(typ,
                                    axes = c("region", "local", "mosaic")) {
        .region_require("ggplot2", "ggalluvial")
        axes <- intersect(axes, names(typ))
        if (length(axes) < 2) stop("Need at least two label columns to draw flows.", call. = FALSE)
        d <- typ[stats::complete.cases(typ[, axes, drop = FALSE]), axes, drop = FALSE]
        # Disambiguate identical codes across axes (e.g. region "1" vs mosaic "1")
        # so flows join correctly and stratum labels stay self-describing.
        for (ax in axes) d[[ax]] <- factor(.region_axis_relabel(d[[ax]], ax))
        d$freq <- 1L
        agg <- stats::aggregate(freq ~ ., data = d, FUN = sum)
        
        mapping <- ggplot2::aes(y = .data$freq)
        for (i in seq_along(axes)) mapping[[paste0("axis", i)]] <- as.name(axes[i])
        
        ggplot2::ggplot(agg, mapping) +
                ggalluvial::geom_alluvium(ggplot2::aes(fill = .data[[axes[1]]]), alpha = 0.6) +
                ggalluvial::geom_stratum(width = 0.28, fill = "grey92", colour = "grey50") +
                ggplot2::geom_text(stat = ggalluvial::StatStratum,
                                   ggplot2::aes(label = ggplot2::after_stat(.data$stratum)), size = 3) +
                ggplot2::scale_x_discrete(limits = axes, expand = c(0.08, 0.08)) +
                ggplot2::scale_fill_manual(values = eco_palette(nlevels(agg[[axes[1]]])),
                                           guide = "none") +
                ggplot2::labs(title = "Typology hierarchy",
                              subtitle = paste(axes, collapse = "  \u2192  "),
                              x = NULL, y = "Mapping Units") +
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
        .region_require("ggplot2")
        for (cc in c(x, y)) if (!cc %in% names(typ))
                stop("Column '", cc, "' not in the typology table.", call. = FALSE)
        ok <- !is.na(typ[[x]]) & !is.na(typ[[y]])
        tab <- table(typ[[x]][ok], typ[[y]][ok])
        if (x == "local"){
                rownames(tab) <- .region_axis_relabel(rownames(tab), axis = "local")
        } else if (y == "local"){
                colnames(tab) <- .region_axis_relabel(colnames(tab), axis = "local")
        }
        prop <- sweep(tab, 1, pmax(rowSums(tab), 1), "/")
        df <- as.data.frame(prop)
        names(df) <- c("xcl", "ycl", "share")
        ari <- .region_ari(typ[[x]], typ[[y]])
        
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
#' @param region,local,mosaic The three pipeline results (any may be NULL).
#' @return A patchwork object.
#' @export
plot_hierarchy_maps <- function(polys, region = NULL, local = NULL, mosaic = NULL) {
        .region_require("patchwork")
        ps <- list()
        if (!is.null(region))  ps$region <- plot_region_map(polys, region) +
                ggplot2::labs(title = "1 \u00b7 Region", subtitle = NULL)
        if (!is.null(local))   ps$local <- plot_local_map(polys, local) +
                ggplot2::labs(title = "2 \u00b7 Local type", subtitle = NULL)
        if (!is.null(mosaic)) ps$mosaic <- plot_mosaic_map(polys, mosaic) +
                ggplot2::labs(title = "3 \u00b7 Mosaic type", subtitle = NULL)
        if (!length(ps)) stop("Supply at least one of `region`, `local`, `mosaic`.", call. = FALSE)
        patchwork::wrap_plots(ps, nrow = 1) +
                patchwork::plot_annotation(
                        title = "The three typology levels in space",
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
        .region_require("ggplot2")
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