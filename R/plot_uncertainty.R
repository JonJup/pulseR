# -----------------------------------------------------------------------------
# Distributions of assignment certainty and margin across typology levels
#
# Neither figure needs geometry: both read only the membership matrices, so
# they run on a headless node straight after fitting.
# -----------------------------------------------------------------------------

#' Distribution of assignment certainty across typology levels
#'
#' One panel per supplied level, showing how strongly units belong to the type
#' they were assigned. Certainty is the unit's membership in its assigned type:
#' own-region \eqn{\Pi} for regions (which can sit below the row maximum where
#' contiguity repair moved a unit), and the row maximum for local, meta- and
#' mosaic types. No geometry is needed.
#'
#' A row of \eqn{k} memberships cannot fall below \eqn{1/k}, so raw certainty
#' is not comparable across levels or fits with different numbers of types.
#' The dotted line marks \eqn{1/k} where \eqn{k} is constant within a level;
#' \code{rescale = TRUE} maps \eqn{[1/k, 1]} onto \eqn{[0, 1]} unit by unit.
#' Local types are conditional on their region, so for them \eqn{k} is the
#' number of local types in the unit's region.
#'
#' Memberships do not mean the same thing at every level (ensemble composition
#' at Level 1, mixture posteriors for local types and \code{method = "dmm"}
#' mosaic types, fuzzifier-dependent geometry for \code{method = "pam"} mosaic
#' types). Compare a level across fits; read differences between levels with
#' that in mind.
#'
#' @param regions,local,meta,mosaic Fitted objects (\code{"pulse_regions"},
#'   \code{"pulse_local_types"}, \code{"pulse_meta_types"},
#'   \code{"mosaic_types"}); any may be \code{NULL}, at least one is required.
#' @param rescale Logical. Map certainty from \eqn{[1/k, 1]} onto \eqn{[0, 1]}.
#'   Default \code{FALSE}.
#' @param drop_trivial Logical. Drop units whose type had no competitor
#'   (\eqn{k = 1}: a one-region fit, or a region with a single local type).
#'   They are certain by construction and pile up at 1. Default \code{TRUE};
#'   the number dropped is reported in the caption.
#' @param type \code{"histogram"} (default, one panel per level) or
#'   \code{"ecdf"} (all levels on one panel).
#' @param bins Number of histogram bins on \eqn{[0, 1]}. Default 50.
#' @return A \code{ggplot} object. The per-unit values are in \code{$data}.
#' @seealso \code{\link{plot_margin}}
#' @examples
#' \dontrun{
#' plot_certainty(regions = reg, local = loc, mosaic = mos)
#' plot_certainty(regions = reg, local = loc, rescale = TRUE, type = "ecdf")
#' }
#' @export
plot_certainty <- function(regions = NULL, local = NULL, meta = NULL, mosaic = NULL,
                           rescale      = FALSE,
                           drop_trivial = TRUE,
                           type         = c("histogram", "ecdf"),
                           bins         = 50L) {
        .need("ggplot2")
        type <- match.arg(type)
        df   <- .uncertainty_df(list(regions = regions, local = local,
                                     meta = meta, mosaic = mosaic), drop_trivial)
        if (rescale) {
                df$value <- ifelse(df$k > 1L,
                                   pmin(pmax((df$certainty - 1 / df$k) / (1 - 1 / df$k), 0), 1),
                                   1)
        } else df$value <- df$certainty

        floor_df <- NULL
        if (!rescale) {
                kk <- tapply(df$k, df$level, function(v) if (length(unique(v)) == 1L) v[1] else NA)
                kk <- kk[!is.na(kk) & kk > 1]
                if (length(kk))
                        floor_df <- data.frame(level = factor(names(kk), levels = levels(df$level)),
                                               x = 1 / as.numeric(kk))
        }
        .plot_unit_distribution(
                df, type, bins, floor_df = floor_df, threshold = NULL,
                title = if (rescale) "Assignment certainty, rescaled to [0, 1]"
                        else "Assignment certainty",
                xlab  = if (rescale) "certainty, (own - 1/k) / (1 - 1/k)"
                        else "certainty (membership in assigned type)",
                subtitle = if (rescale) "solid line = median"
                           else "solid line = median, dotted line = 1/k (no information)")
}


#' Distribution of the membership margin across typology levels
#'
#' The margin is a unit's highest type membership minus its second highest.
#' Near 0, the unit sits between two types; near 1, one type dominates. Unlike
#' certainty, the margin has the same floor (0) at every \eqn{k}, so it is the
#' better quantity for comparing levels or fits with different numbers of
#' types. Local types are compared within the unit's region. No geometry is
#' needed.
#'
#' The dashed line and the per-panel percentage use \code{threshold}, the same
#' quantity \code{plot_regions(what = "transition", margin = )} maps.
#'
#' @inheritParams plot_certainty
#' @param threshold Margin at or below which a unit counts as transitional.
#'   Default 0.1. \code{NULL} drops the line and the percentage.
#' @return A \code{ggplot} object. The per-unit values are in \code{$data}.
#' @seealso \code{\link{plot_certainty}}, \code{\link{plot_regions}}
#' @examples
#' \dontrun{
#' plot_margin(regions = reg, local = loc, meta = met, mosaic = mos)
#' plot_margin(regions = reg, mosaic = mos, threshold = 0.2, type = "ecdf")
#' }
#' @export
plot_margin <- function(regions = NULL, local = NULL, meta = NULL, mosaic = NULL,
                        threshold    = 0.1,
                        drop_trivial = TRUE,
                        type         = c("histogram", "ecdf"),
                        bins         = 50L) {
        .need("ggplot2")
        type <- match.arg(type)
        if (!is.null(threshold))
                stopifnot(is.numeric(threshold), length(threshold) == 1L,
                          threshold >= 0, threshold <= 1)
        df <- .uncertainty_df(list(regions = regions, local = local,
                                   meta = meta, mosaic = mosaic), drop_trivial)
        df$value <- df$margin
        .plot_unit_distribution(
                df, type, bins, floor_df = NULL, threshold = threshold,
                title    = "Membership margin (highest minus second highest)",
                xlab     = "margin",
                subtitle = if (is.null(threshold)) "solid line = median"
                           else sprintf("solid line = median, dashed line = %.2f (transition threshold)",
                                        threshold))
}


# -----------------------------------------------------------------------------
# Internals
# -----------------------------------------------------------------------------

.UNCERTAINTY_LEVELS <- c(regions = "Regions", local = "Local types",
                         meta = "Meta-types", mosaic = "Mosaic types")

.UNCERTAINTY_COLS <- c("Regions" = "#0072B2", "Local types" = "#009E73",
                       "Meta-types" = "#E69F00", "Mosaic types" = "#CC79A7")

#' Highest and second-highest entry per row
#'
#' Processed in row blocks so a wide local-type matrix is never copied whole.
#' Rows with any NA give NA. A one-column matrix gives second = 0.
#' @keywords internal
#' @noRd
.top_two <- function(U, block = 1e7) {
        N <- nrow(U); K <- ncol(U)
        first <- second <- rep(NA_real_, N)
        step <- max(1L, as.integer(block %/% max(K, 1L)))
        for (s in seq(1L, N, by = step)) {
                rows <- s:min(N, s + step - 1L)
                B    <- as.matrix(U[rows, , drop = FALSE])
                ok   <- stats::complete.cases(B)
                if (!any(ok)) next
                B    <- B[ok, , drop = FALSE]
                n    <- nrow(B)
                top  <- max.col(B, ties.method = "first")
                first[rows[ok]] <- B[cbind(seq_len(n), top)]
                if (K > 1L) {
                        B[cbind(seq_len(n), top)] <- -Inf
                        second[rows[ok]] <- B[cbind(seq_len(n), max.col(B, ties.method = "first"))]
                } else second[rows[ok]] <- 0
        }
        list(first = first, second = second)
}

#' Per-unit certainty, margin and effective k for one fitted level
#' @keywords internal
#' @noRd
.level_uncertainty <- function(x, slot) {
        expect <- c(regions = "pulse_regions", local = "pulse_local_types",
                    meta = "pulse_meta_types", mosaic = "mosaic_types")[[slot]]
        if (!inherits(x, expect))
                stop(sprintf("`%s` must be a \"%s\" object, not \"%s\".",
                             slot, expect, class(x)[1]), call. = FALSE)
        U  <- x$memberships
        N  <- nrow(U)
        t2 <- .top_two(U)
        switch(slot,
               regions = {
                       cert <- .confidence_vec(x, "certainty")
                       k    <- rep(ncol(U), N)
               },
               local = {
                       cert <- t2$first
                       G    <- tabulate(match(x$key$region, unique(x$key$region)))
                       k    <- G[match(x$key$region[x$hard_classes], unique(x$key$region))]
               },
               meta = {
                       cert <- x$certainty
                       k    <- rep(ncol(U), N)
               },
               mosaic = {
                       cert <- U[cbind(seq_len(N), x$hard_clusters)]
                       k    <- rep(ncol(U), N)
               })
        data.frame(level     = .UNCERTAINTY_LEVELS[[slot]],
                   certainty = cert,
                   margin    = t2$first - t2$second,
                   k         = as.integer(k),
                   stringsAsFactors = FALSE)
}

#' Stack the supplied levels, drop NA and trivial units, label the panels
#' @keywords internal
#' @noRd
.uncertainty_df <- function(objs, drop_trivial) {
        objs <- objs[!vapply(objs, is.null, logical(1))]
        if (!length(objs))
                stop("Supply at least one of `regions`, `local`, `meta`, `mosaic`.", call. = FALSE)
        stopifnot(is.logical(drop_trivial), length(drop_trivial) == 1L, !is.na(drop_trivial))

        parts <- lapply(names(objs), function(s) .level_uncertainty(objs[[s]], s))
        notes <- character(0)
        parts <- lapply(parts, function(d) {
                lev  <- d$level[1]
                n_na <- sum(is.na(d$certainty) | is.na(d$margin))
                d    <- d[!is.na(d$certainty) & !is.na(d$margin), , drop = FALSE]
                if (n_na) notes <<- c(notes, sprintf("%s: %s units without memberships dropped",
                                                     lev, .fmt_n(n_na)))
                if (drop_trivial) {
                        n_tr <- sum(d$k <= 1L)
                        d    <- d[d$k > 1L, , drop = FALSE]
                        if (n_tr) notes <<- c(notes, sprintf("%s: %s units with a single candidate type dropped",
                                                             lev, .fmt_n(n_tr)))
                }
                if (nrow(d) && all(d$certainty == 1 & d$k > 1L))
                        warning(sprintf("%s memberships are crisp; every unit has certainty 1.", lev),
                                call. = FALSE)
                d
        })
        df <- do.call(rbind, parts)
        if (!nrow(df)) stop("No units left to plot.", call. = FALSE)
        df$level <- factor(df$level, levels = intersect(.UNCERTAINTY_LEVELS, unique(df$level)))
        attr(df, "notes") <- notes
        df
}

#' @keywords internal
#' @noRd
.fmt_n <- function(n) format(n, big.mark = ",", scientific = FALSE, trim = TRUE)

#' @keywords internal
#' @noRd
.k_label <- function(k) {
        r <- range(k)
        if (r[1] == r[2]) sprintf("k = %d", r[1]) else sprintf("k = %d to %d", r[1], r[2])
}

#' Shared histogram / ECDF builder for plot_certainty() and plot_margin()
#' @keywords internal
#' @noRd
.plot_unit_distribution <- function(df, type, bins, floor_df, threshold,
                                    title, xlab, subtitle) {
        stopifnot(is.numeric(bins), length(bins) == 1L, bins >= 2)
        notes <- attr(df, "notes")
        levs  <- levels(df$level)

        stats_df <- do.call(rbind, lapply(levs, function(l) {
                v <- df$value[df$level == l]; k <- df$k[df$level == l]
                lab <- sprintf("%s (%s, n = %s, median %.2f", l, .k_label(k),
                               .fmt_n(length(v)), stats::median(v))
                if (!is.null(threshold))
                        lab <- sprintf("%s, %.1f%% <= %.2f", lab, 100 * mean(v <= threshold), threshold)
                data.frame(level = l, median = stats::median(v), label = paste0(lab, ")"),
                           stringsAsFactors = FALSE)
        }))
        stats_df$level <- factor(stats_df$level, levels = levs)
        labeller <- stats::setNames(stats_df$label, levs)
        caption  <- if (length(notes)) paste(notes, collapse = "\n") else NULL

        if (type == "ecdf") {
                cols <- .UNCERTAINTY_COLS[levs]
                p <- ggplot2::ggplot(df, ggplot2::aes(.data$value, colour = .data$level)) +
                        ggplot2::stat_ecdf(geom = "step", pad = FALSE) +
                        ggplot2::scale_colour_manual(values = cols, labels = labeller, name = NULL) +
                        ggplot2::geom_vline(data = stats_df,
                                            ggplot2::aes(xintercept = .data$median, colour = .data$level),
                                            linewidth = 0.4, show.legend = FALSE) +
                        ggplot2::labs(y = "cumulative share of units") +
                        ggplot2::theme(legend.position = "bottom", legend.direction = "vertical")
        } else {
                p <- ggplot2::ggplot(df, ggplot2::aes(.data$value)) +
                        ggplot2::geom_histogram(binwidth = 1 / bins, boundary = 0,
                                                fill = "#4C72B0", colour = NA) +
                        ggplot2::geom_vline(data = stats_df,
                                            ggplot2::aes(xintercept = .data$median),
                                            colour = "#B2182B", linewidth = 0.5) +
                        ggplot2::facet_wrap(~ level, ncol = 1, scales = "free_y",
                                            labeller = ggplot2::as_labeller(labeller)) +
                        ggplot2::labs(y = "units")
                if (!is.null(floor_df))
                        p <- p + ggplot2::geom_vline(data = floor_df,
                                                     ggplot2::aes(xintercept = .data$x),
                                                     linetype = 3, colour = "grey30")
        }
        if (!is.null(threshold))
                p <- p + ggplot2::geom_vline(xintercept = threshold, linetype = 2, colour = "grey30")

        p + ggplot2::coord_cartesian(xlim = c(0, 1)) +
                ggplot2::labs(title = title, subtitle = subtitle, x = xlab, caption = caption) +
                ggplot2::theme_bw(base_size = 10) +
                ggplot2::theme(strip.text = ggplot2::element_text(hjust = 0),
                               legend.position = if (type == "ecdf") "bottom" else "none",
                               legend.direction = "vertical")
}
