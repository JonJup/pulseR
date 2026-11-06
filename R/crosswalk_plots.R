# -----------------------------------------------------------------------------*
# Figures for typology crosswalks
# -----------------------------------------------------------------------------*

.CW_REL_TAG <- c("1:1" = "=", nested = "<", contains = ">", mixed = "~")

.CW_MAP_COLS <- c(agree_conf    = "#bdbdbd",
                  agree_unc     = "#e6e6e6",
                  disagree_unc  = "#f2a541",
                  disagree_conf = "#5b2a86",
                  none          = "#ffffff")

.CW_MAP_LABS <- c(agree_conf    = "Agree, confident",
                  agree_unc     = "Agree, uncertain",
                  disagree_unc  = "Disagree, uncertain (transition)",
                  disagree_conf = "Disagree, confident (conflict)",
                  none          = "Not compared")

#' Plot a typology crosswalk
#'
#' Figures for the output of \code{\link{crosswalk_types}}. Each answers one
#' question a crosswalk is usually run for:
#' \describe{
#'   \item{\code{"plane"}}{What kind of relationship is this? Each crosswalk is
#'     a point in a unit square. x is the share of the external typology's
#'     information that PULSE does not capture, \eqn{H(E|P)/H(E)}; y is the
#'     share of PULSE information the external typology does not capture,
#'     \eqn{H(P|E)/H(P)}. Bottom-left: same typology. Top-left: PULSE refines
#'     the external typology. Bottom-right: PULSE coarsens it. Top-right:
#'     unrelated. Fill is crisp AMI. Accepts a \code{pulse_crosswalk_set} or a
#'     named list of crosswalks (several external typologies, or one run per
#'     part of the study area). For \code{"within_region"} results the
#'     per-region values are drawn as small grey points.}
#'   \item{\code{"matrix"}}{Which type goes where? Rows are PULSE types,
#'     columns external types, both ordered by the first correspondence
#'     analysis axis so that 1:1 relations form a diagonal and nesting forms
#'     blocks. Fill is the conditional share, square area the joint mass.
#'     Type labels carry their relation: \code{=} 1:1, \code{<} nested,
#'     \code{>} contains, \code{~} mixed.}
#'   \item{\code{"composition"}}{What does each type become? One 100\% bar per
#'     type, split into its \code{top_n} partners and "other", labelled inside
#'     the segments. The lookup view for people who own a typology.}
#'   \item{\code{"map"}}{Where do the typologies disagree, and is it at
#'     uncertain boundaries or in confident cores? Units are classed by
#'     agreement and PULSE certainty; only disagreement is coloured. With
#'     \pkg{patchwork} installed the certainty panel is added on the right.}
#'   \item{\code{"certainty"}}{Match rate by PULSE certainty, with the overall
#'     and chance match rates. A rising curve means disagreement concentrates
#'     where PULSE itself is uncertain.}
#' }
#'
#' @param x A \code{pulse_crosswalk}. For \code{what = "plane"} also a
#'   \code{pulse_crosswalk_set} or a named list of either.
#' @param what Which figure; see Details.
#' @param basis \code{"crisp"} (default) or \code{"fuzzy"} overlap. Fuzzy
#'   shares are membership-weighted shares, not probabilities, and figures say
#'   so in their caption.
#' @param given For \code{"matrix"}: which side the fill is conditioned on.
#'   \code{"auto"} (default) conditions on the typology with more types, so
#'   nesting shows as one dark cell per row or column.
#' @param side For \code{"composition"}: whose types get a bar. Default
#'   \code{"external"}.
#' @param geom For \code{"map"}: an \pkg{sf} object with one row per unit and
#'   an identifier column. Often the same object passed as \code{external}.
#' @param id Identifier column in \code{geom}. Defaults to the \code{id} used in
#'   the crosswalk, then \code{"ID"} or \code{"id"}.
#' @param fill For \code{"map"}: \code{"class"} (default, agreement by
#'   certainty) or \code{"support"} (continuous crosswalk support).
#' @param certainty_cut Rescaled PULSE certainty above which a unit counts as
#'   confident. Certainty is rescaled so that 0 is maximally undecided and 1
#'   is a single type, which keeps the cut comparable across numbers of types.
#'   Default 0.5.
#' @param panel For \code{"map"}: add the certainty panel (needs
#'   \pkg{patchwork}). Default \code{TRUE}.
#' @param dissolve For \code{"map"}: union polygons by class before drawing.
#'   Much faster to draw for large inputs, but the union itself takes time.
#'   Default \code{FALSE}.
#' @param max_types For \code{"matrix"} and \code{"composition"}: types shown
#'   per side. Larger types are kept; in the matrix the rest are merged into one
#'   "other" row or column. Default 40.
#' @param top_n For \code{"composition"}: partners shown per bar. Default 3.
#' @param palette For \code{"composition"}: colours of the partner types, named
#'   by type label or in type order. Pass \code{\link{local_palette}} or
#'   \code{\link{mosaic_palette}} output to match your maps. Default: hues
#'   following the correspondence analysis order, so similar types get similar
#'   hues.
#' @param wording \code{"technical"} (default) or \code{"plain"} axis and
#'   legend text, for audiences without an information theory background.
#' @param label_width Longest type label before truncation. Default 28.
#' @param raster,raster_px Rasterisation control for \code{"map"}, passed to the
#'   package map builder.
#' @param ... Unused.
#' @return A \code{ggplot} object; for \code{what = "map"} with
#'   \code{panel = TRUE} and \pkg{patchwork} installed, a \code{patchwork}.
#' @export
plot_crosswalk <- function(x,
                           what          = c("matrix", "composition", "map",
                                             "certainty", "plane"),
                           basis         = c("crisp", "fuzzy"),
                           given         = c("auto", "pulse", "external"),
                           side          = c("external", "pulse"),
                           geom          = NULL,
                           id            = NULL,
                           fill          = c("class", "support"),
                           certainty_cut = 0.5,
                           panel         = TRUE,
                           dissolve      = FALSE,
                           max_types     = 40L,
                           top_n         = 3L,
                           palette       = NULL,
                           wording       = c("technical", "plain"),
                           label_width   = 28L,
                           raster        = "auto",
                           raster_px     = 2000,
                           ...) {
        .need("ggplot2")
        what    <- match.arg(what)
        basis   <- match.arg(basis)
        given   <- match.arg(given)
        side    <- match.arg(side)
        fill    <- match.arg(fill)
        wording <- match.arg(wording)
        
        if (what == "plane") return(.cw_plot_plane(x, basis, wording))
        if (!inherits(x, "pulse_crosswalk"))
                stop("`x` must be a pulse_crosswalk for what = \"", what, "\". For a set, ",
                     "pick one level, e.g. x$results$local.", call. = FALSE)
        switch(what,
               matrix      = .cw_plot_matrix(x, basis, given, max_types, wording, label_width),
               composition = .cw_plot_composition(x, basis, side, top_n, max_types,
                                                  palette, label_width),
               certainty   = .cw_plot_certainty(x, wording),
               map         = .cw_plot_map(x, geom, id, fill, certainty_cut, panel,
                                          dissolve, wording, raster, raster_px))
}


# -----------------------------------------------------------------------------*
# Shared helpers
# -----------------------------------------------------------------------------*

#' Row and column order from the first correspondence analysis axis
#' @keywords internal
#' @noRd
.cw_ca_order <- function(J) {
        rs <- rowSums(J); cs <- colSums(J)
        out <- list(rows = order(-rs), cols = order(-cs))
        kr <- which(rs > 0); kc <- which(cs > 0)
        if (length(kr) < 2L || length(kc) < 2L) return(out)
        P  <- J[kr, kc, drop = FALSE] / sum(J)
        r  <- rowSums(P); cc <- colSums(P)
        S  <- (P - outer(r, cc)) / sqrt(outer(r, cc))
        sv <- svd(S, nu = 1L, nv = 1L)
        if (!(sv$d[1] > 1e-12)) return(out)
        list(rows = c(kr[order(sv$u[, 1] / sqrt(r))], which(!(rs > 0))),
             cols = c(kc[order(sv$v[, 1] / sqrt(cc))], which(!(cs > 0))))
}

#' @keywords internal
#' @noRd
.cw_trunc <- function(s, w) {
        s <- as.character(s)
        ifelse(nchar(s) > w, paste0(substr(s, 1L, max(w - 3L, 1L)), "..."), s)
}

#' Relation tag per type label, NULL when labels repeat (within-region tables)
#' @keywords internal
#' @noRd
.cw_tags <- function(tab, basis) {
        t <- tab[tab$basis == basis, , drop = FALSE]
        if (anyDuplicated(t$type)) return(NULL)
        stats::setNames(unname(.CW_REL_TAG[t$relation]), t$type)
}

#' @keywords internal
#' @noRd
.cw_provenance <- function(cw, basis) {
        sc <- if (identical(cw$level, "local")) paste0(" (", cw$scope, ")") else ""
        sprintf("%s level%s | %s basis\nweights: %s | N = %s | AMI %s",
                cw$level, sc, basis, cw$settings$weights,
                format(cw$coverage$n_matched, big.mark = ","), cw$settings$ami_method)
}

#' @keywords internal
#' @noRd
.cw_theme <- function(base_size = 10) {
        ggplot2::theme_minimal(base_size = base_size) +
                ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"),
                               plot.title.position = "plot",
                               plot.caption.position = "plot",
                               plot.subtitle = ggplot2::element_text(colour = "grey35"),
                               plot.caption = ggplot2::element_text(colour = "grey35",
                                                                    hjust = 0),
                               panel.grid = ggplot2::element_blank())
}

#' @keywords internal
#' @noRd
.cw_fuzzy_note <- function(basis)
        if (basis == "fuzzy") "Fuzzy basis: membership-weighted shares, not probabilities." else NULL


# -----------------------------------------------------------------------------*
# Matrix
# -----------------------------------------------------------------------------*

#' @keywords internal
#' @noRd
.cw_matrix_data <- function(cw, basis, given, max_types, label_width) {
        J <- cw$overlap[[basis]]
        if (given == "auto") given <- if (nrow(J) >= ncol(J)) "pulse" else "external"
        o <- .cw_ca_order(J)
        pick <- function(ord, mass) {
                big <- order(-mass)[seq_len(min(max_types, sum(mass > 0)))]
                ord[ord %in% big]
        }
        kr <- pick(o$rows, rowSums(J)); kc <- pick(o$cols, colSums(J))
        nr <- nrow(J) - length(kr); nc <- ncol(J) - length(kc)
        tr <- .cw_tags(cw$types$pulse, basis); tc <- .cw_tags(cw$types$external, basis)
        lab <- function(nm, tags) {
                s <- .cw_trunc(nm, label_width)
                if (!is.null(tags)) s <- paste(s, tags[nm])
                make.unique(s, sep = " ")
        }
        J2 <- J[kr, kc, drop = FALSE]
        rl <- lab(rownames(J)[kr], tr); cl <- lab(colnames(J)[kc], tc)
        if (nc > 0) { J2 <- cbind(J2, rowSums(J[kr, -kc, drop = FALSE]))
                      cl <- c(cl, sprintf("other (%d types)", nc)) }
        if (nr > 0) { rest <- J[-kr, , drop = FALSE]
                      oth <- colSums(rest[, kc, drop = FALSE])
                      if (nc > 0) oth <- c(oth, sum(rest[, -kc, drop = FALSE]))
                      J2 <- rbind(J2, oth)
                      rl <- c(rl, sprintf("other (%d types)", nr)) }
        C <- if (given == "pulse") .cw_cond_rows(J2) else t(.cw_cond_rows(t(J2)))
        df <- data.frame(row  = factor(rep(rl, times = length(cl)), levels = rl),
                         col  = factor(rep(cl, each = length(rl)), levels = cl),
                         mass = as.vector(J2) / sum(J),
                         cond = as.vector(C))
        df <- df[df$mass > 0, , drop = FALSE]
        df$side <- sqrt(df$mass / max(df$mass))
        list(df = df, given = given, n_lumped = c(pulse = nr, external = nc))
}

#' @keywords internal
#' @noRd
.cw_plot_matrix <- function(cw, basis, given, max_types, wording, label_width) {
        d  <- .cw_matrix_data(cw, basis, given, max_types, label_width)
        ft <- if (wording == "technical") {
                if (d$given == "pulse") "P(external | PULSE)" else "P(PULSE | external)"
        } else {
                if (d$given == "pulse") "Share of each PULSE type\nin the external type"
                else "Share of each external type\nin the PULSE type"
        }
        cap <- c("Square area: share of the total. Tags: = 1:1, < nested, > contains, ~ mixed.",
                 if (any(d$n_lumped > 0))
                         sprintf("Showing the %d largest types per side; the rest are merged into 'other'.",
                                 max_types),
                 .cw_fuzzy_note(basis))
        ggplot2::ggplot(d$df, ggplot2::aes(.data$col, .data$row)) +
                ggplot2::geom_tile(ggplot2::aes(width = .data$side, height = .data$side,
                                                fill = .data$cond)) +
                ggplot2::scale_fill_gradient(low = "#e8eef5", high = "#08306b",
                                             limits = c(0, 1), name = ft) +
                ggplot2::scale_x_discrete(position = "top") +
                ggplot2::scale_y_discrete(limits = rev(levels(d$df$row))) +
                ggplot2::coord_fixed() +
                ggplot2::labs(x = "External type", y = "PULSE type",
                              title = "Type correspondence",
                              subtitle = .cw_provenance(cw, basis),
                              caption = paste(cap, collapse = "\n")) +
                .cw_theme() +
                ggplot2::theme(axis.text.x.top = ggplot2::element_text(angle = 90, hjust = 0,
                                                                       vjust = 0.5),
                               panel.grid.major = ggplot2::element_line(colour = "grey94",
                                                                        linewidth = 0.3))
}


# -----------------------------------------------------------------------------*
# Composition bars
# -----------------------------------------------------------------------------*

#' @keywords internal
#' @noRd
.cw_composition_data <- function(cw, basis, side, top_n, max_types) {
        J <- cw$overlap[[basis]]
        o <- .cw_ca_order(J)
        if (side == "pulse") {
                bars <- o$rows; partners <- o$cols
        } else {
                J <- t(J); bars <- o$cols; partners <- o$rows
        }
        mass <- rowSums(J)
        bars <- bars[mass[bars] > 0]
        n_all <- length(bars)
        if (n_all > max_types) bars <- bars[bars %in% order(-mass)[seq_len(max_types)]]
        seg <- do.call(rbind, lapply(seq_along(bars), function(b) {
                i  <- bars[b]
                sh <- J[i, ] / mass[i]
                tp <- order(-sh)[seq_len(min(top_n, sum(sh > 0)))]
                s  <- c(sh[tp], other = max(1 - sum(sh[tp]), 0))
                nm <- c(colnames(J)[tp], "other")
                keep <- s > 1e-9
                s <- s[keep]; nm <- nm[keep]
                data.frame(bar = rownames(J)[i], pos = b, partner = nm, share = unname(s),
                           xmax = cumsum(unname(s)), xmin = cumsum(unname(s)) - unname(s),
                           is_other = nm == "other", stringsAsFactors = FALSE)
        }))
        list(seg = seg, bars = rownames(J)[bars], mass = mass[bars],
             partner_order = colnames(J)[partners], n_all = n_all)
}

#' @keywords internal
#' @noRd
.cw_partner_cols <- function(palette, labels) {
        if (is.null(palette)) {
                n <- length(labels)
                return(stats::setNames(grDevices::hcl(h = seq(15, 375, length.out = n + 1L)[seq_len(n)],
                                                      c = 55, l = 68), labels))
        }
        if (!is.null(names(palette)) && all(labels %in% names(palette)))
                return(palette[labels])
        stop("`palette` must be named by type label and cover every partner type.",
             call. = FALSE)
}

#' Black or white text for a fill
#' @keywords internal
#' @noRd
.cw_text_col <- function(cols) {
        L <- grDevices::convertColor(t(grDevices::col2rgb(cols)) / 255, "sRGB", "Luv")[, 1]
        ifelse(L < 55, "white", "grey10")
}

#' @keywords internal
#' @noRd
.cw_plot_composition <- function(cw, basis, side, top_n, max_types, palette, label_width) {
        d <- .cw_composition_data(cw, basis, side, top_n, max_types)
        partner_lab <- if (side == "pulse") colnames(cw$overlap$crisp)
                       else rownames(cw$overlap$crisp)
        ## pulse_palette() names colours "1".."k", not "R1".."Rk": fall back to
        ## position whenever the names do not match but the length does.
        if (!is.null(palette) && !all(partner_lab %in% names(palette)) &&
            length(palette) == length(partner_lab))
                palette <- stats::setNames(palette, partner_lab)
        cols <- c(.cw_partner_cols(palette, d$partner_order), other = "grey88")
        s <- d$seg
        s$fill <- unname(cols[s$partner])
        s$y <- length(d$bars) + 1 - s$pos
        s$txt <- ifelse(s$share >= 0.08 & !s$is_other,
                        paste0(.cw_trunc(s$partner, 12), " ", round(100 * s$share), "%"), "")
        s$txt_col <- .cw_text_col(s$fill)
        tags <- .cw_tags(if (side == "pulse") cw$types$pulse else cw$types$external, basis)
        unit_txt <- if (identical(cw$settings$weights, "unit counts"))
                paste0("n = ", format(round(d$mass), big.mark = ",", trim = TRUE))
        else format(signif(d$mass, 3), big.mark = ",", trim = TRUE)
        ann <- data.frame(y = rev(seq_along(d$bars)),
                          txt = if (is.null(tags)) unit_txt
                                else paste0(unit_txt, "  ", tags[d$bars]))
        whose <- if (side == "pulse") "PULSE" else "external"
        cap <- c(sprintf("Each bar: one %s type, split into its %d largest partners and 'other'.",
                         whose, top_n),
                 "Right: size, then relation (= 1:1, < nested, > contains, ~ mixed).",
                 if (d$n_all > length(d$bars))
                         sprintf("Showing the %d largest of %d types.", length(d$bars), d$n_all),
                 .cw_fuzzy_note(basis))
        ggplot2::ggplot(s) +
                ggplot2::geom_rect(ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax,
                                                ymin = .data$y - 0.4, ymax = .data$y + 0.4,
                                                fill = .data$fill),
                                   colour = "white", linewidth = 0.3) +
                ggplot2::geom_text(ggplot2::aes(x = (.data$xmin + .data$xmax) / 2, y = .data$y,
                                                label = .data$txt, colour = .data$txt_col),
                                   size = 2.6) +
                ggplot2::geom_text(data = ann, ggplot2::aes(x = 1.02, y = .data$y,
                                                            label = .data$txt),
                                   hjust = 0, size = 2.6, colour = "grey30") +
                ggplot2::scale_fill_identity() +
                ggplot2::scale_colour_identity() +
                ggplot2::scale_x_continuous(limits = c(0, 1.3), breaks = seq(0, 1, 0.25),
                                            labels = paste0(seq(0, 100, 25), "%"),
                                            expand = c(0, 0)) +
                ggplot2::scale_y_continuous(breaks = rev(seq_along(d$bars)),
                                            labels = .cw_trunc(d$bars, label_width),
                                            expand = c(0.01, 0.01)) +
                ggplot2::labs(x = NULL, y = paste(if (side == "pulse") "PULSE" else "External",
                                                  "type"),
                              title = sprintf("What each %s type becomes", whose),
                              subtitle = .cw_provenance(cw, basis),
                              caption = paste(cap, collapse = "\n")) +
                .cw_theme()
}


# -----------------------------------------------------------------------------*
# Map and certainty panel
# -----------------------------------------------------------------------------*

#' @keywords internal
#' @noRd
.cw_map_data <- function(cw, geom, id, certainty_cut) {
        if (!inherits(geom, "sf"))
                stop("`geom` must be an sf object with one row per unit.", call. = FALSE)
        if (is.null(id)) id <- cw$settings$id
        if (is.null(id) || !id %in% names(geom)) id <- .id_col(geom)
        if (is.na(id) || !id %in% names(geom))
                stop("No identifier column found in `geom`; set `id`.", call. = FALSE)
        u    <- cw$units
        m    <- match(as.character(geom[[id]]), u$id)
        if (!any(!is.na(m)))
                stop("No `geom` row matches a crosswalk unit; check `id`.", call. = FALSE)
        cert <- .rescale_certainty(u$certainty, nrow(cw$overlap$crisp))
        conf <- cert >= certainty_cut
        cls  <- ifelse(u$match, ifelse(conf, "agree_conf", "agree_unc"),
                       ifelse(conf, "disagree_conf", "disagree_unc"))
        out  <- cls[m]
        out[is.na(m)] <- "none"
        list(class = factor(out, levels = names(.CW_MAP_COLS)), support = u$support[m])
}

#' @keywords internal
#' @noRd
.cw_plot_map <- function(cw, geom, id, fill, certainty_cut, panel, dissolve, wording,
                         raster, raster_px) {
        if (is.null(geom)) stop("`geom` is required for what = \"map\".", call. = FALSE)
        d <- .cw_map_data(cw, geom, id, certainty_cut)
        g <- geom[, attr(geom, "sf_column"), drop = FALSE]
        if (fill == "support") {
                p <- .continuous_map(g, d$support,
                                     title = "Crosswalk support per unit",
                                     limits = c(0, 1), raster = raster, raster_px = raster_px)
        } else {
                if (dissolve) {
                        grp <- split(sf::st_geometry(g), d$class, drop = TRUE)
                        g   <- sf::st_sf(class = factor(names(grp), levels = levels(d$class)),
                                         geometry = do.call(c, lapply(grp, sf::st_union)))
                        cls <- g$class
                } else {
                        cls <- d$class
                        if (nrow(g) > 5e4 && !requireNamespace("terra", quietly = TRUE))
                                message("Drawing ", nrow(g), " polygons. Install 'terra' or ",
                                        "set dissolve = TRUE to speed this up.")
                }
                present <- levels(d$class)[levels(d$class) %in% as.character(d$class)]
                key <- stats::setNames(.CW_MAP_COLS[present], .CW_MAP_LABS[present])
                p <- .map(g, unname(.CW_MAP_COLS[as.character(cls)]),
                          title = "Where the typologies disagree", legend = "none",
                          raster = raster, raster_px = raster_px) +
                        .class_key(key, name = NULL)
        }
        p <- p + ggplot2::labs(subtitle = .cw_provenance(cw, "crisp"),
                               caption = sprintf("Confident: rescaled PULSE certainty >= %.2f.",
                                                 certainty_cut)) +
                ggplot2::theme(legend.position = "bottom",
                               plot.title.position = "plot",
                               plot.caption.position = "plot",
                               plot.caption = ggplot2::element_text(colour = "grey35", hjust = 0),
                               plot.subtitle = ggplot2::element_text(colour = "grey35"))
        if (!panel) return(p)
        if (!requireNamespace("patchwork", quietly = TRUE)) {
                message("Install 'patchwork' to add the certainty panel; returning the map only.")
                return(p)
        }
        p + .cw_plot_certainty(cw, wording) + patchwork::plot_layout(widths = c(3, 1))
}

#' Match rate by rescaled certainty, plus overall and chance rates
#' @keywords internal
#' @noRd
.cw_certainty_data <- function(cw, n_bins = 5L) {
        u    <- cw$units
        cert <- .rescale_certainty(u$certainty, nrow(cw$overlap$crisp))
        br   <- unique(stats::quantile(cert, seq(0, 1, length.out = n_bins + 1L), names = FALSE))
        bin  <- if (length(br) < 3L) rep(1L, nrow(u))
                else cut(cert, br, include.lowest = TRUE, labels = FALSE)
        wm   <- function(v, w) sum(v * w) / max(sum(w), .Machine$double.eps)
        bins <- do.call(rbind, lapply(split(seq_len(nrow(u)), bin), function(i)
                data.frame(certainty = wm(cert[i], u$weight[i]), n_units = length(i),
                           match_rate = wm(u$match[i], u$weight[i]))))
        rownames(bins) <- NULL
        J <- cw$overlap$crisp
        p <- rowSums(J) / sum(J); q <- colSums(J) / sum(J)
        best <- max.col(J, ties.method = "first")
        list(bins = bins, overall = wm(u$match, u$weight),
             chance = sum(p[p > 0] * q[best][p > 0]))
}

#' @keywords internal
#' @noRd
.cw_plot_certainty <- function(cw, wording) {
        d  <- .cw_certainty_data(cw)
        hl <- data.frame(y = c(d$overall, d$chance), what = c("overall", "chance"),
                         lt = c("solid", "dashed"))
        ggplot2::ggplot(d$bins, ggplot2::aes(.data$certainty, .data$match_rate)) +
                ggplot2::geom_hline(data = hl, ggplot2::aes(yintercept = .data$y,
                                                            linetype = .data$lt),
                                    colour = "grey55") +
                ggplot2::geom_text(data = hl, ggplot2::aes(x = 1, y = .data$y, label = .data$what),
                                   hjust = 1, vjust = -0.4, size = 2.8, colour = "grey40") +
                ggplot2::geom_line(colour = "#5b2a86") +
                ggplot2::geom_point(colour = "#5b2a86", size = 2) +
                ggplot2::scale_linetype_identity() +
                ggplot2::scale_x_continuous(limits = c(0, 1)) +
                ggplot2::scale_y_continuous(limits = c(0, 1)) +
                ggplot2::labs(x = if (wording == "technical") "PULSE certainty (rescaled)"
                                  else "How sure PULSE is",
                              y = if (wording == "technical") "Match rate"
                                  else "Share of units where the typologies agree",
                              title = "Agreement by certainty") +
                .cw_theme() +
                ggplot2::theme(panel.grid.major = ggplot2::element_line(colour = "grey92"))
}


# -----------------------------------------------------------------------------*
# Relationship plane
# -----------------------------------------------------------------------------*

#' Points for the relationship plane
#' @keywords internal
#' @noRd
.cw_plane_data <- function(x, basis) {
        items <- if (inherits(x, "pulse_crosswalk")) list(x)
                 else if (inherits(x, "pulse_crosswalk_set")) x$results
                 else if (is.list(x)) x
                 else stop("`x` must be a pulse_crosswalk, a pulse_crosswalk_set or a ",
                           "named list of them.", call. = FALSE)
        if (inherits(x, "pulse_crosswalk")) names(items) <- x$level
        if (is.null(names(items)) || any(!nzchar(names(items))))
                stop("A list for what = \"plane\" must be named.", call. = FALSE)
        flat <- list()
        for (nm in names(items)) {
                it <- items[[nm]]
                if (inherits(it, "pulse_crosswalk_set")) {
                        for (lv in names(it$results)) flat[[paste0(nm, ": ", lv)]] <- it$results[[lv]]
                } else if (inherits(it, "pulse_crosswalk")) {
                        flat[[nm]] <- it
                } else stop("List element `", nm, "` is not a crosswalk.", call. = FALSE)
        }
        coords <- function(get) {
                he <- get("H_ext"); hp <- get("H_pulse")
                c(x = if (he > 0) get("H_ext_given_pulse") / he else NA_real_,
                  y = if (hp > 0) get("H_pulse_given_ext") / hp else NA_real_)
        }
        pts <- do.call(rbind, lapply(names(flat), function(nm) {
                a  <- flat[[nm]]$agreement
                xy <- coords(function(ix) a[[basis]][a$index == ix])
                data.frame(name = nm, x = xy[["x"]], y = xy[["y"]],
                           AMI = a$crisp[a$index == "AMI"],
                           n_types = nrow(flat[[nm]]$overlap$crisp), stringsAsFactors = FALSE)
        }))
        sat <- do.call(rbind, lapply(names(flat), function(nm) {
                b <- flat[[nm]]$by_region
                if (is.null(b)) return(NULL)
                do.call(rbind, lapply(split(b, b$region), function(r) {
                        xy <- coords(function(ix) r[[basis]][r$index == ix])
                        data.frame(name = nm, region = r$region[1], x = xy[["x"]], y = xy[["y"]])
                }))
        }))
        list(points = pts, satellites = sat)
}

#' @keywords internal
#' @noRd
.cw_plot_plane <- function(x, basis, wording) {
        d <- .cw_plane_data(x, basis)
        d$points$hjust <- ifelse(d$points$x > 0.8, 1, ifelse(d$points$x < 0.2, 0, 0.5))
        tech <- wording == "technical"
        corners <- data.frame(
                x = c(0.03, 0.03, 0.97, 0.97), y = c(0.03, 0.97, 0.03, 0.97),
                hjust = c(0, 0, 1, 1), vjust = c(0, 1, 0, 1),
                txt = c("same typology", "PULSE refines", "PULSE coarsens", "unrelated"))
        amin <- min(0, d$points$AMI, na.rm = TRUE)
        p <- ggplot2::ggplot(d$points, ggplot2::aes(.data$x, .data$y)) +
                ggplot2::annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = 1,
                                  fill = NA, colour = "grey75") +
                ggplot2::geom_text(data = corners,
                                   ggplot2::aes(label = .data$txt, hjust = .data$hjust,
                                                vjust = .data$vjust),
                                   colour = "grey55", fontface = "italic", size = 3)
        if (!is.null(d$satellites))
                p <- p + ggplot2::geom_point(data = d$satellites, colour = "grey60", size = 1.1)
        p +
                ggplot2::geom_point(ggplot2::aes(fill = .data$AMI, size = .data$n_types),
                                    shape = 21, colour = "grey15", stroke = 0.4) +
                ggplot2::geom_text(ggplot2::aes(label = .data$name, hjust = .data$hjust),
                                   vjust = -1.3, size = 3.1) +
                ggplot2::scale_fill_viridis_c(option = "cividis", limits = c(amin, 1),
                                              oob = function(v, ...) pmin(pmax(v, amin), 1),
                                              name = if (tech) "AMI (crisp)"
                                                     else "Agreement\nbeyond chance") +
                ggplot2::scale_size(range = c(2.5, 7), name = "Number of\nPULSE types") +
                ggplot2::coord_fixed(xlim = c(-0.02, 1.02), ylim = c(-0.02, 1.08),
                                     expand = FALSE, clip = "off") +
                ggplot2::labs(x = if (tech) "H(external | PULSE) / H(external)"
                                  else "External information PULSE does not capture",
                              y = if (tech) "H(PULSE | external) / H(PULSE)"
                                  else "PULSE information the external typology does not capture",
                              title = "How the typologies relate",
                              subtitle = sprintf("%s basis; fill: crisp AMI%s", basis,
                                                 if (is.null(d$satellites)) ""
                                                 else "; grey points: single regions"),
                              caption = .cw_fuzzy_note(basis)) +
                .cw_theme() +
                ggplot2::theme(panel.grid.major = ggplot2::element_line(colour = "grey94"))
}
