# -----------------------------------------------------------------------------*
# Crosswalk between a pulseR typology and an external typology
# -----------------------------------------------------------------------------*

.CW_CLASSES <- c("pulse_regions", "regions_tuning", "pulse_local_types",
                 "pulse_meta_types", "mosaic_types")

.CW_INDICES <- c("ARI", "AMI", "NMI", "VI", "H_ext_given_pulse",
                 "H_pulse_given_ext", "purity", "inverse_purity", "cramers_v")

#' Crosswalk a pulseR typology against an external typology
#'
#' Cross-tabulates a pulseR typology (any level) against an existing typology
#' defined on the same spatial units, such as a national WFD river typology,
#' and summarises how the two correspond: overall agreement, which types map
#' onto which, and where on the map they disagree.
#'
#' Both typologies are turned into unit-by-type membership matrices (a crisp
#' label becomes a one-hot row) and every summary is derived from the weighted
#' joint overlap \eqn{J = U^\top W V}. The crisp result applies the same
#' formula to the hard labels of both sides; for pulseR these are the labels
#' stored on the object (the consensus partition at level 1), not the argmax
#' of the memberships.
#'
#' @section Interpreting fuzzy results:
#' pulseR memberships describe position in ensemble or signature space; they
#' are not calibrated probabilities. The fuzzy overlap is therefore a
#' membership-weighted overlap, not an expected co-classification. Chance
#' correction has no standard form for soft contingency tables, so ARI and
#' AMI are reported for the crisp basis only. Read crisp AMI as the headline
#' and the fuzzy indices as a check on how much boundary units change it.
#'
#' @section Level 2 scope:
#' Local types are defined within regions. With \code{local_scope = "global"}
#' every region-by-type combination is its own class, so agreement includes
#' whatever the level 1 regions already explain. \code{"within_region"}
#' computes the indices within each region and returns their mass-weighted
#' means, which isolates level 2. To compare pooled river types with an
#' external river typology, pass a \code{pulse_meta_types} object.
#'
#' @param x A \code{pulse_regions}, \code{regions_tuning},
#'   \code{pulse_local_types}, \code{pulse_meta_types} or \code{mosaic_types}
#'   object, or a named list of such objects (for example
#'   \code{list(regions = reg, local = loc, mosaic = mos)}) to crosswalk every
#'   level in one call.
#' @param external An \pkg{sf} object or data frame with one row per unit.
#'   Geometry is ignored; units are matched by \code{id}.
#' @param label Name of the column holding the external type (crisp), or a
#'   character vector of column names holding external memberships (fuzzy).
#'   Fuzzy rows are renormalised to sum to one, with a warning if needed.
#' @param units Unit identifiers of \code{x} in graph vertex order: either a
#'   vector, or the \pkg{igraph} object returned by
#'   \code{polygon_to_network()}, whose \code{polygon_id} vertex attribute is
#'   used. Rows are never matched by position.
#' @param id Name of the identifier column in \code{external}. \code{NULL}
#'   (default) uses \code{"ID"} or \code{"id"}. Identifiers are compared as
#'   character strings.
#' @param weights \code{NULL} (default; every unit counts once) or the name of
#'   a non-negative numeric column in \code{external}, such as catchment area.
#'   Weights enter the overlap tables and information-based indices. ARI and
#'   AMI always use unit counts, because chance correction needs integer counts.
#' @param local_scope \code{"global"} (default) or \code{"within_region"}.
#'   Used for \code{pulse_local_types} only; see the Level 2 scope section.
#' @param nesting Share in \eqn{(0.5, 1]} above which a type counts as lying
#'   within one type of the other typology. Default 0.7.
#' @param ami_max_terms Largest number of hypergeometric terms for the exact
#'   expected mutual information. Above it the asymptotic approximation
#'   \eqn{(K-1)(L-1)/(2N)} is used. The method used is reported in
#'   \code{settings}. Default \code{5e7}.
#'
#' @return For a single object, a list of class \code{"pulse_crosswalk"}:
#'   \describe{
#'     \item{\code{level}, \code{scope}}{Level crosswalked and the local scope.}
#'     \item{\code{overlap}}{K x L matrices (PULSE types in rows): \code{fuzzy}
#'       and \code{crisp} weighted overlap, and unweighted crisp \code{counts}.}
#'     \item{\code{conditional}}{Row and column conditionals of both overlaps:
#'       \code{ext_given_pulse} (what a PULSE type contains) and
#'       \code{pulse_given_ext} (where an external type goes).}
#'     \item{\code{agreement}}{Data frame with columns \code{index},
#'       \code{crisp} and \code{fuzzy}. Indices: ARI, AMI, NMI (arithmetic
#'       normalisation), VI and its two parts \code{H_ext_given_pulse} (zero
#'       when PULSE types nest in external types) and \code{H_pulse_given_ext}
#'       (zero when external types nest in PULSE types), both in nats,
#'       \code{purity}, \code{inverse_purity} and Cramer's V.}
#'     \item{\code{by_region}}{Per-region indices (\code{"within_region"}
#'       only), otherwise \code{NULL}.}
#'     \item{\code{types}}{List with \code{pulse} and \code{external} data
#'       frames, one row per type and \code{basis} (fuzzy or crisp): mass,
#'       share, best match, share in the best match, effective number of
#'       partner types, and relation (\code{"1:1"}, \code{"nested"},
#'       \code{"contains"} or \code{"mixed"}).}
#'     \item{\code{units}}{One row per matched unit: \code{id}, both crisp
#'       labels, \code{match} (the external label equals the crisp best match
#'       of the unit's PULSE type), \code{support} (how strongly the fuzzy
#'       crosswalk predicts the unit's external label), \code{certainty}
#'       (row maximum of the PULSE memberships) and \code{weight}. Join it to
#'       your polygons by \code{id} to map disagreement.}
#'     \item{\code{by_certainty}}{Match rate and mean support by quintile of
#'       \code{certainty}.}
#'     \item{\code{coverage}, \code{settings}}{Matching counts and the
#'       settings used.}
#'   }
#'   For a list input, a \code{"pulse_crosswalk_set"} with \code{results} (one
#'   \code{pulse_crosswalk} per element) and \code{by_level}, a comparison
#'   table. Compare levels on AMI, not NMI: NMI rises with the number of types.
#'
#' @importFrom Matrix sparseMatrix crossprod
#' @export
crosswalk_types <- function(x, external, label, units, id = NULL, weights = NULL,
                            local_scope = c("global", "within_region"),
                            nesting = 0.7, ami_max_terms = 5e7) {
        local_scope <- match.arg(local_scope)
        if (!is.numeric(nesting) || length(nesting) != 1L || !(nesting > 0.5) || nesting > 1)
                stop("`nesting` must be a single number in (0.5, 1].", call. = FALSE)
        
        if (is.list(x) && !inherits(x, .CW_CLASSES))
                return(.crosswalk_set(x, external, label, units, id, weights,
                                      local_scope, nesting, ami_max_terms))
        
        side <- .cw_pulse_side(x)
        if (local_scope == "within_region" && side$level != "local") {
                warning("`local_scope` applies to pulse_local_types only; ignored.",
                        call. = FALSE)
                local_scope <- "global"
        }
        ids <- .cw_unit_ids(units, length(side$row_of_unit))
        ext <- .cw_external(external, label, id, weights)
        
        # --- match units ---------------------------------------------------*
        pos   <- match(ids, ext$id)
        typed <- !is.na(side$row_of_unit)
        keep  <- typed & !is.na(pos)
        in_pulse <- ext$id %in% ids
        coverage <- list(
                n_pulse_units        = length(ids),
                n_typed              = sum(typed),
                n_external_rows      = ext$n_rows,
                n_external_usable    = length(ext$id),
                n_matched            = sum(keep),
                n_typed_unmatched    = sum(typed & is.na(pos)),
                n_external_unmatched = sum(!in_pulse),
                n_no_pulse_type      = sum(!typed & !is.na(pos)),
                share_matched        = sum(keep) / max(sum(typed), 1L))
        if (sum(keep) < 2L)
                stop("Fewer than two units matched. Check that `units` and `external[[id]]` ",
                     "use the same identifiers.", call. = FALSE)
        if (coverage$share_matched < 0.95)
                warning(sprintf(paste0("Only %.1f%% of typed PULSE units matched an external ",
                                       "unit with a usable label; see `coverage`."),
                                100 * coverage$share_matched), call. = FALSE)
        
        U  <- side$U
        r  <- side$row_of_unit[keep]
        e  <- pos[keep]
        n  <- length(r)
        K  <- ncol(U)
        L  <- length(ext$labels)
        w  <- ext$w[e]
        hp <- side$hard[r]
        he <- ext$hard[e]
        
        rs <- rowSums(U)[r]
        if (any(!(rs > 0)))
                stop("Some matched units have all-zero PULSE memberships.", call. = FALSE)
        if (any(abs(rs - 1) > 1e-6))
                warning("PULSE membership rows do not sum to one; renormalised.", call. = FALSE)
        
        # --- overlap tables ------------------------------------------------*
        ## Weights and renormalisation go into the narrow external side, so the
        ## N x K membership matrix is never copied.
        wf <- w / rs
        if (is.null(ext$V)) {
                Se <- Matrix::sparseMatrix(i = r, j = he, x = wf, dims = c(nrow(U), L))
        } else {
                Se <- matrix(0, nrow(U), L)
                Se[r, ] <- ext$V[e, , drop = FALSE] * wf
        }
        Jf <- as.matrix(Matrix::crossprod(U, Se))
        Hp <- Matrix::sparseMatrix(i = seq_len(n), j = hp, x = 1, dims = c(n, K))
        Jc <- as.matrix(Matrix::crossprod(Hp, Matrix::sparseMatrix(i = seq_len(n), j = he,
                                                                   x = w, dims = c(n, L))))
        Cn <- as.matrix(Matrix::crossprod(Hp, Matrix::sparseMatrix(i = seq_len(n), j = he,
                                                                   x = 1, dims = c(n, L))))
        dn <- list(pulse = side$labels, external = ext$labels)
        dimnames(Jf) <- dimnames(Jc) <- dimnames(Cn) <- dn
        
        # --- agreement and type relations ----------------------------------*
        by_region <- NULL
        if (local_scope == "within_region") {
                wr <- .cw_within_region(Jf, Jc, Cn, side$region_of_type, nesting, ami_max_terms)
                agreement <- wr$agreement
                by_region <- wr$by_region
                types     <- wr$types
        } else {
                agreement <- .cw_agreement(Jf, Jc, Cn, ami_max_terms)
                types <- list(
                        pulse    = rbind(.cw_types(Jf, nesting, "fuzzy")$rows,
                                         .cw_types(Jc, nesting, "crisp")$rows),
                        external = rbind(.cw_types(Jf, nesting, "fuzzy")$cols,
                                         .cw_types(Jc, nesting, "crisp")$cols))
        }
        
        # --- unit table ----------------------------------------------------*
        Pf <- .cw_cond_rows(Jf)
        M  <- U %*% Pf                                   # nrow(U) x L, small
        support <- if (is.null(ext$V)) M[cbind(r, he)] / rs
                   else rowSums(M[r, , drop = FALSE] * ext$V[e, , drop = FALSE]) / rs
        certainty <- U[cbind(r, max.col(U, ties.method = "first")[r])] / rs
        best_c <- max.col(Jc, ties.method = "first")
        best_c[rowSums(Jc) <= 0] <- NA_integer_
        units_df <- data.frame(id            = ids[keep],
                               pulse_type    = side$labels[hp],
                               external_type = ext$labels[he],
                               match         = !is.na(best_c[hp]) & best_c[hp] == he,
                               support       = unname(support),
                               certainty     = unname(certainty),
                               weight        = w,
                               stringsAsFactors = FALSE)
        
        out <- list(
                level        = side$level,
                scope        = if (side$level == "local") local_scope else NA_character_,
                overlap      = list(fuzzy = Jf, crisp = Jc, counts = Cn),
                conditional  = list(ext_given_pulse = list(fuzzy = Pf,
                                                           crisp = .cw_cond_rows(Jc)),
                                    pulse_given_ext = list(fuzzy = t(.cw_cond_rows(t(Jf))),
                                                           crisp = t(.cw_cond_rows(t(Jc))))),
                agreement    = agreement,
                by_region    = by_region,
                types        = types,
                units        = units_df,
                by_certainty = .cw_by_certainty(units_df),
                coverage     = coverage,
                settings     = list(label = label, id = ext$id_col,
                                    external_fuzzy = !is.null(ext$V),
                                    weights = if (is.null(weights)) "unit counts" else weights,
                                    nesting = nesting,
                                    ami_method = attr(agreement, "ami_method"),
                                    chance_corrected_on = "unit counts"))
        structure(out, class = "pulse_crosswalk")
}


# -----------------------------------------------------------------------------*
# Inputs
# -----------------------------------------------------------------------------*

#' PULSE side of a crosswalk: memberships, hard labels, unit-to-row map
#' @keywords internal
#' @noRd
.cw_pulse_side <- function(x) {
        if (inherits(x, "regions_tuning")) {
                x <- x$best_result
                if (is.null(x)) stop("`regions_tuning` object carries no `$best_result`.",
                                     call. = FALSE)
        }
        region_of_type <- NULL
        if (inherits(x, "pulse_regions")) {
                U     <- .first_non_null(x$memberships, x$fit$memberships)
                hard  <- .first_non_null(x$hard_clusters, x$fit$hard_clusters)
                level <- "regions"
                pref  <- "R"
        } else if (inherits(x, "pulse_local_types")) {
                U     <- x$memberships
                hard  <- x$hard_classes
                level <- "local"
                region_of_type <- as.integer(x$key$region)
        } else if (inherits(x, "pulse_meta_types")) {
                U     <- x$memberships
                hard  <- x$hard_classes
                level <- "meta"
                pref  <- "meta"
        } else if (inherits(x, "mosaic_types")) {
                U     <- x$memberships
                hard  <- x$hard_clusters
                level <- "mosaic"
                pref  <- "M"
        } else {
                stop("`x` must be one of: ", paste(.CW_CLASSES, collapse = ", "),
                     ", or a named list of them.", call. = FALSE)
        }
        if (is.null(U) || is.null(hard))
                stop("Could not find memberships and hard labels on `x`.", call. = FALSE)
        U    <- as.matrix(U)
        hard <- as.integer(hard)
        K    <- ncol(U)
        if (length(hard) != nrow(U))
                stop("Hard labels and membership rows differ in length.", call. = FALSE)
        if (anyNA(hard) || any(hard < 1L | hard > K))
                stop("Hard labels must index membership columns 1..", K, ".", call. = FALSE)
        
        labels <- if (level == "local") as.character(x$key$name) else colnames(U)
        if (is.null(labels) || length(labels) != K) labels <- paste0(pref, seq_len(K))
        
        if (level == "mosaic") {
                valid <- x$valid
                if (is.null(valid) || sum(valid) != nrow(U))
                        stop("`mosaic_types` object has no `valid` vector matching its ",
                             "membership rows.", call. = FALSE)
                row_of_unit <- rep(NA_integer_, length(valid))
                row_of_unit[valid] <- seq_len(nrow(U))
        } else {
                row_of_unit <- seq_len(nrow(U))
        }
        list(U = U, hard = hard, labels = labels, level = level,
             row_of_unit = row_of_unit, region_of_type = region_of_type)
}

#' Unit identifiers in graph vertex order
#' @keywords internal
#' @noRd
.cw_unit_ids <- function(units, n) {
        if (missing(units) || is.null(units))
                stop("`units` is required: the unit IDs of `x` in graph vertex order, ",
                     "or the graph from polygon_to_network().", call. = FALSE)
        if (inherits(units, "igraph")) {
                ids <- igraph::vertex_attr(units, "polygon_id")
                if (is.null(ids))
                        stop("The graph has no `polygon_id` vertex attribute. Pass the ",
                             "ID vector in graph vertex order instead.", call. = FALSE)
        } else if (is.atomic(units)) {
                ids <- units
        } else {
                stop("`units` must be a vector of IDs or an igraph object.", call. = FALSE)
        }
        if (length(ids) != n)
                stop(sprintf("`units` has %d IDs but `x` has %d units.", length(ids), n),
                     call. = FALSE)
        ids <- as.character(ids)
        if (anyNA(ids)) stop("`units` contains NA IDs.", call. = FALSE)
        if (anyDuplicated(ids)) stop("`units` contains duplicated IDs.", call. = FALSE)
        ids
}

#' External typology: IDs, hard labels, optional memberships, weights
#' @keywords internal
#' @noRd
.cw_external <- function(external, label, id, weights) {
        if (!is.data.frame(external))
                stop("`external` must be an sf object or a data frame.", call. = FALSE)
        df <- if (inherits(external, "sf")) sf::st_drop_geometry(external)
              else as.data.frame(external)
        if (is.null(id)) {
                id <- .id_col(df)
                if (is.na(id)) stop("No `ID` or `id` column in `external`; set `id`.",
                                    call. = FALSE)
        }
        if (!is.character(label) || !length(label))
                stop("`label` must be one or more column names.", call. = FALSE)
        need <- c(id, label, weights)
        miss <- setdiff(need, names(df))
        if (length(miss))
                stop("Column(s) not found in `external`: ", paste(miss, collapse = ", "),
                     call. = FALSE)
        
        ids <- as.character(df[[id]])
        d <- unique(ids[!is.na(ids) & duplicated(ids)])
        if (length(d))
                stop(sprintf("%d duplicated ID(s) in `external`, e.g. %s. One row per unit is required.",
                             length(d), paste(d[seq_len(min(3L, length(d)))], collapse = ", ")), call. = FALSE)
        
        w <- if (is.null(weights)) rep(1, nrow(df)) else df[[weights]]
        if (!is.numeric(w)) stop("`weights` column must be numeric.", call. = FALSE)
        if (any(w < 0, na.rm = TRUE)) stop("`weights` must be non-negative.", call. = FALSE)
        
        V <- NULL
        if (length(label) == 1L) {
                lab <- df[[label]]
                lv  <- if (is.factor(lab)) levels(droplevels(lab))
                       else sort(unique(as.character(lab[!is.na(lab)])))
                hard <- match(as.character(lab), lv)
                ok   <- !is.na(hard)
        } else {
                V <- as.matrix(df[label])
                if (!is.numeric(V)) stop("Fuzzy `label` columns must be numeric.", call. = FALSE)
                if (any(V < 0, na.rm = TRUE))
                        stop("Fuzzy `label` columns must be non-negative.", call. = FALSE)
                lv  <- label
                vs  <- rowSums(V)
                ok  <- stats::complete.cases(V) & vs > 0
                if (any(abs(vs[ok] - 1) > 1e-6))
                        warning("External membership rows do not sum to one; renormalised.",
                                call. = FALSE)
                V[ok, ] <- V[ok, , drop = FALSE] / vs[ok]
                hard <- rep(NA_integer_, nrow(V))
                hard[ok] <- max.col(V[ok, , drop = FALSE], ties.method = "first")
        }
        ok <- ok & !is.na(ids) & !is.na(w)
        list(id = ids[ok], hard = hard[ok],
             V = if (is.null(V)) NULL else V[ok, , drop = FALSE],
             labels = lv, w = as.numeric(w[ok]), n_rows = nrow(df), id_col = id)
}


# -----------------------------------------------------------------------------*
# Indices
# -----------------------------------------------------------------------------*

#' @keywords internal
#' @noRd
.cw_h <- function(p) { p <- p[p > 0]; -sum(p * log(p)) }

#' Row-conditional table; all-zero rows stay zero
#' @keywords internal
#' @noRd
.cw_cond_rows <- function(J) {
        rs <- rowSums(J)
        J / ifelse(rs > 0, rs, 1)
}

#' Information-based indices of an overlap table (weights allowed)
#' @keywords internal
#' @noRd
.cw_info <- function(J) {
        nm  <- .CW_INDICES[-(1:2)]
        tot <- sum(J)
        if (!(tot > 0)) return(stats::setNames(rep(NA_real_, length(nm)), nm))
        P  <- J / tot
        pr <- rowSums(P); pc <- colSums(P)
        Hr <- .cw_h(pr); Hc <- .cw_h(pc); Hj <- .cw_h(P)
        MI <- max(Hr + Hc - Hj, 0)
        kr <- pr > 0; kc <- pc > 0
        q  <- min(sum(kr), sum(kc)) - 1L
        E  <- outer(pr[kr], pc[kc])
        cv <- if (q > 0) sqrt(sum((P[kr, kc, drop = FALSE] - E)^2 / E) / q) else NA_real_
        c(NMI               = if (Hr + Hc > 0) MI / ((Hr + Hc) / 2) else 1,
          VI                = max(2 * Hj - Hr - Hc, 0),
          H_ext_given_pulse = max(Hj - Hr, 0),
          H_pulse_given_ext = max(Hj - Hc, 0),
          purity            = sum(apply(P, 1, max)),
          inverse_purity    = sum(apply(P, 2, max)),
          cramers_v         = cv)
}

#' Adjusted Rand index from a count table
#' @keywords internal
#' @noRd
.cw_ari <- function(C) {
        n <- sum(C)
        if (n < 2) return(NA_real_)
        s  <- sum(choose(C, 2))
        sa <- sum(choose(rowSums(C), 2))
        sb <- sum(choose(colSums(C), 2))
        ex <- sa * sb / choose(n, 2)
        mx <- (sa + sb) / 2
        if (mx - ex == 0) return(1)
        (s - ex) / (mx - ex)
}

#' Exact expected mutual information under the permutation model
#' (Vinh, Epps and Bailey 2010), in nats.
#' @keywords internal
#' @noRd
.cw_emi_exact <- function(a, b, N) {
        lf  <- function(z) lgamma(z + 1)
        lfN <- lf(N)
        emi <- 0
        for (i in seq_along(a)) {
                ai <- a[i]
                for (j in seq_along(b)) {
                        bj <- b[j]
                        lo <- max(1, ai + bj - N); hi <- min(ai, bj)
                        if (lo > hi) next
                        nij <- lo:hi
                        lp  <- lf(ai) + lf(bj) + lf(N - ai) + lf(N - bj) - lfN -
                                lf(nij) - lf(ai - nij) - lf(bj - nij) - lf(N - ai - bj + nij)
                        emi <- emi + sum(nij / N * log(N * nij / (ai * bj)) * exp(lp))
                }
        }
        emi
}

#' ARI and AMI from an unweighted count table
#' @keywords internal
#' @noRd
.cw_chance <- function(C, max_terms) {
        C <- C[rowSums(C) > 0, colSums(C) > 0, drop = FALSE]
        N <- sum(C)
        if (N < 2) return(list(ARI = NA_real_, AMI = NA_real_, method = "none"))
        a  <- rowSums(C); b <- colSums(C)
        Ha <- .cw_h(a / N); Hb <- .cw_h(b / N)
        MI <- max(Ha + Hb - .cw_h(C / N), 0)
        terms  <- sum(pmax(outer(a, b, pmin) - pmax(outer(a, b, "+") - N, 1) + 1, 0))
        method <- if (terms <= max_terms) "exact" else "approx"
        emi <- if (method == "exact") .cw_emi_exact(a, b, N)
               else (length(a) - 1) * (length(b) - 1) / (2 * N)
        den <- (Ha + Hb) / 2 - emi
        ami <- if (length(a) == 1L && length(b) == 1L) 1
               else if (abs(den) < 1e-12) NA_real_
               else (MI - emi) / den
        list(ARI = .cw_ari(C), AMI = ami, method = method)
}

#' Agreement table for one overlap
#' @keywords internal
#' @noRd
.cw_agreement <- function(Jf, Jc, Cn, max_terms) {
        ch  <- .cw_chance(Cn, max_terms)
        out <- data.frame(index = .CW_INDICES,
                          crisp = c(ch$ARI, ch$AMI, unname(.cw_info(Jc))),
                          fuzzy = c(NA_real_, NA_real_, unname(.cw_info(Jf))),
                          stringsAsFactors = FALSE)
        attr(out, "ami_method") <- ch$method
        out
}

#' Per-type relations for both sides of an overlap table
#' @keywords internal
#' @noRd
.cw_types <- function(J, nesting, basis) {
        K <- nrow(J); L <- ncol(J)
        rs <- rowSums(J); cs <- colSums(J); tot <- sum(J)
        Pr <- .cw_cond_rows(J)
        Pc <- t(.cw_cond_rows(t(J)))
        br <- max.col(J, ties.method = "first");    br[!(rs > 0)] <- NA
        bc <- max.col(t(J), ties.method = "first"); bc[!(cs > 0)] <- NA
        pr <- ifelse(is.na(br), NA_real_, Pr[cbind(seq_len(K), ifelse(is.na(br), 1L, br))])
        pc <- ifelse(is.na(bc), NA_real_, Pc[cbind(ifelse(is.na(bc), 1L, bc), seq_len(L))])
        rn <- !is.na(pr) & pr >= nesting           # row type lies within one column type
        cn <- !is.na(pc) & pc >= nesting
        r_cont <- tabulate(bc[cn], K) > 0          # some column type lies within it
        c_cont <- tabulate(br[rn], L) > 0
        bri <- ifelse(is.na(br), 1L, br)
        bci <- ifelse(is.na(bc), 1L, bc)
        r_11 <- rn & cn[bri] & !is.na(bc[bri]) & bc[bri] == seq_len(K)
        c_11 <- cn & rn[bci] & !is.na(br[bci]) & br[bci] == seq_len(L)
        rel <- function(one, nest, cont, mass)
                ifelse(!(mass > 0), NA_character_,
                       ifelse(one, "1:1", ifelse(nest, "nested",
                                                  ifelse(cont, "contains", "mixed"))))
        eff <- function(P) exp(apply(P, 1, .cw_h))
        rows <- data.frame(type = rownames(J), basis = basis, mass = rs,
                           share = rs / tot, best_match = colnames(J)[br],
                           p_best = pr, n_partners = eff(Pr),
                           relation = rel(r_11, rn, r_cont, rs),
                           stringsAsFactors = FALSE, row.names = NULL)
        cols <- data.frame(type = colnames(J), basis = basis, mass = cs,
                           share = cs / tot, best_match = rownames(J)[bc],
                           p_best = pc, n_partners = eff(t(Pc)),
                           relation = rel(c_11, cn, c_cont, cs),
                           stringsAsFactors = FALSE, row.names = NULL)
        rows$n_partners[!(rs > 0)] <- NA; cols$n_partners[!(cs > 0)] <- NA
        list(rows = rows, cols = cols)
}

#' Within-region agreement and relations for local types
#' @keywords internal
#' @noRd
.cw_within_region <- function(Jf, Jc, Cn, region_of_type, nesting, max_terms) {
        regs <- sort(unique(region_of_type))
        per <- lapply(regs, function(g) {
                cols <- which(region_of_type == g)
                a <- .cw_agreement(Jf[cols, , drop = FALSE], Jc[cols, , drop = FALSE],
                                   Cn[cols, , drop = FALSE], max_terms)
                br <- data.frame(region = g, a, n_units = sum(Cn[cols, ]),
                                 mass_crisp = sum(Jc[cols, ]), mass_fuzzy = sum(Jf[cols, ]),
                                 ami_method = attr(a, "ami_method"),
                                 stringsAsFactors = FALSE)
                tf <- .cw_types(Jf[cols, , drop = FALSE], nesting, "fuzzy")
                tc <- .cw_types(Jc[cols, , drop = FALSE], nesting, "crisp")
                list(br = br,
                     pulse = cbind(region = g, rbind(tf$rows, tc$rows)),
                     external = cbind(region = g, rbind(tf$cols, tc$cols)))
        })
        by_region <- do.call(rbind, lapply(per, `[[`, "br"))
        wmean <- function(v, w) { ok <- !is.na(v) & w > 0
                if (!any(ok)) NA_real_ else sum(v[ok] * w[ok]) / sum(w[ok]) }
        agreement <- do.call(rbind, lapply(.CW_INDICES, function(ix) {
                s  <- by_region[by_region$index == ix, ]
                wc <- if (ix %in% c("ARI", "AMI")) s$n_units else s$mass_crisp
                data.frame(index = ix, crisp = wmean(s$crisp, wc),
                           fuzzy = wmean(s$fuzzy, s$mass_fuzzy), stringsAsFactors = FALSE)
        }))
        m <- unique(by_region$ami_method)
        attr(agreement, "ami_method") <- paste(m[m != "none"], collapse = "+")
        ext <- do.call(rbind, lapply(per, `[[`, "external"))
        ext <- ext[ext$mass > 0, , drop = FALSE]
        rownames(ext) <- NULL
        list(agreement = agreement, by_region = by_region,
             types = list(pulse = do.call(rbind, lapply(per, `[[`, "pulse")),
                          external = ext))
}

#' Match rate and support by quintile of PULSE certainty
#' @keywords internal
#' @noRd
.cw_by_certainty <- function(u) {
        br <- unique(stats::quantile(u$certainty, seq(0, 1, 0.2), names = FALSE))
        bin <- if (length(br) < 3L) factor(rep("all", nrow(u)))
               else cut(u$certainty, br, include.lowest = TRUE, dig.lab = 3)
        out <- do.call(rbind, lapply(split(seq_len(nrow(u)), bin), function(i) {
                w <- u$weight[i]
                data.frame(certainty = as.character(bin[i[1]]), n_units = length(i),
                           match_rate = sum(w * u$match[i]) / max(sum(w), .Machine$double.eps),
                           mean_support = sum(w * u$support[i]) / max(sum(w), .Machine$double.eps),
                           stringsAsFactors = FALSE)
        }))
        rownames(out) <- NULL
        out
}


# -----------------------------------------------------------------------------*
# All levels at once
# -----------------------------------------------------------------------------*

#' @keywords internal
#' @noRd
.crosswalk_set <- function(x, external, label, units, id, weights,
                           local_scope, nesting, ami_max_terms) {
        nm <- names(x)
        if (!length(x) || is.null(nm) || any(!nzchar(nm)))
                stop("A list `x` must be named, e.g. list(regions = reg, local = loc).",
                     call. = FALSE)
        res <- lapply(x, function(obj) {
                sc <- if (inherits(obj, "pulse_local_types")) local_scope else "global"
                crosswalk_types(obj, external, label, units, id, weights, sc,
                                nesting, ami_max_terms)
        })
        get <- function(a, ix, b) a[[b]][a$index == ix]
        by_level <- do.call(rbind, lapply(nm, function(k) {
                a <- res[[k]]$agreement
                data.frame(level = k, n_types = nrow(res[[k]]$overlap$crisp),
                           n_units = res[[k]]$coverage$n_matched,
                           ARI = get(a, "ARI", "crisp"), AMI = get(a, "AMI", "crisp"),
                           NMI_crisp = get(a, "NMI", "crisp"), NMI_fuzzy = get(a, "NMI", "fuzzy"),
                           H_ext_given_pulse = get(a, "H_ext_given_pulse", "crisp"),
                           H_pulse_given_ext = get(a, "H_pulse_given_ext", "crisp"),
                           stringsAsFactors = FALSE)
        }))
        structure(list(results = res, by_level = by_level), class = "pulse_crosswalk_set")
}


# -----------------------------------------------------------------------------*
# Print
# -----------------------------------------------------------------------------*

#' @rdname crosswalk_types
#' @param ... Unused.
#' @export
print.pulse_crosswalk <- function(x, ...) {
        cv <- x$coverage
        sc <- if (identical(x$level, "local")) sprintf(" (scope: %s)", x$scope) else ""
        cat(sprintf("<pulse_crosswalk> level: %s%s | %d PULSE x %d external types\n",
                    x$level, sc, nrow(x$overlap$crisp), ncol(x$overlap$crisp)))
        cat(sprintf("Units: %d matched (%.1f%% of typed PULSE units); %d external unmatched; %d without a PULSE type\n",
                    cv$n_matched, 100 * cv$share_matched, cv$n_external_unmatched,
                    cv$n_no_pulse_type))
        cat(sprintf("Weights: %s | AMI: %s\n\n", x$settings$weights, x$settings$ami_method))
        a <- x$agreement
        print(data.frame(index = a$index, crisp = round(a$crisp, 3), fuzzy = round(a$fuzzy, 3)),
              row.names = FALSE)
        tp <- x$types$pulse[x$types$pulse$basis == "crisp", ]
        te <- x$types$external[x$types$external$basis == "crisp", ]
        tab <- function(v) paste(sprintf("%s %d", names(table(v)), table(v)), collapse = ", ")
        cat(sprintf("\nRelations (crisp, nesting >= %.2f)\n  PULSE types:    %s\n  external types: %s\n",
                    x$settings$nesting, tab(tp$relation), tab(te$relation)))
        cat("  'nested' PULSE types refine an external type; 'contains' PULSE types merge external types.\n")
        if (identical(x$level, "local") && identical(x$scope, "global"))
                cat("\nNote: global scope includes the level 1 signal. Use",
                    "local_scope = \"within_region\" to isolate level 2.\n")
        invisible(x)
}

#' @rdname crosswalk_types
#' @export
print.pulse_crosswalk_set <- function(x, ...) {
        cat(sprintf("<pulse_crosswalk_set> %d levels\n\n", length(x$results)))
        b <- x$by_level
        num <- vapply(b, is.double, logical(1))
        b[num] <- lapply(b[num], round, 3)
        print(b, row.names = FALSE)
        cat("\nCompare levels on AMI; NMI rises with the number of types.\n")
        invisible(x)
}
